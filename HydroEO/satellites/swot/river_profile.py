"""Longitudinal river WSE profile extraction from SWOT L2 HR Raster tiles.

Entry point: calculate_river_profile(config, project_dir, global_crs)

Given a user-supplied "chainage" shapefile (points along a river with a
precomputed along-river distance field) and a date range, this:

1. Downloads and preprocesses SWOT L2 HR Raster WSE tiles covering the
   chainage shapefile's extent, reusing
   :func:`HydroEO.satellites.swot.raster.download_raster` with
   ``merge_tiles=False`` so each granule stays as its own clipped,
   quality-masked GeoTIFF (one per acquisition date/pass).
2. Samples WSE at each chainage point for every tile, and runs a
   configurable multi-stage filtering pipeline (see
   :mod:`HydroEO.utils.filters.river_profile_filters`) to produce a cleaned
   longitudinal profile per date.
3. Writes final profile shapefiles + a combined CSV (and, if
   ``keep_intermediates`` is set, every intermediate stage too), optional
   diagnostic plots, and a quality report summarizing how many points each
   stage touched.
"""

from __future__ import annotations

import copy
import logging
import re
import shutil
from pathlib import Path
from typing import Any

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import rowcol
from tqdm import tqdm

from HydroEO.constants import (
    RIVER_PROFILE_DEFAULT_AOI_BUFFER_M,
    RIVER_PROFILE_DEFAULT_FILTERS,
)
from HydroEO.satellites.swot.raster import download_raster
from HydroEO.utils.filters import river_profile_filters as rpf

logger = logging.getLogger(__name__)


def calculate_river_profile(
    config: dict[str, Any],
    project_dir: str,
    global_crs: str = "EPSG:4326",
) -> None:
    """Download SWOT rasters and compute a longitudinal WSE profile.

    Parameters
    ----------
    config:
        The ``river_profile`` section of the project YAML config.
    project_dir:
        Root project directory; outputs are written to
        ``<project_dir>/results/<name>/``.
    global_crs:
        Project-level CRS (from ``gis.global_crs``). Used when the AOI
        search bbox needs a geographic CRS for the CMR granule search.
    """
    name = config.get("name", "river_profile")
    logger.debug(
        "River Profile - name: %s, chainage: %s, temporal range: %s to %s",
        name,
        config["chainage_path"],
        config["startdate"],
        config["enddate"],
    )

    results_dir = Path(project_dir) / "results" / name
    results_dir.mkdir(parents=True, exist_ok=True)

    gdf, x = _load_chainage(config)

    swot_raster_cfg = _build_swot_raster_config(config, name, gdf, global_crs)
    download_raster(config=swot_raster_cfg, project_dir=project_dir, global_crs=global_crs)

    processed_dir = (
        Path(project_dir)
        / "processed"
        / "swot_raster"
        / name
        / swot_raster_cfg["product"]
    )
    wse_files = sorted(processed_dir.glob("*_wse.tif"))
    if not wse_files:
        logger.warning(
            "No processed WSE tiles found in %s; nothing to profile", processed_dir
        )
        return

    variables = config.get("variables", ["wse"])
    geoid_files = sorted(processed_dir.glob("*_geoid.tif")) if "geoid" in variables else []

    filters = _resolve_filters(config.get("filters"))
    orbit_exclusions = config.get("orbit_exclusions") or []
    zero_is_nodata = config.get("zero_is_nodata", True)
    keep_intermediates = config.get("keep_intermediates", False)
    plot_enable = config.get("plot_enable", True)
    plot_dpi = config.get("plot_dpi", 150)
    plot_ylim = config.get("plot_ylim")

    # results/<name> is regenerated fresh each run (raw/processed are left
    # alone above to preserve SWOT download dedup).
    if results_dir.exists():
        shutil.rmtree(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    out_dirs = _make_output_dirs(results_dir, keep_intermediates, plot_enable)

    labels: list[str] = []
    raw_cols, pref_cols, h1_cols, final_cols = {}, {}, {}, {}
    quality_rows: list[dict[str, Any]] = []
    seen_labels: dict[str, int] = {}

    # deferred so interleaved logging doesn't garble the tqdm bar below
    deferred_messages: list[tuple[int, str]] = []

    def _defer(level, message, *args):
        deferred_messages.append((level, message % args if args else message))

    for i, wse_path in enumerate(tqdm(wse_files, desc="Processing WSE tiles"), 1):
        try:
            record = _process_one_date(
                wse_path=wse_path,
                gdf=gdf,
                x=x,
                filters=filters,
                orbit_exclusions=orbit_exclusions,
                zero_is_nodata=zero_is_nodata,
            )
        except Exception as e:
            _defer(
                logging.WARNING, "[WSE %d/%d] ERROR on %s: %s", i, len(wse_files), wse_path.name, e
            )
            continue

        if record is None:
            _defer(
                logging.INFO,
                "[WSE %d/%d] %s | no overlap with chainage -> skipped",
                i,
                len(wse_files),
                wse_path.name,
            )
            continue

        label = _unique_label(record["label"], seen_labels)
        labels.append(label)
        raw_cols[label] = record["y_raw"]
        pref_cols[label] = record["y_pref"]
        h1_cols[label] = record["y_h1"]
        final_cols[label] = record["y_final"]
        record["quality"]["label"] = label
        quality_rows.append(record["quality"])

        _save_profile_shp(
            gdf, record["y_final"], out_dirs["final"], name, label, suffix="profile"
        )
        if keep_intermediates:
            _save_profile_shp(
                gdf, record["y_raw"], out_dirs["raw"], name, label, suffix="profile_raw"
            )
            _save_profile_shp(
                gdf, record["y_pref"], out_dirs["prefilter"], name, label, suffix="profile"
            )
            _save_profile_shp(
                gdf, record["y_h1"], out_dirs["hampel1"], name, label, suffix="profile"
            )

        if plot_enable:
            _plot_per_profile(
                x,
                record,
                label=label,
                out_dir=out_dirs["plots_per"],
                ylim=plot_ylim,
                dpi=plot_dpi,
                river_name=name,
            )

        _defer(
            logging.INFO,
            "[WSE %d/%d] %s | finite raw=%d, final=%d | stage changes: "
            "soft-clamp=%d, hampel1=%d, density-cull=%d, hampel2=%d, "
            "orbit-excluded=%d",
            i,
            len(wse_files),
            wse_path.name,
            record["quality"]["n_finite_raw"],
            record["quality"]["n_finite_final"],
            record["quality"]["n_soft_clamp_changed"],
            record["quality"]["n_hampel1_flagged"],
            record["quality"]["n_density_dropped"],
            record["quality"]["n_hampel2_flagged"],
            record["quality"]["n_orbit_excluded"],
        )

    for level, message in deferred_messages:
        logger.log(level, message)

    if not labels:
        logger.warning("River profile '%s': no dates processed, no outputs written", name)
        return

    _save_combined_csv(x, labels, final_cols, results_dir / f"{name}_profiles_final.csv")
    if keep_intermediates:
        _save_combined_csv(x, labels, raw_cols, results_dir / f"{name}_profiles_raw.csv")
        _save_combined_csv(
            x, labels, pref_cols, results_dir / f"{name}_profiles_prefilter.csv"
        )
        _save_combined_csv(
            x, labels, h1_cols, results_dir / f"{name}_profiles_hampel1.csv"
        )

    if plot_enable:
        _plot_combined(
            x,
            final_cols,
            f"{name} — all dates (final profile)",
            out_dirs["plots_combined"] / f"{name}_combined_final.png",
            ylim=plot_ylim,
            dpi=plot_dpi,
        )

    _write_quality_report(quality_rows, results_dir / "quality_report.csv")

    if geoid_files:
        _process_geoid(geoid_files, gdf, results_dir, name)

    logger.info("River profile '%s': done. %d dates processed.", name, len(labels))


def _resolve_filters(user_filters: dict | None) -> dict[str, dict[str, Any]]:
    """Deep-merge user ``filters`` config over the stage defaults."""
    resolved: dict[str, dict[str, Any]] = copy.deepcopy(RIVER_PROFILE_DEFAULT_FILTERS)
    for stage, overrides in (user_filters or {}).items():
        if stage not in resolved or not isinstance(overrides, dict):
            continue
        resolved[stage].update(overrides)
    return resolved


def _build_swot_raster_config(
    config: dict[str, Any], name: str, gdf: gpd.GeoDataFrame, global_crs: str
) -> dict[str, Any]:
    """Derive an internal swot_raster-shaped config from the chainage AOI."""
    buffer_m = config.get("aoi_buffer_meters", RIVER_PROFILE_DEFAULT_AOI_BUFFER_M)
    gdf_4326 = gdf.to_crs("EPSG:4326")
    if buffer_m:
        utm_crs = gdf_4326.estimate_utm_crs()
        bbox = gdf_4326.to_crs(utm_crs).buffer(buffer_m).to_crs("EPSG:4326").total_bounds
    else:
        # buffering by 0 produces empty (NaN-bounds) geometries
        bbox = gdf_4326.total_bounds

    variables = list(config.get("variables", ["wse"]))
    if "wse" not in variables:
        variables = ["wse"] + variables

    return {
        "aoi": {"name": name, "type": "bbox", "bbox": [float(v) for v in bbox]},
        "product": config.get("product", "SWOT_L2_HR_Raster_D"),
        "startdate": config["startdate"],
        "enddate": config["enddate"],
        "variables": variables,
        "merge_tiles": False,
        "quality_filters": config.get("quality_filters"),
        "granule_filter": config.get("granule_filter"),
    }


def _load_chainage(config: dict[str, Any]) -> tuple[gpd.GeoDataFrame, np.ndarray]:
    """Load the chainage shapefile and resolve the along-river distance array."""
    path = config["chainage_path"]
    gdf = gpd.read_file(path)
    if gdf.crs is None:
        raise ValueError(
            f"Chainage file '{path}' has no CRS defined (missing .prj?). Set the "
            "file's CRS explicitly before using it with HydroEO."
        )

    field = config.get("chainage_field", "cngmeters")
    if field not in gdf.columns:
        raise ValueError(
            f"Chainage field '{field}' not found in '{path}' (columns: "
            f"{list(gdf.columns)}). Set 'river_profile.chainage_field' to the "
            "column holding along-river distance in metres."
        )

    x_raw = gdf[field].to_numpy(dtype=float)
    n_non_finite = int(np.count_nonzero(~np.isfinite(x_raw)))
    if n_non_finite:
        raise ValueError(
            f"Chainage field '{field}' in '{path}' has {n_non_finite} non-finite "
            "(NaN/inf) value(s). Every point needs a valid along-river distance; "
            "fix or drop those rows before using this file with HydroEO."
        )

    if config.get("reverse_chainage", False):
        x = float(np.nanmax(x_raw)) - x_raw
    else:
        x = x_raw

    order = np.argsort(x)
    gdf = gdf.iloc[order].reset_index(drop=True)
    x = x[order]
    gdf[field] = x  # keep in sync with x (may have been reversed above)
    logger.info(
        "Chainage '%s': %d points, range %.1f - %.1f m", field, len(gdf), x.min(), x.max()
    )
    return gdf, x


def _label_from_name(path: Path) -> str:
    """Extract a timestamp label from a SWOT granule filename."""
    m = re.search(r"\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}", path.name)
    if m:
        return m.group(0)
    m2 = re.search(r"\d{8}T\d{6}", path.name)
    if m2:
        return m2.group(0)
    return path.stem


def _unique_label(label: str, seen_labels: dict[str, int]) -> str:
    """Disambiguate labels that collide (e.g. two tiles from the same pass,
    split across a UTM zone boundary, sharing the same start timestamp) so
    one doesn't silently overwrite the other's CSV column/shapefile."""
    count = seen_labels.get(label, 0)
    seen_labels[label] = count + 1
    return label if count == 0 else f"{label}_{count}"


def _read_wse_raster(path: Path, zero_is_nodata: bool):
    with rasterio.open(path) as ds:
        arr = ds.read(1).astype(np.float32)
        mask = np.zeros(arr.shape, bool)
        nodata = ds.nodata
        if nodata is not None and not (isinstance(nodata, float) and np.isnan(nodata)):
            mask |= arr == nodata
        if zero_is_nodata:
            mask |= arr == 0
        if mask.any():
            arr = arr.copy()
            arr[mask] = np.nan
        return arr, ds.transform, ds.crs


def _sample_profile(arr: np.ndarray, transform, xs: np.ndarray, ys: np.ndarray) -> np.ndarray:
    """Sample raster values at point locations (already in the raster's CRS)."""
    rows, cols = rowcol(transform, xs, ys, op=np.round)
    rows = np.asarray(rows, dtype=int)
    cols = np.asarray(cols, dtype=int)
    H, W = arr.shape
    vals = np.full(xs.shape, np.nan, np.float32)
    inside = (rows >= 0) & (rows < H) & (cols >= 0) & (cols < W)
    vals[inside] = arr[rows[inside], cols[inside]]
    return vals


def _apply_orbit_exclusions(
    filename: str, x: np.ndarray, y: np.ndarray, orbit_exclusions: list[dict]
) -> tuple[np.ndarray, int]:
    """
    Mask chainage points for known bad-orbit ranges. For each rule whose
    ``orbit`` substring appears in ``filename``, points with
    ``min_chainage_m <= x <= max_chainage_m`` (either bound optional) are
    set to NaN. Generalizes a per-river manual correction.
    """
    y = y.copy()
    n_excluded = 0
    for rule in orbit_exclusions:
        orbit = str(rule.get("orbit", ""))
        if not orbit or orbit not in filename:
            continue
        lo = rule.get("min_chainage_m")
        hi = rule.get("max_chainage_m")
        exclude = np.ones_like(x, dtype=bool)
        if lo is not None:
            exclude &= x >= lo
        if hi is not None:
            exclude &= x <= hi
        n_excluded += int(np.count_nonzero(exclude & np.isfinite(y)))
        y[exclude] = np.nan
    return y, n_excluded


def _process_one_date(
    wse_path: Path,
    gdf: gpd.GeoDataFrame,
    x: np.ndarray,
    filters: dict,
    orbit_exclusions: list[dict],
    zero_is_nodata: bool,
) -> dict[str, Any] | None:
    arr, transform, raster_crs = _read_wse_raster(wse_path, zero_is_nodata)

    points = gdf.to_crs(raster_crs)
    xs = points.geometry.x.to_numpy()
    ys = points.geometry.y.to_numpy()
    y_raw = _sample_profile(arr, transform, xs, ys)

    n_finite_raw = int(np.count_nonzero(np.isfinite(y_raw)))
    if n_finite_raw == 0:
        return None

    label = _label_from_name(wse_path)

    # left to run through to an all-NaN result rather than returning None,
    # so a fully-excluded date still gets a quality_report row
    y_raw, n_excluded = _apply_orbit_exclusions(wse_path.name, x, y_raw, orbit_exclusions)

    y = y_raw.copy()
    pc = filters["preclip"]
    if pc.get("enabled", True):
        if pc.get("min") is not None:
            y[y < pc["min"]] = np.nan
        if pc.get("max") is not None:
            y[y > pc["max"]] = np.nan
    y_clip = y

    sc = filters["soft_clamp"]
    if sc.get("enabled", True):
        y_pref = rpf.vertical_density_soft_clamp(
            x,
            y_clip,
            bin_width=sc["bin_width_m"],
            y_bin=sc["y_bin_m"],
            fixed_halfw=sc["fixed_halfw_m"],
            mad_factor=sc["mad_factor"],
            min_count=sc["min_count"],
            min_mode_count=sc["min_mode_count"],
            slope_gain=sc["slope_gain"],
            huber_k=sc["huber_k"],
        )
    else:
        y_pref = y_clip.copy()
    both_finite = np.isfinite(y_clip) & np.isfinite(y_pref)
    n_soft = int(np.count_nonzero(both_finite & (np.abs(y_pref - y_clip) > 1e-6)))

    h1 = filters["hampel_1"]
    if h1.get("enabled", True):
        y_h1, ch1 = rpf.hampel_1d_meter(
            x,
            y_pref,
            win_m=h1["win_m"],
            k=h1["sigma"],
            action=h1["action"],
            min_valid=h1["min_valid"],
            replace_mode=h1["replace_mode"],
            huber_k=h1["huber_k"],
            huber_iters=h1["huber_iters"],
        )
    else:
        y_h1, ch1 = y_pref.copy(), np.zeros_like(y_pref, bool)
    n_h1 = int(np.count_nonzero(ch1))

    rq = filters["rolling_quantile"]
    if rq.get("enabled", True):
        y_rq = rpf.rolling_quantile_trend(
            x,
            y_h1,
            win_m=rq["win_m"],
            q=rq["q"],
            min_valid=rq["min_valid"],
            robust_iters=rq["robust_iters"],
            huber_k=rq["huber_k"],
        )
    else:
        y_rq = y_h1.copy()

    dc = filters["density_cull"]
    counts = rpf.rolling_density(x, y_h1, total_win_m=dc["total_win_m"])
    if dc.get("enabled", True):
        usable = np.isfinite(y_rq)
        dens_vals = counts[usable]
        thr = np.percentile(dens_vals, dc["low_pct"]) if dens_vals.size else np.inf
        # strict < avoids culling everything when density is uniform
        low_density = (counts < thr) | (counts < float(dc["abs_min"]))
        low_density = rpf.dilate_mask(low_density, k=int(dc["dilate"]))
        y_after_density = y_rq.copy()
        y_after_density[low_density] = np.nan
        n_dens = int(np.count_nonzero(low_density & np.isfinite(y_rq)))
    else:
        y_after_density = y_rq.copy()
        n_dens = 0

    h2 = filters["hampel_2"]
    if h2.get("enabled", True):
        y_h2, ch2 = rpf.hampel_1d_meter(
            x,
            y_after_density,
            win_m=h2["win_m"],
            k=h2["sigma"],
            action=h2["action"],
            min_valid=h2["min_valid"],
            replace_mode=h2["replace_mode"],
            huber_k=h2["huber_k"],
            huber_iters=h2["huber_iters"],
        )
    else:
        y_h2, ch2 = y_after_density.copy(), np.zeros_like(y_after_density, bool)
    n_h2 = int(np.count_nonzero(ch2))

    sf = filters["spline_fill"]
    if sf.get("enabled", True):
        y_final, y_spline = rpf.lsq_spline_fill(
            x,
            y_h2,
            k=sf["k"],
            counts=counts,
            inner_knots=sf.get("inner_knots"),
            target_spacing=sf["target_spacing_m"],
            weight_scheme=sf["weight_scheme"],
        )
    else:
        y_final, y_spline = y_h2.copy(), y_h2.copy()

    quality = {
        "label": label,
        "file": wse_path.name,
        "n_finite_raw": n_finite_raw,
        "n_orbit_excluded": n_excluded,
        "n_soft_clamp_changed": n_soft,
        "n_hampel1_flagged": n_h1,
        "n_density_dropped": n_dens,
        "n_hampel2_flagged": n_h2,
        "n_finite_final": int(np.count_nonzero(np.isfinite(y_spline))),
    }

    return {
        "label": label,
        "y_raw": y_raw,
        "y_clip": y_clip,
        "y_pref": y_pref,
        "y_h1": y_h1,
        "y_rq": y_rq,
        "y_after_density": y_after_density,
        "y_h2": y_h2,
        "y_final": y_final,
        "y_spline": y_spline,
        "quality": quality,
    }


def _process_geoid(
    geoid_files: list[Path],
    gdf: gpd.GeoDataFrame,
    results_dir: Path,
    name: str,
) -> None:
    """Save per-date geoid profiles. Independent of ``keep_intermediates``:
    geoid is an explicitly requested output variable (via
    ``river_profile.variables``), not an intermediate filtering stage."""
    out_dir = results_dir / "profiles_geoid"
    seen_labels: dict[str, int] = {}
    for f in geoid_files:
        try:
            # geoid can legitimately be 0, unlike WSE
            arr, transform, raster_crs = _read_wse_raster(f, zero_is_nodata=False)
            points = gdf.to_crs(raster_crs)
            vals = _sample_profile(
                arr, transform, points.geometry.x.to_numpy(), points.geometry.y.to_numpy()
            )
            if np.count_nonzero(np.isfinite(vals)) == 0:
                continue
            label = _unique_label(_label_from_name(f), seen_labels)
            _save_profile_shp(
                gdf, vals, out_dir, name, label, suffix="profile_geoid", field_name="GEOID"
            )
        except Exception as e:
            logger.warning("Failed to process geoid tile %s: %s", f.name, e)


def _make_output_dirs(results_dir: Path, keep_intermediates: bool, plot_enable: bool) -> dict:
    dirs = {"final": results_dir / "profiles_final"}
    dirs["final"].mkdir(parents=True, exist_ok=True)
    if keep_intermediates:
        dirs["raw"] = results_dir / "profiles_raw"
        dirs["prefilter"] = results_dir / "profiles_prefilter"
        dirs["hampel1"] = results_dir / "profiles_hampel1"
        for d in ("raw", "prefilter", "hampel1"):
            dirs[d].mkdir(parents=True, exist_ok=True)
    if plot_enable:
        dirs["plots_per"] = results_dir / "plots" / "per_profile"
        dirs["plots_combined"] = results_dir / "plots" / "combined"
        dirs["plots_per"].mkdir(parents=True, exist_ok=True)
        dirs["plots_combined"].mkdir(parents=True, exist_ok=True)
    return dirs


def _save_profile_shp(
    gdf_base: gpd.GeoDataFrame,
    values: np.ndarray,
    out_dir: Path,
    name_prefix: str,
    label: str,
    suffix: str,
    field_name: str = "WSE",
) -> None:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    gdf = gdf_base.copy()
    gdf[field_name] = values
    out_name = f"{name_prefix}_{label}_{suffix}.shp"
    gdf.to_file(out_dir / out_name)


def _save_combined_csv(
    x: np.ndarray, labels: list[str], columns_dict: dict[str, np.ndarray], out_path: Path
) -> None:
    df = pd.DataFrame({"distance_along_river_m": x})
    for lab in labels:
        df[lab] = columns_dict[lab]
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)


def _write_quality_report(quality_rows: list[dict[str, Any]], out_path: Path) -> None:
    df = pd.DataFrame(quality_rows)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)

    totals = df[
        [
            "n_orbit_excluded",
            "n_soft_clamp_changed",
            "n_hampel1_flagged",
            "n_density_dropped",
            "n_hampel2_flagged",
        ]
    ].sum()
    logger.info(
        "River profile quality report (%d dates): orbit-excluded=%d, "
        "soft-clamp changed=%d, Hampel(1) flagged=%d, density-cull dropped=%d, "
        "Hampel(2) flagged=%d. Full per-date report: %s",
        len(df),
        int(totals["n_orbit_excluded"]),
        int(totals["n_soft_clamp_changed"]),
        int(totals["n_hampel1_flagged"]),
        int(totals["n_density_dropped"]),
        int(totals["n_hampel2_flagged"]),
        out_path,
    )


def _plot_per_profile(
    x: np.ndarray,
    record: dict[str, Any],
    label: str,
    out_dir: Path,
    ylim,
    dpi: int,
    river_name: str,
) -> None:
    import matplotlib.pyplot as plt

    plt.switch_backend("Agg")  # always saved to disk, never shown interactively

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(12, 6))
    if ylim is not None:
        plt.ylim(ylim)
    plt.plot(x, record["y_raw"], ".", alpha=0.25, label="RAW (unfiltered)")
    plt.plot(x, record["y_pref"], ".", alpha=0.35, label="Prefilter (soft-clamp)")
    plt.plot(x, record["y_h1"], ".", alpha=0.45, label="Hampel(1)")
    plt.plot(x, record["y_rq"], ".", alpha=0.70, label="Rolling-Quantile (local)")
    plt.plot(x, record["y_after_density"], ".", alpha=0.70, label="RQ + density cull")
    plt.plot(x, record["y_h2"], ".", alpha=0.70, label="Hampel(2)")
    plt.plot(x, record["y_final"], "-", lw=2.0, alpha=0.95, label="Final (spline-filled)")
    plt.plot(x, record["y_spline"], "x", ms=3, alpha=0.5, label="Spline eval (guide)")
    plt.grid(True, alpha=0.3)
    plt.xlabel("Distance along the river [m]")
    plt.ylabel("WSE")
    plt.title(f"WSE profile — {river_name} — {label}")
    plt.legend(ncol=2)
    plt.tight_layout()
    plt.savefig(out_dir / f"{label}.png", dpi=dpi, bbox_inches="tight")
    plt.close()


_COMBINED_PLOT_MAX_LEGEND_ENTRIES = 20
"""Above this many dates, a per-date legend stops being readable (one entry
per line) and starts dominating/squashing the figure, so it's dropped."""


def _plot_combined(
    x: np.ndarray, series_dict: dict[str, np.ndarray], title: str, out_path: Path, ylim, dpi: int
) -> None:
    import matplotlib.pyplot as plt

    plt.switch_backend("Agg")

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n_dates = len(series_dict)
    plt.figure(figsize=(12, 6))
    if ylim is not None:
        plt.ylim(ylim)
    for lab, y in series_dict.items():
        plt.plot(x, y, ".", alpha=0.8, label=lab)
    plt.grid(True, alpha=0.3)
    plt.xlabel("Distance along the river [m]")
    plt.ylabel("WSE")
    plt.title(f"{title} ({n_dates} dates)")
    if n_dates <= _COMBINED_PLOT_MAX_LEGEND_ENTRIES:
        plt.legend(ncol=2, fontsize="small")
    plt.tight_layout()
    plt.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close()
