"""
Per-target run_config persistence, exclusions, and spatial/reach-slope
correction caching -- shared by both reservoirs and rivers.

One YAML file per target ({output}/{id}/run_config.yaml) that is
simultaneously: (a) readable log of decisions made about this
target, (b) what _merge_engine._merge_timeseries
reads to apply those decisions, and (c) hand-editable
for a fully config-driven workflow.
"""

import logging
import os
import datetime
import json

import geopandas as gpd
import pandas as pd
import yaml

from HydroEO.utils import general
from HydroEO.utils.filters import basic_filters

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from HydroEO.project import Project

logger = logging.getLogger(__name__)


def _get_target_ids(prj: "Project", target_type: str):
    """Return the list of target IDs to process, for either 'reservoirs' or 'rivers'."""
    if target_type == "reservoirs":
        return list(prj.reservoirs.download_gdf[prj.reservoirs.id_key])
    if target_type == "rivers":
        return list(prj.rivers.target_ids)
    raise ValueError(f"Unknown target_type: {target_type!r}")


def _target_centroid(prj: "Project", target_type: str, id):
    """
    Return (lat, lon) of a target's own geometry centroid: 
    - the reservoir polygon for target_type="reservoirs"
    - the SWORD node/reach geometry for target_type="rivers"
    computed in a projected (local) CRS for accuracy, then converted back to lat/lon.
    Used as the reference location for apply_distance_penalty/
    apply_spatial_correction. 
    
    Returns (None, None) if the target's geometry can't be found.
    """
    try:
        if target_type == "reservoirs":
            gdf = prj.reservoirs.gdf
            id_key = prj.reservoirs.id_key
        elif target_type == "rivers":
            gdf = prj.rivers.target_features
            id_key = prj.rivers.target_id_col
        else:
            raise ValueError(f"Unknown target_type: {target_type!r}")

        row = gdf.loc[gdf[id_key] == id]
        if len(row) == 0 or row.geometry.isna().all():
            return None, None
        centroid = row.to_crs(prj.local_crs).geometry.centroid.to_crs(prj.global_crs)
        pt = centroid.iloc[0]
        return pt.y, pt.x  # lat, lon
    except Exception as exc:
        logger.warning(
            "Could not compute %s centroid for %s: %s", target_type, id, exc
        )
        return None, None


def _reservoir_centroid(prj: "Project", id):
    """Backward-compatible wrapper, see _target_centroid."""
    return _target_centroid(prj, "reservoirs", id)


def _get_or_fit_spatial_correction_model(
    prj: "Project", target_type: str, id, dense_source_platform="icesat2",
    recalibrate=False, **fit_kwargs,
):
    """
    Load existing spatial correction model if it exists or generates new
    model. Works identically for reservoirs and river targets, see _target_centroid.

    This is only re-fit is requested (recalibrate = True), to avoid 
    changing past corrections silently.

    Returns None if no model exists yet and there isn't enough dense
    source data to fit one (see fit_spatial_correction_model) 
    Equivalent to "no correction available."
    """
    model_path = os.path.join(
        prj.dirs["output"], f"{id}", "spatial_correction_model.json"
    )

    if os.path.exists(model_path) and not recalibrate:
        with open(model_path, "r") as f:
            return json.load(f)

    cleaned_path = os.path.join(
        prj.dirs["output"], f"{id}", "all_cleaned_timeseries.csv"
    )
    if not os.path.exists(cleaned_path):
        logger.info(
            "No cleaned observations yet for %s; cannot fit spatial "
            "correction model.", id,
        )
        return None

    df = pd.read_csv(cleaned_path)
    if "platform" not in df.columns or dense_source_platform not in df["platform"].values:
        logger.info(
            "No %s data available for %s; cannot fit spatial correction "
            "model from it.", dense_source_platform, id,
        )
        return None

    df["date"] = pd.to_datetime(df["date"])
    dense_df = df.loc[df["platform"] == dense_source_platform]

    ref_lat, ref_lon = _target_centroid(prj, target_type, id)
    if ref_lat is None:
        logger.warning(
            "Could not determine %s centroid for %s; cannot fit spatial "
            "correction model.", target_type, id,
        )
        return None

    model = basic_filters.fit_spatial_correction_model(
        dense_df, lat_key="lat", lon_key="lon", height_key="height",
        date_key="date", ref_lat=ref_lat, ref_lon=ref_lon, **fit_kwargs,
    )

    if model is not None:
        general.ifnotmakedirs(os.path.dirname(model_path))
        with open(model_path, "w") as f:
            json.dump(model, f, indent=2)

    return model


def _run_config_path(prj: "Project", id) -> str:
    return os.path.join(prj.dirs["output"], f"{id}", "run_config.yaml")


def _default_run_config(id) -> dict:
    return {
        "target_id": id,
        "last_updated": None,
        "merging_option_overrides": {},
        "exclusions": [],
    }


def _load_run_config(prj: "Project", id) -> dict:
    """Load a target's run_config.yaml, or a fresh default if none exists yet."""
    path = _run_config_path(prj, id)
    if os.path.exists(path):
        with open(path, "r") as f:
            loaded = yaml.safe_load(f)
        if loaded:
            # tolerate a hand-edited file missing a key or two
            defaults = _default_run_config(id)
            defaults.update(loaded)
            return defaults
    return _default_run_config(id)


def _save_run_config(prj: "Project", id, config: dict) -> None:
    config["last_updated"] = datetime.datetime.now().isoformat()
    path = _run_config_path(prj, id)
    general.ifnotmakedirs(os.path.dirname(path))
    with open(path, "w") as f:
        yaml.safe_dump(config, f, sort_keys=False)


def _invalidate_spatial_correction_cache(prj: "Project", id) -> None:
    """
    Delete cached spatial correction model for this target, forcing
    a fresh fit next time use_spatial_correction is used. 
    """
    model_path = os.path.join(prj.dirs["output"], f"{id}", "spatial_correction_model.json")
    if os.path.exists(model_path):
        os.remove(model_path)
        logger.info(
            "Invalidated cached spatial correction model for %s "
            "(exclusions or related options changed).", id,
        )


def _invalidate_reach_slope_correction_cache(prj: "Project", id) -> None:
    """
    Delete any cached reach slope correction model for this target.
    """
    model_path = os.path.join(
        prj.dirs["output"], f"{id}", "reach_slope_correction_model.json"
    )
    if os.path.exists(model_path):
        os.remove(model_path)
        logger.info(
            "Invalidated cached reach slope correction model for %s "
            "(exclusions or related options changed).", id,
        )


def _fit_reach_slope_correction(prj: "Project", target_id, recalibrate: bool = False):
    """
    Fit (or load a persisted) reach-level slope correction from SWOT's
    own directly-measured "slope" field (RiverSP reach product), to
    reference-correct OTHER missions' (ICESat-2/Sentinel-3/6) crossings
    to the reach's geometric midpoint.

    Returns None if no model exists yet and there's no usable SWOT slope
    data to fit one from, "no correction available," not an error.
    """
    model_path = os.path.join(
        prj.dirs["output"], f"{target_id}", "reach_slope_correction_model.json"
    )
    if os.path.exists(model_path) and not recalibrate:
        with open(model_path, "r") as f:
            return json.load(f)

    swot_path = os.path.join(
        prj.dirs["output"], f"{target_id}", "raw_observations", "swot.gpkg"
    )
    if not os.path.exists(swot_path):
        logger.info(
            "No raw SWOT observations for %s; cannot fit reach slope "
            "correction.", target_id,
        )
        return None

    swot_gdf = gpd.read_file(swot_path)
    if "slope" not in swot_gdf.columns:
        logger.warning(
            "No 'slope' field found in SWOT observations for %s -- was "
            "it requested in mission_options['swot']['hydrocron_fields']"
            "['reaches']? Cannot fit reach slope correction.", target_id,
        )
        return None

    valid_slopes = pd.to_numeric(swot_gdf["slope"], errors="coerce").dropna()
    if valid_slopes.empty:
        logger.info(
            "No valid (non-null, numeric) SWOT slope observations for "
            "%s; cannot fit reach slope correction.", target_id,
        )
        return None

    model = {
        "target_id": str(target_id),
        "median_slope": float(valid_slopes.median()),
        "n_observations": int(len(valid_slopes)),
        "fitted_at": datetime.datetime.now().isoformat(),
    }

    general.ifnotmakedirs(os.path.dirname(model_path))
    with open(model_path, "w") as f:
        json.dump(model, f, indent=2)

    return model


def _apply_reach_slope_correction(
    ts_df: pd.DataFrame, prj: "Project", target_id, slope_model: dict,
) -> pd.DataFrame:
    """
    Apply a fitted reach slope correction (see _fit_reach_slope_correction)
    to non-SWOT rows in ts_df.

    NOTE: the sign convention for "correction = slope x distance" 
    has NOT been empirically verified against real data.

    Rows without usable lat/lon (or if the target's geometry can't be
    found) are left uncorrected rather than dropped.
    """
    if "platform" not in ts_df.columns or "lat" not in ts_df.columns or "lon" not in ts_df.columns:
        return ts_df

    target_row = prj.rivers.target_features.loc[
        prj.rivers.target_features[prj.rivers.target_id_col] == target_id
    ]
    if target_row.empty:
        logger.warning(
            "Could not find reach geometry for %s; skipping reach slope "
            "correction.", target_id,
        )
        return ts_df

    reach_geom_global = target_row.geometry.iloc[0]
    local_crs = prj.local_crs
    reach_geom_local = (
        gpd.GeoSeries([reach_geom_global], crs=target_row.crs).to_crs(local_crs).iloc[0]
    )
    reach_midpoint_dist = reach_geom_local.length / 2.0
    median_slope = slope_model["median_slope"]

    mask = (ts_df["platform"] != "swot") & ts_df["lat"].notna() & ts_df["lon"].notna()
    if not mask.any():
        return ts_df

    points_local = gpd.GeoSeries(
        gpd.points_from_xy(ts_df.loc[mask, "lon"], ts_df.loc[mask, "lat"]),
        crs=prj.global_crs,
    ).to_crs(local_crs)

    along_reach_dist = points_local.apply(reach_geom_local.project)
    distance_from_midpoint = along_reach_dist.values - reach_midpoint_dist

    ts_df = ts_df.copy()
    ts_df.loc[mask, "height"] = (
        ts_df.loc[mask, "height"] - median_slope * distance_from_midpoint
    )
    return ts_df


def list_target_observations(prj: "Project", target_type: str, id) -> pd.DataFrame:
    """
    Summarize what observations exist for a target, at (platform, orbit)
    granularity -- the "what could I exclude" view. Meant to be read
    alongside plot_merging's platform-colored progress plots (which show
    WHERE a problem shows up in the actual data), not a replacement for
    looking at the data itself.
    """
    cleaned_path = os.path.join(prj.dirs["output"], f"{id}", "all_cleaned_timeseries.csv")
    if not os.path.exists(cleaned_path):
        logger.warning(
            "No cleaned observations yet for %s -- has create_%s_timeseries() "
            "been run?", id, target_type,
        )
        return pd.DataFrame(columns=["platform", "orbit", "n_points", "date_min", "date_max"])

    df = pd.read_csv(cleaned_path)
    df["date"] = pd.to_datetime(df["date"])
    group_cols = [c for c in ["platform", "orbit"] if c in df.columns]
    summary = (
        df.groupby(group_cols)
        .agg(n_points=("date", "size"), date_min=("date", "min"), date_max=("date", "max"))
        .reset_index()
        .sort_values(group_cols)
        .reset_index(drop=True)
    )
    return summary


def list_exclusions(prj: "Project", target_type: str, id) -> list:
    """Current exclusion rules for a target, from its run_config.yaml."""
    return _load_run_config(prj, id)["exclusions"]


def exclude_from_target(
    prj: "Project", target_type: str, id,
    platform=None, orbit=None, date=None, reason=None,
) -> None:
    """
    Exclude observations from a target's merge, at whatever granularity
    is given -- a whole platform, a specific orbit/pass value, a specific
    date, or any combination (all given fields must match for a row to
    be excluded). Persisted to {output}/{id}/run_config.yaml.

    Applied at the start of _merge_timeseries, before any processing --
    an excluded pass never reaches bias_correct/Kalman/svr_radial at all,
    rather than being fought against downstream.

    Invalidates any cached spatial correction model for this target,
    since it may have been fit using data that's now excluded.

    Examples
    --------
    exclude_from_target(prj, "reservoirs", my_id, platform="S3B")
    exclude_from_target(prj, "reservoirs", my_id, platform="S3B", orbit=1517)
    exclude_from_target(prj, "rivers", my_id, date="2024-03-19")
    """
    if platform is None and orbit is None and date is None:
        raise ValueError(
            "Specify at least one of platform, orbit, or date to exclude."
        )

    config = _load_run_config(prj, id)
    config["exclusions"].append({
        "platform": platform,
        "orbit": orbit,
        "date": str(date) if date is not None else None,
        "reason": reason,
        "added": datetime.datetime.now().isoformat(),
    })
    _save_run_config(prj, id, config)
    _invalidate_spatial_correction_cache(prj, id)
    _invalidate_reach_slope_correction_cache(prj, id)
    logger.info(
        "Added exclusion for %s: platform=%s orbit=%s date=%s (%s)",
        id, platform, orbit, date, reason or "no reason given",
    )


def _exclusion_value_matches(a, b) -> bool:
    """
    Robust equality check for exclusion matching -- tries numeric
    comparison first, so int/float/numeric-string representations of
    the same value all compare correctly (e.g. 1517 == 1517.0 == "1517"
    -- same int/float representation issue confirmed for
    _apply_exclusions' dataframe matching, applied here too for
    consistency), falling back to exact equality for non-numeric values
    (e.g. a string-based orbit identifier) or when either side is None.
    """
    if a is None or b is None:
        return a == b
    try:
        return float(a) == float(b)
    except (TypeError, ValueError):
        return a == b


def remove_exclusion(
    prj: "Project", target_type: str, id,
    index: int = None, platform=None, orbit=None, date=None,
) -> list:
    """
    Remove one or more exclusion rules, either by position (index, from
    list_exclusions()) or by matching criteria -- the same
    platform/orbit/date fields used to add one via exclude_from_target.
    Matching by criteria is usually more convenient than looking up an
    index first: e.g. remove_exclusion(prj, "reservoirs", my_id,
    platform="S3B", orbit=1517) removes exactly the rule that excluded
    that mission+orbit combination (or platform+beam, for ICESat-2 --
    beam values ARE the "orbit" field once concatenated with other
    missions, see PRODUCT_TIMESERIES_KEYS -- there's no separate "beam"
    parameter needed).

    Specify either index, or at least one of platform/orbit/date, not
    both. Criteria matching removes EVERY exclusion rule whose given
    fields match (fields not specified are ignored, not required to be
    None on the stored rule).

    Returns the list of removed rule(s), for confirmation/logging.
    """
    config = _load_run_config(prj, id)
    exclusions = config["exclusions"]
    criteria_given = platform is not None or orbit is not None or date is not None

    if index is not None and criteria_given:
        raise ValueError(
            "Specify either index OR platform/orbit/date criteria, not both."
        )

    if index is not None:
        if index < 0 or index >= len(exclusions):
            raise IndexError(
                f"No exclusion at index {index} for {id}; there are "
                f"{len(exclusions)}. See list_exclusions()."
            )
        removed = [exclusions.pop(index)]
    elif criteria_given:
        date_str = str(date) if date is not None else None
        to_remove = [
            rule for rule in exclusions
            if (platform is None or _exclusion_value_matches(rule.get("platform"), platform))
            and (orbit is None or _exclusion_value_matches(rule.get("orbit"), orbit))
            and (date is None or rule.get("date") == date_str)
        ]
        if not to_remove:
            raise ValueError(
                f"No exclusion found matching platform={platform!r} "
                f"orbit={orbit!r} date={date!r} for {id}. See list_exclusions()."
            )
        for rule in to_remove:
            exclusions.remove(rule)
        removed = to_remove
    else:
        raise ValueError(
            "Specify index, or at least one of platform/orbit/date, to "
            "identify which exclusion(s) to remove."
        )

    _save_run_config(prj, id, config)
    _invalidate_spatial_correction_cache(prj, id)
    _invalidate_reach_slope_correction_cache(prj, id)
    logger.info("Removed %d exclusion(s) for %s: %s", len(removed), id, removed)
    return removed


def set_merging_option(prj: "Project", target_type: str, id, **kwargs) -> None:
    """
    Override one or more merging_options for just this one target,
    persisted the same way as exclusions (highest-priority layer: these
    override prj.reservoirs/rivers.merging_options, which override the
    DEFAULT_*_MERGING_OPTIONS defaults).

    Example: set_merging_option(prj, "reservoirs", my_id, svr_radial_err=0.5)
    """
    config = _load_run_config(prj, id)
    config["merging_option_overrides"].update(kwargs)
    _save_run_config(prj, id, config)
    if "use_spatial_correction" in kwargs or "spatial_correction_dense_source" in kwargs:
        _invalidate_spatial_correction_cache(prj, id)
    if "use_reach_slope_correction" in kwargs:
        _invalidate_reach_slope_correction_cache(prj, id)
    logger.info("Updated merging options for %s: %s", id, kwargs)


def _apply_exclusions(df: pd.DataFrame, exclusions: list) -> pd.DataFrame:
    """
    Filter out rows matching any exclusion rule. Within one rule, every
    specified field (platform/orbit/date) must match for a row to be
    excluded by it; a row is dropped if it matches ANY rule.
    """
    if not exclusions:
        return df

    keep_mask = pd.Series(True, index=df.index)
    for rule in exclusions:
        rule_mask = pd.Series(True, index=df.index)
        if rule.get("platform") is not None:
            rule_mask &= df["platform"] == rule["platform"]
        if rule.get("orbit") is not None:
            if "orbit" in df.columns:
                # Compare numerically when possible, not as strings --
                # a real orbit column commonly gets upcast to float64 by
                # pandas the moment ANY value in it is missing (very
                # common in real satellite data), so a genuine orbit
                # value of 1517 reads as 1517.0 in the dataframe while a
                # YAML-loaded exclusion rule reads it as the plain int
                # 1517. Comparing as strings ("1517.0" vs "1517") then
                # silently matches nothing -- confirmed as a real,
                # reproducible bug, not a hypothetical one. Falls back
                # to string comparison only if the orbit value genuinely
                # isn't numeric (e.g. a string-based identifier).
                try:
                    target_orbit = float(rule["orbit"])
                    rule_mask &= (
                        pd.to_numeric(df["orbit"], errors="coerce") == target_orbit
                    )
                except (TypeError, ValueError):
                    rule_mask &= df["orbit"].astype(str) == str(rule["orbit"])
            else:
                rule_mask &= False
        if rule.get("date") is not None:
            rule_mask &= df["date"].dt.floor("D").astype(str) == str(rule["date"])
        keep_mask &= ~rule_mask

    return df.loc[keep_mask].reset_index(drop=True)


