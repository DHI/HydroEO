"""SWOT raster download and preprocessing.

Entry point: download_raster(config, project_dir, credentials)
"""

from __future__ import annotations

import logging
import re
import shutil
import tempfile
from pathlib import Path
from typing import Any

import geopandas as gpd
import numpy as np
import rasterio
import xarray as xr
import rioxarray  # noqa: F401 - required for xarray .rio accessor
from pyproj import CRS
from rasterio.merge import merge
from rasterio.warp import calculate_default_transform, reproject, Resampling
from shapely.geometry import box, mapping
from tqdm import tqdm

from HydroEO import FLOAT32_NODATA_VALUE
from HydroEO.satellites.swot._download import (
    _download_files,
    _filter_new,
    _login,
    _search,
)

logger = logging.getLogger(__name__)

DEFAULT_VARIABLES = [
    "wse",
    "wse_uncert",
    "wse_qual",
    "height_cor_xover",
    "geoid",
    "n_wse_pix",
    "n_other_pix",
    "layover_impact",
]

# Always extracted: used for quality masking and primary output
REQUIRED_VARIABLES = {"wse", "wse_uncert", "layover_impact"}


def download_raster(
    config: dict[str, Any],
    project_dir: str,
    credentials: tuple[str | None, str | None],
    global_crs: str = "EPSG:4326",
) -> None:
    """Download, preprocess, and merge SWOT raster data for an AOI.

    Parameters
    ----------
    config:
        The ``swot_raster`` section of the project YAML config.
    project_dir:
        Root project directory; outputs are written to
        ``<project_dir>/swot_raster/<aoi_name>/``.
    credentials:
        ``(earthdata_username, earthdata_password)`` tuple.
    global_crs:
        Project-level CRS (from ``gis.global_crs``). Used as the default
        output CRS for the merge phase unless ``target_crs`` is set in
        ``config``.
    """
    logger.debug(
        "SWOT Raster Download - AOI: %s, Product: %s, Temporal range: %s to %s",
        config["aoi"]["name"],
        config["product"],
        config["startdate"],
        config["enddate"],
    )

    aoi_name = config["aoi"]["name"]
    product = config["product"]
    base_dir = Path(project_dir) / "swot_raster" / aoi_name
    raw_dir = base_dir / "raw" / product
    processed_dir = base_dir / "processed" / product
    raw_dir.mkdir(parents=True, exist_ok=True)
    processed_dir.mkdir(parents=True, exist_ok=True)

    log_path = raw_dir / "downloaded.log"
    if not log_path.exists():
        log_path.touch()

    with open(log_path, "r") as log_file:
        processed_granules = {line.rstrip() for line in log_file}

    downloaded_files = _download_granules(
        config, raw_dir, processed_granules, credentials
    )

    # Always preprocess if there are NC files in raw_dir (new or from previous runs)
    nc_files_in_raw = list(raw_dir.glob("*.nc"))
    if nc_files_in_raw or downloaded_files:
        _preprocess_granules(config, raw_dir, processed_dir, log_path)
    else:
        logger.info("No new granules to download and no unprocessed netCDF files found")

    processed_files = list(processed_dir.glob("*.tif"))
    if processed_files:
        logger.info(
            "Found %d processed TIF files, proceeding with merge",
            len(processed_files),
        )
        _merge_and_reproject_granules(config, processed_dir, global_crs)
    else:
        logger.info("No processed TIF files found, skipping merge phase")

    logger.debug("SWOT raster download and processing complete")


def _download_granules(
    config: dict,
    raw_dir: Path,
    processed_granules: set[str],
    credentials: tuple[str | None, str | None],
) -> list[str]:
    """Download SWOT granules matching query, skipping already processed ones."""
    logger.debug("=== SWOT Raster Download Phase ===")

    _login(credentials)

    aoi_config = config["aoi"]
    if aoi_config["type"] == "bbox":
        bbox = aoi_config["bbox"]
        bounds = (bbox[0], bbox[1], bbox[2], bbox[3])
        logger.info("AOI Type: bbox, Bounds: %s", bounds)
    else:
        logger.error(
            "AOI type '%s' not yet implemented for download", aoi_config["type"]
        )
        return []

    startdate = config["startdate"]
    enddate = config["enddate"]
    time_start = f"{startdate[0]:04d}-{startdate[1]:02d}-{startdate[2]:02d} 00:00:00"
    time_end = f"{enddate[0]:04d}-{enddate[1]:02d}-{enddate[2]:02d} 23:59:59"
    logger.debug("Temporal range: %s to %s", time_start, time_end)

    granule_filter = config.get("granule_filter", None)
    if granule_filter:
        logger.debug("Granule filter: %s", granule_filter)

    product = config["product"]
    logger.info(
        "Searching for SWOT product '%s' in AOI %s", product, aoi_config["name"]
    )

    search_params = {
        "short_name": product,
        "bounding_box": bounds,
        "temporal": (time_start, time_end),
        "granule_name": granule_filter,
    }
    swot_results = _search(**search_params)
    logger.debug("Found %d granules matching query", len(swot_results))

    if not swot_results:
        logger.warning("No SWOT granules found matching the query")
        return []

    new_results = _filter_new(swot_results, processed_granules)
    logger.debug(
        "%d new granules to download (already processed: %d)",
        len(new_results),
        len(processed_granules),
    )

    if not new_results:
        logger.info("All granules already processed, skipping download")
        return []

    logger.debug("Downloading %d granules to %s", len(new_results), raw_dir)
    downloaded_files = _download_files(new_results, str(raw_dir))
    logger.debug("Successfully downloaded %d files", len(downloaded_files))
    return downloaded_files or []


def _preprocess_granules(
    config: dict,
    raw_dir: Path,
    processed_dir: Path,
    log_path: Path,
) -> None:
    """Extract layers from netCDF files, clip to AOI, and clean up raw files.

    Processes only NC files that haven't already produced output TIFFs.
    """
    logger.debug("=== SWOT Raster Preprocessing Phase ===")

    aoi_config = config["aoi"]

    variables = list(
        REQUIRED_VARIABLES | set(config.get("variables", DEFAULT_VARIABLES))
    )

    aoi_gdf = None
    if aoi_config["type"] in ["shapefile", "geopackage"]:
        aoi_path = aoi_config.get("path")
        if aoi_path and Path(aoi_path).exists():
            try:
                aoi_gdf = gpd.read_file(aoi_path)
                logger.debug("Loaded AOI from %s", aoi_path)
            except Exception as e:
                logger.warning("Could not load AOI file: %s", e)
    elif aoi_config["type"] == "bbox":
        bbox = aoi_config["bbox"]
        aoi_gdf = gpd.GeoDataFrame(geometry=[box(*bbox)], crs="EPSG:4326")
        logger.debug("Created AOI GeoDataFrame from bbox %s", bbox)

    all_nc_files = list(raw_dir.glob("*.nc"))

    # Filter: only process NC files that haven't produced output TIFFs yet
    nc_files = []
    already_processed = []
    for nc_path in all_nc_files:
        # Check if any output TIFFs exist for this NC file
        existing_tiffs = list(processed_dir.glob(f"{nc_path.stem}_*.tif"))
        if not existing_tiffs:
            nc_files.append(nc_path)
        else:
            already_processed.append(nc_path.name)

    if not nc_files:
        logger.info(
            "No unprocessed netCDF files found (all %d already processed)",
            len(all_nc_files),
        )
        return

    logger.debug(
        "Found %d unprocessed netCDF files to process (skipping %d already processed)",
        len(nc_files),
        len(already_processed),
    )

    deferred_warnings = []

    def _defer_warning(message, *args):
        if args:
            deferred_warnings.append(message % args)
        else:
            deferred_warnings.append(message)

    for nc_path in tqdm(nc_files, desc="Processing netCDF files"):
        try:
            ds = xr.open_dataset(nc_path)
        except Exception as e:
            _defer_warning("Cannot open %s: %s", nc_path.name, e)
            continue

        native_crs = _detect_crs(ds, nc_path)
        if native_crs is None:
            _defer_warning("Could not detect CRS for %s, skipping", nc_path.name)
            continue

        try:
            qf = config.get("quality_filters") or {}
            max_wse_uncert = qf.get("max_wse_uncert", 0.3)
            max_layover_impact = qf.get("max_layover_impact", 0.3)
            mask = (ds["wse_uncert"] < max_wse_uncert) & (
                ds["layover_impact"] < max_layover_impact
            )
        except (KeyError, TypeError):
            _defer_warning(
                "Missing required quality fields in %s, skipping", nc_path.name
            )
            continue

        if aoi_gdf is not None:
            try:
                aoi_reproj = aoi_gdf.to_crs(native_crs)
                xmin, ymin, xmax, ymax = aoi_reproj.geometry.total_bounds
                wse_filt = ds["wse"].where(mask)
                wse_filt.rio.write_crs(native_crs.to_wkt(), inplace=True)
                left, bottom, right, top = wse_filt.rio.bounds()

                if xmax < left or xmin > right or ymax < bottom or ymin > top:
                    continue
            except Exception as e:
                _defer_warning("Bounds check failed for %s: %s", nc_path.name, e)

        for var in tqdm(variables, desc=f"Variables in {nc_path.stem}", leave=False):
            if var not in ds:
                continue

            try:
                da = ds[var].where(mask)
                da.rio.write_crs(native_crs.to_wkt(), inplace=True)

                if aoi_gdf is not None:
                    try:
                        aoi_reproj = aoi_gdf.to_crs(native_crs)
                        geom = [mapping(aoi_reproj.unary_union)]
                        clipped = da.rio.clip(
                            geom, aoi_reproj.crs, drop=True, all_touched=False
                        )
                    except Exception as e:
                        logger.debug("Clip failed for '%s': %s", var, e)
                        clipped = da
                else:
                    clipped = da

                out_tif = processed_dir / f"{nc_path.stem}_{var}.tif"
                clipped.rio.to_raster(str(out_tif), driver="GTiff")
            except Exception as e:
                _defer_warning(
                    "Processing failed for '%s' in %s: %s", var, nc_path.name, e
                )

        ds.close()

    for warning in deferred_warnings:
        logger.debug(warning)

    nc_names = [f.stem for f in nc_files]
    with open(log_path, "a") as log_file:
        for name in nc_names:
            log_file.write(name + "\n")

    logger.debug("Cleaning up raw netCDF files")
    deferred_warnings = []
    for nc_path in tqdm(nc_files, desc="Cleaning up"):
        try:
            nc_path.unlink()
        except Exception as e:
            _defer_warning("Failed to delete %s: %s", nc_path.name, e)

    for warning in deferred_warnings:
        logger.debug(warning)


def _detect_crs(ds: xr.Dataset, nc_path: Path) -> CRS | None:
    """Detect CRS from netCDF metadata or filename."""
    variables = DEFAULT_VARIABLES

    for var in ("crs", "spatial_ref", "transverse_mercator"):
        if var in ds:
            try:
                crs_wkt = ds[var].attrs.get("crs_wkt") or ds[var].attrs.get(
                    "spatial_ref"
                )
                if crs_wkt:
                    return CRS.from_wkt(crs_wkt)
                epsg = ds[var].attrs.get("epsg_code") or ds[var].attrs.get("EPSG")
                if epsg:
                    return CRS.from_epsg(int(str(epsg).replace("EPSG:", "")))
            except Exception:
                pass

    for var in variables:
        if var in ds:
            try:
                detected = ds[var].rio.crs
                if detected:
                    return CRS.from_user_input(detected)
            except Exception:
                pass

    match = re.search(r"UTM(\d{1,2})([A-Z])", nc_path.stem, re.IGNORECASE)
    if match:
        zone = int(match.group(1))
        band = match.group(2).upper()
        hemisphere = "north" if band >= "N" else "south"
        return CRS.from_dict(
            {"proj": "utm", "zone": zone, "south": hemisphere == "south"}
        )

    return None


def _merge_and_reproject_granules(
    config: dict, processed_dir: Path, global_crs: str = "EPSG:4326"
) -> None:
    """Organize processed TIFFs by date and merge/reproject tiles."""
    logger.debug("=== SWOT Raster Merge and Reproject Phase ===")

    if not config.get("merge_tiles", True):
        logger.debug("merge_tiles is disabled in config, skipping merge phase")
        return

    # target_crs in swot_raster config is an optional override; falls back to global_crs
    crs_to_use = config.get("target_crs", global_crs)

    try:
        epsg_code = int(crs_to_use.replace("EPSG:", ""))
    except (ValueError, AttributeError):
        logger.error("Invalid CRS format: %s", crs_to_use)
        return

    merged_dir = processed_dir.parent.parent / "merged"
    merged_dir.mkdir(parents=True, exist_ok=True)

    tif_files = list(processed_dir.glob("*.tif"))
    if not tif_files:
        logger.warning("No TIF files found in %s", processed_dir)
        return

    logger.info("Found %d TIF files to merge", len(tif_files))

    file_groups: dict[tuple, list] = {}
    for tif_path in tif_files:
        match = re.search(r"(\d{8})T", tif_path.stem)
        if not match:
            logger.warning("Could not extract date from %s", tif_path.name)
            continue

        date_str = match.group(1)
        var_name = tif_path.stem.split("_")[-1]
        key = (var_name, date_str)
        file_groups.setdefault(key, []).append(tif_path)

    logger.debug("Grouped into %d variable-date combinations", len(file_groups))

    for (var_name, date_str), tif_list in tqdm(
        file_groups.items(), desc="Merging tiles"
    ):
        out_name = f"{date_str}_{var_name}_merged.tif"
        dst_path = merged_dir / out_name
        if len(tif_list) == 1:
            _reproject_raster(tif_list[0], dst_path, epsg_code, var_name)
        else:
            _merge_rasters(tif_list, dst_path, epsg_code, var_name)

    logger.info("Merge and reproject complete. Results in %s", merged_dir)


def _reproject_raster(src_path: Path, dst_path: Path, epsg: int, var_name: str) -> None:
    """Reproject a single raster to target EPSG."""
    try:
        with rasterio.open(src_path) as src:
            dst_crs = f"EPSG:{epsg}"
            transform, width, height = calculate_default_transform(
                src.crs, dst_crs, src.width, src.height, *src.bounds
            )
            profile = src.profile.copy()
            profile.update(
                {
                    "crs": dst_crs,
                    "transform": transform,
                    "width": width,
                    "height": height,
                }
            )
            nodata = profile.get("nodata", src.nodata)
            profile["nodata"] = nodata

            with rasterio.open(dst_path, "w", **profile) as dst:
                for b in range(1, src.count + 1):
                    reproject(
                        source=rasterio.band(src, b),
                        destination=rasterio.band(dst, b),
                        src_transform=src.transform,
                        src_crs=src.crs,
                        dst_transform=transform,
                        dst_crs=dst_crs,
                        resampling=_resampling_for(var_name),
                    )
    except Exception as e:
        logger.warning("Failed to reproject %s: %s", src_path.name, e)


def _merge_rasters(
    tif_list: list[Path], dst_path: Path, epsg: int, var_name: str
) -> None:
    """Merge multiple rasters and reproject to target EPSG."""
    tmpdir = tempfile.mkdtemp(prefix="swot_merge_")
    reproj_paths = []

    try:
        for i, tif_path in enumerate(tif_list):
            tmp_out = Path(tmpdir) / f"reproj_{i}.tif"
            _reproject_raster(tif_path, tmp_out, epsg, var_name)
            reproj_paths.append(tmp_out)

        datasets = [rasterio.open(p) for p in reproj_paths]
        try:
            mosaic, out_trans = merge(datasets)
            out_meta = datasets[0].meta.copy()
            out_meta.update(
                {
                    "driver": "GTiff",
                    "count": mosaic.shape[0],
                    "width": mosaic.shape[2],
                    "height": mosaic.shape[1],
                    "transform": out_trans,
                    "crs": f"EPSG:{epsg}",
                    "nodata": FLOAT32_NODATA_VALUE,
                }
            )
            mosaic[mosaic == 0] = FLOAT32_NODATA_VALUE
            mosaic[mosaic == np.nan] = FLOAT32_NODATA_VALUE
            with rasterio.open(dst_path, "w", **out_meta) as dest:
                dest.write(mosaic)
        finally:
            for ds in datasets:
                ds.close()

    except Exception as e:
        logger.warning("Failed to merge rasters: %s", e)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def _resampling_for(var_name: str) -> Resampling:
    """Choose resampling method based on variable type."""
    var_lower = var_name.lower()
    continuous_exact = {"wse", "wse_uncert", "geoid", "height_cor_xover"}
    if var_lower in continuous_exact:
        return Resampling.bilinear
    return Resampling.nearest
