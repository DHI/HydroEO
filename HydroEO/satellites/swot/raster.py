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

import earthaccess
import geopandas as gpd
import rasterio
import xarray as xr
import rioxarray  # noqa: F401 - required for xarray .rio accessor
from pyproj import CRS
from rasterio.merge import merge
from rasterio.warp import calculate_default_transform, reproject, Resampling
from shapely.geometry import mapping
from tqdm import tqdm

logger = logging.getLogger(__name__)


def download_raster(
    config: dict[str, Any],
    project_dir: str,
    credentials: tuple[str | None, str | None],
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
    """
    logger.info(
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
    if downloaded_files:
        _preprocess_granules(config, raw_dir, processed_dir, log_path, downloaded_files)
    else:
        logger.info("No new granules to download")

    processed_files = list(processed_dir.glob("*.tif"))
    if processed_files:
        logger.info(
            "Found %d processed TIF files, proceeding with merge",
            len(processed_files),
        )
        _merge_and_reproject_granules(config, processed_dir)
    else:
        logger.info("No processed TIF files found, skipping merge phase")

    logger.info("SWOT raster download and processing complete")


def _download_granules(
    config: dict,
    raw_dir: Path,
    processed_granules: set[str],
    credentials: tuple[str | None, str | None],
) -> list[str]:
    """Download SWOT granules matching query, skipping already processed ones."""
    logger.info("=== SWOT Raster Download Phase ===")

    username, password = credentials
    if username and password:
        try:
            earthaccess.login(strategy="environment", persist=True)
        except Exception:
            earthaccess.login()
    else:
        logger.warning("No Earthdata credentials provided, attempting anonymous login")
        earthaccess.login()

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
    logger.info("Temporal range: %s to %s", time_start, time_end)

    granule_filter = config.get("granule_filter", None)
    if granule_filter:
        logger.info("Granule filter: %s", granule_filter)

    product = config["product"]
    logger.info(
        "Searching for SWOT product '%s' in AOI %s", product, aoi_config["name"]
    )
    try:
        swot_results = earthaccess.search_data(
            short_name=product,
            bounding_box=bounds,
            temporal=(time_start, time_end),
            granule_name=granule_filter,
        )
        logger.info("Found %d granules matching query", len(swot_results))
    except Exception as e:
        logger.error("Error searching for SWOT data: %s", e)
        return []

    if not swot_results:
        logger.warning("No SWOT granules found matching the query")
        return []

    granule_ids = [result["umm"]["GranuleUR"] for result in swot_results]
    new_granules = [gid for gid in granule_ids if gid not in processed_granules]
    logger.info(
        "%d new granules to download (already processed: %d)",
        len(new_granules),
        len(processed_granules),
    )

    if not new_granules:
        logger.info("All granules already processed, skipping download")
        return []

    new_results = [r for r in swot_results if r["umm"]["GranuleUR"] in new_granules]
    logger.info("Downloading %d granules to %s", len(new_results), raw_dir)
    try:
        downloaded_files = earthaccess.download(new_results, str(raw_dir))
        logger.info("Successfully downloaded %d files", len(downloaded_files))
        return downloaded_files or []
    except Exception as e:
        logger.error("Error downloading SWOT data: %s", e)
        return []


def _preprocess_granules(
    config: dict,
    raw_dir: Path,
    processed_dir: Path,
    log_path: Path,
    downloaded_files: list[str],
) -> None:
    """Extract layers from netCDF files, clip to AOI, and clean up raw files."""
    logger.info("=== SWOT Raster Preprocessing Phase ===")

    aoi_config = config["aoi"]

    variables = [
        "wse",
        "wse_uncert",
        "wse_qual",
        "height_cor_xover",
        "geoid",
        "n_wse_pix",
        "n_other_pix",
        "layover_impact",
    ]

    aoi_gdf = None
    if aoi_config["type"] in ["shapefile", "geopackage"]:
        aoi_path = aoi_config.get("path")
        if aoi_path and Path(aoi_path).exists():
            try:
                aoi_gdf = gpd.read_file(aoi_path)
                logger.info("Loaded AOI from %s", aoi_path)
            except Exception as e:
                logger.warning("Could not load AOI file: %s", e)

    nc_files = list(raw_dir.glob("*.nc"))
    logger.info("Found %d netCDF files to process", len(nc_files))

    for nc_path in tqdm(nc_files, desc="Processing netCDF files"):
        try:
            ds = xr.open_dataset(nc_path)
        except Exception as e:
            logger.warning("Cannot open %s: %s", nc_path.name, e)
            continue

        native_crs = _detect_crs(ds, nc_path)
        if native_crs is None:
            logger.warning("Could not detect CRS for %s, skipping", nc_path.name)
            continue

        try:
            mask = (ds["wse_uncert"] < 0.3) & (ds["layover_impact"] < 0.3)
        except (KeyError, TypeError):
            logger.warning(
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
                logger.warning("Bounds check failed for %s: %s", nc_path.name, e)

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
                logger.warning(
                    "Processing failed for '%s' in %s: %s", var, nc_path.name, e
                )

        ds.close()

    nc_names = [f.stem for f in nc_files]
    with open(log_path, "a") as log_file:
        for name in nc_names:
            log_file.write(name + "\n")

    logger.info("Cleaning up raw netCDF files")
    for nc_path in tqdm(nc_files, desc="Cleaning up"):
        try:
            nc_path.unlink()
        except Exception as e:
            logger.warning("Failed to delete %s: %s", nc_path.name, e)


def _detect_crs(ds: xr.Dataset, nc_path: Path) -> CRS | None:
    """Detect CRS from netCDF metadata or filename."""
    variables = [
        "wse",
        "wse_uncert",
        "wse_qual",
        "height_cor_xover",
        "geoid",
        "n_wse_pix",
        "n_other_pix",
        "layover_impact",
    ]

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


def _merge_and_reproject_granules(config: dict, processed_dir: Path) -> None:
    """Organize processed TIFFs by date and merge/reproject tiles."""
    logger.info("=== SWOT Raster Merge and Reproject Phase ===")

    target_crs = config.get("target_crs", None)
    if not target_crs:
        logger.warning("No target_crs specified in config, skipping merge phase")
        return

    try:
        epsg_code = int(target_crs.replace("EPSG:", ""))
    except (ValueError, AttributeError):
        logger.error("Invalid target_crs format: %s", target_crs)
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

    logger.info("Grouped into %d variable-date combinations", len(file_groups))

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
                }
            )
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
