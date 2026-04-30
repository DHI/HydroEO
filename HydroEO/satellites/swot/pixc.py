"""SWOT Pixel Cloud (PIXC) download, preprocessing, and rasterization.

Entry point: download_pixc(config, project_dir, credentials)
"""

from __future__ import annotations

import datetime
import logging
import re
from pathlib import Path
from typing import Any

import earthaccess
import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
import xarray as xr
from pyproj import CRS
from rasterio.transform import from_origin
from scipy.stats import binned_statistic_2d
from shapely.geometry import box
from tqdm import tqdm

logger = logging.getLogger(__name__)

# Classification enum for SWOT PIXC data
CLASS_MAP = {
    "land": 1,
    "land_near_water": 2,
    "water_near_land": 3,
    "open_water": 4,
    "dark_water": 5,
    "low_coh_water_near_land": 6,
    "open_low_coh_water": 7,
}

DEFAULT_FIELDS = ["heightEGM"]
DEFAULT_STAT_METHOD = "median"
DEFAULT_GRID_RESOLUTION = 100  # meters


def download_pixc(
    config: dict[str, Any],
    project_dir: str,
    credentials: tuple[str | None, str | None],
) -> None:
    """Download, preprocess, and rasterize SWOT Pixel Cloud data for an AOI.

    Parameters
    ----------
    config:
        The ``swot_pixc`` section of the project YAML config.
    project_dir:
        Root project directory; outputs are written to
        ``<project_dir>/swot_pixc/<aoi_name>/``.
    credentials:
        ``(earthdata_username, earthdata_password)`` tuple.
    """
    logger.debug(
        "SWOT Pixel Cloud Download - AOI: %s, Product: %s, Temporal range: %s to %s",
        config["aoi"]["name"],
        config.get("product", "SWOT_L2_HR_PIXC_2.0"),
        config["startdate"],
        config["enddate"],
    )

    aoi_name = config["aoi"]["name"]
    product = config.get("product", "SWOT_L2_HR_PIXC_2.0")
    base_dir = Path(project_dir) / "swot_pixc" / aoi_name
    raw_dir = base_dir / "raw" / product
    trimmed_dir = base_dir / "trimmed"
    raster_dir = base_dir / "raster"
    raw_dir.mkdir(parents=True, exist_ok=True)
    trimmed_dir.mkdir(parents=True, exist_ok=True)
    raster_dir.mkdir(parents=True, exist_ok=True)

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
        _preprocess_granules(config, raw_dir, trimmed_dir, log_path)
    else:
        logger.info("No new granules to download and no unprocessed netCDF files found")

    trimmed_files = list(trimmed_dir.glob("*_trimmed.geojson"))
    if trimmed_files:
        logger.debug(
            "Found %d trimmed GeoJSON files, proceeding with rasterization",
            len(trimmed_files),
        )
        _rasterize_granules(config, trimmed_dir, raster_dir)
    else:
        logger.info("No trimmed GeoJSON files found, skipping rasterization phase")

    logger.debug("SWOT Pixel Cloud download and processing complete")


def _download_granules(
    config: dict,
    raw_dir: Path,
    processed_granules: set[str],
    credentials: tuple[str | None, str | None],
) -> list[str]:
    """Download SWOT PIXC granules matching query, skipping already processed ones."""
    logger.debug("=== SWOT Pixel Cloud Download Phase ===")

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
        logger.debug("AOI Type: bbox, Bounds: %s", bounds)
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

    product = config.get("product", "SWOT_L2_HR_PIXC_2.0")
    logger.info(
        "Searching for SWOT PIXC product '%s' in AOI %s", product, aoi_config["name"]
    )
    try:
        swot_results = earthaccess.search_data(
            short_name=product,
            bounding_box=bounds,
            temporal=(time_start, time_end),
        )
        logger.debug("Found %d granules matching query", len(swot_results))
    except Exception as e:
        logger.error("Error searching for SWOT data: %s", e)
        return []

    if not swot_results:
        logger.warning("No SWOT PIXC granules found matching the query")
        return []

    granule_ids = [result["umm"]["GranuleUR"] for result in swot_results]
    new_granules = [gid for gid in granule_ids if gid not in processed_granules]
    logger.debug(
        "%d new granules to download (already processed: %d)",
        len(new_granules),
        len(processed_granules),
    )

    if not new_granules:
        logger.info("All granules already processed, skipping download")
        return []

    new_results = [r for r in swot_results if r["umm"]["GranuleUR"] in new_granules]
    logger.debug("Downloading %d granules to %s", len(new_results), raw_dir)
    try:
        downloaded_files = earthaccess.download(
            new_results, str(raw_dir), show_progress=True
        )
        logger.debug("Successfully downloaded %d files", len(downloaded_files))
        return downloaded_files or []
    except Exception as e:
        logger.error("Error downloading SWOT data: %s", e)
        return []


def _preprocess_granules(
    config: dict,
    raw_dir: Path,
    trimmed_dir: Path,
    log_path: Path,
) -> None:
    """Extract pixel cloud data, filter by class, clip to AOI, and save as GeoJSON.

    Skips granules that have already produced trimmed output.
    """
    logger.debug("=== SWOT Pixel Cloud Preprocessing Phase ===")

    aoi_config = config["aoi"]
    classes = config.get("classes", ["open_water", "water_near_land"])

    # Convert class names to integer flags
    class_flags = [CLASS_MAP.get(c) for c in classes if c in CLASS_MAP]
    if not class_flags:
        logger.warning("No valid classes found in config, defaulting to open_water")
        class_flags = [CLASS_MAP["open_water"]]
    logger.debug("Filtering by classes: %s (flags: %s)", classes, class_flags)

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

    # Filter: only process NC files that haven't produced trimmed output yet
    nc_files = []
    already_processed = []
    for nc_path in all_nc_files:
        # Infer trimmed filename
        trimmed_path = trimmed_dir / (nc_path.stem + "_trimmed.geojson")
        if not trimmed_path.exists():
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
            # Open pixel_cloud group from PIXC netCDF
            try:
                ds = xr.open_dataset(
                    str(nc_path), group="pixel_cloud", engine="netcdf4"
                )  # TODO: previously: open_mfdataset. Install dask to use it.
            except (ValueError, RuntimeError) as e:
                # Fall back to h5netcdf if netcdf4 fails
                if "chunk manager" in str(e) or "netcdf4" in str(e):
                    logger.debug("Trying h5netcdf engine for %s", nc_path.name)
                    ds = xr.open_dataset(
                        str(nc_path),
                        group="pixel_cloud",
                        engine="h5netcdf",
                    )
                else:
                    raise

            # Set CRS to WGS84 for point data
            ds = ds.rio.write_crs("EPSG:4326", inplace=True)

            # Flatten and filter by classification
            classification_flat = ds.classification.values.ravel()
            class_condition = np.isin(classification_flat, class_flags)

            if not np.any(class_condition):
                logger.debug("No points matching classification in %s", nc_path.name)
                continue

            # Extract and filter all required fields
            lon_flat = ds.longitude.values.ravel()[class_condition]
            lat_flat = ds.latitude.values.ravel()[class_condition]
            height = ds.height.values.ravel()[class_condition]
            geoid = ds.geoid.values.ravel()[class_condition]
            solid_earth_tide = ds.solid_earth_tide.values.ravel()[class_condition]
            load_tide = ds.load_tide_fes.values.ravel()[class_condition]
            pole_tide = ds.pole_tide.values.ravel()[class_condition]

            # Optional fields (may not exist in all PIXC versions)
            water_frac = None
            phase_noise_std = None
            dheight_dphase = None
            sig0 = None
            classification = classification_flat[class_condition]

            try:
                water_frac = ds.water_frac.values.ravel()[class_condition]
            except (KeyError, AttributeError):
                pass
            try:
                phase_noise_std = ds.phase_noise_std.values.ravel()[class_condition]
            except (KeyError, AttributeError):
                pass
            try:
                dheight_dphase = ds.dheight_dphase.values.ravel()[class_condition]
            except (KeyError, AttributeError):
                pass
            try:
                sig0 = ds.sig0.values.ravel()[class_condition]
            except (KeyError, AttributeError):
                pass

            # Compute heightEGM: height - geoid - solid_earth_tide - load_tide - pole_tide
            heightEGM = height - geoid - solid_earth_tide - load_tide - pole_tide

            # Build GeoDataFrame
            data = {
                "height": height,
                "heightEGM": heightEGM,
                "geoid": geoid,
                "solid_earth_tide": solid_earth_tide,
                "load_tide": load_tide,
                "pole_tide": pole_tide,
                "classification": classification,
            }

            if water_frac is not None:
                data["water_frac"] = water_frac
            if phase_noise_std is not None:
                data["phase_noise_std"] = phase_noise_std
            if dheight_dphase is not None:
                data["dheight_dphase"] = dheight_dphase
            if sig0 is not None:
                data["sig0"] = sig0

            gdf = gpd.GeoDataFrame(
                pd.DataFrame(data),
                geometry=gpd.points_from_xy(lon_flat, lat_flat),
                crs="EPSG:4326",
            )

            # Clip to AOI if provided
            if aoi_gdf is not None and len(gdf) > 0:
                try:
                    gdf = gpd.clip(gdf, aoi_gdf)
                except Exception as e:
                    _defer_warning("Clip failed for %s: %s", nc_path.name, e)

            # Skip if resulting GeoDataFrame is empty or very small
            if len(gdf) == 0:
                logger.debug(
                    "No points remaining after filtering/clipping in %s",
                    nc_path.name,
                )
                continue

            # Save to GeoJSON
            trimmed_path = trimmed_dir / (nc_path.stem + "_trimmed.geojson")
            try:
                gdf.to_file(str(trimmed_path), driver="GeoJSON")
                logger.debug("Wrote trimmed GeoJSON: %s", trimmed_path.name)
            except Exception as e:
                _defer_warning("Failed to write trimmed GeoJSON: %s", e)
                continue

        except Exception as e:
            _defer_warning("Cannot process %s: %s", nc_path.name, e)
            continue

    # Log processed NC stems
    nc_names = [f.stem for f in nc_files]
    with open(log_path, "a") as log_file:
        for name in nc_names:
            log_file.write(name + "\n")

    for warning in deferred_warnings:
        logger.debug(warning)

    logger.debug("Cleaning up raw netCDF files")
    deferred_warnings = []
    for nc_path in tqdm(nc_files, desc="Cleaning up"):
        try:
            nc_path.unlink()
        except Exception as e:
            _defer_warning("Failed to delete %s: %s", nc_path.name, e)

    for warning in deferred_warnings:
        logger.debug(warning)


def _rasterize_granules(
    config: dict,
    trimmed_dir: Path,
    raster_dir: Path,
) -> None:
    """Grid trimmed PIXC GeoJSON point data to raster TIFFs.

    Groups points by date, rasterizes to a regular grid, and outputs GeoTIFFs.
    """
    logger.debug("=== SWOT Pixel Cloud Rasterization Phase ===")

    fields = config.get("fields", DEFAULT_FIELDS)
    stat_method = config.get("stat_method", DEFAULT_STAT_METHOD)
    grid_resolution = config.get("grid_resolution", DEFAULT_GRID_RESOLUTION)

    # Get target CRS, auto-detecting if not specified
    target_crs = config.get("target_crs")
    if not target_crs:
        aoi_bbox = config["aoi"].get("bbox")
        if aoi_bbox:
            # Estimate UTM from bbox center
            center_lon = (aoi_bbox[0] + aoi_bbox[2]) / 2
            center_lat = (aoi_bbox[1] + aoi_bbox[3]) / 2
            center_gdf = gpd.GeoDataFrame(
                geometry=[box(center_lon, center_lat, center_lon, center_lat)],
                crs="EPSG:4326",
            )
            utm_crs = center_gdf.estimate_utm_crs()
            target_crs = utm_crs.to_string() if utm_crs else "EPSG:32633"
        else:
            target_crs = "EPSG:32633"  # Default to UTM 33N
        logger.debug("Auto-detected target CRS: %s", target_crs)

    try:
        target_crs_obj = CRS.from_string(target_crs)
    except Exception as e:
        logger.error("Invalid target CRS %s: %s", target_crs, e)
        return

    trimmed_files = sorted(trimmed_dir.glob("*_trimmed.geojson"))
    if not trimmed_files:
        logger.warning("No trimmed GeoJSON files found in %s", trimmed_dir)
        return

    logger.debug("Loading %d trimmed GeoJSON files", len(trimmed_files))

    # Load all GeoJSONs and group by date
    date_map = {}  # date_str -> list of trimmed_file paths

    # Regex to find datetime in format YYYYMMDDTHHMMSS (year must be 19xx or 20xx)
    datetime_pattern = re.compile(r"((?:19|20)\d{6}T\d{6})")

    for trimmed_file in tqdm(trimmed_files, desc="Loading GeoJSON files"):
        try:
            # Parse date from filename using regex to find YYYYMMDDTHHMMSS pattern
            match = datetime_pattern.search(trimmed_file.stem)
            if match:
                date_str = match.group(1)  # YYYYMMDDTHHMMSS
                try:
                    date_obj = datetime.datetime.strptime(date_str, "%Y%m%dT%H%M%S")
                    date_key = date_obj.strftime("%Y%m%d")  # YYYYMMDD for grouping
                except ValueError:
                    logger.debug("Could not parse date from %s", trimmed_file.name)
                    continue
            else:
                logger.debug("Could not extract date from %s", trimmed_file.name)
                continue

            if date_key not in date_map:
                date_map[date_key] = []
            date_map[date_key].append(trimmed_file)

        except Exception as e:
            logger.warning("Error processing %s: %s", trimmed_file.name, e)
            continue

    if not date_map:
        logger.warning("No valid dates found in GeoJSON files")
        return

    logger.debug("Found %d unique dates", len(date_map))

    # Process each date
    for date_key in tqdm(sorted(date_map.keys()), desc="Rasterizing by date"):
        try:
            # Load all GeoJSON files for this date
            gdf_list = []
            for gj_file in date_map[date_key]:
                try:
                    gdf = gpd.read_file(gj_file)
                    gdf_list.append(gdf)
                except Exception as e:
                    logger.warning("Could not read %s: %s", gj_file.name, e)
                    continue

            if not gdf_list:
                logger.debug("No valid GeoJSON data for date %s", date_key)
                continue

            # Concatenate all points for this date
            gdf_date = pd.concat(gdf_list, ignore_index=True)
            gdf_date = gpd.GeoDataFrame(gdf_date, geometry="geometry", crs="EPSG:4326")

            if len(gdf_date) == 0:
                logger.debug("No points after concatenation for date %s", date_key)
                continue

            # Reproject to target CRS
            gdf_date = gdf_date.to_crs(target_crs_obj)

            # Get bounds in projected coordinates
            bounds = gdf_date.total_bounds  # (minx, miny, maxx, maxy)
            xmin, ymin, xmax, ymax = bounds

            # Create grid cells
            nx = int(np.ceil((xmax - xmin) / grid_resolution))
            ny = int(np.ceil((ymax - ymin) / grid_resolution))

            x_edges = xmin + grid_resolution * np.arange(nx + 1)
            y_edges = ymin + grid_resolution * np.arange(ny + 1)

            logger.debug(
                "Rasterizing %d points for %s: grid %dx%d at %d m resolution",
                len(gdf_date),
                date_key,
                nx,
                ny,
                grid_resolution,
            )

            # Extract x, y coords in projected space
            x_coords = np.array([geom.x for geom in gdf_date.geometry])
            y_coords = np.array([geom.y for geom in gdf_date.geometry])

            # Rasterize each field
            for field in fields:
                if field not in gdf_date.columns:
                    logger.warning(
                        "Field '%s' not found in GeoJSON for date %s", field, date_key
                    )
                    continue

                try:
                    z_values = gdf_date[field].values

                    # Use binned_statistic_2d to grid the data
                    stat_raster, _, _, _ = binned_statistic_2d(
                        y_coords,
                        x_coords,
                        z_values,
                        statistic=stat_method,
                        bins=[y_edges, x_edges],
                    )

                    # Flip vertically to match rasterio convention
                    stat_raster = np.flipud(stat_raster)

                    # Create transform
                    transform = from_origin(
                        xmin, ymax, grid_resolution, grid_resolution
                    )

                    # Write GeoTIFF
                    output_file = raster_dir / f"{date_key}_{field}.tif"
                    with rasterio.open(
                        output_file,
                        "w",
                        driver="GTiff",
                        height=stat_raster.shape[0],
                        width=stat_raster.shape[1],
                        count=1,
                        dtype=stat_raster.dtype,
                        crs=target_crs_obj,
                        transform=transform,
                        nodata=np.nan,
                    ) as dst:
                        dst.write(stat_raster, 1)

                    logger.debug("Wrote GeoTIFF: %s", output_file.name)

                except Exception as e:
                    logger.warning(
                        "Failed to rasterize field '%s' for date %s: %s",
                        field,
                        date_key,
                        e,
                    )
                    continue

        except Exception as e:
            logger.warning("Error processing date %s: %s", date_key, e)
            continue

    logger.debug("Rasterization complete")
