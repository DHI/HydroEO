import os
import shutil
import zipfile
import logging

import py_hydroweb

import pandas as pd
import geopandas as gpd
import shapely

from HydroEO.utils import general

logger = logging.getLogger(__name__)

help_message = """
Download products from your hydroweb.next projects (https://hydroweb.next.theia-land.fr) using the py-hydroweb lib (https://pypi.org/project/py-hydroweb/)
This script is an example tuned for your last hydroweb.next project but feel free to adapt it for future requests.
Follow these steps:
1. If not already done, install py-hydroweb latest version using `pip install -U py-hydroweb` or `conda install py-hydroweb` (WARNING: python >= 3.8 is required)
2a. Generate an API-Key from hydroweb.next portal in your user settings
2b. Carefully store your API-Key (2 options):
- either in an environment variable `export HYDROWEB_API_KEY="<your_key_here>"`
- or in below script by replacing <your_key_here>
3. You can change download directory by adding an `output_folder` parameter when calling `submit_and_download_zip` (see below). By default, current path is used.
4. You are all set, run this script `python download_script.py`

For more documentation about how to use the py-hydroweb lib, please refer to https://pypi.org/project/py-hydroweb/.
"""


def download_PLD(
    download_dir: str,
    bounds: list,
    raw_pld_path: str = None,
    keep_raw: bool = True,
    continent_codes: list = None,
    ):
    """Download and subset SWOT Prior Lake Database.

    Parameters
    ----------
    download_dir : str
        Directory where PLD files will be stored (typically {main_dir}/aux/PLD/)
    bounds : list
        Geographic bounds [lon_min, lat_min, lon_max, lat_max] for subsetting
    raw_pld_path : str, optional
        Path to existing PLD zip file or extracted folder. If provided, skips download.
    keep_raw : bool, default True
        If False, delete raw zip and temp extraction folder after subset is created.
    continent_codes : list[str], optional
        Continent-tile suffixes (e.g. ["AS", "AU"]) identifying which
        full-schema PLD files (containing res_id, not just lake_id +
        geometry) to use for backfilling res_id only. 
        Inspect your downloaded PLD folder to see which ones are actually present. 
        If omitted, all full-schema tiles found are used.
     """
    # create download directory if needed
    general.ifnotmakedirs(download_dir)

    # Define paths early to enable skip logic
    zipped_name = "PLD_temp.zip"
    downloaded_zip_path = os.path.join(download_dir, zipped_name)
    unzipped_dir = os.path.join(download_dir, "PLD_temp")
    extracted_dir = os.path.join(
        unzipped_dir, "SWOT_PRIOR_LAKE_DATABASE", "SWOT_PRIOR_LAKE_DATABASE"
    )

    # Track whether the zip was downloaded by HydroEO (vs user-provided)
    hydroweb_managed_zip = True

    def _dir_has_pld_files(d):
        return any(f.endswith((".sqlite", ".gpkg")) for f in os.listdir(d))

     # Check if extracted files already exist
    extracted_files_exist = os.path.isdir(extracted_dir) and _dir_has_pld_files(
        extracted_dir)

    # Handle user-provided raw_pld_path
    if raw_pld_path is not None and os.path.exists(raw_pld_path):
        if os.path.isfile(raw_pld_path) and raw_pld_path.endswith(".zip"):
            # User provided a zip file
            logger.info("Using provided PLD zip file: %s", raw_pld_path)
            downloaded_zip_path = raw_pld_path
            hydroweb_managed_zip = False
        elif os.path.isdir(raw_pld_path):
            # User provided a directory
            logger.info("Using provided PLD directory: %s", raw_pld_path)
            # Check if it contains PLD files (.sqlite or .gpkg) directly
            if _dir_has_pld_files(raw_pld_path):
                extracted_dir = raw_pld_path
                extracted_files_exist = True
            else:
                # Check for nested SWOT_PRIOR_LAKE_DATABASE structure
                nested_path = os.path.join(
                    raw_pld_path, "SWOT_PRIOR_LAKE_DATABASE", "SWOT_PRIOR_LAKE_DATABASE"
                )
                if os.path.isdir(nested_path) and _dir_has_pld_files(nested_path):
                    extracted_dir = nested_path
                    extracted_files_exist = True
                    unzipped_dir = raw_pld_path  # track parent for deletion logic
                else:
                    logger.warning(
                        "Provided directory does not contain .sqlite or "
                        ".gpkg PLD files: %s",
                        raw_pld_path,
                    )
                    return
            hydroweb_managed_zip = False
    # Download only if zip file doesn't exist and extracted files don't exist
    elif os.path.isfile(downloaded_zip_path):
        logger.info(
            "Zip file already exists at %s, skipping download", downloaded_zip_path
        )
    elif extracted_files_exist:
        logger.info(
            "Extracted files already exist in %s, skipping download", extracted_dir
        )
    else:
        logger.info("Downloading SWOT Prior Lake Database")
        # Create a client
        #  - either using the API-Key environment variable
        client: py_hydroweb.Client = py_hydroweb.Client(
            "https://hydroweb.next.theia-land.fr/api"
        )

        # Initiate a new download basket (input the name you want here)
        basket: py_hydroweb.DownloadBasket = py_hydroweb.DownloadBasket("pld_download")

        # Add collections in our basket
        basket.add_collection("SWOT_PRIOR_LAKE_DATABASE", bbox=bounds)

        # Do download (input the archive name you want here, and optionally an output folder)
        downloaded_zip_path = client.submit_and_download_zip(
            basket, zip_filename=zipped_name, output_folder=download_dir
        )

    # Extract only if extracted files don't already exist
    if not extracted_files_exist:
        if os.path.isfile(downloaded_zip_path):
            logger.info("Extracting zip file")
            with zipfile.ZipFile(downloaded_zip_path, "r") as zip_ref:
                zip_ref.extractall(unzipped_dir)
        else:
            logger.warning(
                "No zip file found at %s and no extracted files at %s",
                downloaded_zip_path,
                extracted_dir,
            )
            return
    else:
        logger.info("Extracted files already exist, skipping extraction")

    # Ensure extracted_dir is set (in case raw_pld_path was a directory)
    if not os.path.isdir(extracted_dir):
        extracted_dir = os.path.join(
            unzipped_dir, "SWOT_PRIOR_LAKE_DATABASE", "SWOT_PRIOR_LAKE_DATABASE"
        )

    # Get list of downloaded files
    downloaded_files = os.listdir(extracted_dir)

    def _read_pld_file(filepath, is_sqlite):
        if is_sqlite:
            # SWOT PLD .sqlite tiles store lake_id as the feature ID
            g = gpd.read_file(filepath, layer="lake", fid_as_index=True)
            g.columns = [c.lower() for c in g.columns]
            g.index.name = "lake_id"
            g = g.reset_index()
        else:
            # .gpkg PLD files already carry lake_id as a normal column.
            g = gpd.read_file(filepath)
            g.columns = [c.lower() for c in g.columns]
        return g

    logger.info("Merging products and removing temporary files")
    light_files = [
        f
        for f in downloaded_files
        if "_light" in f.lower() and f.lower().endswith((".gpkg", ".sqlite"))
    ]
    all_tile_files = [
        f
        for f in downloaded_files
        if f not in light_files and (f.endswith(".sqlite") or f.endswith(".gpkg"))
    ]
    bbox_poly = shapely.Polygon.from_bounds(*bounds)

    if light_files:
        if len(light_files) > 1:
            logger.warning(
                "Multiple '_light' PLD files found (%s) -- using all of "
                "them, but normally there should be exactly one.",
                ", ".join(light_files),
            )
        coverage_list = [
            _read_pld_file(os.path.join(extracted_dir, f), f.lower().endswith(".sqlite"))
            for f in light_files
        ]
        coverage_gdf = pd.concat(coverage_list, ignore_index=True)
        if "res_id" not in coverage_gdf.columns:
            coverage_gdf["res_id"] = None
        coverage_gdf = coverage_gdf[["lake_id", "res_id", "geometry"]]

        # Find lake using intersection with PLD
        coverage_gdf = coverage_gdf.loc[
            coverage_gdf.intersects(bbox_poly)
        ].reset_index(drop=True)

        backfill_tile_files = all_tile_files
        if continent_codes:
            backfill_tile_files = [
                f
                for f in all_tile_files
                if any(code.lower() in f.lower() for code in continent_codes)
            ]
            unmatched_codes = [
                code
                for code in continent_codes
                if not any(code.lower() in f.lower() for f in all_tile_files)
            ]
            if unmatched_codes:
                logger.warning(
                    "hydroweb.continent_codes %s requested for PLD res_id "
                    "backfill, but no matching tile file was found among "
                    "the downloaded PLD data (%s). Re-download the PLD "
                    "if you expect a tile for these regions.",
                    unmatched_codes,
                    ", ".join(downloaded_files),
                )

        if backfill_tile_files and len(coverage_gdf) > 0:
            res_id_map = {}
            for file in backfill_tile_files:
                filepath = os.path.join(extracted_dir, file)
                logger.info("found %s (for res_id backfill)", filepath)
                tile_gdf = _read_pld_file(filepath, file.endswith(".sqlite"))
                if "res_id" in tile_gdf.columns:
                    res_id_map.update(dict(zip(tile_gdf["lake_id"], tile_gdf["res_id"])))
            if res_id_map:
                coverage_gdf["res_id"] = coverage_gdf["lake_id"].map(
                    res_id_map
                ).combine_first(coverage_gdf["res_id"])

        still_missing = int(coverage_gdf["res_id"].isna().sum())
        if still_missing > 0:
            logger.warning(
                "%d of %d PLD lake(s) within this project's bounds could "
                "not have res_id backfilled from the available "
                "tile file(s) (%s). Set hydroweb.continent_codes in "
                "your config to the correct region code(s) "
                "and/or re-download the PLD to include the missing tile for backfilling res_id.",
                still_missing,
                len(coverage_gdf),
                ", ".join(all_tile_files) if all_tile_files else "none present",
            )
    else:
        logger.warning(
            "No global '_light' PLD file found among the downloaded PLD "
            "data (%s) -- falling back to per-continent/region "
            "tile files. If this PLD download "
            "Check missing_res_id and redownload PLD if data is unexpectedly missing. " \
            "If you expect a global '_light' file, re-download the PLD to include it.",
            ", ".join(downloaded_files),
        )
        coverage_list = [
            _read_pld_file(os.path.join(extracted_dir, f), f.endswith(".sqlite"))
            for f in all_tile_files
        ]
        if coverage_list:
            coverage_gdf = pd.concat(coverage_list, ignore_index=True)
        else:
            coverage_gdf = gpd.GeoDataFrame(
                {"lake_id": [], "res_id": [], "geometry": []}, crs="EPSG:4326"
            )
        if "res_id" not in coverage_gdf.columns:
            coverage_gdf["res_id"] = None
        coverage_gdf = coverage_gdf[["lake_id", "res_id", "geometry"]]
        coverage_gdf = coverage_gdf.loc[
            coverage_gdf.intersects(bbox_poly)
        ].reset_index(drop=True)

    gdf = coverage_gdf

    gdf["res_id"] = gdf["res_id"].astype("float64")


    # save the concatenated dataframe as GPKG
    export_path = os.path.join(download_dir, "PLD_subset.gpkg")
    gdf.to_file(export_path, driver="GPKG")
    logger.info("Merged data saved to: %s", export_path)

    # Cleanup raw files if requested
    if not keep_raw:
        # Only delete the zip if HydroEO downloaded it (not user-provided)
        if hydroweb_managed_zip and os.path.isfile(downloaded_zip_path):
            logger.info("Deleting raw zip file: %s", downloaded_zip_path)
            os.remove(downloaded_zip_path)
        # Always delete the temp extraction folder (HydroEO always manages it)
        if os.path.isdir(unzipped_dir):
            logger.info("Deleting temp extraction directory: %s", unzipped_dir)
            shutil.rmtree(unzipped_dir)

    # we also remove the .downloads cache that is left from hydroweb
    cache_path = os.path.join(download_dir, ".downloaded")
    if os.path.exists(cache_path):
        shutil.rmtree(cache_path)
