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


def download_PLD(download_dir: str, file_name: str, bounds: list):
    # create download directory if needed
    general.ifnotmakedirs(download_dir)

    # Define paths early to enable skip logic
    zipped_name = "PLD_temp.zip"
    downloaded_zip_path = os.path.join(download_dir, zipped_name)
    unzipped_dir = os.path.join(download_dir, "PLD_temp")
    extracted_dir = os.path.join(
        unzipped_dir, "SWOT_PRIOR_LAKE_DATABASE", "SWOT_PRIOR_LAKE_DATABASE"
    )

    # Check if extracted files already exist
    extracted_files_exist = os.path.isdir(extracted_dir) and any(
        f.endswith(".sqlite") for f in os.listdir(extracted_dir)
    )

    # Download only if zip file doesn't exist and extracted files don't exist
    if os.path.isfile(downloaded_zip_path):
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

    # Get list of downloaded files
    downloaded_files = os.listdir(extracted_dir)

    # now clean up and merge lake datafiles
    logger.info("Merging products and removing temporary files")
    gdf_list = list()
    for file in downloaded_files:
        if file.endswith(".sqlite"):
            filepath = os.path.join(extracted_dir, file)
            logger.info("found %s", filepath)

            # Read layer directly as GeoDataFrame and normalize column names
            gdf = gpd.read_file(filepath, layer="lake")
            gdf.columns = [c.lower() for c in gdf.columns]
            # Select basin_id (unique lake identifier) and rename to lake_id
            gdf = gdf[["basin_id", "res_id", "geometry"]].rename(
                columns={"basin_id": "lake_id"}
            )

            # Filter to bounds and append
            gdf = gdf.loc[gdf.within(shapely.Polygon.from_bounds(*bounds))]
            gdf_list.append(gdf)

    # once we have processed all files, concatenate them
    gdf = pd.concat(gdf_list).reset_index(drop=True)

    # save the concatenated dataframe in the output folder
    export_path = os.path.join(download_dir, file_name)
    gdf.to_file(export_path)
    logger.info("Merged data saved to: %s", export_path)

    # once we have processed the files within the download directory, remove the folder
    # TODO
    # os.remove(downloaded_zip_path)
    # shutil.rmtree(unzipped_dir)

    # we also remove the .downloads cache that is left from dag
    cache_path = os.path.join(download_dir, ".downloaded")
    if os.path.exists(cache_path):
        shutil.rmtree(cache_path)
