import os
import shutil
import zipfile

import py_hydroweb
import sqlite3

import pandas as pd
import geopandas as gpd
import shapely

from HydroEO.utils import general

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
    zipped_name = "PLD_temp.zip"
    downloaded_zip_path = client.submit_and_download_zip(
        basket, zip_filename=zipped_name, output_folder=download_dir
    )

    # unzip files
    with zipfile.ZipFile(downloaded_zip_path, "r") as zip_ref:
        unzipped_dir = os.path.join(download_dir, "PLD_temp")
        zip_ref.extractall(unzipped_dir)

    extracted_dir = os.path.join(
        unzipped_dir, "SWOT_PRIOR_LAKE_DATABASE", "SWOT_PRIOR_LAKE_DATABASE"
    )
    downloaded_files = os.listdir(extracted_dir)

    # now clean up and merge lake datafiles
    print("Merging products and removing temporary files")
    gdf_list = list()
    for file in downloaded_files:
        # read any sqlite files within the downloads folders and extract data
        if file.endswith(".sqlite"):
            print("found %s" % os.path.join(extracted_dir, file))

            # Establish sql connection
            con = sqlite3.connect(os.path.join(extracted_dir, file))
            con.enable_load_extension(True)
            con.execute('SELECT load_extension("mod_spatialite")')

            # query the needed columns, make into df and close connection
            query = "SELECT lake_id, res_id, AsBinary(GEOMETRY) as 'geometry' FROM lake"
            df = pd.read_sql_query(query, con)
            con.close()

            # transfer df to gdf
            geometry = [shapely.wkb.loads(binary_rep) for binary_rep in df.geometry]
            gdf = gpd.GeoDataFrame(
                df[["lake_id", "res_id"]], geometry=geometry, crs=4326
            )

            # filter loaded data to bounds
            gdf = gdf.loc[gdf.within(shapely.Polygon.from_bounds(*bounds))]

            # add the remaining shapes to the list for concatenation
            gdf_list.append(gdf)

    # once we have processed the files within the download directory, remove the folder
    shutil.rmtree(downloaded_zip_path)
    shutil.rmtree(unzipped_dir)

    # we also remove the .downloads cache that is left from dag
    cache_path = os.path.join(download_dir, ".downloaded")
    if os.path.exists(cache_path):
        shutil.rmtree(cache_path)

    # once we have processed all files, concatenate them
    gdf = pd.concat(gdf_list).reset_index(drop=True)

    # save the concatenated dataframe in the output folder
    export_path = os.path.join(download_dir, file_name)
    gdf.to_file(export_path)
    print(f"Merged data saved to: {export_path}")
