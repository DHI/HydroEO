import os
import shutil

from eodag import EODataAccessGateway

import sqlite3

import pandas as pd
import geopandas as gpd
import shapely

from HydroEO.utils import general

help_message = """
Download products from your Hydroweb.next project (https://hydroweb.next.theia-land.fr) using EODAG (https://github.com/CS-SI/eodag)
This script is an example tuned for your last Hydroweb.next project but feel free to adapt it for future requests.
Follow these steps:
1. If not already done, install EODAG and packaging latest version using `pip install -U eodag packaging` or `conda update eodag packaging`
2a. Generate an API-Key from Hydroweb.next portal in your user settings
2b. Carefully store your API-Key
- either in your eodag configuration file (usually ~/.config/eodag/eodag.yml, automatically generated the first time you use eodag) in auth/credentials/apikey="PLEASE_CHANGE_ME"
- or in an environment variable `export EODAG__HYDROWEB_NEXT__AUTH__CREDENTIALS__APIKEY="PLEASE_CHANGE_ME"`
3. You can change download directory by modifying the variable path_out. By default, current path is used.
4. You are all set, run this script `python download_SWOT_Prior_Lake_Database.py`

For more information, please refer to EODAG Documentation https://eodag.readthedocs.io
"""


def download_PLD(download_dir: str, file_name: str, bounds: list):
    # create download directory if needed
    general.ifnotmakedirs(download_dir)

    # Set timeout to 30s
    os.environ["EODAG__HYDROWEB_NEXT__SEARCH__TIMEOUT"] = (
        "30"  # TODO: maybe just set in query parameters
    )

    dag = EODataAccessGateway()

    query_args = {
        "productType": "SWOT_PRIOR_LAKE_DATABASE",
        "geom": bounds,
        "items_per_page": 2000,  # Default search criteria when iterating over collection pages
    }

    # Iterate over all pages to find all products
    print("Searching products")
    search_results = []
    for i, page_results in enumerate(dag.search_iter_page(**query_args)):
        search_results.extend(page_results)

    # This command actually downloads the matching products
    print("Downloading products")
    downloaded_dirs = dag.download_all(search_results, output_dir=download_dir)

    # now clean up and merge lake datafiles
    print("Merging products and removing temporary files")
    gdf_list = list()
    for dir in downloaded_dirs:
        for file in os.listdir(dir):
            # read any sqlite files within the downloads folders and extract data
            if file.endswith(".sqlite"):
                print("found %s" % os.path.join(dir, file))

                # Establish sql connection
                con = sqlite3.connect(os.path.join(dir, file))
                con.enable_load_extension(True)
                con.execute('SELECT load_extension("mod_spatialite")')

                # query the needed columns, make into df and close connection
                query = (
                    "SELECT lake_id, res_id, AsBinary(GEOMETRY) as 'geometry' FROM lake"
                )
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
        shutil.rmtree(dir)

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
