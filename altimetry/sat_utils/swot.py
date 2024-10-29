import shapely
import geopandas as gpd
import glob
from pathlib import Path
import numpy as np
import pandas as pd
import os
import shutil
import zipfile
import earthaccess

import datetime

from altimetry import geometry, utils

def query(aoi: list, startdate: datetime.date, enddate: datetime.date, earthdata_credentials: tuple, product: str  ='SWOT_L2_HR_LakeSP_2.0') -> object:

    # format coordinates and extract bounds
    aoi = geometry.format_coord_list(aoi)

    # login and authenticate earthacess
    auth = earthaccess.login()

    # define query parameters
    params = {
        'short_name' : product,
        'temporal' : (startdate, enddate),
        'bounding_box' : shapely.Polygon(aoi).bounds,
    }

    # make the search
    results = earthaccess.search_data(**params)

    # make list of file names for later reference
    # result_files = list()
    # for result in results:
    #     with (earthaccess.open([result])[0] as file):
    #         result_files.append(file)

    return results


def download(result, download_directory: str):

    files = earthaccess.download(result, download_directory)

    return files


def subset_by_id(files:list, ids:list):

    for file in files:

        # extract file properties
        file_dir = os.path.dirname(file)
        file_name = os.path.basename(file)

        # make a temporary directory
        temp_dir = os.path.join(file_dir, '.temp')
        utils.ifnotmakedirs(temp_dir)

        # unzip file
        with zipfile.ZipFile(file, 'r') as zip_ref:
            zip_ref.extractall(temp_dir)

        # open shape file
        shp_file = file_name.split('.')[0]+'.shp'
        gdf = gpd.read_file(os.path.join(temp_dir, shp_file))

        # remove water bodies that have no data
        gdf = gdf.loc[gdf.obs_id != 'no_data'].reset_index(drop=True)

        # extract entries that are in id list
        gdf = gdf.loc[np.in1d(gdf.lake_id.astype(int).values, ids)]

        # save file
        export_file = 'sub_'+shp_file
        gdf.to_file(os.path.join(file_dir, export_file))

        # clean up original file and temp directory
        os.remove(file)
        shutil.rmtree(temp_dir)