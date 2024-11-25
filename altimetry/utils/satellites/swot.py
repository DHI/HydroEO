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

from tqdm import tqdm

import datetime

from altimetry.utils import utils, geometry


def query(
    aoi: list,
    startdate: datetime.date,
    enddate: datetime.date,
    earthdata_credentials: tuple,
    product: str = "SWOT_L2_HR_LakeSP_2.0",
) -> object:
    # format coordinates and extract bounds
    aoi = geometry.format_coord_list(aoi)

    # login and authenticate earthacess
    auth = earthaccess.login()

    # define query parameters
    params = {
        "short_name": product,
        "temporal": (startdate, enddate),
        "bounding_box": shapely.Polygon(aoi).bounds,
    }

    # make the search
    results = earthaccess.search_data(**params)

    return results


def download(results, download_directory: str):
    # Check if we have a progress log file in this directory, if not make it
    log_path = os.path.join(download_directory, "downloaded.log")
    if not os.path.exists(log_path):
        with open(log_path, "w") as log:
            pass

    # open the log file with reading and writing access
    with open(log_path, "r") as log:
        # first read all of the downloaded ids
        downloaded_ids = [line.rstrip() for line in log]

    to_download = list()
    for result in results:
        # check file name to see if its downloaded and download if needed
        file_name = result.data_links()[0].split("/")[-1].split(".")[0]
        if file_name not in downloaded_ids:
            to_download.append(result)

    print(f"{len(results)-len(to_download)} files shown as downloaded in log")
    print(f"{len(to_download)} will be downloaded")
    if to_download:
        files = earthaccess.download(to_download, download_directory)

        with open(log_path, "a") as log:
            for file in files:
                log.write(file.split("\\")[-1].split(".zip")[0] + "\n")

        return files

    else:
        return []


def subset_by_id(files: list, ids: list):
    for file in tqdm(files, desc="Subsetting files"):
        # extract file properties
        file_dir = os.path.dirname(file)
        file_name = os.path.basename(file)

        # make a temporary directory
        temp_dir = os.path.join(file_dir, ".temp")
        utils.ifnotmakedirs(temp_dir)

        # unzip file
        with zipfile.ZipFile(file, "r") as zip_ref:
            zip_ref.extractall(temp_dir)

        # open shape file
        shp_file = file_name.split(".")[0] + ".shp"
        gdf = gpd.read_file(os.path.join(temp_dir, shp_file))

        # remove water bodies that have no data
        gdf = gdf.loc[gdf.obs_id != "no_data"].reset_index(drop=True)

        # extract entries that are in id list
        gdf = gdf.loc[np.in1d(gdf.lake_id.astype(int).values, ids)]

        if len(gdf) > 0:
            # save file
            export_file = "sub_" + shp_file
            gdf.to_file(os.path.join(file_dir, export_file))

        # clean up original file and temp directory
        os.remove(file)
        shutil.rmtree(temp_dir)


def merge_shps(dir):
    gdf_list = list()
    for file in os.listdir(dir):
        if file.endswith(".shp"):
            gdf_list.append(gpd.read_file(os.path.join(dir, file)))

    # return the combined gdf
    gdf = pd.concat(gdf_list).reset_index(drop=True)

    return gdf


def extract_observations(src_dir, dst_dir, dst_file_name, features, id_key):
    # load in combined observations from individual files in download directory
    data_gdf = merge_shps(src_dir)

    # now loop through the ids in the features gdf to extract the observations from the main one
    for i in tqdm(features.index, desc="Extracting SWOT Lake SP product"):
        dl_id = str(features.loc[i, id_key])
        lake_id = str(int(features.loc[i, "prior_lake_id"]))

        # filter observations to keep only the ones associated with this lake/reservoir
        observations = (
            data_gdf.loc[data_gdf.lake_id.astype(int).astype(str) == lake_id]
            .reset_index(drop=True)
            .sort_values(by="time")
        )

        # if we have observations for this reservoir export it
        if len(observations) > 0:
            observations["platform"] = "swot"
            observations["product"] = "SWOT_L2_HR_LakeSP_2.0"
            observations["height"] = observations.wse
            observations["date"] = pd.to_datetime(observations.time_str)

            dst_sub_dir = os.path.join(dst_dir, f"{dl_id}", "raw_observations")
            utils.ifnotmakedirs(dst_sub_dir)
            dst_path = os.path.join(dst_sub_dir, dst_file_name)

            observations.to_file(dst_path)


def get_latest_obs_date(data_dir):
    dates = list()

    for dir in os.listdir(data_dir):
        shp_path = os.path.join(data_dir, dir, "raw_observations", "swot.shp")

        if os.path.exists(shp_path):
            gdf = gpd.read_file(shp_path)
            dates.append(max(gdf.date.values).astype(datetime.date))

    return max(dates)
