import datetime
import logging
import os
import shutil
import zipfile

import geopandas as gpd
import numpy as np
import pandas as pd
from tqdm import tqdm

from HydroEO.utils import general

logger = logging.getLogger(__name__)


def subset_by_id(files: list, ids: list):
    for file in tqdm(files, desc="Subsetting files"):
        # extract file properties
        file_dir = os.path.dirname(file)
        file_name = os.path.basename(file)

        # make a temporary directory
        temp_dir = os.path.join(file_dir, ".temp")
        general.ifnotmakedirs(temp_dir)

        # unzip file
        try:
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

        except Exception:
            logger.warning("Unable to unzip: %s", file, exc_info=True)


def merge_shps(dir):
    gdf_list = list()
    for file in os.listdir(dir):
        if file.endswith(".shp"):
            gdf_list.append(gpd.read_file(os.path.join(dir, file)))

    if not gdf_list:
        return None

    # return the combined gdf
    gdf = pd.concat(gdf_list).reset_index(drop=True)

    return gdf


def extract_observations(
    src_dir,
    dst_dir,
    dst_file_name,
    features,
    id_key,
    exclude_obs_id_values=None,
    product_name="SWOT_L2_HR_LakeSP_D",
):
    # load in combined observations from individual files in download directory
    data_gdf = merge_shps(src_dir)
    if data_gdf is None:
        return []
    excluded_obs_ids = set(exclude_obs_id_values or ["no_data"])

    empty_ids = []
    # now loop through the ids in the features gdf to extract the observations from the main one
    for _, feat in tqdm(
        features.iterrows(), total=len(features), desc="Extracting SWOT Lake SP product"
    ):
        dl_id = str(feat[id_key])
        if not np.isnan(feat["prior_lake_id"]):
            lake_id = str(int(feat["prior_lake_id"]))

            # filter observations to keep only the ones associated with this lake/reservoir
            observations = (
                data_gdf.loc[data_gdf.lake_id.astype(int).astype(str) == lake_id]
                .reset_index(drop=True)
                .sort_values(by="time")
            )

            if "obs_id" in observations.columns and excluded_obs_ids:
                observations = observations.loc[
                    ~observations["obs_id"].astype(str).isin(excluded_obs_ids)
                ].reset_index(drop=True)

            # if we have observations for this reservoir export it
            if len(observations) > 0:
                observations["platform"] = "swot"
                observations["product"] = product_name
                observations["height"] = observations.wse
                observations["date"] = pd.to_datetime(observations.time_str)
                observations["orbit"] = (
                    observations.lake_id
                )  # TODO: edit to SWOT equivalent, ask PASE

                dst_sub_dir = os.path.join(dst_dir, f"{dl_id}", "raw_observations")
                general.ifnotmakedirs(dst_sub_dir)
                dst_path = os.path.join(dst_sub_dir, dst_file_name)

                observations.to_file(dst_path)
            else:
                empty_ids.append(dl_id)

    return empty_ids


def get_latest_obs_date(data_dir):
    dates = list()

    for dir in os.listdir(data_dir):
        shp_path = os.path.join(data_dir, dir, "raw_observations", "swot.shp")

        if os.path.exists(shp_path):
            gdf = gpd.read_file(shp_path)
            dates.append(max(gdf.date.values).astype(datetime.date))

    return max(dates)
