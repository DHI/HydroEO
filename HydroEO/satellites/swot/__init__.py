"""SWOT mission public API."""

import datetime

import earthaccess

from HydroEO.satellites.swot import download as _download
from HydroEO.satellites.swot import preprocess as _preprocess

SWOT_LAKE_SHORT_NAME = _download.SWOT_LAKE_SHORT_NAME


def query(
    aoi: list,
    startdate: datetime.date,
    enddate: datetime.date,
    product: str = SWOT_LAKE_SHORT_NAME,
) -> object:
    return _download.query(
        aoi=aoi,
        startdate=startdate,
        enddate=enddate,
        product=product,
        earthaccess_client=earthaccess,
    )


def download(results, download_directory: str):
    return _download.download(
        results=results,
        download_directory=download_directory,
        earthaccess_client=earthaccess,
    )


def subset_by_id(files: list, ids: list):
    return _preprocess.subset_by_id(files=files, ids=ids)


def merge_shps(dir):
    return _preprocess.merge_shps(dir=dir)


def extract_observations(
    src_dir,
    dst_dir,
    dst_file_name,
    features,
    id_key,
    exclude_obs_id_values=None,
):
    return _preprocess.extract_observations(
        src_dir=src_dir,
        dst_dir=dst_dir,
        dst_file_name=dst_file_name,
        features=features,
        id_key=id_key,
        exclude_obs_id_values=exclude_obs_id_values,
        product_name=SWOT_LAKE_SHORT_NAME,
    )


def get_latest_obs_date(data_dir):
    return _preprocess.get_latest_obs_date(data_dir=data_dir)


__all__ = [
    "SWOT_LAKE_SHORT_NAME",
    "query",
    "download",
    "subset_by_id",
    "merge_shps",
    "extract_observations",
    "get_latest_obs_date",
]
