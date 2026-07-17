"""Sentinel-3/6 mission public API."""

import logging
import datetime

from HydroEO.downloaders import creodias
from HydroEO.satellites.sentinel import download as _download
from HydroEO.satellites.sentinel import preprocess as _preprocess
from HydroEO.utils import geometry

logger = logging.getLogger(__name__)


def query(
    aoi: list,
    startdate: datetime.date,
    enddate: datetime.date,
    product: str = "S3",
    creodias_credentials: tuple = None,
) -> object:
    return _download.query(
        aoi=aoi,
        startdate=startdate,
        enddate=enddate,
        product=product,
        creodias_client=creodias,
        format_coord_list=geometry.format_coord_list,
    )


def download(
    ids: list,
    download_directory: str,
    creodias_credentials: tuple,
    token: str = None,
    session_start_time: str = None,
    threads: int = 1,
):
    return _download.download(
        ids=ids,
        download_directory=download_directory,
        creodias_credentials=creodias_credentials,
        creodias_client=creodias,
        token=token,
        session_start_time=session_start_time,
        threads=threads,
    )


S6_HR_SHORT_NAMES = _download.S6_HR_SHORT_NAMES


def query_earthdata(
    aoi: list,
    startdate: datetime.date,
    enddate: datetime.date,
    latency: str = "NTC",
    short_name: str = None,
) -> list:
    """
    Query PO.DAAC/EarthData for the Sentinel-6 HR product (not available
    via CREODIAS). Reuses earthaccess, same as SWOT.
    """
    return _download.query_earthdata(
        aoi=aoi,
        startdate=startdate,
        enddate=enddate,
        format_coord_list=geometry.format_coord_list,
        latency=latency,
        short_name=short_name,
    )


def download_earthdata(results: list, download_directory: str):
    return _download.download_earthdata(
        results=results,
        download_directory=download_directory,
    )


def subset(
    aoi: list,
    download_dir: str,
    dest_dir: str,
    file_id: str = "enhanced_measurement.nc",
    product: str = "S3",
    show_progress=False,
):
    return _preprocess.subset(
        aoi=aoi,
        download_dir=download_dir,
        dest_dir=dest_dir,
        file_id=file_id,
        product=product,
        show_progress=show_progress,
    )


def extract_observations(
    src_dir, dst_path, features, sigma0_max=1e5, processed_log_path=None, overwrite=False
):
    return _preprocess.extract_observations(
        src_dir=src_dir,
        dst_path=dst_path,
        features=features,
        sigma0_max=sigma0_max,
        processed_log_path=processed_log_path,
        overwrite=overwrite,
    )


def get_latest_obs_date(data_dir, product):
    return _preprocess.get_latest_obs_date(data_dir=data_dir, product=product)


__all__ = [
    "query",
    "download",
    "query_earthdata",
    "download_earthdata",
    "S6_HR_SHORT_NAMES",
    "subset",
    "extract_observations",
    "get_latest_obs_date",
]