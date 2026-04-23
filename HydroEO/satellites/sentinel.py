import logging
import datetime

from HydroEO.downloaders import creodias
from HydroEO.satellites import sentinel_download, sentinel_preprocess
from HydroEO.utils import geometry


logger = logging.getLogger(__name__)


def query(
    aoi: list,
    startdate: datetime.date,
    enddate: datetime.date,
    product: str = "S3",
    creodias_credentials: tuple = None,
) -> object:
    return sentinel_download.query(
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
    return sentinel_download.download(
        ids=ids,
        download_directory=download_directory,
        creodias_credentials=creodias_credentials,
        creodias_client=creodias,
        token=token,
        session_start_time=session_start_time,
        threads=threads,
    )


def subset(
    aoi: list,
    download_dir: str,
    dest_dir: str,
    file_id: str = "enhanced_measurement.nc",
    product: str = "S3",
    show_progress=False,
):
    """Subset sentinel-3/6 netCDF files to area of interest."""
    return sentinel_preprocess.subset(
        aoi=aoi,
        download_dir=download_dir,
        dest_dir=dest_dir,
        file_id=file_id,
        product=product,
        show_progress=show_progress,
    )


def extract_observations(src_dir, dst_path, features, sigma0_max=1e5):
    """Extract Sentinel observations to shapefile."""
    return sentinel_preprocess.extract_observations(
        src_dir=src_dir,
        dst_path=dst_path,
        features=features,
        sigma0_max=sigma0_max,
    )


def get_latest_obs_date(data_dir, product):
    """Get latest observation date from Sentinel product shapefile."""
    return sentinel_preprocess.get_latest_obs_date(data_dir=data_dir, product=product)
