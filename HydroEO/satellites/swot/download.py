import contextlib
import datetime
import logging
import os
import warnings

import earthaccess
import shapely

from HydroEO.utils import geometry

logger = logging.getLogger(__name__)


# Baseline D short name (supersedes Version C / 2.0)
SWOT_LAKE_SHORT_NAME = "SWOT_L2_HR_LakeSP_D"


@contextlib.contextmanager
def _suppress_granule_size_warning():
    """Suppress known earthaccess DataGranule.size deprecation warning."""
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=r"As of version 1\.0, `DataGranule\.size` will be accessed as an attribute",
            category=FutureWarning,
            module=r"earthaccess\.(results|store)",
        )
        yield


def _login(credentials=None):
    """Authenticate with earthaccess, handling credentials if provided.

    Parameters
    ----------
    credentials : tuple(str, str) | None
        (username, password) tuple, or None for anonymous login.
    """
    if credentials:
        username, password = credentials
        if username and password:
            try:
                earthaccess.login(strategy="environment", persist=True)
            except Exception:
                earthaccess.login()
        else:
            logger.warning(
                "No Earthdata credentials provided, attempting anonymous login"
            )
            earthaccess.login()
    else:
        earthaccess.login()


def _search(**params):
    """Search earthaccess for granules matching query parameters.

    Parameters
    ----------
    **params
        Query parameters (short_name, temporal, bounding_box, granule_name, etc.).

    Returns
    -------
    list
        Matching granule results, or empty list on error.
    """
    try:
        with _suppress_granule_size_warning():
            results = earthaccess.search_data(**params)
        return results
    except Exception as e:
        logger.error("Error searching for data: %s", e)
        return []


def _filter_new(results, processed_granules):
    """Filter results to keep only granules not yet processed.

    Parameters
    ----------
    results : list
        Granule result objects from earthaccess.search_data().
    processed_granules : set[str]
        GranuleUR values already processed, to exclude.

    Returns
    -------
    list
        Granule result objects not in processed_granules.
    """
    granule_ids = [result["umm"]["GranuleUR"] for result in results]
    new_granule_ids = [gid for gid in granule_ids if gid not in processed_granules]
    return [r for r in results if r["umm"]["GranuleUR"] in new_granule_ids]


def _download_files(results, directory):
    """Download granule files to a local directory.

    Parameters
    ----------
    results : list
        Granule result objects to download.
    directory : str
        Local directory to write files to.

    Returns
    -------
    list
        Paths to successfully downloaded files, or empty list on error.
    """
    if not results:
        return []

    try:
        with _suppress_granule_size_warning():
            files = earthaccess.download(results, directory, show_progress=True)
        return files or []
    except Exception as e:
        logger.error("Error downloading files: %s", e)
        return []


def query(
    aoi: list,
    startdate: datetime.date,
    enddate: datetime.date,
    product: str = SWOT_LAKE_SHORT_NAME,
) -> object:
    # format coordinates and extract bounds
    aoi = geometry.format_coord_list(aoi)

    # login and authenticate earthaccess
    _login()

    # define query parameters
    params = {
        "short_name": product,
        "temporal": (startdate, enddate),
        "bounding_box": shapely.Polygon(aoi).bounds,
    }

    results = _search(**params)
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

    logger.info("%s files shown as downloaded in log", len(results) - len(to_download))
    logger.info("%s files will be downloaded", len(to_download))
    if to_download:
        files = _download_files(to_download, download_directory)

        with open(log_path, "a") as log:
            for file in files:
                file_name = str(file).replace("\\", "/").split("/")[-1]
                log.write(file_name.split(".zip")[0] + "\n")

        return files

    return []
