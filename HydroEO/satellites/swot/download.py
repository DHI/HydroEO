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


def query(
    aoi: list,
    startdate: datetime.date,
    enddate: datetime.date,
    product: str = SWOT_LAKE_SHORT_NAME,
    earthaccess_client=earthaccess,
) -> object:
    # format coordinates and extract bounds
    aoi = geometry.format_coord_list(aoi)

    # login and authenticate earthacess
    earthaccess_client.login()

    # define query parameters
    params = {
        "short_name": product,
        "temporal": (startdate, enddate),
        "bounding_box": shapely.Polygon(aoi).bounds,
    }

    # Silence a known earthaccess deprecation warning until upstream migrates
    # DataGranule.size() to DataGranule.size attribute access.
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=r"As of version 1\.0, `DataGranule\.size` will be accessed as an attribute",
            category=FutureWarning,
            module=r"earthaccess\.results",
        )
        results = earthaccess_client.search_data(**params)

    return results


def download(results, download_directory: str, earthaccess_client=earthaccess):
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
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message=r"As of version 1\.0, `DataGranule\.size` will be accessed as an attribute",
                category=FutureWarning,
                module=r"earthaccess\.(results|store)",
            )
            files = earthaccess_client.download(
                to_download, download_directory, show_progress=True
            )

        with open(log_path, "a") as log:
            for file in files:
                file_name = str(file).replace("\\", "/").split("/")[-1]
                log.write(file_name.split(".zip")[0] + "\n")

        return files

    return []
