import datetime
import logging
import os

import shapely

logger = logging.getLogger(__name__)


def query(
    aoi: list,
    startdate: datetime.date,
    enddate: datetime.date,
    product: str,
    creodias_client,
    format_coord_list,
) -> list[str]:
    # ensure product label is in all caps to standardize
    product = product.upper()

    # normalize coordinate list
    aoi = format_coord_list(aoi)

    # build query parameters based on requested product
    if product == "S3":
        params = {
            "collection": "SENTINEL-3",
            "start_date": startdate,
            "end_date": enddate,
            "geometry": shapely.Polygon(aoi).wkt,
            "productType": "SR_2_LAN_HY",
            "instrumentShortName": "SRAL",
        }
    elif product == "S6":
        params = {
            "collection": "SENTINEL-6",
            "start_date": startdate,
            "end_date": enddate,
            "geometry": shapely.Polygon(aoi).wkt,
            "productType": "P4_2__LR_____",
            "instrumentShortName": "P4",
        }
    else:
        raise ValueError(
            f'"{product}" is unrecognized as a valid sentinel product for query'
        )

    results = creodias_client.query(**params)
    return [result["Id"] for result in results.values()]


def download(
    ids: list,
    download_directory: str,
    creodias_credentials: tuple,
    creodias_client,
    token: str = None,
    session_start_time: str = None,
    threads: int = 1,
):
    # Check if we have a progress log file in this directory, if not make it
    log_path = os.path.join(download_directory, "downloaded.log")
    if not os.path.exists(log_path):
        with open(log_path, "w") as log:
            pass

    # open the log file with reading and writing access
    with open(log_path, "r") as log:
        downloaded_ids = [line.rstrip() for line in log]
        already_downloaded = [uid for uid in ids if uid in downloaded_ids]
        ids_to_download = [uid for uid in ids if uid not in downloaded_ids]

    logger.debug("%s products returned by query.", len(ids))
    logger.debug("%s products already downloaded.", len(already_downloaded))
    logger.debug("Downloading %s files.", len(ids_to_download))

    return creodias_client.download_list(
        ids_to_download,
        username=creodias_credentials[0],
        password=creodias_credentials[1],
        token=token,
        session_start_time=session_start_time,
        outdir=download_directory,
        log_file=log_path,
        threads=threads,
    )
