import contextlib
import datetime
import logging
import os
import warnings

import earthaccess
import shapely

logger = logging.getLogger(__name__)

# Sentinel-6 MF HR product short names on PO.DAAC, by latency tier. NTC
# (non-time-critical, ~60 day latency) is the final, most complete
# reprocessed product -- the right default for a research/monitoring
# pipeline that favors completeness over speed. STC (~36hr latency)
# trades some completeness for lower latency if that's ever needed.
# G01 is the current reprocessing baseline as of this writing (Dec 2025) --
# confirmed directly against PO.DAAC's own dataset pages, not guessed.
# NRT is not included here: no NRT short name for the HR product was
# directly confirmed, so it's deliberately left out rather than guessed --
# add it once verified, if needed.
S6_HR_SHORT_NAMES = {
    "NTC": "JASON_CS_S6A_L2_ALT_HR_STD_OST_NTC_G01",
    "STC": "JASON_CS_S6A_L2_ALT_HR_STD_OST_STC_F",
}


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


def _earthdata_login():
    earthaccess.login()


def query_earthdata(
    aoi: list,
    startdate: datetime.date,
    enddate: datetime.date,
    format_coord_list,
    latency: str = "NTC",
    short_name: str = None,
) -> list:
    """
    Query PO.DAAC/EarthData for Sentinel-6 HR granules via earthaccess --
    High Rate (HR) product (20Hz Ku-band data, not just 1Hz) is only
    available via PO.DAAC/EarthData, not CDSE.

    Parameters
    ----------
    aoi : list
        Corner coordinates of the area of interest.
    startdate, enddate : datetime.date
    format_coord_list : callable
        Same coordinate-normalizing function used by query() (CREODIAS
        path) -- passed in rather than imported directly to avoid a new
        import-time dependency here.
    latency : {"NTC", "STC"}, optional
        Which latency tier's short name to use if short_name isn't given
        explicitly. Default "NTC" (final, ~60-day-latency reprocessed
        product). See S6_HR_SHORT_NAMES.
    short_name : str, optional
        Explicit PO.DAAC collection short name, overriding `latency` --
        use this to pin an exact reprocessing baseline (e.g. an older
        F08 short name) rather than relying on whatever S6_HR_SHORT_NAMES
        currently points to.

    Returns
    -------
    list
        earthaccess granule result objects
    """
    if short_name is None:
        if latency not in S6_HR_SHORT_NAMES:
            raise ValueError(
                f"Unknown latency tier {latency!r}; expected one of "
                f"{list(S6_HR_SHORT_NAMES)} or an explicit short_name."
            )
        short_name = S6_HR_SHORT_NAMES[latency]

    aoi = format_coord_list(aoi)

    _earthdata_login()

    params = {
        "short_name": short_name,
        "temporal": (startdate, enddate),
        "bounding_box": shapely.Polygon(aoi).bounds,
    }

    try:
        with _suppress_granule_size_warning():
            results = earthaccess.search_data(**params)
        return results
    except Exception as exc:
        logger.error("Error searching EarthData for Sentinel-6 HR: %s", exc)
        return []


def download_earthdata(results: list, download_directory: str) -> list:
    """
    Download Sentinel-6 HR granules found via query_earthdata. Mirrors
    SWOT's download() (see HydroEO.satellites.swot._download) exactly,
    including the same downloaded.log progress tracking to avoid
    re-downloading files across runs.
    """
    log_path = os.path.join(download_directory, "downloaded.log")
    if not os.path.exists(log_path):
        with open(log_path, "w") as log:
            pass

    with open(log_path, "r") as log:
        downloaded_ids = [line.rstrip() for line in log]

    to_download = []
    for result in results:
        file_name = result.data_links()[0].split("/")[-1].split(".")[0]
        if file_name not in downloaded_ids:
            to_download.append(result)

    logger.info(
        "%s files shown as downloaded in log", len(results) - len(to_download)
    )
    logger.info("%s files will be downloaded", len(to_download))

    if not to_download:
        return []

    try:
        with _suppress_granule_size_warning():
            files = earthaccess.download(to_download, download_directory)
    except Exception as exc:
        logger.error("Error downloading Sentinel-6 HR files: %s", exc)
        return []

    with open(log_path, "a") as log:
        for file in files or []:
            file_name = str(file).replace("\\", "/").split("/")[-1]
            log.write(file_name.split(".nc")[0] + "\n")

    return files or []


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