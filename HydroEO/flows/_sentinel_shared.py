"""
Sentinel-3/6 download logic shared by both reservoirs and rivers.
"""

import logging

from HydroEO.satellites import sentinel
from HydroEO.utils import general

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from HydroEO.project import Project

logger = logging.getLogger(__name__)


def _sentinel6_use_earthdata(prj: "Project") -> bool:
    """
    Whether Sentinel-6 should be downloaded from PO.DAAC/EarthData (HR
    product, 20Hz Ku-band) rather than CREODIAS (LR product only, see
    query()'s productType="P4_2__LR_____"). Set via
    mission_options["sentinel6"]["source"] = "earthdata" in config.
    """
    return (
        prj.mission_options.get("sentinel6", {})
        .get("source", "creodias")
        .lower()
        == "earthdata"
    )


def _download_sentinel_for_target(
    prj: "Project", mission: str, product: str, coords, download_dir,
    startdate, enddate, sentinel_creds, session_token, session_start_time,
) -> tuple:
    """
    Download + subset Sentinel-3/6 data for one target's AOI (a
    reservoir polygon, or a river waterbody corridor's envelope) --

    For Sentinel-6, if _sentinel6_use_earthdata(prj) is True, uses
    PO.DAAC/EarthData (see sentinel.query_earthdata/download_earthdata)
    to get the HR product in .nc format instead of CREODIAS's LR-only product
    in SAFE-zip format.

    Returns (session_token, session_start_time).
    """
    dir_key = mission
    use_earthdata = mission == "sentinel6" and _sentinel6_use_earthdata(prj)

    logger.info(
        "Searching for Sentinel-%s (%s) from %s to %s",
        product,
        "PO.DAAC/EarthData HR" if use_earthdata else "CREODIAS",
        startdate,
        enddate,
    )

    if use_earthdata:
        s6_opts = prj.mission_options.get("sentinel6", {})
        results = sentinel.query_earthdata(
            aoi=coords,
            startdate=startdate,
            enddate=enddate,
            latency=s6_opts.get("latency", "NTC"),
            short_name=s6_opts.get("short_name"),
        )
        sentinel.download_earthdata(results, download_directory=download_dir)
        # EarthData granules arrive flat already -- no SAFE-zip to unzip.
    else:
        ids = sentinel.query(
            aoi=coords,
            startdate=startdate,
            enddate=enddate,
            product=product,
            creodias_credentials=sentinel_creds,
        )

        session_token, session_start_time = sentinel.download(
            ids,
            download_directory=download_dir,
            creodias_credentials=sentinel_creds,
            token=session_token,
            session_start_time=session_start_time,
            threads=prj.mission_options.get(dir_key, {}).get("download_threads", 1),
        )

        general.unzip_dir_files_with_ext(
            download_dir, download_dir, ".nc", show_progress=True
        )

    sentinel.subset(
        aoi=coords,
        download_dir=download_dir,
        dest_dir=download_dir,
        file_id=prj.mission_options.get(dir_key, {}).get(
            "subset_file_id", "enhanced_measurement.nc"
        ),
        product=product,
        show_progress=True,
    )

    return session_token, session_start_time


