"""Sentinel-3/6 download logic shared by both reservoirs and rivers.

_sentinel6_use_earthdata and _download_sentinel_for_target are used by
both _reservoir_download.py and _river_download.py (the CREODIAS/EarthData
branching only needs to exist in one place); neither is itself the target
of a sibling patch in the test suite, so this module has no test-imposed
co-location constraint of its own.
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
    shared by both _download_reservoirs_sentinel and
    _download_rivers_sentinel so the CREODIAS/EarthData branching logic
    only needs to exist in one place.

    For Sentinel-6, if _sentinel6_use_earthdata(prj) is True, uses
    PO.DAAC/EarthData (see sentinel.query_earthdata/download_earthdata)
    to get the HR product instead of CREODIAS's LR-only product.
    EarthData files arrive flat (no SAFE-zip directory), so the unzip
    step is skipped for that path -- subset() already handles both flat
    and zipped-folder inputs (see sentinel/preprocess.py's file
    discovery, extended for this).

    Returns (session_token, session_start_time) -- unchanged from what
    was passed in when using the EarthData path, since that mechanism
    (CREODIAS session reuse) doesn't apply to it.
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


