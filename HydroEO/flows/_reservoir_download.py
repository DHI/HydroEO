"""Reservoirs: multi-mission download orchestration.

download_reservoirs and its three per-mission workers
(_download_reservoirs_swot, _download_reservoirs_icesat2,
_download_reservoirs_sentinel) are tested together via
patch.object(flows, "_name") in tests/unit/test_flows.py -- keep them in
this one module.
"""

import logging
import os
import datetime

from HydroEO.satellites import swot, icesat2
from HydroEO.utils import general
from ._river_common import _simplify_to_one_polygon
from ._sentinel_shared import _sentinel6_use_earthdata, _download_sentinel_for_target

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from HydroEO.project import Project

logger = logging.getLogger(__name__)


def download_reservoirs(prj: "Project") -> None:
    """Download altimetry data for all configured missions (reservoirs mode).

    Parameters
    ----------
    prj : Project
        Project instance with download configuration
    """
    for mission in prj.to_download:
        if mission == "swot":
            _download_reservoirs_swot(prj)
        elif mission == "icesat2":
            _download_reservoirs_icesat2(prj)
        elif mission in ["sentinel3", "sentinel6"]:
            _download_reservoirs_sentinel(prj, mission)
        else:
            logger.warning("Skipping unsupported mission in download: %s", mission)


def _download_reservoirs_swot(prj: "Project") -> None:
    """Download SWOT Lake SP data for reservoirs."""
    download_dir = prj.dirs["swot"]
    general.ifnotmakedirs(download_dir)

    startdate = prj.startdates["swot"]
    enddate = prj.enddates["swot"]

    if isinstance(startdate, list):
        startdate = datetime.date(*startdate)
    if isinstance(enddate, list):
        enddate = datetime.date(*enddate)

    coords = [
        (x, y)
        for x, y in prj.reservoirs.download_gdf.unary_union.envelope.exterior.coords
    ]

    logger.info(
        "Searching for %s for aoi from %s to %s",
        swot.SWOT_LAKE_SHORT_NAME,
        startdate,
        enddate,
    )
    results = swot.query(aoi=coords, startdate=startdate, enddate=enddate)

    # Filter to prior-lake granules
    logger.info("%s products returned from query", len(results))
    to_download = []
    for result in results:
        link = result.data_links()[0]
        filename = link.split("/")[-1].lower()
        if "_prior_" in filename or "prior" in filename.split("_"):
            to_download.append(result)
    logger.info("%s prior-lake granules selected for download", len(to_download))

    _ = swot.download(to_download, download_directory=download_dir)

    files_in_dir = [
        os.path.join(download_dir, f)
        for f in os.listdir(download_dir)
        if f.endswith(".zip")
    ]
    swot.subset_by_id(
        files_in_dir, prj.reservoirs.download_gdf["prior_lake_id"].astype(int).values
    )


def _download_reservoirs_icesat2(prj: "Project") -> None:
    """Download ICESat-2 ATL13 data for reservoirs."""
    for i in prj.reservoirs.download_gdf.index:
        id = prj.reservoirs.download_gdf.loc[i, prj.reservoirs.id_key]
        logger.info("Downloading data for id %s", id)

        geom = _simplify_to_one_polygon(prj.reservoirs.download_gdf.loc[i, "geometry"])
        coords = list(geom.exterior.coords)

        parquet_dir = os.path.join(prj.dirs["icesat2_processed"], rf"{id}")
        general.ifnotmakedirs(parquet_dir)

        startdate = prj.startdates["icesat2"]
        enddate = prj.enddates["icesat2"]

        if isinstance(startdate, list):
            startdate = datetime.date(*startdate)
        if isinstance(enddate, list):
            enddate = datetime.date(*enddate)

        logger.info(
            "Searching for Icesat2 ATL13 for aoi from %s to %s", startdate, enddate
        )
        try:
            _ = icesat2.query(
                aoi=coords,
                startdate=startdate,
                enddate=enddate,
                download_directory=parquet_dir,
                atl13_options=prj.mission_options.get("icesat2", {}).get("atl13", {}),
                atl13_fields=prj.mission_options.get("icesat2", {}).get("atl13_fields")
                or None,
            )
        except Exception as exc:
            logger.warning("ICESat-2 download skipped for %s: %s", id, exc)


def _download_reservoirs_sentinel(prj: "Project", mission: str) -> None:
    """Download Sentinel-3 or Sentinel-6 data for reservoirs."""
    product = "S3" if mission == "sentinel3" else "S6"

    session_token = None
    session_start_time = None

    # EarthData (Sentinel-6 HR) needs no CREODIAS credentials at all --
    # only require them if we're actually going to use CREODIAS. But it
    # does need its OWN upfront check -- without it, earthaccess.login()
    # silently falls through to interactive prompting when nothing else
    # is configured, which hangs in a non-interactive run instead of
    # failing clearly (see Project._require_earthdata_credentials).
    sentinel_creds = None
    use_earthdata_s6 = mission == "sentinel6" and _sentinel6_use_earthdata(prj)
    if use_earthdata_s6:
        prj._require_earthdata_credentials()
    else:
        sentinel_creds = prj._require_creodias_credentials()

    for i in prj.reservoirs.download_gdf.index:
        id = prj.reservoirs.download_gdf.loc[i, prj.reservoirs.id_key]
        logger.info("Downloading data for id %s", id)

        coords = [
            (x, y)
            for x, y in prj.reservoirs.download_gdf.loc[
                i, "geometry"
            ].envelope.exterior.coords
        ]

        download_dir = os.path.join(prj.dirs[mission], rf"{id}")
        general.ifnotmakedirs(download_dir)

        startdate = prj.startdates[mission]
        enddate = prj.enddates[mission]

        if isinstance(startdate, list):
            startdate = datetime.date(*startdate)
        if isinstance(enddate, list):
            enddate = datetime.date(*enddate)

        session_token, session_start_time = _download_sentinel_for_target(
            prj, mission, product, coords, download_dir,
            startdate, enddate, sentinel_creds, session_token, session_start_time,
        )


# ============================================================================
# RIVERS: Download
# ============================================================================


