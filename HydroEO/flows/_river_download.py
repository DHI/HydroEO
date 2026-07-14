"""Rivers: multi-mission download orchestration.

download_rivers and _download_swot_hydrocron_timeseries are tested
together via patch.object(flows, "_name") in tests/unit/test_flows.py --
keep them in this one module. _download_rivers_icesat2/_sentinel and
_get_latest_hydrocron_obs_date are not themselves patched as siblings of
anything, so they're free to live here too (this mirrors
_reservoir_download.py's structure).
"""

import logging
import os
import datetime

import pandas as pd

from HydroEO.satellites import icesat2, sentinel
from HydroEO.utils import general
from ._river_common import _group_river_targets_by_waterbody, _iter_geometry_pieces
from ._sentinel_shared import _sentinel6_use_earthdata, _download_sentinel_for_target
from ._summaries import _river_target_corridor

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from HydroEO.project import Project

logger = logging.getLogger(__name__)


def download_rivers(prj: "Project") -> None:
    """Download altimetry data for all configured missions (rivers mode).

    SWOT uses the Hydrocron timeseries API directly (per node/reach, no
    clustering needed -- see _download_swot_hydrocron_timeseries).
    ICESat-2/Sentinel-3/6 download raw observations over a buffered
    corridor around each waterbody's SWORD targets (see
    _river_target_corridor); associating individual points with a
    specific target happens later, during extraction.

    Parameters
    ----------
    prj : Project
        Project instance with rivers configuration
    """
    if not hasattr(prj, "rivers"):
        return

    if "swot" not in prj.to_download and not any(
        m in prj.to_download for m in ("icesat2", "sentinel3", "sentinel6")
    ):
        logger.warning(
            "Rivers are configured but no mission is enabled for download. "
            "Add a top-level mission section (e.g. 'swot:', 'icesat2:', "
            "'sentinel3:'/'sentinel6:') with 'download: true' to actually "
            "download river observations. SWORD itself will still have "
            "been prepared by initialize(), which is why you may see "
            "SWORD files but no timeseries data."
        )
        return

    if "swot" in prj.to_download:
        startdate = prj.startdates.get("swot")
        enddate = prj.enddates.get("swot")

        if isinstance(startdate, list):
            startdate = datetime.date(*startdate)
        if isinstance(enddate, list):
            enddate = datetime.date(*enddate)

        _download_swot_hydrocron_timeseries(prj, startdate, enddate)

    if "icesat2" in prj.to_download:
        _download_rivers_icesat2(prj)

    if "sentinel3" in prj.to_download:
        _download_rivers_sentinel(prj, "sentinel3")

    if "sentinel6" in prj.to_download:
        _download_rivers_sentinel(prj, "sentinel6")


def _download_swot_hydrocron_timeseries(prj: "Project", startdate, enddate) -> None:
    """Download SWOT Hydrocron timeseries for river targets."""
    from urllib import parse, request as url_request
    import json

    HYDROCRON_TIMESERIES_URL = (
        "https://soto.podaac.earthdatacloud.nasa.gov/hydrocron/v1/timeseries"
    )

    # Determine feature config
    if prj.rivers.target_id_col == "node_id":
        feature = "Node"
        quality_column = "node_q"
        fields = (
            prj.mission_options.get("swot", {})
            .get("hydrocron_fields", {})
            .get("nodes", [])
        )
        max_q = (
            prj.mission_options.get("swot", {})
            .get("quality_filters", {})
            .get("nodes", {})
            .get("max_q", 2)
        )
    else:
        feature = "Reach"
        quality_column = "reach_q"
        fields = (
            prj.mission_options.get("swot", {})
            .get("hydrocron_fields", {})
            .get("reaches", [])
        )
        max_q = (
            prj.mission_options.get("swot", {})
            .get("quality_filters", {})
            .get("reaches", {})
            .get("max_q", 2)
        )

    # Group targets by waterbody
    waterbody_groups = _group_river_targets_by_waterbody(prj)
    id_label = "nodes" if prj.rivers.target_id_col == "node_id" else "reaches"

    summary = {
        "requested": 0,
        "successful": 0,
        "failed": 0,
        "empty_after_filter": 0,
    }

    for wb_id, target_ids in waterbody_groups.items():
        summary["requested"] = summary["requested"] + len(target_ids)
        output_path = os.path.join(
            prj.dirs["swot"], str(wb_id), f"{id_label}_timeseries.csv"
        )
        general.ifnotmakedirs(os.path.dirname(output_path))

        wb_startdate = startdate
        latest_obs = _get_latest_hydrocron_obs_date(output_path)
        if latest_obs is not None:
            wb_startdate = latest_obs

        deferred_warnings = []

        def _defer_warning(message, *args):
            if args:
                deferred_warnings.append(message % args)
            else:
                deferred_warnings.append(message)

        frames = []
        for target_id in tqdm(target_ids, desc="Downloading hydrocron data"):
            try:
                query_params = {
                    "feature": feature,
                    "feature_id": str(target_id),
                    "start_time": wb_startdate.strftime("%Y-%m-%dT%H:%M:%SZ"),
                    "end_time": enddate.strftime("%Y-%m-%dT%H:%M:%SZ"),
                    "output": "csv",
                    "fields": ",".join(fields),
                }
                request_url = (
                    f"{HYDROCRON_TIMESERIES_URL}?{parse.urlencode(query_params)}"
                )

                with url_request.urlopen(request_url) as response:
                    status_code = getattr(response, "status", response.getcode())
                    payload = json.loads(response.read().decode("utf-8"))
            except Exception as exc:
                _defer_warning(
                    "Hydrocron request failed for %s %s: %s",
                    prj.rivers.target_id_col,
                    target_id,
                    exc,
                )
                summary["failed"] += 1
                continue

            csv_payload = (
                payload.get("results", {}).get("csv")
                if isinstance(payload, dict)
                else None
            )
            if status_code != 200 or not csv_payload:
                _defer_warning(
                    "Hydrocron returned status %s for %s %s",
                    status_code,
                    prj.rivers.target_id_col,
                    target_id,
                )
                summary["failed"] += 1
                continue

            try:
                df = pd.read_csv(StringIO(csv_payload))
            except Exception as exc:
                _defer_warning(
                    "Failed to parse CSV for %s %s: %s",
                    prj.rivers.target_id_col,
                    target_id,
                    exc,
                )
                summary["failed"] += 1
                continue

            if df.empty or quality_column not in df.columns:
                if df.empty:
                    _defer_warning(
                        "Hydrocron returned no data for %s %s",
                        prj.rivers.target_id_col,
                        target_id,
                    )
                else:
                    _defer_warning(
                        "Quality column %s not in Hydrocron response for %s %s",
                        quality_column,
                        prj.rivers.target_id_col,
                        target_id,
                    )
                continue

            df = df[df[quality_column] <= max_q]
            if df.empty:
                _defer_warning(
                    "All Hydrocron observations filtered for %s %s (quality > %s)",
                    prj.rivers.target_id_col,
                    target_id,
                    max_q,
                )
                summary["empty_after_filter"] += 1
                continue

            frames.append(df)
            summary["successful"] += 1

        if frames:
            combined = pd.concat(frames, ignore_index=True)
            combined.to_csv(output_path, index=False)

        for warning in deferred_warnings:
            logger.debug(warning)

    logger.info(
        "Hydrocron download complete: %s requested, %s successful, %s failed, %s empty after filtering. See file logs for more info.",
        summary["requested"],
        summary["successful"],
        summary["failed"],
        summary["empty_after_filter"],
    )


def _download_rivers_icesat2(prj: "Project") -> None:
    """Download ICESat-2 ATL13 data for river waterbody groups.

    Mirrors _download_reservoirs_icesat2, but queries over a buffered
    corridor around each waterbody's SWORD targets (see
    _river_target_corridor) rather than a single reservoir polygon. If
    a waterbody's corridor comes out as disconnected pieces (a
    MultiPolygon), queries each piece separately (see
    _iter_geometry_pieces) rather than only the first -- unlike
    reservoirs, a river waterbody's targets can legitimately be
    disjoint (e.g. separate reaches far apart), so collapsing to one
    query would either miss coverage or require an artificially large
    combined shape.
    """
    waterbody_groups = _group_river_targets_by_waterbody(prj)
    explicit_buffer = getattr(prj.rivers, "extraction_buffer_meters", None)
    width_buffer_factor = getattr(prj.rivers, "width_buffer_factor", 1.05)

    startdate = prj.startdates["icesat2"]
    enddate = prj.enddates["icesat2"]
    if isinstance(startdate, list):
        startdate = datetime.date(*startdate)
    if isinstance(enddate, list):
        enddate = datetime.date(*enddate)

    for wb_id, target_ids in waterbody_groups.items():
        logger.info("Downloading ICESat-2 data for waterbody %s", wb_id)

        corridor_gdf = _river_target_corridor(
            prj, target_ids, buffer_meters=explicit_buffer,
            width_buffer_factor=width_buffer_factor,
        )
        if corridor_gdf is None:
            logger.warning(
                "No SWORD geometry found for waterbody %s; skipping "
                "ICESat-2 download.", wb_id,
            )
            continue

        pieces = _iter_geometry_pieces(corridor_gdf.geometry.iloc[0])
        parquet_dir = os.path.join(prj.dirs["icesat2_processed"], f"{wb_id}")
        general.ifnotmakedirs(parquet_dir)

        for piece_idx, geom in enumerate(pieces):
            coords = list(geom.exterior.coords)

            logger.info(
                "Searching for Icesat2 ATL13 for waterbody %s (piece %d/%d) "
                "from %s to %s", wb_id, piece_idx + 1, len(pieces), startdate, enddate,
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
                logger.warning(
                    "ICESat-2 download skipped for waterbody %s (piece %d/%d): %s",
                    wb_id, piece_idx + 1, len(pieces), exc,
                )


def _download_rivers_sentinel(prj: "Project", mission: str) -> None:
    """Download Sentinel-3 or Sentinel-6 data for river waterbody groups.

    Mirrors _download_reservoirs_sentinel (both now share
    _download_sentinel_for_target, including the CREODIAS/EarthData
    branching for Sentinel-6). NOTE: sentinel.query/query_earthdata take
    a bounding box (envelope), not the exact corridor polygon -- for a
    long or winding river corridor this can query/download a
    meaningfully larger area than the actual buffered corridor. This is
    an existing limitation inherited from the reservoir path (where it
    matters far less, since a reservoir's envelope is close to its
    actual extent), not something new introduced here -- worth
    revisiting if it turns out to matter in practice for a large or
    winding waterbody.
    """
    product = "S3" if mission == "sentinel3" else "S6"

    sentinel_creds = None
    use_earthdata_s6 = mission == "sentinel6" and _sentinel6_use_earthdata(prj)
    if use_earthdata_s6:
        prj._require_earthdata_credentials()
    else:
        sentinel_creds = prj._require_creodias_credentials()

    waterbody_groups = _group_river_targets_by_waterbody(prj)
    explicit_buffer = getattr(prj.rivers, "extraction_buffer_meters", None)
    width_buffer_factor = getattr(prj.rivers, "width_buffer_factor", 1.05)

    startdate = prj.startdates[mission]
    enddate = prj.enddates[mission]
    if isinstance(startdate, list):
        startdate = datetime.date(*startdate)
    if isinstance(enddate, list):
        enddate = datetime.date(*enddate)

    session_token = None
    session_start_time = None

    for wb_id, target_ids in waterbody_groups.items():
        logger.info("Downloading data for waterbody %s", wb_id)

        corridor_gdf = _river_target_corridor(
            prj, target_ids, buffer_meters=explicit_buffer,
            width_buffer_factor=width_buffer_factor,
        )
        if corridor_gdf is None:
            logger.warning(
                "No SWORD geometry found for waterbody %s; skipping "
                "Sentinel-%s download.", wb_id, product,
            )
            continue

        # Envelope each disconnected piece separately rather than the
        # whole (possibly MultiPolygon) corridor at once -- sentinel's
        # API only accepts a bounding box, so one envelope covering
        # widely separated pieces could be far larger than any of them
        # individually.
        pieces = _iter_geometry_pieces(corridor_gdf.geometry.iloc[0])

        download_dir = os.path.join(prj.dirs[mission], f"{wb_id}")
        general.ifnotmakedirs(download_dir)

        for piece_idx, geom in enumerate(pieces):
            coords = [(x, y) for x, y in geom.envelope.exterior.coords]

            logger.info(
                "Searching for Sentinel-%s for waterbody %s (piece %d/%d) "
                "from %s to %s", product, wb_id, piece_idx + 1, len(pieces),
                startdate, enddate,
            )
            session_token, session_start_time = _download_sentinel_for_target(
                prj, mission, product, coords, download_dir,
                startdate, enddate, sentinel_creds, session_token, session_start_time,
            )


def _get_latest_hydrocron_obs_date(output_path) -> datetime.date:
    """Get latest observation date from existing Hydrocron output."""
    if not os.path.exists(output_path):
        return None

    try:
        existing = pd.read_csv(output_path)
    except Exception as exc:
        logger.warning(
            "Unable to read existing Hydrocron output %s: %s", output_path, exc
        )
        return None

    if "time_str" not in existing.columns or existing.empty:
        return None

    timestamps = pd.to_datetime(existing["time_str"], errors="coerce", utc=True)
    timestamps = timestamps.dropna()
    if timestamps.empty:
        return None

    latest_obs = timestamps.max().to_pydatetime()
    return datetime.date(latest_obs.year, latest_obs.month, latest_obs.day)


# ============================================================================
# RIVERS: Timeseries Processing (extraction)
# ============================================================================


