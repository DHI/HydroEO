"""Standalone flow functions for HydroEO pipelines.

These functions implement the core download, processing, and visualization logic
previously embedded in Reservoirs and Rivers classes. They operate on Project state
and external data, with no direct method dependencies.
"""

import logging
import os
import datetime
import json
import yaml
from io import StringIO
from typing import TYPE_CHECKING

import geopandas as gpd
import mikeio
import pandas as pd
from tqdm import tqdm

from HydroEO.satellites import swot, icesat2, sentinel
from HydroEO.utils import general, timeseries
from HydroEO.utils.filters import basic_filters
from HydroEO.downloaders import hydroweb
from HydroEO import plotting

if TYPE_CHECKING:
    from HydroEO.project import Project

logger = logging.getLogger(__name__)


# ============================================================================
# RESERVOIRS: Initialization
# ============================================================================


def initialize_reservoirs(prj: "Project") -> None:
    """Initialize PLD matching for reservoirs mode.

    Downloads PLD database, matches it to input reservoirs, and stores
    prior_lake_id values on prj.reservoirs.gdf.

    Parameters
    ----------
    prj : Project
        Project instance with reservoirs config/state populated
    """
    if not hasattr(prj, "reservoirs"):
        return

    if "swot" not in prj.to_download and "swot" not in prj.to_process:
        return

    # Download PLD if needed
    _download_pld(prj)

    # Match reservoirs to PLD
    _assign_pld_id(prj)

    # Export flags for missing priors
    _flag_missing_priors(prj)

    # Set download geometry (for reservoirs, same as input boundaries)
    prj.reservoirs.download_gdf = prj.reservoirs.gdf


def _download_pld(prj: "Project") -> None:
    """Download PLD database to project directory."""
    pld_path = prj.dirs["pld"]

    if os.path.exists(pld_path):
        logger.info("PLD located")
        return

    logger.info("Downloading PLD")
    download_dir = os.path.dirname(pld_path)
    bounds = list(prj.reservoirs.gdf.unary_union.bounds)
    raw_pld_path = prj.dirs.get("pld_raw")

    # Determine if raw_pld_path is inside project main_dir
    keep_raw = getattr(prj, "keep_raw_pld", False)
    effective_keep_raw = keep_raw
    if raw_pld_path is not None and os.path.exists(raw_pld_path):
        if not os.path.abspath(raw_pld_path).startswith(
            os.path.abspath(prj.dirs["main"])
        ):
            logger.warning(
                "raw_pld_path '%s' is outside project folder '%s'. "
                "Skipping deletion of raw PLD files to preserve external data.",
                raw_pld_path,
                prj.dirs["main"],
            )
            effective_keep_raw = True

    hydroweb.download_PLD(
        download_dir=download_dir,
        bounds=bounds,
        raw_pld_path=raw_pld_path,
        keep_raw=effective_keep_raw,
    )


def _assign_pld_id(prj: "Project") -> None:
    """Spatial join reservoirs with PLD to assign prior_lake_id."""
    pld = gpd.read_file(prj.dirs["pld"])

    pld = pld.rename(
        columns={"lake_id": "prior_lake_id", "res_id": "prior_res_id"}
        )
    joined_gdf = gpd.sjoin_nearest(
        prj.reservoirs.gdf.to_crs(prj.local_crs),
        pld.to_crs(prj.local_crs),
        how="left",
        max_distance=prj.mission_options.get("swot", {}).get(
            "pld_match_max_distance_m", 100
        ),
        distance_col="dist_to_pld",
    )
    joined_gdf = joined_gdf.to_crs(prj.reservoirs.gdf.crs)

    joined_gdf.loc[joined_gdf.prior_lake_id.isnull(), "prior_lake_id"] = -9999

    prj.reservoirs.gdf = joined_gdf


def _flag_missing_priors(prj: "Project") -> None:
    """Export geopackages of reservoirs present/missing in PLD to aux/PLD folder."""
    gdf = prj.reservoirs.gdf
    present = gdf.loc[gdf.prior_lake_id > 0].reset_index(drop=True)
    missing = gdf.loc[gdf.prior_lake_id < 0].reset_index(drop=True)

    # Output to aux/PLD/ folder
    pld_dir = os.path.dirname(prj.dirs["pld"])
    present_path = os.path.join(pld_dir, "present_in_pld.gpkg")
    missing_path = os.path.join(pld_dir, "missing_in_pld.gpkg")

    present.to_file(present_path, driver="GPKG")
    missing.to_file(missing_path, driver="GPKG")

    logger.info(
        "Out of the %s reservoirs, %s are present and %s are missing from the PLD.",
        len(gdf),
        len(present),
        len(missing),
    )


# ============================================================================
# RIVERS: Initialization
# ============================================================================


def initialize_rivers(prj: "Project") -> None:
    """Initialize SWORD target IDs for rivers mode.

    Parameters
    ----------
    prj : Project
        Project instance with rivers config/state populated
    """
    if not hasattr(prj, "rivers"):
        return

    if prj.rivers.input_mode == "aoi_path":
        _prepare_rivers_from_sword(prj)

    id_label = "node" if prj.rivers.target_id_col == "node_id" else "reach"
    logger.info(
        "Found river %s %s ids",
        len(prj.rivers.target_ids),
        id_label,
    )
    logger.debug(
        "Found river %s ids: %s",
        id_label,
        ", ".join(str(target_id) for target_id in prj.rivers.target_ids),
    )


def _prepare_rivers_from_sword(prj: "Project") -> None:
    """Prepare SWORD target features by spatial intersection with AOI or from saved subset.

    If SWORD_subset.gpkg exists, reads from it directly (skips download and spatial operations).
    Otherwise, ensures SWORD database, performs spatial intersection with AOI, saves subset.
    """
    subset_path = prj.dirs.get("sword_subset")

    # Gate 1: Check if subset already exists
    if subset_path and os.path.exists(subset_path):
        logger.info("SWORD subset located at %s", subset_path)
        subset = gpd.read_file(subset_path)
    else:
        # Gate 2: Ensure SWORD database and perform spatial intersection
        _ensure_sword_database(prj)

        gpkg_name = (
            f"{prj.rivers.continent_key}_sword_{prj.rivers.feature_type}_v17b.gpkg"
        )
        gpkg_path = os.path.join(prj.dirs["sword"], gpkg_name)

        if not os.path.exists(gpkg_path):
            raise FileNotFoundError(f"Expected SWORD file not found: {gpkg_path}")

        sword_gdf = gpd.read_file(gpkg_path)

        # Buffer AOI if requested
        aoi_local = prj.rivers.aoi_gdf.to_crs(prj.local_crs).copy()
        if prj.rivers.buffer_meters and prj.rivers.buffer_meters > 0:
            aoi_local["geometry"] = aoi_local.geometry.buffer(prj.rivers.buffer_meters)

        # Intersect with SWORD
        sword_local = sword_gdf.to_crs(prj.local_crs)
        subset = sword_local.loc[sword_local.intersects(aoi_local.unary_union)].copy()

        if prj.rivers.id_key not in prj.rivers.aoi_gdf.columns:
            raise KeyError(
                f"Expected AOI column '{prj.rivers.id_key}' missing from river input file"
            )

        aoi_join = aoi_local[[prj.rivers.id_key, "geometry"]].copy()
        subset = gpd.sjoin(
            subset,
            aoi_join,
            how="inner",
            predicate="intersects",
        ).drop(columns=["index_right"], errors="ignore")
        subset = subset.drop_duplicates().to_crs(prj.rivers.aoi_gdf.crs)

        # Save subset to disk
        if subset_path:
            general.ifnotmakedirs(os.path.dirname(subset_path))
            subset.to_file(subset_path, driver="GPKG")
            logger.info("SWORD subset saved to %s", subset_path)

        # Cleanup: delete gpkg folder if keep_raw_sword is False
        if not prj.keep_raw_sword:
            try:
                import shutil

                sword_dir = prj.dirs.get("sword")
                if sword_dir and os.path.isdir(sword_dir):
                    shutil.rmtree(sword_dir)
                    logger.info("Deleted SWORD gpkg folder (keep_raw_sword=False)")
            except Exception as e:
                logger.warning("Failed to delete SWORD gpkg folder: %s", e)
        else:
            logger.info("Kept SWORD gpkg folder (keep_raw_sword=True)")

    # Extract target IDs from subset
    source_id_col = "node_id" if prj.rivers.feature_type == "nodes" else "reach_id"
    if source_id_col not in subset.columns:
        raise KeyError(f"Expected SWORD column '{source_id_col}' missing from subset")

    prj.rivers.target_features = subset
    prj.rivers.target_id_col = source_id_col
    prj.rivers.target_ids = [int(value) for value in subset[source_id_col]]


def _ensure_sword_database(prj: "Project") -> None:
    """Ensure SWORD database is available locally.

    Checks if full SWORD database (GPKGs) already exists in prj.dirs["sword"].
    If not, handles three scenarios:
    1. User-provided zip: extract to {main_dir}/aux/SWORD/
    2. User-provided directory: use it directly
    3. Auto-download from Zenodo: download and extract to {main_dir}/aux/SWORD/

    Respects keep_raw_sword config to optionally delete downloaded zip.
    """
    from urllib import request as url_request
    import zipfile

    SWORD_V17B_ZIP_URL = (
        "https://zenodo.org/records/15299138/files/SWORD_v17b_gpkg.zip?download=1"
    )

    sword_dir = prj.dirs["sword"]

    # Check if SWORD database already exists
    if os.path.isdir(sword_dir):
        # Check if directory contains any SWORD GPKGs
        gpkg_files = [f for f in os.listdir(sword_dir) if f.endswith("_v17b.gpkg")]
        if gpkg_files:
            logger.info("SWORD database located at %s", sword_dir)
            return

    logger.info("SWORD database not found. Preparing it now.")

    # User-provided raw_sword_path
    if "sword_raw" in prj.dirs:
        raw_path = prj.dirs["sword_raw"]

        # Case 1: User provided a zip file
        if raw_path.lower().endswith(".zip") and os.path.isfile(raw_path):
            logger.info("Using user-provided SWORD zip: %s", raw_path)
            general.ifnotmakedirs(os.path.dirname(sword_dir))
            with zipfile.ZipFile(raw_path, "r") as zip_ref:
                zip_ref.extractall(os.path.dirname(sword_dir))
            logger.info("SWORD extracted to %s", sword_dir)
            return

        # Case 2: User provided a directory
        elif os.path.isdir(raw_path):
            logger.info("Using user-provided SWORD directory: %s", raw_path)
            # Check if GPKGs are in raw_path/gpkg/ or directly in raw_path
            gpkg_subdir = os.path.join(raw_path, "gpkg")
            if os.path.isdir(gpkg_subdir):
                prj.dirs["sword"] = gpkg_subdir
                logger.info("SWORD database found in %s", gpkg_subdir)
            else:
                prj.dirs["sword"] = raw_path
                logger.info("SWORD database found in %s", raw_path)
            return

    # Case 3: Auto-download from Zenodo
    logger.info("Downloading SWORD v17b from Zenodo...")
    general.ifnotmakedirs(os.path.dirname(sword_dir))

    zip_path = os.path.join(os.path.dirname(sword_dir), "SWORD_v17b_gpkg.zip")
    url_request.urlretrieve(SWORD_V17B_ZIP_URL, zip_path)
    logger.info("Downloaded SWORD v17b to %s", zip_path)

    logger.info("Extracting SWORD v17b...")
    with zipfile.ZipFile(zip_path, "r") as zip_ref:
        zip_ref.extractall(os.path.dirname(sword_dir))

    logger.info("SWORD extracted to %s", sword_dir)

    # Cleanup: delete zip if keep_raw_sword is False
    if not prj.keep_raw_sword:
        try:
            os.remove(zip_path)
            logger.info("Deleted raw SWORD zip file (keep_raw_sword=False)")
        except Exception as e:
            logger.warning("Failed to delete SWORD zip %s: %s", zip_path, e)
    else:
        logger.info("Kept raw SWORD zip file at %s (keep_raw_sword=True)", zip_path)


# ============================================================================
# RESERVOIRS: Download
# ============================================================================


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


def _group_river_targets_by_waterbody(prj: "Project") -> dict:
    """Return {waterbody_id: [target_id, ...]} grouping."""
    if prj.rivers.target_features is not None and len(prj.rivers.target_features) > 0:
        groups: dict = {}
        seen_target_ids: set = set()
        for _, row in prj.rivers.target_features.iterrows():
            target_id = int(row[prj.rivers.target_id_col])
            if target_id in seen_target_ids:
                continue
            seen_target_ids.add(target_id)
            wb_id = str(row[prj.rivers.id_key])
            groups.setdefault(wb_id, []).append(target_id)
        return groups

    if prj.rivers.configured_id:
        return {str(prj.rivers.configured_id): list(prj.rivers.target_ids)}

    raise ValueError(
        "Unable to group river targets by waterbody. "
        "Configure rivers.id or provide rivers.aoi_path with rivers.id_key."
    )


def _river_target_corridor(
    prj: "Project", target_ids, buffer_meters=None, width_buffer_factor=1.05,
):
    """
    Build one buffered, dissolved corridor polygon covering the given
    river targets (nodes or reaches), for use as the spatial AOI when
    downloading/extracting ICESat-2 and Sentinel-3/6 observations.

    This is deliberately a SEPARATE buffer distance from
    prj.rivers.buffer_meters (used earlier to decide which SWORD
    targets intersect the user's AOI at all) -- that question ("is this
    target in scope") and this one ("how far from the centerline could
    real river water still be, for a raw altimetry point to plausibly
    belong to this target") are different, and conflating them risks
    the same "one parameter doing two jobs badly" issue found elsewhere
    in this pipeline.

    Parameters
    ----------
    buffer_meters : float or None, optional
        Explicit, uniform buffer distance (meters), applied to every
        target regardless of its actual width. If None (default), uses
        each target's own SWORD "width" attribute instead: buffer
        distance = (width / 2) * width_buffer_factor. This is HALF the
        width, not the full width -- buffering a line expands it
        symmetrically by the given distance on EACH side, so a buffer
        of width/2 gives a corridor whose TOTAL span is approximately
        width * width_buffer_factor, matching the river's actual extent
        plus a margin, rather than doubling it. Falls back to
        _river_extraction_buffer_meters() (a flat scalar) if no usable
        "width" column is found -- e.g. if your SWORD data names it
        differently than assumed here, this degrades gracefully with a
        log message rather than failing.
    width_buffer_factor : float, optional
        Margin applied on top of each target's own width when using the
        width-based default. Default 1.05 -- 5% wider than the river's
        actual channel width. Only used when buffer_meters is None.

    NOTE: "width" is the expected SWORD column name per the standard
    SWORD data dictionary -- this has NOT been verified against a real
    downloaded SWORD file in this session (no sample data was
    available), unlike most other assumptions in this codebase. Check
    your actual target_features columns if width-based buffering
    doesn't seem to be kicking in.

    Returns a single-row GeoDataFrame in prj.global_crs (matching what
    icesat2.extract_observations/sentinel.extract_observations expect
    for their `features` argument, same as reservoirs), or None if no
    matching SWORD geometry is found for target_ids. Note the returned
    geometry may be a MultiPolygon if targets form disconnected pieces
    (e.g. separate reaches far enough apart that their buffers never
    touch) -- see _iter_geometry_pieces for how downloads handle this.
    """
    features = prj.rivers.target_features
    subset = features.loc[features[prj.rivers.target_id_col].isin(target_ids)]
    if subset.empty:
        return None

    local = subset.to_crs(prj.local_crs)

    if buffer_meters is not None:
        distances = buffer_meters
    elif "width" in local.columns and local["width"].notna().any():
        fallback_width = local["width"].median()
        distances = (local["width"].fillna(fallback_width) / 2) * width_buffer_factor
    else:
        logger.info(
            "No 'width' column found in SWORD target_features for this "
            "waterbody -- falling back to a flat extraction buffer "
            "instead of width-based sizing. Check your SWORD data's "
            "actual column names if this is unexpected."
        )
        distances = _river_extraction_buffer_meters(prj)

    buffered = local.buffer(distances)
    corridor = buffered.unary_union
    corridor_gdf = gpd.GeoDataFrame(
        geometry=[corridor], crs=prj.local_crs
    ).to_crs(prj.global_crs)
    return corridor_gdf


def _iter_geometry_pieces(geom):
    """
    Yield each individual polygon from a geometry: every part of a
    MultiPolygon, or the geometry itself for a plain Polygon. Used for
    river downloads, where disconnected corridor pieces should each get
    their own query rather than only querying the first piece (silently
    dropping coverage of the rest) or merging them into one shape that
    would also cover the (possibly large, irrelevant) gap between them.
    """
    if hasattr(geom, "geoms"):
        return list(geom.geoms)
    return [geom]


def _simplify_to_one_polygon(geom):
    """
    Collapse a MultiPolygon into a single encompassing polygon via
    convex hull. Used for reservoirs: unlike rivers, a reservoir is
    treated as one target regardless of how many disconnected parts its
    input polygon has, so a single combined query is preferred over
    splitting into several separate ones. Convex hull guarantees full
    coverage of every part, at the cost of also covering some in-between
    area that may not be real water -- an accepted tradeoff for treating
    one reservoir as one query rather than several.
    """
    if hasattr(geom, "geoms"):
        return geom.convex_hull
    return geom


def _river_extraction_buffer_meters(prj: "Project") -> float:
    """
    Resolve the extraction-corridor buffer distance: prefer an explicit
    prj.rivers.extraction_buffer_meters if set, else fall back to the
    SWORD-intersection prj.rivers.buffer_meters, else a conservative
    default. Kept as its own small function since this fallback chain
    is used by both download and extraction.
    """
    explicit = getattr(prj.rivers, "extraction_buffer_meters", None)
    if explicit:
        return explicit
    if prj.rivers.buffer_meters:
        return prj.rivers.buffer_meters
    return 500.0


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


def _assign_points_to_river_targets(
    points, targets, target_id_col, max_distance_meters, local_crs
):
    """
    Assign each point in `points` to its nearest feature in `targets`
    (SWORD node or reach geometries, whichever prj.rivers.target_id_col
    is configured for), dropping points farther than max_distance_meters
    from any target.

    Uses gpd.sjoin_nearest rather than a custom NearestNeighbors/DBSCAN
    approach -- it handles point-to-line matching natively (needed for
    reaches, not just nodes), and max_distance is expressed directly in
    real distance units once both inputs are reprojected to local_crs.

    Returns points (unprojected, original CRS) with target_id_col and a
    _dist_to_target_m column added; rows with no target within range are
    dropped entirely.
    """
    points_local = points.to_crs(local_crs)
    targets_local = targets[[target_id_col, "geometry"]].to_crs(local_crs)

    joined = gpd.sjoin_nearest(
        points_local,
        targets_local,
        how="inner",
        max_distance=max_distance_meters,
        distance_col="_dist_to_target_m",
    )

    result = points.loc[joined.index].copy()
    result[target_id_col] = joined[target_id_col].values
    result["_dist_to_target_m"] = joined["_dist_to_target_m"].values
    return result


def _extract_rivers_timeseries(prj: "Project", overwrite: bool = False) -> None:
    """Extract timeseries observations from raw downloaded files, for rivers.

    SWOT still needs a (lightweight) extraction step here: Hydrocron
    already returns a per-node/reach timeseries directly, but grouped
    per WATERBODY (one CSV covering every target in that waterbody) --
    see _extract_rivers_swot_observations for splitting that into the
    same per-target file structure ICESat-2/Sentinel-3/6 use, so the
    shared clean/merge pipeline can treat every mission identically.

    Parameters
    ----------
    overwrite : bool, optional
        Same semantics as _extract_reservoirs_timeseries: if False
        (default), any target whose output .gpkg already exists is
        skipped rather than re-extracted.
    """
    if "icesat2" in prj.to_process:
        _extract_rivers_icesat2_observations(prj, overwrite=overwrite)

    if "sentinel3" in prj.to_process:
        _extract_rivers_sentinel_observations(
            prj, "sentinel3", "S3", overwrite=overwrite
        )

    if "sentinel6" in prj.to_process:
        _extract_rivers_sentinel_observations(
            prj, "sentinel6", "S6", overwrite=overwrite
        )

    if "swot" in prj.to_process:
        _extract_rivers_swot_observations(prj, overwrite=overwrite)


def _extract_rivers_icesat2_observations(prj: "Project", overwrite: bool = False) -> None:
    """Extract ICESat-2 ATL13 observations for each river target.

    Unlike reservoirs (one polygon = one target), a river waterbody's
    raw download covers many targets at once. This extracts once per
    waterbody -- reusing icesat2.extract_observations exactly as
    reservoirs use it, with the buffered corridor (see
    _river_target_corridor) as the spatial filter instead of a single
    reservoir polygon -- then assigns each surviving point to its
    nearest target via sjoin_nearest, and splits the result into the
    same per-target {output}/{target_id}/raw_observations/icesat2.gpkg
    structure reservoirs already use, so everything downstream
    (clean/merge) can treat a river target exactly like a reservoir.
    """
    waterbody_groups = _group_river_targets_by_waterbody(prj)
    explicit_buffer = getattr(prj.rivers, "extraction_buffer_meters", None)
    width_buffer_factor = getattr(prj.rivers, "width_buffer_factor", 1.05)
    max_assign_dist = (
        getattr(prj.rivers, "max_node_assignment_meters", None)
        or _river_extraction_buffer_meters(prj)
    )

    tmp_dir = os.path.join(prj.dirs["output"], "_tmp_river_extraction")

    for wb_id, target_ids in waterbody_groups.items():
        parquet_dir = os.path.join(prj.dirs["icesat2_processed"], f"{wb_id}")
        if not os.path.exists(os.path.join(parquet_dir, "atl13.parquet")):
            continue

        if not overwrite:
            remaining = [
                t
                for t in target_ids
                if not os.path.exists(
                    os.path.join(
                        prj.dirs["output"], f"{t}", "raw_observations", "icesat2.gpkg"
                    )
                )
            ]
            if not remaining:
                continue
            target_ids = remaining

        corridor_gdf = _river_target_corridor(
            prj, target_ids, buffer_meters=explicit_buffer,
            width_buffer_factor=width_buffer_factor,
        )
        if corridor_gdf is None:
            continue

        general.ifnotmakedirs(tmp_dir)
        tmp_dst = os.path.join(tmp_dir, f"{wb_id}_icesat2.gpkg")

        try:
            icesat2.extract_observations(
                src_dir=parquet_dir,
                dst_path=tmp_dst,
                features=corridor_gdf,
                atl13_fields=prj.mission_options.get("icesat2", {}).get("atl13_fields"),
                track_keys=prj.mission_options.get("icesat2", {}).get("track_keys"),
            )
        except Exception as exc:
            logger.warning(
                "Failed to extract ICESat-2 for waterbody %s: %s", wb_id, exc
            )
            continue

        if not os.path.exists(tmp_dst):
            continue

        points = gpd.read_file(tmp_dst)
        os.remove(tmp_dst)
        if points.empty:
            continue

        targets = prj.rivers.target_features.loc[
            prj.rivers.target_features[prj.rivers.target_id_col].isin(target_ids)
        ]
        assigned = _assign_points_to_river_targets(
            points, targets, prj.rivers.target_id_col, max_assign_dist, prj.local_crs
        )
        if assigned.empty:
            logger.warning(
                "No ICESat-2 points within %sm of any target for waterbody %s",
                max_assign_dist, wb_id,
            )
            continue

        for target_id, group in assigned.groupby(prj.rivers.target_id_col):
            dst_dir = os.path.join(
                prj.dirs["output"], f"{target_id}", "raw_observations"
            )
            general.ifnotmakedirs(dst_dir)
            dst_path = os.path.join(dst_dir, "icesat2.gpkg")
            group.drop(
                columns=["index_right", "_dist_to_target_m"], errors="ignore"
            ).to_file(dst_path, driver="GPKG")


def _extract_rivers_sentinel_observations(
    prj: "Project", mission_key: str, product: str, overwrite: bool = False
) -> None:
    """Extract Sentinel-3 or Sentinel-6 observations for each river target.

    Same per-waterbody-then-split approach as
    _extract_rivers_icesat2_observations, plus a water-only filter: a
    buffered river corridor is much looser than a reservoir polygon (it
    genuinely includes riverbank, fields, vegetation alongside the
    channel), and unlike ICESat-2, Sentinel-3/6 have no built-in water
    classification. sigma0_min filters this out as a self-contained
    post-processing step here (rather than modifying
    sentinel.extract_observations itself, whose internals haven't been
    verified) -- water gives a strong, consistent specular radar
    return; land gives a weaker, noisier one. Needs empirical tuning
    against real river data, same as every other threshold in this
    pipeline -- the default here (0.0, i.e. no-op) is a safe starting
    point, not a verified value; set mission_options[mission_key]
    ['sigma0_min'] once you have real data to check it against.
    """
    waterbody_groups = _group_river_targets_by_waterbody(prj)
    explicit_buffer = getattr(prj.rivers, "extraction_buffer_meters", None)
    width_buffer_factor = getattr(prj.rivers, "width_buffer_factor", 1.05)
    max_assign_dist = (
        getattr(prj.rivers, "max_node_assignment_meters", None)
        or _river_extraction_buffer_meters(prj)
    )
    sigma0_min = prj.mission_options.get(mission_key, {}).get("sigma0_min", 0.0)

    tmp_dir = os.path.join(prj.dirs["output"], "_tmp_river_extraction")

    for wb_id, target_ids in waterbody_groups.items():
        download_dir = os.path.join(prj.dirs[mission_key], f"{wb_id}")
        if not os.path.exists(download_dir):
            continue

        if not overwrite:
            remaining = [
                t
                for t in target_ids
                if not os.path.exists(
                    os.path.join(
                        prj.dirs["output"],
                        f"{t}",
                        "raw_observations",
                        f"{mission_key}.gpkg",
                    )
                )
            ]
            if not remaining:
                continue
            target_ids = remaining

        corridor_gdf = _river_target_corridor(
            prj, target_ids, buffer_meters=explicit_buffer,
            width_buffer_factor=width_buffer_factor,
        )
        if corridor_gdf is None:
            continue

        general.ifnotmakedirs(tmp_dir)
        tmp_dst = os.path.join(tmp_dir, f"{wb_id}_{mission_key}.gpkg")

        try:
            sentinel.extract_observations(
                src_dir=download_dir,
                dst_path=tmp_dst,
                features=corridor_gdf,
                sigma0_max=prj.mission_options.get(mission_key, {}).get(
                    "sigma0_max", 1e5
                ),
            )
        except Exception as exc:
            logger.warning(
                "Failed to extract %s for waterbody %s: %s", mission_key, wb_id, exc
            )
            continue

        if not os.path.exists(tmp_dst):
            continue

        points = gpd.read_file(tmp_dst)
        os.remove(tmp_dst)
        if points.empty:
            continue

        if "sigma0" in points.columns and sigma0_min:
            before = len(points)
            points = points.loc[points["sigma0"] >= sigma0_min].reset_index(drop=True)
            logger.info(
                "%s waterbody %s: sigma0_min=%s kept %d/%d points",
                mission_key, wb_id, sigma0_min, len(points), before,
            )
        if points.empty:
            continue

        targets = prj.rivers.target_features.loc[
            prj.rivers.target_features[prj.rivers.target_id_col].isin(target_ids)
        ]
        assigned = _assign_points_to_river_targets(
            points, targets, prj.rivers.target_id_col, max_assign_dist, prj.local_crs
        )
        if assigned.empty:
            logger.warning(
                "No %s points within %sm of any target for waterbody %s",
                mission_key, max_assign_dist, wb_id,
            )
            continue

        for target_id, group in assigned.groupby(prj.rivers.target_id_col):
            dst_dir = os.path.join(
                prj.dirs["output"], f"{target_id}", "raw_observations"
            )
            general.ifnotmakedirs(dst_dir)
            dst_path = os.path.join(dst_dir, f"{mission_key}.gpkg")
            group.drop(
                columns=["index_right", "_dist_to_target_m"], errors="ignore"
            ).to_file(dst_path, driver="GPKG")


def _extract_rivers_swot_observations(prj: "Project", overwrite: bool = False) -> None:
    """
    Split Hydrocron's per-waterbody timeseries CSV into the same
    per-target {output}/{target_id}/raw_observations/swot.gpkg structure
    every other mission uses, so the shared clean/merge pipeline can
    treat SWOT identically to ICESat-2/Sentinel-3/6 for rivers.

    Unlike LakeSP for reservoirs, Hydrocron's own CSV doesn't include
    per-observation coordinates in the default field lists (see
    rivers.yaml) -- but nothing downstream actually needs per-observation
    lat/lon for SWOT (see PRODUCT_TIMESERIES_KEYS: no lat_key/lon_key for
    "swot"), so this attaches the target's own SWORD geometry as a
    constant placeholder purely so the file can be saved/read as .gpkg
    like every other mission's output -- the geometry's actual value is
    never used downstream, only the height/date/platform/orbit columns.

    Quality filtering (max_q) is already applied at download time (see
    _download_swot_hydrocron_timeseries), so it isn't repeated here.
    """
    waterbody_groups = _group_river_targets_by_waterbody(prj)
    id_label = "nodes" if prj.rivers.target_id_col == "node_id" else "reaches"

    for wb_id, target_ids in waterbody_groups.items():
        src_path = os.path.join(
            prj.dirs["swot"], str(wb_id), f"{id_label}_timeseries.csv"
        )
        if not os.path.exists(src_path):
            continue

        if not overwrite:
            remaining = [
                t
                for t in target_ids
                if not os.path.exists(
                    os.path.join(
                        prj.dirs["output"], f"{t}", "raw_observations", "swot.gpkg"
                    )
                )
            ]
            if not remaining:
                continue
            target_ids = remaining

        try:
            df = pd.read_csv(src_path)
        except Exception as exc:
            logger.warning(
                "Failed to read Hydrocron CSV for waterbody %s: %s", wb_id, exc
            )
            continue

        if df.empty or prj.rivers.target_id_col not in df.columns:
            continue

        df = df.loc[df[prj.rivers.target_id_col].isin(target_ids)].copy()
        if df.empty:
            continue

        df["height"] = df["wse"]
        df["date"] = pd.to_datetime(df["time_str"])
        df["platform"] = "swot"
        df["product"] = f"SWOT_Hydrocron_{id_label}"
        # Matches the reservoir SWOT convention (orbit = lake_id, constant
        # per target) -- there's no meaningful "which persistent track"
        # concept distinct from the target itself for Hydrocron data.
        df["orbit"] = df[prj.rivers.target_id_col]

        for target_id, group in df.groupby(prj.rivers.target_id_col):
            target_row = prj.rivers.target_features.loc[
                prj.rivers.target_features[prj.rivers.target_id_col] == target_id
            ]
            if target_row.empty:
                continue

            group = group.copy()
            group["geometry"] = target_row.geometry.iloc[0]
            gdf = gpd.GeoDataFrame(
                group, geometry="geometry", crs=prj.rivers.target_features.crs
            )

            dst_dir = os.path.join(
                prj.dirs["output"], f"{target_id}", "raw_observations"
            )
            general.ifnotmakedirs(dst_dir)
            gdf.to_file(os.path.join(dst_dir, "swot.gpkg"), driver="GPKG")


# ============================================================================
# RESERVOIRS: Timeseries Processing
# ============================================================================


def create_reservoirs_timeseries(prj: "Project") -> None:
    """Extract, clean, and merge timeseries for reservoirs.

    Parameters
    ----------
    prj : Project
        Project instance with reservoirs configuration
    """
    if not hasattr(prj, "reservoirs"):
        return

    # Extract raw observations from downloaded files (skips reservoirs whose
    # gpkg already exists, unless prj.reservoirs.overwrite_extraction=True --
    # confirmed this re-read/re-extraction was the dominant real-world cost,
    # far more than anything in clean()/merge())
    _extract_reservoirs_timeseries(
        prj, overwrite=getattr(prj.reservoirs, "overwrite_extraction", False)
    )

    # Clean observations with filters
    _clean_reservoirs_timeseries(prj)

    # Export to dfs0 if enabled
    if getattr(prj.reservoirs, "export_to_dfs0", False):
        _export_cleaned_to_dfs0(prj)

    # Merge multi-mission timeseries
    _merge_reservoirs_timeseries(prj)


def create_rivers_timeseries(prj: "Project") -> None:
    """Extract, clean, and merge timeseries for river targets (nodes/reaches).

    Mirrors create_reservoirs_timeseries. Not yet done: an export_to_dfs0
    equivalent for rivers, since _export_cleaned_to_dfs0 currently
    iterates prj.reservoirs.download_gdf specifically -- left out here
    rather than silently generalizing something not explicitly asked
    for yet.

    Parameters
    ----------
    prj : Project
        Project instance with rivers configuration
    """
    if not hasattr(prj, "rivers"):
        return

    _extract_rivers_timeseries(
        prj, overwrite=getattr(prj.rivers, "overwrite_extraction", False)
    )

    _clean_rivers_timeseries(prj)

    _merge_rivers_timeseries(prj)


def _extract_reservoirs_timeseries(prj: "Project", overwrite: bool = False) -> None:
    """Extract timeseries observations from raw downloaded files.

    Parameters
    ----------
    overwrite : bool, optional
        If False (default), any reservoir/mission whose output .gpkg
        already exists is skipped entirely rather than re-read and
        re-extracted. Confirmed on real data that this re-extraction --
        not clean()/merge() -- was the dominant cost in real end-to-end
        runs (orders of magnitude larger than the merge pipeline itself).
        Set True to force re-extraction (e.g. new raw downloads arrived).
    """
    if "icesat2" in prj.to_process:
        _extract_icesat2_observations(prj, overwrite=overwrite)

    if "sentinel3" in prj.to_process:
        _extract_sentinel_observations(prj, "sentinel3", "S3", overwrite=overwrite)

    if "sentinel6" in prj.to_process:
        _extract_sentinel_observations(prj, "sentinel6", "S6", overwrite=overwrite)

    if "swot" in prj.to_process:
        _extract_swot_observations(prj, overwrite=overwrite)


def _extract_icesat2_observations(prj: "Project", overwrite: bool = False) -> None:
    """Extract ICESat-2 ATL13 observations for each reservoir."""
    available_ids = [
        id
        for id in prj.reservoirs.download_gdf[prj.reservoirs.id_key]
        if os.path.exists(
            os.path.join(prj.dirs["icesat2_processed"], f"{id}", "atl13.parquet")
        )
    ]
    if not available_ids:
        logger.warning("No ICESat-2 downloads found; skipping timeseries extraction.")
        return

    if not overwrite:
        skip_count = 0
        remaining_ids = []
        for id in available_ids:
            dst_path = os.path.join(
                prj.dirs["output"], f"{id}", "raw_observations", "icesat2.gpkg"
            )
            if os.path.exists(dst_path):
                skip_count += 1
            else:
                remaining_ids.append(id)
        if skip_count:
            logger.info(
                "ICESat-2 extraction: skipping %d reservoir(s) with an existing "
                "icesat2.gpkg (pass overwrite=True to force re-extraction).",
                skip_count,
            )
        available_ids = remaining_ids
        if not available_ids:
            return

    empty_ids = []
    for id in tqdm(available_ids, desc="Extracting ICESat-2 ATL13 product"):
        sub_gdf = prj.reservoirs.download_gdf.loc[
            prj.reservoirs.download_gdf[prj.reservoirs.id_key] == id
        ]
        download_dir = os.path.join(prj.dirs["icesat2_processed"], f"{id}")
        dst_dir = os.path.join(prj.dirs["output"], f"{id}", "raw_observations")
        general.ifnotmakedirs(dst_dir)
        dst_path = os.path.join(dst_dir, "icesat2.gpkg")

        try:
            icesat2.extract_observations(
                src_dir=download_dir,
                dst_path=dst_path,
                features=sub_gdf,
                atl13_fields=prj.mission_options.get("icesat2", {}).get("atl13_fields"),
                track_keys=prj.mission_options.get("icesat2", {}).get("track_keys"),
            )
        except Exception as exc:
            logger.warning("Failed to extract ICESat-2 for %s: %s", id, exc)

        if not os.path.exists(dst_path):
            empty_ids.append(id)

    if empty_ids:
        logger.warning(
            "ICESat-2 timeseries empty for: %s (no observations passed the spatial filter or the download returned no data)",
            ", ".join(str(i) for i in empty_ids),
        )


def _extract_sentinel_observations(
    prj: "Project", mission_key: str, product: str, overwrite: bool = False
) -> None:
    """Extract Sentinel-3 or Sentinel-6 observations for each reservoir."""
    available_ids = [
        id
        for id in prj.reservoirs.download_gdf[prj.reservoirs.id_key]
        if os.path.exists(os.path.join(prj.dirs[mission_key], f"{id}"))
    ]
    if not available_ids:
        logger.warning(
            "No %s downloads found; skipping timeseries extraction.", mission_key
        )
        return

    if not overwrite:
        skip_count = 0
        remaining_ids = []
        for id in available_ids:
            dst_path = os.path.join(
                prj.dirs["output"], f"{id}", "raw_observations", f"{mission_key}.gpkg"
            )
            if os.path.exists(dst_path):
                skip_count += 1
            else:
                remaining_ids.append(id)
        if skip_count:
            logger.info(
                "%s extraction: skipping %d reservoir(s) with an existing "
                "%s.gpkg (pass overwrite=True to force re-extraction).",
                mission_key, skip_count, mission_key,
            )
        available_ids = remaining_ids
        if not available_ids:
            return

    empty_ids = []
    for id in tqdm(available_ids, desc=f"Extracting Sentinel-{product} product"):
        sub_gdf = prj.reservoirs.download_gdf.loc[
            prj.reservoirs.download_gdf[prj.reservoirs.id_key] == id
        ]
        download_dir = os.path.join(prj.dirs[mission_key], f"{id}")
        dst_dir = os.path.join(prj.dirs["output"], f"{id}", "raw_observations")
        general.ifnotmakedirs(dst_dir)
        dst_path = os.path.join(dst_dir, f"{mission_key}.gpkg")

        try:
            sentinel.extract_observations(
                src_dir=download_dir,
                dst_path=dst_path,
                features=sub_gdf,
                sigma0_max=prj.mission_options.get(mission_key, {}).get(
                    "sigma0_max", 1e5
                ),
            )
        except Exception as exc:
            logger.warning("Failed to extract %s for %s: %s", mission_key, id, exc)

        if not os.path.exists(dst_path):
            empty_ids.append(id)

    if empty_ids:
        logger.warning(
            "Sentinel-%s timeseries empty for: %s (no observations passed the spatial filter or the download returned no data)",
            product,
            ", ".join(str(i) for i in empty_ids),
        )


def _extract_swot_observations(prj: "Project", overwrite: bool = False) -> None:
    """Extract SWOT Lake SP observations for all reservoirs."""
    download_dir = prj.dirs["swot"]
    if not os.path.exists(download_dir):
        logger.warning("No SWOT downloads found; skipping timeseries extraction.")
        return

    features = prj.reservoirs.download_gdf
    id_key = prj.reservoirs.id_key

    if not overwrite:
        def _has_output(id):
            return os.path.exists(
                os.path.join(prj.dirs["output"], f"{id}", "raw_observations", "swot.gpkg")
            )
        all_ids = features[id_key].tolist()
        remaining_ids = [i for i in all_ids if not _has_output(i)]
        skip_count = len(all_ids) - len(remaining_ids)
        if skip_count:
            logger.info(
                "SWOT extraction: skipping %d reservoir(s) with an existing "
                "swot.gpkg (pass overwrite=True to force re-extraction).",
                skip_count,
            )
        if not remaining_ids:
            return
        features = features.loc[features[id_key].isin(remaining_ids)]

    empty_ids = swot.extract_observations(
        src_dir=download_dir,
        dst_dir=prj.dirs["output"],
        dst_file_name="swot.gpkg",
        features=features,
        id_key=id_key,
        exclude_obs_id_values=prj.mission_options.get("swot", {}).get(
            "exclude_obs_id_values", ["no_data"]
        ),
    )
    if empty_ids:
        logger.warning(
            "SWOT timeseries empty for: %s (no observations matched the prior lake ID or all were excluded)",
            ", ".join(str(i) for i in empty_ids),
        )


# Per-product mapping from generic Timeseries key attributes to the actual
# column names each mission's extractor writes. Sentinel-3 shares
# Sentinel-6's extractor/schema (same sentinel.extract_observations
# function, see _extract_sentinel_observations).
#
# ****************************************************************************
# TERMINOLOGY TRAP -- read before touching orbit_key/pass_key for Sentinel:
# The raw Sentinel-3/6 data has TWO similarly-named but opposite-meaning
# columns:
#   - "orbit": the absolute revolution counter. Unique on every single
#     crossing, never repeats. USELESS as orbit_key (bias_correct needs a
#     persistent identifier to accumulate overlap against -- grouping by
#     something that's different every time means every "source" has
#     exactly 1 observation and nothing can ever be calibrated: this is
#     exactly the bug that caused every S3A/S3B track to be dropped as
#     unanchored in practice).
#   - "pass": the satellite-engineering term for the STABLE, REPEATING
#     ground track number (same value every ~27-day repeat cycle for
#     S3A/S3B). This is what orbit_key actually needs.
# Confusingly, our own framework's `pass_key` means the OPPOSITE thing (one
# specific, one-time crossing -- e.g. file_name) from what "pass" means in
# the satellite data itself (the repeating track). Do not be tempted to
# point pass_key at the raw "pass" column -- file_name is correct there.
# ****************************************************************************
#
# ICESat-2's orbit_key is "beam" (the persistent ground track/virtual
# station) -- cycle_number only matters as an ingredient of the compound
# "pass" column built at extraction time (see
# HydroEO.satellites.icesat2.preprocess.extract_observations). SWOT's
# LakeSP product is already one integrated WSE per crossing with its own
# formal uncertainty (wse_u), so it needs neither lat/lon nor pass_key --
# see preset_error_key, and daily_mad_error's handling of it.
PRODUCT_TIMESERIES_KEYS = {
    "sentinel3": dict(
        lat_key="lat", lon_key="lon", pass_key="file_name",
        platform_key="platform", orbit_key="relative_orbit",
    ),
    # NOTE: sentinel6 still uses "pass" as orbit_key -- NOT verified to be
    # unstable the way it was for sentinel3 (confirmed empirically: on
    # real data, "pass" was unique-per-crossing for every S3A/S3B visit,
    # i.e. not stable at all, while "relative_orbit" genuinely repeated
    # across multiple visits -- e.g. S3B crossed via 2 distinct stable
    # configurations, with real biases of -0.14m and +0.22m that a
    # platform-only grouping was averaging into one misleading +0.04m).
    # Sentinel-6 may have the same "pass" instability and may also have
    # its own "relative_orbit"-equivalent column, but this hasn't been
    # checked against real S6 data -- don't assume the same fix applies
    # without verifying first.
    "sentinel6": dict(
        lat_key="lat", lon_key="lon", pass_key="file_name",
        platform_key="platform", orbit_key="pass",
    ),
    "icesat2": dict(
        lat_key="lat", lon_key="lon", pass_key="pass",
        platform_key="platform", orbit_key="beam",
    ),
    "swot": dict(
        platform_key="platform", orbit_key="orbit", preset_error_key="wse_u",
    ),
}

# Default .merge() tuning for reservoirs, mirroring the shape of
# processing_options (a project-level dict of pipeline parameters) but
# applied once per reservoir rather than per-product, since .merge() runs
# on the already-combined multi-product timeseries. Override via
# prj.merging_options in project config; falls back to these reservoir-
# appropriate defaults if that attribute isn't set.
DEFAULT_RESERVOIR_MERGING_OPTIONS = {
    "window_km": 1.5,
    "svr_linear_err": 0.1,
    "svr_linear_epsilon": 0.1,
    # Both updated from DAHITI's lake-tuned defaults (err=0.1, gamma=
    # 0.0000438) based on real reservoir data validated this session --
    # the lake-tuned gamma implied a ~151-day smoothing lengthscale, far
    # too coarse for a reservoir with real multi-week transitions (see
    # the svr_radial oversmoothing discussion). err=1.0 (river-like,
    # rather than the stricter lake value) and gamma x50 (~21-day
    # lengthscale instead of ~151 days) let the trend actually track
    # real fast changes instead of rejecting them as if they were noise.
    "svr_radial_err": 1.0,
    "svr_radial_rbf_c": 10000,
    "svr_radial_gamma": 0.0000438 * 50,
    "svr_radial_epsilon": 0.1,
    # Confirmed on real data across two reservoirs: revisit sparsity varies a
    # lot (e.g. one reservoir had icesat2/S3A/S3B visiting only 7/14/13
    # distinct days all year). At "10D"/3, sparse sources can fail to ever
    # find 3 overlapping bins and get dropped as unanchored ENTIRELY (not
    # just trimmed) -- confirmed: this silently dropped 2 of 3 missions
    # (icesat2, S3B) for one real reservoir. "20D"/1 recovered all of it.
    # Widening is monotonically safe against data loss (a wider window can
    # only find equal-or-more overlapping bins, never fewer) -- the
    # tradeoff is a very wide bin could blur real water-level change within
    # the window into the bias estimate; 20D is a modest widening, not an
    # extreme one.
    "bias_time_bin": "20D",
    "bias_min_overlap": 1,
    # Confirmed empirically on real data: "platform_orbit" (using
    # orbit_key -- now sentinel3's verified-stable "relative_orbit"
    # column, see PRODUCT_TIMESERIES_KEYS) reveals genuine within-platform
    # bias heterogeneity that "platform" alone was masking. One real
    # reservoir's S3B crosses via two distinct, independently stable
    # configurations (5 days on one, 8 on the other) with biases of
    # -0.14m and +0.22m respectively -- "platform" grouping averaged
    # these into one misleading +0.04m. Same pattern for ICESat-2's
    # beams (orbit_key="beam"): per-beam biases ranged 0.06-0.18m under
    # "platform_orbit", collapsed to one number under "platform". Total
    # kept-row count was IDENTICAL either way on the reservoir tested
    # (3182/4901) -- this is a precision gain, not a data-loss risk, at
    # least for sentinel3/icesat2. NOTE: sentinel6 still uses "pass" as
    # orbit_key (unverified whether it's stable or has a
    # relative_orbit-equivalent -- see PRODUCT_TIMESERIES_KEYS) -- if
    # it's actually unstable like sentinel3's old "pass" mapping was,
    # "platform_orbit" could fragment sentinel6 into single-crossing
    # sources. Recheck against real sentinel6 data before trusting this
    # default for a project relying heavily on sentinel6.
    "bias_group_by": "platform_orbit",
    # Not a spatial correction -- just flags (and records in
    # ts.bias_correct_diagnostics) when a source's observations are
    # centered far from the anchor's, since for a large/elongated
    # reservoir some of the estimated bias could be real spatial signal.
    # Worth a closer look per-reservoir if this fires, not an error.
    "bias_centroid_warn_km": 5.0,
    # Off by default -- inflates Kalman input error by distance from the
    # reservoir polygon's own centroid, addressing crossings that may be
    # hydraulically unrepresentative (e.g. far upstream, subject to real
    # slope bias) even when ADM alone reports them as highly precise. Set
    # to a real value (m of extra error per km of distance) to enable --
    # the right scale depends on the true magnitude of upstream slope bias
    # for your reservoirs, which needs empirical tuning, not a guessed
    # default.
    "distance_penalty_scale_per_km": None,
    # Off by default -- a genuine height correction (not just error
    # inflation) using a spatial deviation model fit once from a dense
    # source (default ICESat-2) and persisted to disk per reservoir (see
    # _get_or_fit_spatial_correction_model) so past corrections don't
    # shift retroactively as new data arrives. Turn on once you've
    # confirmed (as we did empirically) that the target reservoir shows a
    # real, day-to-day-consistent spatial deviation pattern -- fitting
    # requires several qualifying dense-source days (see
    # fit_spatial_correction_model's min_days), and silently does nothing
    # if there isn't enough dense-source data yet.
    "use_spatial_correction": False,
    "spatial_correction_dense_source": "icesat2",
}


def _clean_timeseries(prj: "Project", target_type: str) -> None:
    """Apply quality filters to extracted timeseries, for either
    reservoirs or river targets (nodes/reaches)."""
    target_ids = _get_target_ids(prj, target_type)
    ids_with_raw = [
        id
        for id in target_ids
        if os.path.exists(os.path.join(prj.dirs["output"], f"{id}", "raw_observations"))
    ]
    if not ids_with_raw:
        logger.warning(
            "No raw observations found for any %s; skipping timeseries cleaning.",
            target_type,
        )
        return

    for id in tqdm(ids_with_raw, desc=f"Cleaning product timeseries ({target_type})"):
        for product in prj.to_process:
            df = _load_product_timeseries(
                os.path.join(prj.dirs["output"], f"{id}", "raw_observations"),
                ".gpkg",
                [product],
                lambda path: gpd.read_file(path).drop(columns=["geometry"]),
            )
            if df is not None:
                product_options = prj.processing_options.get(
                    product,
                    {
                        "processing_filters": ["elevation", "MAD"],
                        "elevation_min_m": 0.0,
                        "elevation_max_m": 8000.0,
                        "mad_threshold": 5.0,
                    },
                )

                ts = timeseries.Timeseries(
                    df, date_key="date", height_key="height",
                    **PRODUCT_TIMESERIES_KEYS.get(product, {}),
                )

                ts.clean(
                    product_options.get("processing_filters", ["elevation", "MAD"]),
                    filter_params={
                        "elevation_min_m": product_options.get("elevation_min_m", 0.0),
                        "elevation_max_m": product_options.get(
                            "elevation_max_m", 8000.0
                        ),
                        "mad_threshold": product_options.get("mad_threshold", 5.0),
                    },
                )

                export_dir = os.path.join(
                    prj.dirs["output"], f"{id}", "cleaned_observations"
                )
                general.ifnotmakedirs(export_dir)
                ts.export_csv(os.path.join(export_dir, f"{product}.csv"))


def _clean_reservoirs_timeseries(prj: "Project") -> None:
    """Apply quality filters to extracted reservoir timeseries."""
    _clean_timeseries(prj, "reservoirs")


def _clean_rivers_timeseries(prj: "Project") -> None:
    """Apply quality filters to extracted river timeseries."""
    _clean_timeseries(prj, "rivers")


def _load_product_timeseries(data_dir, ext, products, reader_fn):
    """Load files of given extension from directory, optionally filtered to products."""
    if not os.path.exists(data_dir):
        return None
    df_list = []
    for file in os.listdir(data_dir):
        if file.endswith(ext):
            if not products or file.split(".")[0] in products:
                try:
                    df_list.append(reader_fn(os.path.join(data_dir, file)))
                except Exception as exc:
                    logger.warning(
                        "Failed to load %s from %s: %s",
                        file,
                        data_dir,
                        exc,
                    )
    return pd.concat(df_list) if df_list else None


def _get_target_ids(prj: "Project", target_type: str):
    """Return the list of target IDs to process, for either 'reservoirs' or 'rivers'."""
    if target_type == "reservoirs":
        return list(prj.reservoirs.download_gdf[prj.reservoirs.id_key])
    if target_type == "rivers":
        return list(prj.rivers.target_ids)
    raise ValueError(f"Unknown target_type: {target_type!r}")


def _target_centroid(prj: "Project", target_type: str, id):
    """
    Return (lat, lon) of a target's own geometry centroid -- the
    reservoir polygon for target_type="reservoirs", or the SWORD
    node/reach geometry for target_type="rivers" -- computed in a
    projected (local) CRS for accuracy, then converted back to lat/lon.
    Used as the reference location for apply_distance_penalty/
    apply_spatial_correction. Returns (None, None) if the target's
    geometry can't be found, so callers can treat that as "skip" rather
    than fail.
    """
    try:
        if target_type == "reservoirs":
            gdf = prj.reservoirs.gdf
            id_key = prj.reservoirs.id_key
        elif target_type == "rivers":
            gdf = prj.rivers.target_features
            id_key = prj.rivers.target_id_col
        else:
            raise ValueError(f"Unknown target_type: {target_type!r}")

        row = gdf.loc[gdf[id_key] == id]
        if len(row) == 0 or row.geometry.isna().all():
            return None, None
        centroid = row.to_crs(prj.local_crs).geometry.centroid.to_crs(prj.global_crs)
        pt = centroid.iloc[0]
        return pt.y, pt.x  # lat, lon
    except Exception as exc:
        logger.warning(
            "Could not compute %s centroid for %s: %s", target_type, id, exc
        )
        return None, None


def _reservoir_centroid(prj: "Project", id):
    """Backward-compatible wrapper -- see _target_centroid."""
    return _target_centroid(prj, "reservoirs", id)


def _get_or_fit_spatial_correction_model(
    prj: "Project", target_type: str, id, dense_source_platform="icesat2",
    recalibrate=False, **fit_kwargs,
):
    """
    Load a persisted spatial correction model for this target if one
    exists, or fit a fresh one and persist it. Works identically for
    reservoirs and river targets -- see _target_centroid.

    This is deliberately NOT re-fit automatically every run: doing so
    would make past corrections shift retroactively every time new
    dense-source data arrives, since the fitted slope would change.
    Pass recalibrate=True to explicitly force a re-fit (e.g. as a
    deliberate, occasional recalibration step) -- not something that
    should happen as a silent side effect of routine reprocessing.

    Returns None if no model exists yet and there isn't enough dense
    source data to fit one (see fit_spatial_correction_model) -- callers
    should treat this the same as "no correction available."
    """
    model_path = os.path.join(
        prj.dirs["output"], f"{id}", "spatial_correction_model.json"
    )

    if os.path.exists(model_path) and not recalibrate:
        with open(model_path, "r") as f:
            return json.load(f)

    cleaned_path = os.path.join(
        prj.dirs["output"], f"{id}", "all_cleaned_timeseries.csv"
    )
    if not os.path.exists(cleaned_path):
        logger.info(
            "No cleaned observations yet for %s; cannot fit spatial "
            "correction model.", id,
        )
        return None

    df = pd.read_csv(cleaned_path)
    if "platform" not in df.columns or dense_source_platform not in df["platform"].values:
        logger.info(
            "No %s data available for %s; cannot fit spatial correction "
            "model from it.", dense_source_platform, id,
        )
        return None

    df["date"] = pd.to_datetime(df["date"])
    dense_df = df.loc[df["platform"] == dense_source_platform]

    ref_lat, ref_lon = _target_centroid(prj, target_type, id)
    if ref_lat is None:
        logger.warning(
            "Could not determine %s centroid for %s; cannot fit spatial "
            "correction model.", target_type, id,
        )
        return None

    model = basic_filters.fit_spatial_correction_model(
        dense_df, lat_key="lat", lon_key="lon", height_key="height",
        date_key="date", ref_lat=ref_lat, ref_lon=ref_lon, **fit_kwargs,
    )

    if model is not None:
        general.ifnotmakedirs(os.path.dirname(model_path))
        with open(model_path, "w") as f:
            json.dump(model, f, indent=2)

    return model


# ============================================================================
# Per-target run config: exclusions + per-target merging option overrides
# ============================================================================
#
# One YAML file per target ({output}/{id}/run_config.yaml) that is
# simultaneously: (a) a human-readable log of decisions made about this
# target, (b) the actual source of truth _merge_timeseries reads to apply
# those decisions, and (c) something a user can hand-edit directly for a
# fully config-driven workflow. Interactive functions below
# (exclude_from_target, set_merging_option, ...) read-modify-write this
# same file, so a decision made once in a notebook session is exactly the
# same artifact you'd edit by hand or check into version control -- there
# is no separate "notebook state" to keep in sync with "the config".


def _run_config_path(prj: "Project", id) -> str:
    return os.path.join(prj.dirs["output"], f"{id}", "run_config.yaml")


def _default_run_config(id) -> dict:
    return {
        "target_id": id,
        "last_updated": None,
        "merging_option_overrides": {},
        "exclusions": [],
    }


def _load_run_config(prj: "Project", id) -> dict:
    """Load a target's run_config.yaml, or a fresh default if none exists yet."""
    path = _run_config_path(prj, id)
    if os.path.exists(path):
        with open(path, "r") as f:
            loaded = yaml.safe_load(f)
        if loaded:
            # tolerate a hand-edited file missing a key or two
            defaults = _default_run_config(id)
            defaults.update(loaded)
            return defaults
    return _default_run_config(id)


def _save_run_config(prj: "Project", id, config: dict) -> None:
    config["last_updated"] = datetime.datetime.now().isoformat()
    path = _run_config_path(prj, id)
    general.ifnotmakedirs(os.path.dirname(path))
    with open(path, "w") as f:
        yaml.safe_dump(config, f, sort_keys=False)


def _invalidate_spatial_correction_cache(prj: "Project", id) -> None:
    """
    Delete any cached spatial correction model for this target, forcing
    a fresh fit next time use_spatial_correction is used. Called whenever
    exclusions or spatial-correction-relevant options change -- the
    model may have been fit using observations that are no longer
    included, and this is exactly the kind of deliberate, explicit
    trigger (not routine reprocessing) that recalibration is meant for --
    see _get_or_fit_spatial_correction_model.
    """
    model_path = os.path.join(prj.dirs["output"], f"{id}", "spatial_correction_model.json")
    if os.path.exists(model_path):
        os.remove(model_path)
        logger.info(
            "Invalidated cached spatial correction model for %s "
            "(exclusions or related options changed).", id,
        )


def _invalidate_reach_slope_correction_cache(prj: "Project", id) -> None:
    """
    Delete any cached reach slope correction model for this target,
    forcing a fresh fit next time use_reach_slope_correction is used.
    Called whenever exclusions change -- an exclusion could target SWOT
    observations specifically, which is exactly what this model is fit
    from (see _fit_reach_slope_correction), so a cached model could
    otherwise silently keep reflecting now-excluded SWOT slope values.
    """
    model_path = os.path.join(
        prj.dirs["output"], f"{id}", "reach_slope_correction_model.json"
    )
    if os.path.exists(model_path):
        os.remove(model_path)
        logger.info(
            "Invalidated cached reach slope correction model for %s "
            "(exclusions or related options changed).", id,
        )


def _fit_reach_slope_correction(prj: "Project", target_id, recalibrate: bool = False):
    """
    Fit (or load a persisted) reach-level slope correction from SWOT's
    own directly-measured "slope" field (RiverSP reach product), used to
    reference-correct OTHER missions' (ICESat-2/Sentinel-3/6) crossings
    to what they'd read at the reach's geometric midpoint.

    ONLY meaningful for reaches (prj.rivers.target_id_col == "reach_id")
    -- a node is a single ~200m-spaced point, not a ~10km segment with
    its own along-reach slope in the same sense. Callers must not invoke
    this for node-mode projects.

    Rationale: SWOT's reach-level WSE is an aggregate over the reach's
    ~50 constituent, roughly-evenly-spaced nodes, not a value evaluated
    at one specific point -- for an evenly-sampled linear profile, the
    mean equals the value at the mean position, so this is treated as
    approximately midpoint-referenced. This is an evidence-based
    inference from the RiverSP processing chain, NOT a fact directly
    confirmed in SWOT's product documentation (which does not explicitly
    state a reference point) -- validate against real Hydrocron
    node-vs-reach output for a known reach before trusting this deeply.

    Uses the MEDIAN of all available SWOT slope observations for this
    reach as a single, persistent correction -- not a per-date-specific
    one -- consistent with this pipeline's existing preference (see
    fit_spatial_correction_model) for a stable, once-fit value over a
    per-observation one, and avoiding the complexity/fragility of
    matching a specific SWOT overpass date to each individual non-SWOT
    observation's date.

    Persisted to {output}/{target_id}/reach_slope_correction_model.json
    -- fit once, loaded thereafter, only refit on explicit
    recalibrate=True -- so past corrections don't shift retroactively
    as new SWOT data arrives, same reasoning as the spatial correction
    model's caching.

    Returns None if no model exists yet and there's no usable SWOT slope
    data to fit one from -- callers should treat this as "no correction
    available," not an error.
    """
    model_path = os.path.join(
        prj.dirs["output"], f"{target_id}", "reach_slope_correction_model.json"
    )
    if os.path.exists(model_path) and not recalibrate:
        with open(model_path, "r") as f:
            return json.load(f)

    swot_path = os.path.join(
        prj.dirs["output"], f"{target_id}", "raw_observations", "swot.gpkg"
    )
    if not os.path.exists(swot_path):
        logger.info(
            "No raw SWOT observations for %s; cannot fit reach slope "
            "correction.", target_id,
        )
        return None

    swot_gdf = gpd.read_file(swot_path)
    if "slope" not in swot_gdf.columns:
        logger.warning(
            "No 'slope' field found in SWOT observations for %s -- was "
            "it requested in mission_options['swot']['hydrocron_fields']"
            "['reaches']? Cannot fit reach slope correction.", target_id,
        )
        return None

    valid_slopes = pd.to_numeric(swot_gdf["slope"], errors="coerce").dropna()
    if valid_slopes.empty:
        logger.info(
            "No valid (non-null, numeric) SWOT slope observations for "
            "%s; cannot fit reach slope correction.", target_id,
        )
        return None

    model = {
        "target_id": str(target_id),
        "median_slope": float(valid_slopes.median()),
        "n_observations": int(len(valid_slopes)),
        "fitted_at": datetime.datetime.now().isoformat(),
    }

    general.ifnotmakedirs(os.path.dirname(model_path))
    with open(model_path, "w") as f:
        json.dump(model, f, indent=2)

    return model


def _apply_reach_slope_correction(
    ts_df: pd.DataFrame, prj: "Project", target_id, slope_model: dict,
) -> pd.DataFrame:
    """
    Apply a fitted reach slope correction (see _fit_reach_slope_correction)
    to non-SWOT rows in ts_df -- adjusts "height" to what each row would
    read at the reach's geometric midpoint, using its along-reach
    projected position and the reach's persistent median slope. SWOT's
    own rows are left untouched (already assumed midpoint-referenced --
    see _fit_reach_slope_correction's docstring for the reasoning and
    its caveats).

    NOTE: the sign convention for "correction = slope x distance" here
    has NOT been empirically verified against real data in this
    session -- confirm it actually reduces cross-mission scatter for a
    real reach (not increases it) before trusting this in production;
    flip the sign if it doesn't.

    Rows without usable lat/lon (or if the target's geometry can't be
    found) are left uncorrected rather than dropped.
    """
    if "platform" not in ts_df.columns or "lat" not in ts_df.columns or "lon" not in ts_df.columns:
        return ts_df

    target_row = prj.rivers.target_features.loc[
        prj.rivers.target_features[prj.rivers.target_id_col] == target_id
    ]
    if target_row.empty:
        logger.warning(
            "Could not find reach geometry for %s; skipping reach slope "
            "correction.", target_id,
        )
        return ts_df

    reach_geom_global = target_row.geometry.iloc[0]
    local_crs = prj.local_crs
    reach_geom_local = (
        gpd.GeoSeries([reach_geom_global], crs=target_row.crs).to_crs(local_crs).iloc[0]
    )
    reach_midpoint_dist = reach_geom_local.length / 2.0
    median_slope = slope_model["median_slope"]

    mask = (ts_df["platform"] != "swot") & ts_df["lat"].notna() & ts_df["lon"].notna()
    if not mask.any():
        return ts_df

    points_local = gpd.GeoSeries(
        gpd.points_from_xy(ts_df.loc[mask, "lon"], ts_df.loc[mask, "lat"]),
        crs=prj.global_crs,
    ).to_crs(local_crs)

    along_reach_dist = points_local.apply(reach_geom_local.project)
    distance_from_midpoint = along_reach_dist.values - reach_midpoint_dist

    ts_df = ts_df.copy()
    ts_df.loc[mask, "height"] = (
        ts_df.loc[mask, "height"] - median_slope * distance_from_midpoint
    )
    return ts_df


def list_target_observations(prj: "Project", target_type: str, id) -> pd.DataFrame:
    """
    Summarize what observations exist for a target, at (platform, orbit)
    granularity -- the "what could I exclude" view. Meant to be read
    alongside plot_merging's platform-colored progress plots (which show
    WHERE a problem shows up in the actual data), not a replacement for
    looking at the data itself.
    """
    cleaned_path = os.path.join(prj.dirs["output"], f"{id}", "all_cleaned_timeseries.csv")
    if not os.path.exists(cleaned_path):
        logger.warning(
            "No cleaned observations yet for %s -- has create_%s_timeseries() "
            "been run?", id, target_type,
        )
        return pd.DataFrame(columns=["platform", "orbit", "n_points", "date_min", "date_max"])

    df = pd.read_csv(cleaned_path)
    df["date"] = pd.to_datetime(df["date"])
    group_cols = [c for c in ["platform", "orbit"] if c in df.columns]
    summary = (
        df.groupby(group_cols)
        .agg(n_points=("date", "size"), date_min=("date", "min"), date_max=("date", "max"))
        .reset_index()
        .sort_values(group_cols)
        .reset_index(drop=True)
    )
    return summary


def list_exclusions(prj: "Project", target_type: str, id) -> list:
    """Current exclusion rules for a target, from its run_config.yaml."""
    return _load_run_config(prj, id)["exclusions"]


def exclude_from_target(
    prj: "Project", target_type: str, id,
    platform=None, orbit=None, date=None, reason=None,
) -> None:
    """
    Exclude observations from a target's merge, at whatever granularity
    is given -- a whole platform, a specific orbit/pass value, a specific
    date, or any combination (all given fields must match for a row to
    be excluded). Persisted to {output}/{id}/run_config.yaml.

    Applied at the start of _merge_timeseries, before any processing --
    an excluded pass never reaches bias_correct/Kalman/svr_radial at all,
    rather than being fought against downstream.

    Invalidates any cached spatial correction model for this target,
    since it may have been fit using data that's now excluded.

    Examples
    --------
    exclude_from_target(prj, "reservoirs", my_id, platform="S3B")
    exclude_from_target(prj, "reservoirs", my_id, platform="S3B", orbit=1517)
    exclude_from_target(prj, "rivers", my_id, date="2024-03-19")
    """
    if platform is None and orbit is None and date is None:
        raise ValueError(
            "Specify at least one of platform, orbit, or date to exclude."
        )

    config = _load_run_config(prj, id)
    config["exclusions"].append({
        "platform": platform,
        "orbit": orbit,
        "date": str(date) if date is not None else None,
        "reason": reason,
        "added": datetime.datetime.now().isoformat(),
    })
    _save_run_config(prj, id, config)
    _invalidate_spatial_correction_cache(prj, id)
    _invalidate_reach_slope_correction_cache(prj, id)
    logger.info(
        "Added exclusion for %s: platform=%s orbit=%s date=%s (%s)",
        id, platform, orbit, date, reason or "no reason given",
    )


def _exclusion_value_matches(a, b) -> bool:
    """
    Robust equality check for exclusion matching -- tries numeric
    comparison first, so int/float/numeric-string representations of
    the same value all compare correctly (e.g. 1517 == 1517.0 == "1517"
    -- same int/float representation issue confirmed for
    _apply_exclusions' dataframe matching, applied here too for
    consistency), falling back to exact equality for non-numeric values
    (e.g. a string-based orbit identifier) or when either side is None.
    """
    if a is None or b is None:
        return a == b
    try:
        return float(a) == float(b)
    except (TypeError, ValueError):
        return a == b


def remove_exclusion(
    prj: "Project", target_type: str, id,
    index: int = None, platform=None, orbit=None, date=None,
) -> list:
    """
    Remove one or more exclusion rules, either by position (index, from
    list_exclusions()) or by matching criteria -- the same
    platform/orbit/date fields used to add one via exclude_from_target.
    Matching by criteria is usually more convenient than looking up an
    index first: e.g. remove_exclusion(prj, "reservoirs", my_id,
    platform="S3B", orbit=1517) removes exactly the rule that excluded
    that mission+orbit combination (or platform+beam, for ICESat-2 --
    beam values ARE the "orbit" field once concatenated with other
    missions, see PRODUCT_TIMESERIES_KEYS -- there's no separate "beam"
    parameter needed).

    Specify either index, or at least one of platform/orbit/date, not
    both. Criteria matching removes EVERY exclusion rule whose given
    fields match (fields not specified are ignored, not required to be
    None on the stored rule).

    Returns the list of removed rule(s), for confirmation/logging.
    """
    config = _load_run_config(prj, id)
    exclusions = config["exclusions"]
    criteria_given = platform is not None or orbit is not None or date is not None

    if index is not None and criteria_given:
        raise ValueError(
            "Specify either index OR platform/orbit/date criteria, not both."
        )

    if index is not None:
        if index < 0 or index >= len(exclusions):
            raise IndexError(
                f"No exclusion at index {index} for {id}; there are "
                f"{len(exclusions)}. See list_exclusions()."
            )
        removed = [exclusions.pop(index)]
    elif criteria_given:
        date_str = str(date) if date is not None else None
        to_remove = [
            rule for rule in exclusions
            if (platform is None or _exclusion_value_matches(rule.get("platform"), platform))
            and (orbit is None or _exclusion_value_matches(rule.get("orbit"), orbit))
            and (date is None or rule.get("date") == date_str)
        ]
        if not to_remove:
            raise ValueError(
                f"No exclusion found matching platform={platform!r} "
                f"orbit={orbit!r} date={date!r} for {id}. See list_exclusions()."
            )
        for rule in to_remove:
            exclusions.remove(rule)
        removed = to_remove
    else:
        raise ValueError(
            "Specify index, or at least one of platform/orbit/date, to "
            "identify which exclusion(s) to remove."
        )

    _save_run_config(prj, id, config)
    _invalidate_spatial_correction_cache(prj, id)
    _invalidate_reach_slope_correction_cache(prj, id)
    logger.info("Removed %d exclusion(s) for %s: %s", len(removed), id, removed)
    return removed


def set_merging_option(prj: "Project", target_type: str, id, **kwargs) -> None:
    """
    Override one or more merging_options for just this one target,
    persisted the same way as exclusions (highest-priority layer: these
    override prj.reservoirs/rivers.merging_options, which override the
    DEFAULT_*_MERGING_OPTIONS defaults).

    Example: set_merging_option(prj, "reservoirs", my_id, svr_radial_err=0.5)
    """
    config = _load_run_config(prj, id)
    config["merging_option_overrides"].update(kwargs)
    _save_run_config(prj, id, config)
    if "use_spatial_correction" in kwargs or "spatial_correction_dense_source" in kwargs:
        _invalidate_spatial_correction_cache(prj, id)
    if "use_reach_slope_correction" in kwargs:
        _invalidate_reach_slope_correction_cache(prj, id)
    logger.info("Updated merging options for %s: %s", id, kwargs)


def _apply_exclusions(df: pd.DataFrame, exclusions: list) -> pd.DataFrame:
    """
    Filter out rows matching any exclusion rule. Within one rule, every
    specified field (platform/orbit/date) must match for a row to be
    excluded by it; a row is dropped if it matches ANY rule.
    """
    if not exclusions:
        return df

    keep_mask = pd.Series(True, index=df.index)
    for rule in exclusions:
        rule_mask = pd.Series(True, index=df.index)
        if rule.get("platform") is not None:
            rule_mask &= df["platform"] == rule["platform"]
        if rule.get("orbit") is not None:
            if "orbit" in df.columns:
                # Compare numerically when possible, not as strings --
                # a real orbit column commonly gets upcast to float64 by
                # pandas the moment ANY value in it is missing (very
                # common in real satellite data), so a genuine orbit
                # value of 1517 reads as 1517.0 in the dataframe while a
                # YAML-loaded exclusion rule reads it as the plain int
                # 1517. Comparing as strings ("1517.0" vs "1517") then
                # silently matches nothing -- confirmed as a real,
                # reproducible bug, not a hypothetical one. Falls back
                # to string comparison only if the orbit value genuinely
                # isn't numeric (e.g. a string-based identifier).
                try:
                    target_orbit = float(rule["orbit"])
                    rule_mask &= (
                        pd.to_numeric(df["orbit"], errors="coerce") == target_orbit
                    )
                except (TypeError, ValueError):
                    rule_mask &= df["orbit"].astype(str) == str(rule["orbit"])
            else:
                rule_mask &= False
        if rule.get("date") is not None:
            rule_mask &= df["date"].dt.floor("D").astype(str) == str(rule["date"])
        keep_mask &= ~rule_mask

    return df.loc[keep_mask].reset_index(drop=True)


DEFAULT_RIVER_MERGING_OPTIONS = {
    # Mostly a starting point copied from the reservoir defaults and NOT
    # independently validated against real river data the way the
    # reservoir defaults were validated this session -- river dynamics
    # differ genuinely (e.g. a real, expected along-reach gradient), so
    # do not assume the rest of these are correct without checking.
    # svr_radial_err/gamma below ARE an explicit exception (set directly,
    # not copied): gamma x100 (~15-day lengthscale, vs DAHITI's ~151-day
    # lake value) and err=1.0, matching the same oversmoothing reasoning
    # as the reservoir defaults, just with a shorter lengthscale given
    # rivers can change faster still.
    "window_km": 1.5,
    "svr_linear_err": 0.1,
    "svr_linear_epsilon": 0.1,
    "svr_radial_err": 1.0,
    "svr_radial_rbf_c": 10000,
    "svr_radial_gamma": 0.0000438 * 100,
    "svr_radial_epsilon": 0.1,
    "bias_time_bin": "20D",
    "bias_min_overlap": 1,
    # Same reasoning/evidence as the reservoir default (see
    # DEFAULT_RESERVOIR_MERGING_OPTIONS) for switching from "platform" to
    # "platform_orbit" -- but this is carried over, not independently
    # verified against real river data. A single river target (node/reach)
    # is a much smaller footprint than a reservoir, so it's genuinely
    # unclear whether the same within-platform configuration split
    # (e.g. S3B's two distinct crossing geometries) would even occur at
    # this scale -- check real per-target bias diagnostics once river
    # data exists before trusting this.
    "bias_group_by": "platform_orbit",
    "bias_centroid_warn_km": 5.0,
    # Off by default, same reasoning as reservoirs. NOTE: an earlier
    # version of this comment claimed a river target's footprint is
    # "much smaller than a reservoir" -- that's wrong for reaches
    # specifically (confirmed ~10km typical length, comparable to or
    # larger than many reservoirs), so distance_penalty/spatial
    # correction may matter just as much for reaches as for reservoirs.
    # It remains true that these tools address spread WITHIN one
    # target's own crossing footprint, never the natural gradient
    # BETWEEN different targets, which should never be "corrected away".
    "distance_penalty_scale_per_km": None,
    "use_spatial_correction": False,
    "spatial_correction_dense_source": "icesat2",
    # Off by default. ONLY meaningful when
    # prj.rivers.target_id_col == "reach_id" -- reference-corrects
    # non-SWOT crossings (ICESat-2/Sentinel-3/6) to what they'd read at
    # the reach's geometric midpoint, using SWOT's own directly-measured
    # "slope" field (see _fit_reach_slope_correction/
    # _apply_reach_slope_correction). Requires "slope" to be present in
    # mission_options["swot"]["hydrocron_fields"]["reaches"]. The
    # midpoint-referenced assumption for SWOT's own reach WSE is an
    # evidence-based inference from the RiverSP processing chain, not a
    # fact directly confirmed in SWOT's documentation -- and the sign of
    # the correction has not been empirically verified against real
    # data in this session. Validate both before trusting this in
    # production.
    "use_reach_slope_correction": False,
}


def _merge_timeseries(prj: "Project", target_type: str) -> None:
    """Merge multi-mission timeseries into combined datasets, for either
    reservoirs or river targets (nodes/reaches)."""
    target_ids = _get_target_ids(prj, target_type)
    ids_with_cleaned = [
        id
        for id in target_ids
        if os.path.exists(
            os.path.join(prj.dirs["output"], f"{id}", "cleaned_observations")
        )
    ]
    if not ids_with_cleaned:
        logger.warning(
            "No cleaned observations found for any %s; skipping timeseries merging.",
            target_type,
        )
        return

    default_options = (
        DEFAULT_RESERVOIR_MERGING_OPTIONS
        if target_type == "reservoirs"
        else DEFAULT_RIVER_MERGING_OPTIONS
    )
    # prj.reservoirs.merging_options / prj.rivers.merging_options are the
    # intended per-target-type override locations (set from the
    # respective YAML config sections); fall back to the older shared
    # prj.merging_options for backward compatibility if the per-type one
    # isn't set yet.
    target_owner = prj.reservoirs if target_type == "reservoirs" else prj.rivers
    overrides = getattr(target_owner, "merging_options", None)
    if overrides is None:
        overrides = getattr(prj, "merging_options", None)

    for id in tqdm(ids_with_cleaned, desc=f"Merging product timeseries ({target_type})"):
        ts_list = []
        for product in prj.to_process:
            df = _load_product_timeseries(
                os.path.join(prj.dirs["output"], f"{id}", "cleaned_observations"),
                ".csv",
                [product],
                pd.read_csv,
            )
            if df is not None:
                df["date"] = pd.to_datetime(
                    df.date, format="mixed", utc=True
                ).dt.tz_convert(None)
                df = df.sort_values(by="date")
                ts_list.append(
                    timeseries.Timeseries(
                        df, date_key="date", height_key="height",
                        **PRODUCT_TIMESERIES_KEYS.get(product, {}),
                    )
                )

        if len(ts_list) > 0:
            ts = timeseries.concat(ts_list)

            data_dir = os.path.join(prj.dirs["output"], f"{id}")
            general.ifnotmakedirs(data_dir)

            run_config = _load_run_config(prj, id)

            merging_options = dict(default_options)
            merging_options.update(overrides or {})
            # Per-target overrides (from notebook calls to
            # set_merging_option, or hand-edited in run_config.yaml) take
            # priority over project-wide settings for just this target.
            merging_options.update(run_config.get("merging_option_overrides", {}))
            distance_penalty_scale = merging_options.pop(
                "distance_penalty_scale_per_km", None
            )
            use_spatial_correction = merging_options.pop(
                "use_spatial_correction", False
            )
            spatial_correction_dense_source = merging_options.pop(
                "spatial_correction_dense_source", "icesat2"
            )
            # Off by default. Only meaningful for reach-mode river
            # projects (a node is a single ~200m point, not a ~10km
            # segment with its own along-reach slope) -- see
            # _fit_reach_slope_correction for the full reasoning and the
            # sign-convention caveat that should be checked against real
            # data before trusting this in production.
            use_reach_slope_correction = merging_options.pop(
                "use_reach_slope_correction", False
            )

            # Apply exclusions BEFORE exporting all_cleaned_timeseries.csv
            # (not just before merge processing) -- this file is meant to
            # reflect what's actually being worked with, and writing it
            # before exclusions were applied meant it always showed
            # excluded data regardless of how many times you re-ran,
            # which looked exactly like a stale file from an old run but
            # was actually happening on every single run. The full,
            # pre-exclusion record is still available per-mission in
            # cleaned_observations/{product}.csv (written earlier, in
            # _clean_timeseries, before any exclusion is applied) -- so
            # nothing is lost by making this file reflect exclusions.
            exclusions = run_config.get("exclusions", [])
            if exclusions:
                before = len(ts.df)
                ts.df = _apply_exclusions(ts.df, exclusions)
                logger.info(
                    "%s: %d exclusion rule(s) applied, %d/%d observations kept.",
                    id, len(exclusions), len(ts.df), before,
                )

            if (
                use_reach_slope_correction
                and target_type == "rivers"
                and getattr(prj.rivers, "target_id_col", None) == "reach_id"
            ):
                slope_model = _fit_reach_slope_correction(prj, id)
                if slope_model is not None:
                    ts.df = _apply_reach_slope_correction(ts.df, prj, id, slope_model)
                    logger.info(
                        "%s: applied reach slope correction (median_slope=%.6g "
                        "from %d SWOT observations).",
                        id, slope_model["median_slope"], slope_model["n_observations"],
                    )
                else:
                    logger.info(
                        "%s: use_reach_slope_correction is enabled but no "
                        "usable SWOT slope data was found; skipping "
                        "correction for this target.", id,
                    )

            ts.export_csv(os.path.join(data_dir, "all_cleaned_timeseries.csv"))

            ref_lat, ref_lon = _target_centroid(prj, target_type, id)

            spatial_correction_model = None
            if use_spatial_correction:
                spatial_correction_model = _get_or_fit_spatial_correction_model(
                    prj, target_type, id,
                    dense_source_platform=spatial_correction_dense_source,
                )

            ts = ts.merge(
                save_progress=True,
                dir=os.path.join(data_dir, "merged_progress"),
                ref_lat=ref_lat,
                ref_lon=ref_lon,
                distance_penalty_scale_per_km=distance_penalty_scale,
                spatial_correction_model=spatial_correction_model,
                **merging_options,
            )
            ts.export_csv(os.path.join(data_dir, "merged_timeseries.csv"))


def _merge_reservoirs_timeseries(prj: "Project") -> None:
    """Merge multi-mission timeseries into combined datasets, for reservoirs."""
    _merge_timeseries(prj, "reservoirs")


def _merge_rivers_timeseries(prj: "Project") -> None:
    """Merge multi-mission timeseries into combined datasets, for river targets."""
    _merge_timeseries(prj, "rivers")


# ============================================================================
# RESERVOIRS: Summaries & Visualization
# ============================================================================


def generate_reservoirs_summaries(
    prj: "Project", show: bool = False, save: bool = True
) -> None:
    """Generate per-reservoir plotting summaries.

    Parameters
    ----------
    prj : Project
        Project instance with reservoirs configuration
    show : bool
        Whether to display plots interactively
    save : bool
        Whether to save plots to disk
    """
    if not hasattr(prj, "reservoirs"):
        return

    for reservoir_id in prj.reservoirs.download_gdf[prj.reservoirs.id_key]:
        plotting.plot_crossings(
            gdf=prj.reservoirs.gdf,
            id_key=prj.reservoirs.id_key,
            reservoir_id=reservoir_id,
            output_dir=prj.dirs["output"],
            reservoir_type=prj.reservoirs.type,
            show=show,
            save=save,
        )

        plotting.plot_cleaning(
            reservoir_id=reservoir_id,
            output_dir=prj.dirs["output"],
            get_unfiltered_fn=lambda id, products: _load_product_timeseries(
                os.path.join(prj.dirs["output"], f"{id}", "raw_observations"),
                ".gpkg",
                products,
                lambda path: gpd.read_file(path).drop(columns=["geometry"]),
            ),
            get_cleaned_fn=lambda id, products: _load_and_parse_cleaned_timeseries(
                prj, id, products
            ),
            get_merged_fn=lambda id: _load_merged_timeseries(prj, id),
            reservoir_type=prj.reservoirs.type,
            show=show,
            save=save,
            products=getattr(prj, "to_process", None),
        )

        plotting.plot_merging(
            reservoir_id=reservoir_id,
            output_dir=prj.dirs["output"],
            reservoir_type=prj.reservoirs.type,
            show=show,
            save=save,
        )


def _load_and_parse_cleaned_timeseries(prj, id, products):
    """Load cleaned observations and parse dates."""
    df = _load_product_timeseries(
        os.path.join(prj.dirs["output"], f"{id}", "cleaned_observations"),
        ".csv",
        products,
        pd.read_csv,
    )
    if df is not None:
        df["date"] = pd.to_datetime(df.date, format="mixed", utc=True).dt.tz_convert(
            None
        )
        df = df.sort_values(by="date")
    return df


def _load_merged_timeseries(prj, id):
    """Load merged timeseries if it exists."""
    data_path = os.path.join(prj.dirs["output"], f"{id}", "merged_timeseries.csv")

    if os.path.exists(data_path):
        df = pd.read_csv(data_path)
        df["date"] = pd.to_datetime(df.date)
        df = df.sort_values(by="date")
        return df
    else:
        logger.warning(
            "%s does not exist, be sure to merge product timeseries first!",
            data_path,
        )
        return None


# ============================================================================
# RIVERS: Summaries & Visualization
# ============================================================================


def _project_num_months(prj: "Project") -> int:
    """
    Approximate number of months spanned by the project's configured
    date range -- used as a minimum-observation-count threshold for
    plotting (see _has_enough_observations_to_plot). Falls back to 1 if
    the project-level dates aren't resolvable for some reason.
    """
    project_cfg = prj.config.get("project", {})
    start = project_cfg.get("startdate")
    end = project_cfg.get("enddate")
    if not start or not end:
        return 1

    start_date = datetime.date(*start) if isinstance(start, list) else start
    end_date = datetime.date(*end) if isinstance(end, list) else end
    months = (
        (end_date.year - start_date.year) * 12
        + (end_date.month - start_date.month)
        + 1
    )
    return max(months, 1)


def _has_enough_observations_to_plot(prj: "Project", target_id, min_months: int) -> bool:
    """
    Whether a target has enough merged observations to be worth
    plotting -- more than min_months (the project's date range in
    months) or more than 2, whichever is larger. A reach/reservoir with
    only 1-2 points produces a plot that adds noise without telling you
    anything.
    """
    df = _load_merged_timeseries(prj, target_id)
    if df is None:
        return False
    threshold = max(min_months, 2)
    return len(df) > threshold


def generate_rivers_summaries(
    prj: "Project", show: bool = False, save: bool = True
) -> None:
    """Generate per-river plotting summaries.

    Parameters
    ----------
    prj : Project
        Project instance with rivers configuration
    show : bool
        Whether to display plots interactively
    save : bool
        Whether to save plots to disk
    """
    if not hasattr(prj, "rivers"):
        return

    waterbody_groups = _group_river_targets_by_waterbody(prj)
    min_months = _project_num_months(prj)

    for wb_id, target_ids in waterbody_groups.items():
        # Only plot targets with enough observations to be worth looking
        # at -- applies to all three plot types (map, time series, merge
        # progress) so a target excluded from one isn't confusingly still
        # shown in another.
        plottable_ids = [
            t for t in target_ids
            if _has_enough_observations_to_plot(prj, t, min_months)
        ]
        if not plottable_ids:
            logger.info(
                "Skipping plots for waterbody %s -- no targets with more "
                "than %d observations.", wb_id, max(min_months, 2),
            )
            continue

        # Compute the actual extraction corridor (same buffer resolution
        # used for real extraction, see _river_target_corridor) so the
        # shaded area shown is exactly what extraction uses, not an
        # approximation -- lets you visually judge whether width-based
        # buffering produced a reasonable corridor for this waterbody
        # (e.g. a lake-flagged reach whose SWORD width reflects a much
        # wider lake extent) without needing external knowledge of the
        # real river geometry.
        corridor_gdf = _river_target_corridor(
            prj, plottable_ids,
            buffer_meters=getattr(prj.rivers, "extraction_buffer_meters", None),
            width_buffer_factor=getattr(prj.rivers, "width_buffer_factor", 1.05),
        )
        corridor_geometry = corridor_gdf.geometry.iloc[0] if corridor_gdf is not None else None

        plotting.plot_river_crossings(
            prj, wb_id, plottable_ids, prj.dirs["output"], show=show, save=save,
            corridor_geometry=corridor_geometry,
        )

        plotting.plot_river_data(
            prj, wb_id, plottable_ids, prj.dirs["output"],
            get_merged_fn=lambda id: _load_merged_timeseries(prj, id),
            show=show, save=save,
        )

        for target_id in plottable_ids:
            plotting.plot_merging(
                reservoir_id=target_id,
                output_dir=prj.dirs["output"],
                reservoir_type="river",
                show=show,
                save=save,
            )


# ============================================================================
# MIKEIO
# ============================================================================


def _export_cleaned_to_dfs0(prj: "Project") -> None:
    """Export cleaned timeseries observations to dfs0 format.

    Parameters
    ----------
    prj : Project
        Project instance with reservoirs configuration

    Notes
    -----
    Exports height data from cleaned CSV observations to dfs0 format
    for each product in each reservoir's cleaned_observations folder.
    """
    ids_with_cleaned = [
        id
        for id in prj.reservoirs.download_gdf[prj.reservoirs.id_key]
        if os.path.exists(
            os.path.join(prj.dirs["output"], f"{id}", "cleaned_observations")
        )
    ]

    if not ids_with_cleaned:
        logger.warning(
            "No cleaned observations found for any reservoir; skipping dfs0 export."
        )
        return

    for id in ids_with_cleaned:
        cleaned_dir = os.path.join(prj.dirs["output"], f"{id}", "cleaned_observations")

        for product in tqdm(
            prj.to_process, desc="Exporting cleaned observations to dfs0"
        ):
            csv_path = os.path.join(cleaned_dir, f"{product}.csv")

            if not os.path.exists(csv_path):
                continue

            try:
                # Load and prepare data
                df = pd.read_csv(csv_path)

                # Parse and set datetime index
                df["date"] = pd.to_datetime(
                    df.date, format="mixed", utc=True
                ).dt.tz_convert(None)
                df = df.set_index("date")
                df = df.sort_index()

                # Remove rows with NaN heights
                df = df.dropna(subset=["height"])

                if len(df) == 0:
                    logger.warning(
                        "No valid height data for %s in %s after cleaning",
                        product,
                        id,
                    )
                    continue

                # Create Dataset and export to dfs0
                items = {"height": mikeio.ItemInfo(mikeio.EUMType.Water_Level)}
                ds = mikeio.from_pandas(df[["height"]], items=items)

                # Write dfs0 file
                dfs0_path = os.path.join(cleaned_dir, f"{product}.dfs0")
                ds.to_dfs(dfs0_path)

                logger.debug("Exported dfs0: %s", dfs0_path)

            except Exception as exc:
                logger.error(
                    "Failed to export %s to dfs0 for %s: %s",
                    product,
                    id,
                    exc,
                )