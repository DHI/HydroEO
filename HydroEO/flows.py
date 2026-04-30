"""Standalone flow functions for HydroEO pipelines.

These functions implement the core download, processing, and visualization logic
previously embedded in Reservoirs and Rivers classes. They operate on Project state
and external data, with no direct method dependencies.
"""

import logging
import os
import datetime
from io import StringIO
from typing import TYPE_CHECKING

import geopandas as gpd
import pandas as pd
from tqdm import tqdm

from HydroEO.satellites import swot, icesat2, sentinel
from HydroEO.utils import general, timeseries
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
    file_name = os.path.basename(pld_path)
    bounds = list(prj.reservoirs.gdf.unary_union.bounds)
    hydroweb.download_PLD(download_dir=download_dir, file_name=file_name, bounds=bounds)


def _assign_pld_id(prj: "Project") -> None:
    """Spatial join reservoirs with PLD to assign prior_lake_id."""
    pld = gpd.read_file(prj.dirs["pld"])

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

    joined_gdf = joined_gdf.rename(
        columns={"lake_id": "prior_lake_id", "res_id": "prior_res_id"}
    )

    joined_gdf.loc[joined_gdf.prior_lake_id.isnull(), "prior_lake_id"] = -9999

    prj.reservoirs.gdf = joined_gdf


def _flag_missing_priors(prj: "Project") -> None:
    """Export shapefiles of reservoirs present/missing in PLD."""
    gdf = prj.reservoirs.gdf
    present = gdf.loc[gdf.prior_lake_id > 0].reset_index(drop=True)
    missing = gdf.loc[gdf.prior_lake_id < 0].reset_index(drop=True)

    shp_field_map = {
        "index_right": "idx_right",
        "prior_lake_id": "prior_lake",
        "prior_res_id": "prior_res",
        "dist_to_pld": "dist_pld",
    }
    present_out = present.rename(columns=shp_field_map)
    missing_out = missing.rename(columns=shp_field_map)

    present_out.to_file(os.path.join(prj.dirs["output"], "present_in_pld.shp"))
    missing_out.to_file(os.path.join(prj.dirs["output"], "missing_in_pld.shp"))

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
        "Initialized river %s ids: %s",
        id_label,
        ", ".join(str(target_id) for target_id in prj.rivers.target_ids),
    )


def _prepare_rivers_from_sword(prj: "Project") -> None:
    """Prepare SWORD target features by spatial intersection with AOI."""
    sword_dir = _ensure_sword_database(prj)
    gpkg_name = f"{prj.rivers.continent_key}_sword_{prj.rivers.feature_type}_v17b.gpkg"
    gpkg_path = os.path.join(sword_dir, gpkg_name)

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

    source_id_col = "node_id" if prj.rivers.feature_type == "nodes" else "reach_id"
    if source_id_col not in subset.columns:
        raise KeyError(
            f"Expected SWORD column '{source_id_col}' missing from {gpkg_name}"
        )

    prj.rivers.target_features = subset
    prj.rivers.target_id_col = source_id_col
    prj.rivers.target_ids = [int(value) for value in subset[source_id_col]]


def _ensure_sword_database(prj: "Project") -> str:
    """Ensure SWORD database is downloaded and extracted."""
    from urllib import request as url_request
    import zipfile

    SWORD_V17B_ZIP_URL = (
        "https://zenodo.org/records/15299138/files/SWORD_v17b_gpkg.zip?download=1"
    )

    sword_dir = os.path.join(prj.dirs["main"], "SWORD_v17b_gpkg", "gpkg")
    prj.dirs["sword"] = sword_dir

    if os.path.isdir(sword_dir):
        return sword_dir

    logger.warning(
        "SWORD_v17b database not found in %s. Downloading and extracting it now.",
        sword_dir,
    )
    general.ifnotmakedirs(sword_dir)

    zip_path = os.path.join(prj.dirs["main"], "SWORD_v17b_gpkg.zip")
    url_request.urlretrieve(SWORD_V17B_ZIP_URL, zip_path)

    with zipfile.ZipFile(zip_path, "r") as zip_ref:
        zip_ref.extractall(sword_dir)

    os.remove(zip_path)
    return sword_dir


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
    download_dir = os.path.join(prj.dirs["swot"], "reservoirs")
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

        geom = prj.reservoirs.download_gdf.loc[i, "geometry"]
        if hasattr(geom, "geoms"):
            geom = geom.geoms[0]
        coords = list(geom.exterior.coords)

        download_dir = os.path.join(prj.dirs["icesat2"], "reservoirs", rf"{id}")
        general.ifnotmakedirs(download_dir)

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
                download_directory=download_dir,
                atl13_options=prj.mission_options.get("icesat2", {}).get("atl13", {}),
                atl13_fields=prj.mission_options.get("icesat2", {}).get("atl13_fields")
                or None,
            )
        except Exception as exc:
            logger.warning("ICESat-2 download skipped for %s: %s", id, exc)


def _download_reservoirs_sentinel(prj: "Project", mission: str) -> None:
    """Download Sentinel-3 or Sentinel-6 data for reservoirs."""
    product = "S3" if mission == "sentinel3" else "S6"
    dir_key = mission

    session_token = None
    session_start_time = None

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

        download_dir = os.path.join(prj.dirs[dir_key], "reservoirs", rf"{id}")
        general.ifnotmakedirs(download_dir)

        startdate = prj.startdates[mission]
        enddate = prj.enddates[mission]

        if isinstance(startdate, list):
            startdate = datetime.date(*startdate)
        if isinstance(enddate, list):
            enddate = datetime.date(*enddate)

        logger.info(
            "Searching for Sentinel-%s for aoi from %s to %s",
            product,
            startdate,
            enddate,
        )
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


# ============================================================================
# RIVERS: Download
# ============================================================================


def download_rivers(prj: "Project") -> None:
    """Download SWOT Hydrocron data for rivers.

    Parameters
    ----------
    prj : Project
        Project instance with rivers configuration
    """
    if not hasattr(prj, "rivers"):
        return

    if "swot" not in prj.to_download:
        return

    startdate = prj.startdates.get("swot")
    enddate = prj.enddates.get("swot")

    if isinstance(startdate, list):
        startdate = datetime.date(*startdate)
    if isinstance(enddate, list):
        enddate = datetime.date(*enddate)

    _download_swot_hydrocron_timeseries(prj, startdate, enddate)


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

    summary = {
        "requested": len(waterbody_groups),
        "successful": 0,
        "failed": 0,
        "empty_after_filter": 0,
    }

    for wb_id, target_ids in waterbody_groups.items():
        output_path = os.path.join(
            prj.dirs["swot"], "rivers", str(wb_id), "timeseries.csv"
        )
        general.ifnotmakedirs(os.path.dirname(output_path))

        wb_startdate = startdate
        latest_obs = _get_latest_hydrocron_obs_date(output_path)
        if latest_obs is not None:
            wb_startdate = latest_obs

        frames = []
        for target_id in target_ids:
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
                logger.warning(
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
                logger.warning(
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
                logger.warning(
                    "Failed to parse CSV for %s %s: %s",
                    prj.rivers.target_id_col,
                    target_id,
                    exc,
                )
                summary["failed"] += 1
                continue

            if df.empty or quality_column not in df.columns:
                if df.empty:
                    logger.info(
                        "Hydrocron returned no data for %s %s",
                        prj.rivers.target_id_col,
                        target_id,
                    )
                else:
                    logger.warning(
                        "Quality column %s not in Hydrocron response for %s %s",
                        quality_column,
                        prj.rivers.target_id_col,
                        target_id,
                    )
                continue

            df = df[df[quality_column] <= max_q]
            if df.empty:
                logger.info(
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

    logger.info(
        "Hydrocron download complete: %s requested, %s successful, %s failed, %s empty after filtering",
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

    # Extract raw observations from downloaded files
    _extract_reservoirs_timeseries(prj)

    # Clean observations with filters
    _clean_reservoirs_timeseries(prj)

    # Merge multi-mission timeseries
    _merge_reservoirs_timeseries(prj)


def _extract_reservoirs_timeseries(prj: "Project") -> None:
    """Extract timeseries observations from raw downloaded files."""
    if "icesat2" in prj.to_process:
        _extract_icesat2_observations(prj)

    if "sentinel3" in prj.to_process:
        _extract_sentinel_observations(prj, "sentinel3", "S3")

    if "sentinel6" in prj.to_process:
        _extract_sentinel_observations(prj, "sentinel6", "S6")

    if "swot" in prj.to_process:
        _extract_swot_observations(prj)


def _extract_icesat2_observations(prj: "Project") -> None:
    """Extract ICESat-2 ATL13 observations for each reservoir."""
    available_ids = [
        id
        for id in prj.reservoirs.download_gdf[prj.reservoirs.id_key]
        if os.path.exists(
            os.path.join(prj.dirs["icesat2"], "reservoirs", f"{id}", "atl13.parquet")
        )
    ]
    if not available_ids:
        logger.warning("No ICESat-2 downloads found; skipping timeseries extraction.")
        return

    empty_ids = []
    for id in tqdm(available_ids, desc="Extracting ICESat-2 ATL13 product"):
        sub_gdf = prj.reservoirs.download_gdf.loc[
            prj.reservoirs.download_gdf[prj.reservoirs.id_key] == id
        ]
        download_dir = os.path.join(prj.dirs["icesat2"], "reservoirs", f"{id}")
        dst_dir = os.path.join(prj.dirs["output"], f"{id}", "raw_observations")
        general.ifnotmakedirs(dst_dir)
        dst_path = os.path.join(dst_dir, "icesat2.shp")

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
    prj: "Project", mission_key: str, product: str
) -> None:
    """Extract Sentinel-3 or Sentinel-6 observations for each reservoir."""
    available_ids = [
        id
        for id in prj.reservoirs.download_gdf[prj.reservoirs.id_key]
        if os.path.exists(os.path.join(prj.dirs[mission_key], "reservoirs", f"{id}"))
    ]
    if not available_ids:
        logger.warning(
            "No %s downloads found; skipping timeseries extraction.", mission_key
        )
        return

    empty_ids = []
    for id in tqdm(available_ids, desc=f"Extracting Sentinel-{product} product"):
        sub_gdf = prj.reservoirs.download_gdf.loc[
            prj.reservoirs.download_gdf[prj.reservoirs.id_key] == id
        ]
        download_dir = os.path.join(prj.dirs[mission_key], "reservoirs", f"{id}")
        dst_dir = os.path.join(prj.dirs["output"], f"{id}", "raw_observations")
        general.ifnotmakedirs(dst_dir)
        dst_path = os.path.join(dst_dir, f"{mission_key}.shp")

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


def _extract_swot_observations(prj: "Project") -> None:
    """Extract SWOT Lake SP observations for all reservoirs."""
    download_dir = os.path.join(prj.dirs["swot"], "reservoirs")
    if not os.path.exists(download_dir):
        logger.warning("No SWOT downloads found; skipping timeseries extraction.")
        return

    empty_ids = swot.extract_observations(
        src_dir=download_dir,
        dst_dir=prj.dirs["output"],
        dst_file_name="swot.shp",
        features=prj.reservoirs.download_gdf,
        id_key=prj.reservoirs.id_key,
        exclude_obs_id_values=prj.mission_options.get("swot", {}).get(
            "exclude_obs_id_values", ["no_data"]
        ),
    )
    if empty_ids:
        logger.warning(
            "SWOT timeseries empty for: %s (no observations matched the prior lake ID or all were excluded)",
            ", ".join(str(i) for i in empty_ids),
        )


def _clean_reservoirs_timeseries(prj: "Project") -> None:
    """Apply quality filters to extracted timeseries."""
    ids_with_raw = [
        id
        for id in prj.reservoirs.download_gdf[prj.reservoirs.id_key]
        if os.path.exists(os.path.join(prj.dirs["output"], f"{id}", "raw_observations"))
    ]
    if not ids_with_raw:
        logger.warning(
            "No raw observations found for any reservoir; skipping timeseries cleaning."
        )
        return

    for id in tqdm(ids_with_raw, desc="Cleaning product timeseries"):
        for product in prj.to_process:
            df = _load_product_timeseries(
                os.path.join(prj.dirs["output"], f"{id}", "raw_observations"),
                ".shp",
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

                ts = timeseries.Timeseries(df, date_key="date", height_key="height")

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


def _merge_reservoirs_timeseries(prj: "Project") -> None:
    """Merge multi-mission timeseries into combined datasets."""
    ids_with_cleaned = [
        id
        for id in prj.reservoirs.download_gdf[prj.reservoirs.id_key]
        if os.path.exists(
            os.path.join(prj.dirs["output"], f"{id}", "cleaned_observations")
        )
    ]
    if not ids_with_cleaned:
        logger.warning(
            "No cleaned observations found for any reservoir; skipping timeseries merging."
        )
        return

    for id in tqdm(ids_with_cleaned, desc="Merging product timeseries"):
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
                    timeseries.Timeseries(df, date_key="date", height_key="height")
                )

        if len(ts_list) > 0:
            ts = timeseries.concat(ts_list)

            data_dir = os.path.join(prj.dirs["output"], f"{id}")
            general.ifnotmakedirs(data_dir)
            ts.export_csv(os.path.join(data_dir, "all_cleaned_timeseries.csv"))

            ts = ts.merge(
                save_progress=True,
                dir=os.path.join(data_dir, "merged_progress"),
            )
            ts.export_csv(os.path.join(data_dir, "merged_timeseries.csv"))


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
        logger.info("Summarizing crossings")
        plotting.plot_crossings(
            gdf=prj.reservoirs.gdf,
            id_key=prj.reservoirs.id_key,
            reservoir_id=reservoir_id,
            output_dir=prj.dirs["output"],
            reservoir_type=prj.reservoirs.type,
            show=show,
            save=save,
        )

        logger.info("Summarizing cleaning results")
        plotting.plot_cleaning(
            reservoir_id=reservoir_id,
            output_dir=prj.dirs["output"],
            get_unfiltered_fn=lambda id, products: _load_product_timeseries(
                os.path.join(prj.dirs["output"], f"{id}", "raw_observations"),
                ".shp",
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

        logger.info("Summarizing merged results")
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
