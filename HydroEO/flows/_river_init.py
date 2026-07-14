"""Rivers: SWORD database initialization/subsetting.

initialize_rivers -> _prepare_rivers_from_sword -> _ensure_sword_database
are tested together via patch.object(flows, "_name") in
tests/unit/test_flows.py -- keep them in this one module.
"""

import logging
import os

import geopandas as gpd

from HydroEO.utils import general

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from HydroEO.project import Project

logger = logging.getLogger(__name__)


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


