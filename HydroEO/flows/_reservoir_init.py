"""
Reservoirs: PLD (Prior Lake Database) initialization.
"""

import logging
import os

import geopandas as gpd

from HydroEO.downloaders import hydroweb

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from HydroEO.project import Project

logger = logging.getLogger(__name__)


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
    bounds = list(prj.reservoirs.gdf.union_all().bounds)
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
        continent_codes=getattr(prj, "pld_continent_codes", None),
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
    joined_gdf = joined_gdf.drop(columns=["index_right", "index_left"], errors="ignore")

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
    if len(missing) > 0:
        logger.warning(
            "%d reservoir(s) not matched to any Prior Lake Database (PLD) lake "
            "within pld_match_max_distance_m: %s (see aux/PLD/missing_in_pld.gpkg). "
            "SWOT Lake SP cannot report observations for these. "
            "Consider increasing pld_match_max_distance_m if it's a "
            "near miss.",
            len(missing),
            ", ".join(str(v) for v in missing[prj.reservoirs.id_key]),
        )