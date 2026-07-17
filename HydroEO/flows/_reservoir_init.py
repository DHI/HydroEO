"""
Reservoirs: PLD (Prior Lake Database) initialization.
"""

import logging
import os

import geopandas as gpd
import pandas as pd

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
    """Spatial join reservoirs with PLD to assign prior_lake_id.

    Prefers the PLD lake with the LARGEST overlapping area over merely the
    nearest one.
    """
    pld = gpd.read_file(prj.dirs["pld"])
    pld = pld.rename(columns={"lake_id": "prior_lake_id", "res_id": "prior_res_id"})
    if "prior_res_id" not in pld.columns:
        pld["prior_res_id"] = None

    id_key = prj.reservoirs.id_key
    reservoirs_local = prj.reservoirs.gdf.to_crs(prj.local_crs)
    pld_local = pld.to_crs(prj.local_crs)
    reservoir_areas = reservoirs_local.set_index(id_key).geometry.area
    min_overlap_pct = prj.mission_options.get("swot", {}).get(
        "pld_match_min_overlap_pct", 10.0
    )

    # Match by largest overlapping area, computed via a true geometric
    # intersection (not just "does it intersect" or "how far apart are the
    # boundaries").
    overlap = gpd.overlay(
        reservoirs_local[[id_key, "geometry"]],
        pld_local[["prior_lake_id", "prior_res_id", "geometry"]],
        how="intersection",
        keep_geom_type=False,
    )
    if len(overlap) > 0:
        overlap["dist_to_pld"] = 0.0
        overlap["pld_match_method"] = "overlap"
        overlap["_overlap_area"] = overlap.geometry.area
        matches = (
            overlap.sort_values("_overlap_area", ascending=False)
            .drop_duplicates(subset=id_key, keep="first")
        )

        # Overlap area as a percentage of the RESERVOIR's own area.
        matches["_reservoir_area"] = matches[id_key].map(reservoir_areas)
        matches["_overlap_pct"] = (
            matches["_overlap_area"] / matches["_reservoir_area"] * 100
        )
        low_overlap = matches.loc[matches["_overlap_pct"] < min_overlap_pct]
        if len(low_overlap) > 0:
            logger.warning(
                "%d reservoir(s) matched a PLD lake covering less than "
                "%.0f%% of the reservoir's own area: %s -- still using "
                "this match (it's the best candidate available), but "
                "verify it's correct rather than a small, unrelated lake "
                "that happens to clip the reservoir's edge.",
                len(low_overlap),
                min_overlap_pct,
                ", ".join(
                    f"{row[id_key]} ({row['_overlap_pct']:.1f}%)"
                    for _, row in low_overlap.iterrows()
                ),
            )

        matches = matches.drop(
            columns=["_overlap_area", "_overlap_pct", "_reservoir_area", "geometry"]
        )
    else:
        matches = pd.DataFrame(
            columns=[id_key, "prior_lake_id", "prior_res_id", "dist_to_pld", "pld_match_method"]
        )

    joined_gdf = prj.reservoirs.gdf.merge(matches, on=id_key, how="left")

    # Real PLD files may store lake_id as text rather than integer (e.g.
    # the '_light' product) -- coerce here, once, rather than leaving every
    # downstream numeric comparison (_flag_missing_priors's "> 0" / "< 0",
    # and the equivalent check during SWOT extraction) vulnerable to a
    # TypeError comparing str > int.
    joined_gdf["prior_lake_id"] = pd.to_numeric(
        joined_gdf["prior_lake_id"], errors="coerce"
    )
    joined_gdf.loc[joined_gdf.prior_lake_id.isnull(), "prior_lake_id"] = -9999
    joined_gdf.loc[joined_gdf["prior_lake_id"] == -9999, "pld_match_method"] = "unmatched"

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
            "%d reservoir(s) have no geometrically overlapping Prior Lake "
            "Database (PLD) lake: %s (see aux/PLD/missing_in_pld.gpkg). SWOT "
            "Lake SP cannot report observations for these regardless of "
            "download success. "
            "Check the geometry in missing_in_pld.gpkg against your source "
            "reservoir polygon before proceeding, or exclude them from this "
            "project if they're SWOT-only.",
            len(missing),
            ", ".join(str(v) for v in missing[prj.reservoirs.id_key]),
        )
