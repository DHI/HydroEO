"""
Rivers: extraction + clean + merge orchestration.
"""

import logging
import os

import geopandas as gpd
import pandas as pd

from HydroEO.satellites import icesat2, sentinel
from HydroEO.utils import general
from ._river_common import (
    _group_river_targets_by_waterbody,
    _river_extraction_buffer_meters,
    _assign_points_to_river_targets,
)
from ._summaries import _river_target_corridor
from ._clean_engine import _clean_timeseries
from ._merge_engine import _merge_timeseries

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from HydroEO.project import Project

logger = logging.getLogger(__name__)


def _extract_rivers_timeseries(prj: "Project", overwrite: bool = False) -> None:
    """Extract timeseries observations from raw downloaded files, for rivers.

    overwrite:  Skip extraction if the target's output .gpkg already exists 
                (False, default) or re-extract anyway (True). 
                This is the same semantics as _extract_reservoirs_timeseries.
    Lightweight extraction always required for SWOT to match other missions 
    in post-processing.

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
    """Extract ICESat-2 ATL13 observations for each river target (virtual station).
    Uses buffered river corridor and assigns to nearest target. 
    Clean and merge treats river target (VS) exactly like a reservoir.

    Parameters
    ----------
    overwrite : bool, optional
        Skip extraction if the target's output .gpkg already exists
        (False, default) or re-extract anyway (True). This is the same
        semantics as _extract_reservoirs_timeseries.
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

    Option to use sigma0 as a quality filter (sigma0_min) 
    is available in the mission's YAML config section.

    Parameters
    ----------
    mission_key : str
        "sentinel3" or "sentinel6"
    product : str
        "S3" or "S6" (used for output file naming)
    overwrite : bool, optional
        Skip extraction if the target's output .gpkg already exists
        (False, default) or re-extract anyway (True). This is the same
        semantics as _extract_reservoirs_timeseries.
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

    Uses SWORD geometry as placeholder in lieu of "lat/lon".

    Parameters
    ----------
    overwrite : bool, optional
        Skip extraction if the target's output .gpkg already exists
        (False, default) or re-extract anyway (True). This is the same
        semantics as _extract_reservoirs_timeseries.
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


def create_rivers_timeseries(prj: "Project") -> None:
    """Extract, clean, and merge timeseries for river targets (nodes/reaches).

    Mirrors create_reservoirs_timeseries. 
    NOTE: Extraction to dsf0 missing.

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


def _clean_rivers_timeseries(prj: "Project") -> None:
    """Apply quality filters to extracted river timeseries."""
    _clean_timeseries(prj, "rivers")


def _merge_rivers_timeseries(prj: "Project") -> None:
    """Merge multi-mission timeseries into combined datasets, for river targets."""
    _merge_timeseries(prj, "rivers")

