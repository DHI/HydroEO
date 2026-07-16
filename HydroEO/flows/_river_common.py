"""
River geometry helpers shared by download, extraction, and summaries.
"""

import logging

import geopandas as gpd

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from HydroEO.project import Project

logger = logging.getLogger(__name__)


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


def _iter_geometry_pieces(geom):
    """
    Yield each individual polygon from a geometry: every part of a
    MultiPolygon, or the geometry itself for a plain Polygon.
    """
    if hasattr(geom, "geoms"):
        return list(geom.geoms)
    return [geom]


def _simplify_to_one_polygon(geom):
    """
    Collapse a MultiPolygon into a single encompassing polygon via
    convex hull.
    """
    if hasattr(geom, "geoms"):
        return geom.convex_hull
    return geom


def _river_extraction_buffer_meters(prj: "Project") -> float:
    """
    Resolve the extraction-corridor buffer distance: use
    prj.rivers.extraction_buffer_meters if set, else fall back to the
    SWORD-intersection prj.rivers.buffer_meters, else a conservative
    default (500 m).
    """
    explicit = getattr(prj.rivers, "extraction_buffer_meters", None)
    if explicit:
        return explicit
    if prj.rivers.buffer_meters:
        return prj.rivers.buffer_meters
    return 500.0


def _assign_points_to_river_targets(
    points, targets, target_id_col, max_distance_meters, local_crs
):
    """
    Assign each altimetry point (`points`) to its nearest feature in `targets`
    (SWORD node or reach geometries, whichever prj.rivers.target_id_col
    is configured for), dropping points farther than max_distance_meters
    from any target.

    Return points (unprojected, original CRS) with target_id_col and a
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


