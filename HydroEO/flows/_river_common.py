"""River geometry helpers shared by download, extraction, and summaries.

None of these are themselves the target of a sibling patch in the test
suite (only _river_target_corridor is, which is why that one function
lives in _summaries.py instead -- see generate_rivers_summaries's test).
"""

import logging

import geopandas as gpd

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


