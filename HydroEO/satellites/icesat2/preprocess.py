from __future__ import annotations

import datetime
import os

import geopandas as gpd
import pandas as pd


def extract_observations(
    src_dir,
    dst_path,
    features,
    atl13_fields=None,  # kept for API compat; no-op in SlideRule path
    track_keys=None,  # kept for API compat; no-op in SlideRule path
):
    """Extract ICESat-2 ATL13 observations from cached parquet and save shapefile."""
    parquet_path = os.path.join(src_dir, "atl13.parquet")
    if not os.path.exists(parquet_path):
        return

    gdf = gpd.read_parquet(parquet_path)
    if gdf.empty:
        return

    # Filter observations to those inside the reservoir geometry.
    gdf = gdf.loc[gdf.within(features.unary_union)].reset_index(drop=True)
    if len(gdf) == 0:
        return

    gdf["platform"] = "icesat2"
    gdf["product"] = "ATL13"

    # Ensure CRS matches features; convert if necessary.
    if gdf.crs is None:
        gdf = gdf.set_crs(features.crs)
    elif gdf.crs != features.crs:
        gdf = gdf.to_crs(features.crs)

    # Strip timezone info so shapefiles can represent the date column.
    if "date" in gdf.columns:
        gdf["date"] = pd.to_datetime(gdf["date"]).dt.tz_localize(None)

    gdf = _to_shapefile_safe_columns(gdf)
    gdf.to_file(dst_path)


def _to_shapefile_safe_columns(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Truncate column names to 10 chars (ESRI Shapefile field name limit)."""
    rename_map: dict[str, str] = {}
    used: set[str] = set()

    for col in gdf.columns:
        if col == "geometry":
            continue

        base = str(col)[:10]
        candidate = base
        suffix = 1
        while candidate in used:
            suffix_str = str(suffix)
            candidate = f"{base[: 10 - len(suffix_str)]}{suffix_str}"
            suffix += 1

        used.add(candidate)
        if candidate != col:
            rename_map[col] = candidate

    return gdf.rename(columns=rename_map)


def get_latest_obs_date(data_dir):
    shp_path = os.path.join(data_dir, "icesat2.shp")
    if os.path.exists(shp_path):
        gdf = gpd.read_file(shp_path)
        last_obs_date = max(gdf.date.values).astype(datetime.date)
        return last_obs_date
