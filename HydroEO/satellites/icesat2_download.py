from __future__ import annotations

import datetime
import logging
import os

import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon

from HydroEO.utils import general

logger = logging.getLogger(__name__)

# Column mapping: SlideRule atl13x default output -> HydroEO schema
SR_DEFAULT_COLUMN_MAP: dict[str, str] = {
    "ht_ortho": "height",  # orthometric water surface height (EGM2008, server-side)
    "cycle": "cycle_number",
    "gt": "beam",
    # rgt is kept as-is (no rename needed)
}

# Ancillary beam-group fields fetchable via parms["atl13_fields"]
SR_ANCILLARY_COLUMN_MAP: dict[str, str] = {
    "inland_water_body_id": "wb_id",
    "inland_water_body_size": "wb_size",
    "inland_water_body_type": "wb_type",
    "segment_slope_trk_bdy": "wb_slope",
    "err_ht_water_surf": "height_err",
    "segment_quality": "quality_seg",
    "segment_geoid": "geoid_track",
    "segment_geoid_free2mean": "geoid_corr_track",
    "segment_dem_ht": "dem",
    "segment_near_sat_fract": "sat_frac_track",
}

# Valid ancillary field names for config validation.
SR_ATL13_VALID_ANCILLARY_FIELDS: set[str] = set(SR_ANCILLARY_COLUMN_MAP)

# Kept for backward API compatibility — no ancillary fields are fetched by default.
ATL13_DEFAULT_FIELDS: list[str] = []

ATL13_CORE_FIELDS = ["height", "lat", "lon", "date"]


def query(
    aoi: list,
    startdate: datetime.date,
    enddate: datetime.date,
    download_directory: str | None,
    atl13_options: dict | None = None,
    atl13_fields: list[str] | None = None,
    sliderule_client=None,
) -> gpd.GeoDataFrame:
    """Download ATL13 water surface elevations via SlideRule's atl13x endpoint."""
    from HydroEO.system import HydroEODownloadError

    if sliderule_client is None:
        raise ValueError("sliderule_client must be provided")

    poly = Polygon(aoi)
    centroid = poly.representative_point()

    atl13_sub: dict = {"coord": {"lon": centroid.x, "lat": centroid.y}}
    atl13_sub.update(atl13_options or {})

    parms: dict = {
        "atl13": atl13_sub,  # required — water-body AMS lookup
        "t0": startdate.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "t1": enddate.strftime("%Y-%m-%dT%H:%M:%SZ"),
    }

    if atl13_fields:
        parms["atl13_fields"] = atl13_fields

    try:
        logger.info(
            "Submitting ICESat-2 ATL13 request to SlideRule; this may take time."
        )
        gdf = sliderule_client.run("atl13x", parms)
    except FileNotFoundError as exc:
        raise HydroEODownloadError(
            f"SlideRule lookup failed for centroid "
            f"(lon={centroid.x:.4f}, lat={centroid.y:.4f})."
        ) from exc
    except Exception as exc:
        raise HydroEODownloadError(f"SlideRule atl13x request failed: {exc}") from exc

    if gdf is None or gdf.empty:
        raise HydroEODownloadError(
            f"SlideRule atl13x returned no data for centroid "
            f"(lon={centroid.x:.4f}, lat={centroid.y:.4f}). "
            "This usually means the reservoir is not registered in the AMS "
            "database (GRWL / HydroLAKES). Verify that the water body appears "
            "in HydroLAKES and that the polygon centroid falls inside it."
        )

    rename_map = {k: v for k, v in SR_DEFAULT_COLUMN_MAP.items() if k in gdf.columns}
    if atl13_fields:
        rename_map.update(
            {k: v for k, v in SR_ANCILLARY_COLUMN_MAP.items() if k in gdf.columns}
        )
    gdf = gdf.rename(columns=rename_map)

    gdf["lat"] = gdf.geometry.y
    gdf["lon"] = gdf.geometry.x
    gdf["date"] = gdf.index
    gdf = gdf.reset_index(drop=True)

    start_dt = pd.Timestamp(startdate).tz_localize("UTC")
    end_dt = pd.Timestamp(enddate).tz_localize("UTC")
    date_col = pd.to_datetime(gdf["date"], utc=True)
    gdf = gdf.loc[(date_col >= start_dt) & (date_col <= end_dt)].reset_index(drop=True)

    if download_directory is not None:
        general.ifnotmakedirs(download_directory)
        gdf.to_parquet(os.path.join(download_directory, "atl13.parquet"))

    return gdf
