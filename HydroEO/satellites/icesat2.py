"""Compatibility facade for ICESat-2 ATL13 workflows.

This module preserves the historical public API while delegating to explicit
flow surfaces split by responsibility.
"""

from __future__ import annotations

import datetime

import geopandas as gpd
import sliderule

from HydroEO.satellites import icesat2_download, icesat2_preprocess

# Initialise SlideRule client once at module load time.
# sliderule.init() only sets the server URL — it does not open a network
# connection — so it is safe to call at import time.
sliderule.init("slideruleearth.io")

SR_DEFAULT_COLUMN_MAP = icesat2_download.SR_DEFAULT_COLUMN_MAP
SR_ANCILLARY_COLUMN_MAP = icesat2_download.SR_ANCILLARY_COLUMN_MAP
SR_ATL13_VALID_ANCILLARY_FIELDS = icesat2_download.SR_ATL13_VALID_ANCILLARY_FIELDS
ATL13_DEFAULT_FIELDS = icesat2_download.ATL13_DEFAULT_FIELDS
ATL13_CORE_FIELDS = icesat2_download.ATL13_CORE_FIELDS


def query(
    aoi: list,
    startdate: datetime.date,
    enddate: datetime.date,
    download_directory: str | None,
    atl13_options: dict | None = None,
    atl13_fields: list[str] | None = None,
) -> gpd.GeoDataFrame:
    return icesat2_download.query(
        aoi=aoi,
        startdate=startdate,
        enddate=enddate,
        download_directory=download_directory,
        atl13_options=atl13_options,
        atl13_fields=atl13_fields,
        sliderule_client=sliderule,
    )


def extract_observations(
    src_dir,
    dst_path,
    features,
    atl13_fields=None,  # kept for API compat; no-op in SlideRule path
    track_keys=None,  # kept for API compat; no-op in SlideRule path
):
    return icesat2_preprocess.extract_observations(
        src_dir=src_dir,
        dst_path=dst_path,
        features=features,
        atl13_fields=atl13_fields,
        track_keys=track_keys,
    )


def _to_shapefile_safe_columns(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    return icesat2_preprocess._to_shapefile_safe_columns(gdf)


def get_latest_obs_date(data_dir):
    return icesat2_preprocess.get_latest_obs_date(data_dir=data_dir)


__all__ = [
    "SR_DEFAULT_COLUMN_MAP",
    "SR_ANCILLARY_COLUMN_MAP",
    "SR_ATL13_VALID_ANCILLARY_FIELDS",
    "ATL13_DEFAULT_FIELDS",
    "ATL13_CORE_FIELDS",
    "query",
    "extract_observations",
    "_to_shapefile_safe_columns",
    "get_latest_obs_date",
]
