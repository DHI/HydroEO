"""ICESat-2 mission public API."""

from __future__ import annotations

import datetime
import logging

import geopandas as gpd
import sliderule

from HydroEO.satellites.icesat2 import download as _download
from HydroEO.satellites.icesat2 import preprocess as _preprocess

# Initialise SlideRule client once at module load time.
# sliderule.init() only sets the server URL — it does not open a network
# connection — so it is safe to call at import time.
sliderule.init("slideruleearth.io")
sliderule.set_verbose(True, loglevel=logging.WARNING)

SR_DEFAULT_COLUMN_MAP = _download.SR_DEFAULT_COLUMN_MAP
SR_ANCILLARY_COLUMN_MAP = _download.SR_ANCILLARY_COLUMN_MAP
SR_ATL13_VALID_ANCILLARY_FIELDS = _download.SR_ATL13_VALID_ANCILLARY_FIELDS
ATL13_DEFAULT_FIELDS = _download.ATL13_DEFAULT_FIELDS
ATL13_CORE_FIELDS = _download.ATL13_CORE_FIELDS


def query(
    aoi: list,
    startdate: datetime.date,
    enddate: datetime.date,
    download_directory: str | None,
    atl13_options: dict | None = None,
    atl13_fields: list[str] | None = None,
) -> gpd.GeoDataFrame:
    return _download.query(
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
    atl13_fields=None,
    track_keys=None,
):
    return _preprocess.extract_observations(
        src_dir=src_dir,
        dst_path=dst_path,
        features=features,
        atl13_fields=atl13_fields,
        track_keys=track_keys,
    )


def _to_shapefile_safe_columns(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    return _preprocess._to_shapefile_safe_columns(gdf)


def get_latest_obs_date(data_dir):
    return _preprocess.get_latest_obs_date(data_dir=data_dir)


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
