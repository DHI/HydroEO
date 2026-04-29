"""HydroEO constants and default configuration values.

Centralized location for all constants, defaults, and configuration parameters
organized by satellite mission and functional domain.
"""

from HydroEO.satellites.icesat2 import ATL13_DEFAULT_FIELDS

# ============================================================================
# Version
# ============================================================================

__version__ = "0.1.0"

# ============================================================================
# Data Constants
# ============================================================================

FLOAT32_NODATA_VALUE = -99999.0
"""NODATA sentinel value for float32 rasters"""

# ============================================================================
# ICESat-2 Configuration
# ============================================================================

ICESAT2_DEFAULT_FIELDS = ATL13_DEFAULT_FIELDS
ICESAT2_REQUIRED_FIELDS: list[str] = []
"""SlideRule's atl13x always returns core fields (height, lat/lon, date, beam, rgt, 
cycle_number) — no forced field merging is needed."""

ICESAT2_SUPPORTED_TRACK_KEYS = ["gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"]

# ============================================================================
# SWOT Configuration
# ============================================================================

SWOT_DEFAULT_HYDROCRON_FIELDS = {
    "nodes": [
        "node_id",
        "node_q",
        "reach_id",
        "time_str",
        "wse",
        "wse_u",
        "p_wse",
        "geoid_hght",
        "sword_version",
        "solid_tide",
        "load_tidef",
        "pole_tide",
        "width",
        "width_u",
        "p_width",
        "xovr_cal_q",
        "rdr_sig0",
        "xovr_cal_c",
        "dark_frac",
    ],
    "reaches": [
        "reach_id",
        "reach_q",
        "time_str",
        "wse",
        "wse_u",
        "slope",
        "slope_u",
        "slope2",
        "slope2_u",
        "width",
        "width_u",
        "geoid_hght",
        "solid_tide",
        "load_tidef",
        "pole_tide",
        "p_wse",
        "p_width",
    ],
}

SWOT_DEFAULT_QUALITY_FILTERS = {
    "nodes": {"max_q": 2},
    "reaches": {"max_q": 2},
}

# ============================================================================
# Timeseries & Processing
# ============================================================================

SUPPORTED_CLEAN_FILTERS = ["elevation", "MAD", "daily_mean", "hampel", "rolling_median"]

# ============================================================================
# Mission Defaults
# ============================================================================

MISSION_DEFAULTS = {
    "swot": {
        "download": False,
        "process": False,
        "pld_match_max_distance_m": 100.0,
        "exclude_obs_id_values": ["no_data"],
        "hydrocron_fields": SWOT_DEFAULT_HYDROCRON_FIELDS,
        "quality_filters": SWOT_DEFAULT_QUALITY_FILTERS,
        "processing_filters": ["elevation", "MAD"],
        "elevation_min_m": 0.0,
        "elevation_max_m": 8000.0,
        "mad_threshold": 5.0,
    },
    "icesat2": {
        "download": False,
        "process": False,
        "atl13_fields": ICESAT2_DEFAULT_FIELDS,
        "atl13": {"pass_invalid": False, "beams": [], "spots": []},
        "track_keys": ICESAT2_SUPPORTED_TRACK_KEYS,
        "processing_filters": ["elevation", "MAD"],
        "elevation_min_m": 0.0,
        "elevation_max_m": 8000.0,
        "mad_threshold": 5.0,
    },
    "sentinel3": {
        "download": False,
        "process": False,
        "subset_file_id": "enhanced_measurement.nc",
        "sigma0_max": 1e5,
        "download_threads": 1,
        "processing_filters": ["elevation", "MAD"],
        "elevation_min_m": 0.0,
        "elevation_max_m": 8000.0,
        "mad_threshold": 5.0,
    },
    "sentinel6": {
        "download": False,
        "process": False,
        "subset_file_id": "enhanced_measurement.nc",
        "sigma0_max": 1e5,
        "download_threads": 1,
        "processing_filters": ["elevation", "MAD"],
        "elevation_min_m": 0.0,
        "elevation_max_m": 8000.0,
        "mad_threshold": 5.0,
    },
}
