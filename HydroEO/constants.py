"""HydroEO constants and default configuration values.

Centralized location for all constants, defaults, and configuration parameters
organized by satellite mission and functional domain.
"""

from typing import Any

from HydroEO.satellites.icesat2 import ATL13_DEFAULT_FIELDS

# ============================================================================
# Version
# ============================================================================

__version__ = "0.2.4"

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

# Project.__sat_init and HydroEO.flows._clean_engine both
# reference these constants
DEFAULT_PROCESSING_FILTERS = ["elevation", "MAD"]
DEFAULT_ELEVATION_MIN_M = 0.0
DEFAULT_ELEVATION_MAX_M = 8000.0
DEFAULT_MAD_THRESHOLD = 5.0

# ============================================================================
# Mission Defaults
# ============================================================================

MISSION_DEFAULTS = {
    "swot": {
        "download": False,
        "process": False,
        # Minimum PLD-lake overlap as a percentage (0-100) of the
        # reservoir's own area. Previously used pld_match_max_distance_m.
        # NOTE: Note backward compatible.
        "pld_match_min_overlap_pct": 10.0,
        "exclude_obs_id_values": ["no_data"],
        "hydrocron_fields": SWOT_DEFAULT_HYDROCRON_FIELDS,
        "quality_filters": SWOT_DEFAULT_QUALITY_FILTERS,
        "processing_filters": DEFAULT_PROCESSING_FILTERS,
        "elevation_min_m": DEFAULT_ELEVATION_MIN_M,
        "elevation_max_m": DEFAULT_ELEVATION_MAX_M,
        "mad_threshold": DEFAULT_MAD_THRESHOLD,
    },
    "icesat2": {
        "download": False,
        "process": False,
        "atl13_fields": ICESAT2_DEFAULT_FIELDS,
        "atl13": {"pass_invalid": False, "beams": [], "spots": []},
        "track_keys": ICESAT2_SUPPORTED_TRACK_KEYS,
        "processing_filters": DEFAULT_PROCESSING_FILTERS,
        "elevation_min_m": DEFAULT_ELEVATION_MIN_M,
        "elevation_max_m": DEFAULT_ELEVATION_MAX_M,
        "mad_threshold": DEFAULT_MAD_THRESHOLD,
    },
    "sentinel3": {
        "download": False,
        "process": False,
        "subset_file_id": "enhanced_measurement.nc",
        "sigma0_max": 1e5,
        "download_threads": 1,
        "processing_filters": DEFAULT_PROCESSING_FILTERS,
        "elevation_min_m": DEFAULT_ELEVATION_MIN_M,
        "elevation_max_m": DEFAULT_ELEVATION_MAX_M,
        "mad_threshold": DEFAULT_MAD_THRESHOLD,
    },
    "sentinel6": {
        "download": False,
        "process": False,
        "subset_file_id": "enhanced_measurement.nc",
        "sigma0_max": 1e5,
        "download_threads": 1,
        "processing_filters": DEFAULT_PROCESSING_FILTERS,
        "elevation_min_m": DEFAULT_ELEVATION_MIN_M,
        "elevation_max_m": DEFAULT_ELEVATION_MAX_M,
        "mad_threshold": DEFAULT_MAD_THRESHOLD,
    },
}

# ============================================================================
# River Profile (SWOT longitudinal WSE profile extraction)
# ============================================================================

RIVER_PROFILE_DEFAULT_AOI_BUFFER_M = 2_000.0

RIVER_PROFILE_DEFAULT_FILTERS: dict[str, dict[str, Any]] = {
    "preclip": {"enabled": True, "min": -5.0, "max": 10.0},
    "soft_clamp": {
        "enabled": True,
        "bin_width_m": 18_000.0,
        "y_bin_m": 0.25,
        "fixed_halfw_m": 0.40,
        "mad_factor": 1.8,
        "min_count": 4,
        "min_mode_count": 2,
        "slope_gain": 0.35,
        "huber_k": 1.5,
    },
    "hampel_1": {
        "enabled": True,
        "win_m": 20_000.0,
        "sigma": 1.6,
        "min_valid": 8,
        "action": "mask",
        "replace_mode": "trend",
        "huber_k": 1.6,
        "huber_iters": 2,
    },
    "rolling_quantile": {
        "enabled": True,
        "win_m": 12_000.0,
        "q": 0.50,
        "min_valid": 8,
        "robust_iters": 2,
        "huber_k": 1.5,
    },
    "density_cull": {
        "enabled": True,
        "total_win_m": 5_000.0,
        "low_pct": 3.0,
        "abs_min": 20,
        "dilate": 0,
    },
    "hampel_2": {
        "enabled": True,
        "win_m": 20_000.0,
        "sigma": 1.6,
        "min_valid": 8,
        "action": "mask",
        "replace_mode": "trend",
        "huber_k": 1.6,
        "huber_iters": 2,
    },
    "spline_fill": {
        "enabled": True,
        "k": 3,
        "inner_knots": None,
        "target_spacing_m": 18_000.0,
        "weight_scheme": "by_density",
    },
}
"""Default parameters for each stage of the river-profile filtering pipeline.
See configs/river_profile.md for what each stage does. Deep-merged with the
user's `river_profile.filters` config in
HydroEO.satellites.swot.river_profile._resolve_filters.
"""