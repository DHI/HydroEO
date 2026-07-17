"""
Shared constants for the merge/clean pipeline (both reservoirs and rivers).
"""

# Default cleaning-filter parameters, used by _clean_engine.py as the
# fallback whenever a product's processing_options doesn't specify these
# explicitly. Re-exported from HydroEO.constant.
from HydroEO.constants import (
    DEFAULT_PROCESSING_FILTERS as DEFAULT_PROCESSING_FILTERS,
    DEFAULT_ELEVATION_MIN_M as DEFAULT_ELEVATION_MIN_M,
    DEFAULT_ELEVATION_MAX_M as DEFAULT_ELEVATION_MAX_M,
    DEFAULT_MAD_THRESHOLD as DEFAULT_MAD_THRESHOLD,
)


PRODUCT_TIMESERIES_KEYS = {
    "sentinel3": dict(
        lat_key="lat", lon_key="lon", pass_key="file_name",
        platform_key="platform", orbit_key="relative_orbit",
    ),
    # NOTE: sentinel6 still uses "pass" as orbit_key -- NOT verified
    "sentinel6": dict(
        lat_key="lat", lon_key="lon", pass_key="file_name",
        platform_key="platform", orbit_key="pass",
    ),
    "icesat2": dict(
        lat_key="lat", lon_key="lon", pass_key="pass",
        platform_key="platform", orbit_key="beam",
    ),
    "swot": dict(
        platform_key="platform", orbit_key="orbit", preset_error_key="wse_u",
    ),
}

# Default .merge() tuning for reservoirs, mirroring the shape of
# processing_options (a project-level dict of pipeline parameters) but
# applied once per reservoir rather than per-product, since .merge() runs
# on the already-combined multi-product timeseries. Override via
# prj.merging_options in project config; falls back to these reservoir-
# appropriate defaults if that attribute isn't set.
DEFAULT_RESERVOIR_MERGING_OPTIONS = {
    "window_km": 1.5,
    "svr_linear_err": 0.1,
    "svr_linear_epsilon": 0.1,
    "svr_radial_err": 1.0, # DAHITI default err=0.1 too strict for reservoirs with real multi-week transitions
    "svr_radial_rbf_c": 10000,
    "svr_radial_gamma": 0.0000438 * 50, # DAHITI default gamma=0.0000438 (~151 days) too smooth for reservoirs with real multi-week transitions
    "svr_radial_epsilon": 0.1,
    "bias_time_bin": "20D", #Based on examples - sparse sources can fail and get dropped silently but in order of magnitude of S3 return period
    "bias_min_overlap": 1,
    "bias_group_by": "platform_orbit", # Should represent individual observation groups (platform and orbit/beam), see note for S6
    "bias_centroid_warn_km": 5.0, #Warning threshold for large/elongated reservoirs 
    "distance_penalty_scale_per_km": None, # Off by default: inflate Kalman input error by distance from the reservoir centroid
    "use_spatial_correction": False, # Off by default: apply a spatial correction model fit from a dense spatial source (ICESat-2) to other missions' observations, if enough qualifying dense-source days exist
    "spatial_correction_dense_source": "icesat2",
}


DEFAULT_RIVER_MERGING_OPTIONS = {
    # Mostly a starting point copied from the reservoir defaults
    # NOT independently validated against real river data.
    "window_km": 1.5,
    "svr_linear_err": 0.1,
    "svr_linear_epsilon": 0.1,
    "svr_radial_err": 1.0, #Matches reservoir default, DAHITI default err=0.1 too strict for rivers with real multi-week transitions
    "svr_radial_rbf_c": 10000,
    "svr_radial_gamma": 0.0000438 * 100, #gamma x100 (~15-day lengthscale, vs DAHITI's ~151-day lake value)
    "svr_radial_epsilon": 0.1,
    "bias_time_bin": "20D",
    "bias_min_overlap": 1,
    "bias_group_by": "platform_orbit", # Less likely for river targets but useful flag
    "bias_centroid_warn_km": 5.0, # Relevant for reaches spanning multiple km 
    "distance_penalty_scale_per_km": None,
    "use_spatial_correction": False,
    "spatial_correction_dense_source": "icesat2",
    "use_reach_slope_correction": False, #Off by default, meaningful to correct non-SWOT, uses SWOT slope. Not validated.
}


