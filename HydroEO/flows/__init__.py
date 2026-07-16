"""Standalone flow functions for HydroEO pipelines.

These functions implement the core download, processing, and visualization
logic previously embedded in Reservoirs and Rivers classes. They operate on
Project state and external data, with no direct method dependencies.

Splits single flows.py file by concern (init/download/extract-clean-merge/run-config/summaries)
Every name below (public and private) is re-exported here so that
`from HydroEO import flows; flows.<name>` keeps working exactly as it did
against the old single-file module, including for tests that patch
private helpers via `patch.object(flows, "_name")`.
"""
import mikeio  # noqa: F401 -- re-exported so `flows.mikeio` resolves (see _reservoir_pipeline.py)

from HydroEO import plotting  # noqa: F401 -- re-exported so `flows.plotting`/`HydroEO.flows.plotting` resolves

from ._reservoir_init import (
    _assign_pld_id as _assign_pld_id,
    _download_pld as _download_pld,
    _flag_missing_priors as _flag_missing_priors,
    initialize_reservoirs as initialize_reservoirs,
)
from ._river_init import (
    _ensure_sword_database as _ensure_sword_database,
    _prepare_rivers_from_sword as _prepare_rivers_from_sword,
    initialize_rivers as initialize_rivers,
)
from ._sentinel_shared import (
    _download_sentinel_for_target as _download_sentinel_for_target,
    _sentinel6_use_earthdata as _sentinel6_use_earthdata,
)
from ._reservoir_download import (
    _download_reservoirs_icesat2 as _download_reservoirs_icesat2,
    _download_reservoirs_sentinel as _download_reservoirs_sentinel,
    _download_reservoirs_swot as _download_reservoirs_swot,
    download_reservoirs as download_reservoirs,
)
from ._river_download import (
    _download_rivers_icesat2 as _download_rivers_icesat2,
    _download_rivers_sentinel as _download_rivers_sentinel,
    _download_swot_hydrocron_timeseries as _download_swot_hydrocron_timeseries,
    _get_latest_hydrocron_obs_date as _get_latest_hydrocron_obs_date,
    download_rivers as download_rivers,
)
from ._river_common import (
    _assign_points_to_river_targets as _assign_points_to_river_targets,
    _group_river_targets_by_waterbody as _group_river_targets_by_waterbody,
    _iter_geometry_pieces as _iter_geometry_pieces,
    _river_extraction_buffer_meters as _river_extraction_buffer_meters,
    _simplify_to_one_polygon as _simplify_to_one_polygon,
)
from ._river_pipeline import (
    _clean_rivers_timeseries as _clean_rivers_timeseries,
    _extract_rivers_icesat2_observations as _extract_rivers_icesat2_observations,
    _extract_rivers_sentinel_observations as _extract_rivers_sentinel_observations,
    _extract_rivers_swot_observations as _extract_rivers_swot_observations,
    _extract_rivers_timeseries as _extract_rivers_timeseries,
    _merge_rivers_timeseries as _merge_rivers_timeseries,
    create_rivers_timeseries as create_rivers_timeseries,
)
from ._reservoir_pipeline import (
    _clean_reservoirs_timeseries as _clean_reservoirs_timeseries,
    _export_cleaned_to_dfs0 as _export_cleaned_to_dfs0,
    _extract_icesat2_observations as _extract_icesat2_observations,
    _extract_reservoirs_timeseries as _extract_reservoirs_timeseries,
    _extract_sentinel_observations as _extract_sentinel_observations,
    _extract_swot_observations as _extract_swot_observations,
    _merge_reservoirs_timeseries as _merge_reservoirs_timeseries,
    create_reservoirs_timeseries as create_reservoirs_timeseries,
)
from ._constants import (
    DEFAULT_RESERVOIR_MERGING_OPTIONS as DEFAULT_RESERVOIR_MERGING_OPTIONS,
    DEFAULT_RIVER_MERGING_OPTIONS as DEFAULT_RIVER_MERGING_OPTIONS,
    PRODUCT_TIMESERIES_KEYS as PRODUCT_TIMESERIES_KEYS,
)
from ._clean_engine import (
    _clean_timeseries as _clean_timeseries,
)
from ._run_config import (
    _apply_exclusions as _apply_exclusions,
    _apply_reach_slope_correction as _apply_reach_slope_correction,
    _default_run_config as _default_run_config,
    _exclusion_value_matches as _exclusion_value_matches,
    _fit_reach_slope_correction as _fit_reach_slope_correction,
    _get_or_fit_spatial_correction_model as _get_or_fit_spatial_correction_model,
    _get_target_ids as _get_target_ids,
    _invalidate_reach_slope_correction_cache as _invalidate_reach_slope_correction_cache,
    _invalidate_spatial_correction_cache as _invalidate_spatial_correction_cache,
    _load_run_config as _load_run_config,
    _reservoir_centroid as _reservoir_centroid,
    _run_config_path as _run_config_path,
    _save_run_config as _save_run_config,
    _target_centroid as _target_centroid,
    exclude_from_target as exclude_from_target,
    list_exclusions as list_exclusions,
    list_target_observations as list_target_observations,
    remove_exclusion as remove_exclusion,
    set_merging_option as set_merging_option,
)
from ._merge_engine import (
    _merge_timeseries as _merge_timeseries,
)
from ._summaries import (
    _has_enough_observations_to_plot as _has_enough_observations_to_plot,
    _load_and_parse_cleaned_timeseries as _load_and_parse_cleaned_timeseries,
    _load_merged_timeseries as _load_merged_timeseries,
    _load_product_timeseries as _load_product_timeseries,
    _project_num_months as _project_num_months,
    _river_target_corridor as _river_target_corridor,
    generate_reservoirs_summaries as generate_reservoirs_summaries,
    generate_rivers_summaries as generate_rivers_summaries,
)

__all__ = [
    # Reservoirs
    "initialize_reservoirs",
    "download_reservoirs",
    "create_reservoirs_timeseries",
    "generate_reservoirs_summaries",
    # Rivers
    "initialize_rivers",
    "download_rivers",
    "create_rivers_timeseries",
    "generate_rivers_summaries",
    # Per-target exclusion/merging-option API (used by Project)
    "list_target_observations",
    "exclude_from_target",
    "list_exclusions",
    "remove_exclusion",
    "set_merging_option",
]
