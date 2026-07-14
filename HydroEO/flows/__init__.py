"""Standalone flow functions for HydroEO pipelines.

These functions implement the core download, processing, and visualization
logic previously embedded in Reservoirs and Rivers classes. They operate on
Project state and external data, with no direct method dependencies.

This package replaces what used to be a single ~3,100-line flows.py file,
split by concern (init/download/extract-clean-merge/run-config/summaries)
-- see each submodule's docstring for why its particular grouping was
chosen. Every name below (public and private) is re-exported here so that
`from HydroEO import flows; flows.<name>` keeps working exactly as it did
against the old single-file module, including for tests that patch
private helpers via `patch.object(flows, "_name")`.
"""

import mikeio  # noqa: F401 -- re-exported so `flows.mikeio` resolves (see _reservoir_pipeline.py)

from HydroEO import plotting  # noqa: F401 -- re-exported so `flows.plotting`/`HydroEO.flows.plotting` resolves

from ._reservoir_init import (
    _assign_pld_id,
    _download_pld,
    _flag_missing_priors,
    initialize_reservoirs,
)
from ._river_init import (
    _ensure_sword_database,
    _prepare_rivers_from_sword,
    initialize_rivers,
)
from ._sentinel_shared import (
    _download_sentinel_for_target,
    _sentinel6_use_earthdata,
)
from ._reservoir_download import (
    _download_reservoirs_icesat2,
    _download_reservoirs_sentinel,
    _download_reservoirs_swot,
    download_reservoirs,
)
from ._river_download import (
    _download_rivers_icesat2,
    _download_rivers_sentinel,
    _download_swot_hydrocron_timeseries,
    _get_latest_hydrocron_obs_date,
    download_rivers,
)
from ._river_common import (
    _assign_points_to_river_targets,
    _group_river_targets_by_waterbody,
    _iter_geometry_pieces,
    _river_extraction_buffer_meters,
    _simplify_to_one_polygon,
)
from ._river_pipeline import (
    _clean_rivers_timeseries,
    _extract_rivers_icesat2_observations,
    _extract_rivers_sentinel_observations,
    _extract_rivers_swot_observations,
    _extract_rivers_timeseries,
    _merge_rivers_timeseries,
    create_rivers_timeseries,
)
from ._reservoir_pipeline import (
    _clean_reservoirs_timeseries,
    _export_cleaned_to_dfs0,
    _extract_icesat2_observations,
    _extract_reservoirs_timeseries,
    _extract_sentinel_observations,
    _extract_swot_observations,
    _merge_reservoirs_timeseries,
    create_reservoirs_timeseries,
)
from ._constants import (
    DEFAULT_RESERVOIR_MERGING_OPTIONS,
    DEFAULT_RIVER_MERGING_OPTIONS,
    PRODUCT_TIMESERIES_KEYS,
)
from ._clean_engine import (
    _clean_timeseries,
)
from ._run_config import (
    _apply_exclusions,
    _apply_reach_slope_correction,
    _default_run_config,
    _exclusion_value_matches,
    _fit_reach_slope_correction,
    _get_or_fit_spatial_correction_model,
    _get_target_ids,
    _invalidate_reach_slope_correction_cache,
    _invalidate_spatial_correction_cache,
    _load_run_config,
    _reservoir_centroid,
    _run_config_path,
    _save_run_config,
    _target_centroid,
    exclude_from_target,
    list_exclusions,
    list_target_observations,
    remove_exclusion,
    set_merging_option,
)
from ._merge_engine import (
    _merge_timeseries,
)
from ._summaries import (
    _has_enough_observations_to_plot,
    _load_and_parse_cleaned_timeseries,
    _load_merged_timeseries,
    _load_product_timeseries,
    _project_num_months,
    _river_target_corridor,
    generate_reservoirs_summaries,
    generate_rivers_summaries,
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
