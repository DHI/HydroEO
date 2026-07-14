"""Shared merge engine used by both reservoirs and rivers.

_merge_timeseries applies the merge()/Kalman/svr_radial pipeline
generically for either target_type -- called by
_merge_reservoirs_timeseries (in _reservoir_pipeline.py) and
_merge_rivers_timeseries (in _river_pipeline.py). Not itself patched as a
sibling of either wrapper in the test suite, so it's free to live in its
own module.
"""

import logging
import os

import pandas as pd
from tqdm import tqdm

from HydroEO.utils import general, timeseries
from ._constants import (
    PRODUCT_TIMESERIES_KEYS,
    DEFAULT_RESERVOIR_MERGING_OPTIONS,
    DEFAULT_RIVER_MERGING_OPTIONS,
)
from ._run_config import (
    _get_target_ids,
    _target_centroid,
    _load_run_config,
    _apply_exclusions,
    _get_or_fit_spatial_correction_model,
    _fit_reach_slope_correction,
    _apply_reach_slope_correction,
)
from ._summaries import _load_product_timeseries

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from HydroEO.project import Project

logger = logging.getLogger(__name__)


def _merge_timeseries(prj: "Project", target_type: str) -> None:
    """Merge multi-mission timeseries into combined datasets, for either
    reservoirs or river targets (nodes/reaches)."""
    target_ids = _get_target_ids(prj, target_type)
    ids_with_cleaned = [
        id
        for id in target_ids
        if os.path.exists(
            os.path.join(prj.dirs["output"], f"{id}", "cleaned_observations")
        )
    ]
    if not ids_with_cleaned:
        logger.warning(
            "No cleaned observations found for any %s; skipping timeseries merging.",
            target_type,
        )
        return

    default_options = (
        DEFAULT_RESERVOIR_MERGING_OPTIONS
        if target_type == "reservoirs"
        else DEFAULT_RIVER_MERGING_OPTIONS
    )
    # prj.reservoirs.merging_options / prj.rivers.merging_options are the
    # intended per-target-type override locations (set from the
    # respective YAML config sections); fall back to the older shared
    # prj.merging_options for backward compatibility if the per-type one
    # isn't set yet.
    target_owner = prj.reservoirs if target_type == "reservoirs" else prj.rivers
    overrides = getattr(target_owner, "merging_options", None)
    if overrides is None:
        overrides = getattr(prj, "merging_options", None)

    for id in tqdm(ids_with_cleaned, desc=f"Merging product timeseries ({target_type})"):
        ts_list = []
        for product in prj.to_process:
            df = _load_product_timeseries(
                os.path.join(prj.dirs["output"], f"{id}", "cleaned_observations"),
                ".csv",
                [product],
                pd.read_csv,
            )
            if df is not None:
                df["date"] = pd.to_datetime(
                    df.date, format="mixed", utc=True
                ).dt.tz_convert(None)
                df = df.sort_values(by="date")
                ts_list.append(
                    timeseries.Timeseries(
                        df, date_key="date", height_key="height",
                        **PRODUCT_TIMESERIES_KEYS.get(product, {}),
                    )
                )

        if len(ts_list) > 0:
            ts = timeseries.concat(ts_list)

            data_dir = os.path.join(prj.dirs["output"], f"{id}")
            general.ifnotmakedirs(data_dir)

            run_config = _load_run_config(prj, id)

            merging_options = dict(default_options)
            merging_options.update(overrides or {})
            # Per-target overrides (from notebook calls to
            # set_merging_option, or hand-edited in run_config.yaml) take
            # priority over project-wide settings for just this target.
            merging_options.update(run_config.get("merging_option_overrides", {}))
            distance_penalty_scale = merging_options.pop(
                "distance_penalty_scale_per_km", None
            )
            use_spatial_correction = merging_options.pop(
                "use_spatial_correction", False
            )
            spatial_correction_dense_source = merging_options.pop(
                "spatial_correction_dense_source", "icesat2"
            )
            # Off by default. Only meaningful for reach-mode river
            # projects (a node is a single ~200m point, not a ~10km
            # segment with its own along-reach slope) -- see
            # _fit_reach_slope_correction for the full reasoning and the
            # sign-convention caveat that should be checked against real
            # data before trusting this in production.
            use_reach_slope_correction = merging_options.pop(
                "use_reach_slope_correction", False
            )

            # Apply exclusions BEFORE exporting all_cleaned_timeseries.csv
            # (not just before merge processing) -- this file is meant to
            # reflect what's actually being worked with, and writing it
            # before exclusions were applied meant it always showed
            # excluded data regardless of how many times you re-ran,
            # which looked exactly like a stale file from an old run but
            # was actually happening on every single run. The full,
            # pre-exclusion record is still available per-mission in
            # cleaned_observations/{product}.csv (written earlier, in
            # _clean_timeseries, before any exclusion is applied) -- so
            # nothing is lost by making this file reflect exclusions.
            exclusions = run_config.get("exclusions", [])
            if exclusions:
                before = len(ts.df)
                ts.df = _apply_exclusions(ts.df, exclusions)
                logger.info(
                    "%s: %d exclusion rule(s) applied, %d/%d observations kept.",
                    id, len(exclusions), len(ts.df), before,
                )

            if (
                use_reach_slope_correction
                and target_type == "rivers"
                and getattr(prj.rivers, "target_id_col", None) == "reach_id"
            ):
                slope_model = _fit_reach_slope_correction(prj, id)
                if slope_model is not None:
                    ts.df = _apply_reach_slope_correction(ts.df, prj, id, slope_model)
                    logger.info(
                        "%s: applied reach slope correction (median_slope=%.6g "
                        "from %d SWOT observations).",
                        id, slope_model["median_slope"], slope_model["n_observations"],
                    )
                else:
                    logger.info(
                        "%s: use_reach_slope_correction is enabled but no "
                        "usable SWOT slope data was found; skipping "
                        "correction for this target.", id,
                    )

            ts.export_csv(os.path.join(data_dir, "all_cleaned_timeseries.csv"))

            ref_lat, ref_lon = _target_centroid(prj, target_type, id)

            spatial_correction_model = None
            if use_spatial_correction:
                spatial_correction_model = _get_or_fit_spatial_correction_model(
                    prj, target_type, id,
                    dense_source_platform=spatial_correction_dense_source,
                )

            ts = ts.merge(
                save_progress=True,
                dir=os.path.join(data_dir, "merged_progress"),
                ref_lat=ref_lat,
                ref_lon=ref_lon,
                distance_penalty_scale_per_km=distance_penalty_scale,
                spatial_correction_model=spatial_correction_model,
                **merging_options,
            )
            ts.export_csv(os.path.join(data_dir, "merged_timeseries.csv"))


