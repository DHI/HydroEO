"""
Shared timeseries cleaning engine for both reservoirs and rivers.
"""

import logging
import os

import geopandas as gpd
from tqdm import tqdm

from HydroEO.utils import general, timeseries
from ._constants import (
    PRODUCT_TIMESERIES_KEYS,
    DEFAULT_PROCESSING_FILTERS,
    DEFAULT_ELEVATION_MIN_M,
    DEFAULT_ELEVATION_MAX_M,
    DEFAULT_MAD_THRESHOLD,
)
from ._run_config import _get_target_ids
from ._summaries import _load_product_timeseries

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from HydroEO.project import Project

logger = logging.getLogger(__name__)


def _clean_timeseries(prj: "Project", target_type: str) -> None:
    """Apply quality filters to extracted timeseries, for either
    reservoirs or river targets (nodes/reaches)."""
    target_ids = _get_target_ids(prj, target_type)
    ids_with_raw = [
        id
        for id in target_ids
        if os.path.exists(os.path.join(prj.dirs["output"], f"{id}", "raw_observations"))
    ]
    if not ids_with_raw:
        logger.warning(
            "No raw observations found for any %s; skipping timeseries cleaning.",
            target_type,
        )
        return

    for id in tqdm(ids_with_raw, desc=f"Cleaning product timeseries ({target_type})"):
        for product in prj.to_process:
            df = _load_product_timeseries(
                os.path.join(prj.dirs["output"], f"{id}", "raw_observations"),
                ".gpkg",
                [product],
                lambda path: gpd.read_file(path).drop(columns=["geometry"]),
            )
            if df is not None:
                product_options = prj.processing_options.get(
                    product,
                    {
                        "processing_filters": ["elevation", "MAD"],
                        "elevation_min_m": 0.0,
                        "elevation_max_m": 8000.0,
                        "mad_threshold": 5.0,
                    },
                )

                ts = timeseries.Timeseries(
                    df, date_key="date", height_key="height",
                    **PRODUCT_TIMESERIES_KEYS.get(product, {}),
                )

                ts.clean(
                    product_options.get("processing_filters", ["elevation", "MAD"]),
                    filter_params={
                        "elevation_min_m": product_options.get("elevation_min_m", 0.0),
                        "elevation_max_m": product_options.get(
                            "elevation_max_m", 8000.0
                        ),
                        "mad_threshold": product_options.get("mad_threshold", 5.0),
                    },
                )

                export_dir = os.path.join(
                    prj.dirs["output"], f"{id}", "cleaned_observations"
                )
                general.ifnotmakedirs(export_dir)
                ts.export_csv(os.path.join(export_dir, f"{product}.csv"))
