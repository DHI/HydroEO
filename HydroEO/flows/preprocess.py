"""Timeseries preprocess flow orchestration."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from HydroEO.waterbody import WaterBody


@dataclass
class PreprocessFlow:
    """Run extract/clean/merge in sequence for selected missions.

    Parameters
    ----------
    reservoirs : WaterBody
        WaterBody instance (Reservoirs or Rivers) with extract/clean/merge methods.
    to_process : list[str]
        Product names to process (e.g., ['swot', 'icesat2', 'sentinel3']).
    processing_options : dict[str, dict[str, Any]]
        Per-product filter options mapping product name to filter config.
    """

    reservoirs: WaterBody
    to_process: list[str]
    processing_options: dict[str, dict[str, Any]]

    def run(self) -> None:
        self.reservoirs.extract_product_timeseries(self.to_process)
        self.reservoirs.clean_product_timeseries(
            products=self.to_process,
            filter_options_by_product=self.processing_options,
        )
        self.reservoirs.merge_product_timeseries(products=self.to_process)
