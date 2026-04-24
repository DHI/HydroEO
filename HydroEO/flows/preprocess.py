"""Timeseries preprocess flow orchestration."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class PreprocessFlow:
    """Run extract/clean/merge in sequence for selected missions."""

    reservoirs: object
    to_process: list[str]
    processing_options: dict

    def run(self) -> None:
        self.reservoirs.extract_product_timeseries(self.to_process)
        self.reservoirs.clean_product_timeseries(
            products=self.to_process,
            filter_options_by_product=self.processing_options,
        )
        self.reservoirs.merge_product_timeseries(products=self.to_process)
