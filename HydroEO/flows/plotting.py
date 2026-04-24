"""Plotting flow orchestration."""

from __future__ import annotations

from dataclasses import dataclass
import logging

logger = logging.getLogger(__name__)


@dataclass
class PlottingFlow:
    """Generate per-reservoir plotting summaries."""

    reservoirs: object

    def run(self, show: bool = False, save: bool = True) -> None:
        for reservoir_id in self.reservoirs.download_gdf[self.reservoirs.id_key]:
            logger.info("Summarizing crossings")
            self.reservoirs.summarize_crossings_by_id(
                reservoir_id,
                show=show,
                save=save,
            )

            logger.info("Summarizing cleaning results")
            self.reservoirs.summarize_cleaning_by_id(
                reservoir_id,
                show=show,
                save=save,
            )

            logger.info("Summarizing merged results")
            self.reservoirs.summarize_merging_by_id(
                reservoir_id,
                show=show,
                save=save,
            )
