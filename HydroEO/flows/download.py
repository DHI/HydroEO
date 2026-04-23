"""Download flow orchestration.

This module contains a non-breaking orchestration layer that delegates to
existing Reservoirs download APIs.
"""

from __future__ import annotations

from dataclasses import dataclass
import logging
from typing import Callable

logger = logging.getLogger(__name__)


@dataclass
class ReservoirDownloadFlow:
    """Run mission downloads while preserving existing behavior."""

    reservoirs: object
    to_download: list[str]
    startdates: dict[str, list[int]]
    enddates: dict[str, list[int]]
    earthdata_credentials: tuple[str | None, str | None]
    creodias_credentials_provider: Callable[[], tuple[str, str]]

    def run(
        self,
        update_existing: bool = False,
        enddate_overrides: dict[str, list[int]] | None = None,
    ) -> None:
        enddate_overrides = enddate_overrides or {}

        for mission in self.to_download:
            if mission == "swot":
                self.reservoirs.download_altimetry(
                    product="SWOT_LAKE",
                    startdate=self.startdates["swot"],
                    enddate=enddate_overrides.get("swot", self.enddates["swot"]),
                    update_existing=update_existing,
                )
                continue

            if mission == "icesat2":
                self.reservoirs.download_altimetry(
                    product="ATL13",
                    startdate=self.startdates["icesat2"],
                    enddate=enddate_overrides.get("icesat2", self.enddates["icesat2"]),
                    update_existing=update_existing,
                    credentials=self.earthdata_credentials,
                )
                continue

            if mission in ["sentinel3", "sentinel6"]:
                product = "S3" if mission == "sentinel3" else "S6"
                sentinel_creds = self.creodias_credentials_provider()
                self.reservoirs.download_altimetry(
                    product=product,
                    startdate=self.startdates[mission],
                    enddate=enddate_overrides.get(mission, self.enddates[mission]),
                    credentials=sentinel_creds,
                    update_existing=update_existing,
                )
                continue

            logger.warning("Skipping unsupported mission in download flow: %s", mission)
