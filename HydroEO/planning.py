"""Dry-run planning: what a run would download, without downloading anything.

Only metadata searches (NASA CMR) are made. Plans cover the two scenarios a UI
needs before the user commits to a long download: ``river_profile`` and
``reservoirs``.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

from HydroEO.satellites.swot._download import SWOT_LAKE_SHORT_NAME, _login, _search

if TYPE_CHECKING:
    from HydroEO.project import Project

logger = logging.getLogger(__name__)

# More granules than this is slow on a small machine; surfaced as a warning.
LARGE_GRANULE_COUNT = 500


@dataclass
class PlanItem:
    """One dataset the run would fetch."""

    label: str
    granules: int | None  # None = not countable before the run (per-target queries)
    size_mb: float | None = None
    note: str = ""


@dataclass
class Plan:
    """What a run would do, plus anything the user should know first."""

    scenario: str
    items: list[PlanItem] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    bbox: list[float] | None = None  # lon/lat: west, south, east, north
    targets: int | None = None

    @property
    def total_granules(self) -> int:
        return sum(i.granules or 0 for i in self.items)

    @property
    def total_size_mb(self) -> float:
        return sum(i.size_mb or 0.0 for i in self.items)


def _granule_size_mb(results: list) -> float | None:
    """Total size of CMR results in MB; ``None`` if CMR does not report sizes."""
    try:
        # earthaccess >= 1.0 makes ``size`` an attribute; older versions a method.
        sizes = [r.size() if callable(r.size) else r.size for r in results]
        return float(sum(sizes))
    except Exception:
        return None


def _date_str(date: list[int], end: bool) -> str:
    return f"{date[0]:04d}-{date[1]:02d}-{date[2]:02d} " + ("23:59:59" if end else "00:00:00")


def _count_swot(
    short_name: str, bbox: tuple[float, ...], startdate: list[int], enddate: list[int],
    granule_filter: str | None = None,
) -> tuple[int, float | None]:
    params: dict[str, Any] = {
        "short_name": short_name,
        "bounding_box": tuple(bbox),
        "temporal": (_date_str(startdate, False), _date_str(enddate, True)),
    }
    if granule_filter:
        params["granule_name"] = granule_filter
    results = _search(**params)
    return len(results), _granule_size_mb(results)


def _finish(plan: Plan) -> Plan:
    if plan.total_granules > LARGE_GRANULE_COUNT:
        plan.warnings.append(
            f"{plan.total_granules} granules will be downloaded; consider a shorter period or smaller area."
        )
    return plan


def plan_river_profile(config: dict[str, Any], global_crs: str = "EPSG:4326") -> Plan:
    """Plan a river-profile run from its ``river_profile`` config section.

    Raises the same errors as the run itself for an invalid chainage file.
    """
    from HydroEO.satellites.swot.river_profile import _build_swot_raster_config, _load_chainage

    name = config.get("name", "river_profile")
    gdf, _ = _load_chainage(config)
    raster_cfg = _build_swot_raster_config(config, name, gdf, global_crs)
    bbox = raster_cfg["aoi"]["bbox"]

    _login()
    count, size = _count_swot(
        raster_cfg["product"], tuple(bbox), raster_cfg["startdate"], raster_cfg["enddate"],
        raster_cfg.get("granule_filter"),
    )
    plan = Plan("river_profile", bbox=bbox, targets=len(gdf))
    plan.items.append(PlanItem(raster_cfg["product"], count, size))
    if count == 0:
        plan.warnings.append("No SWOT raster granules found for this reach and period.")
    return _finish(plan)


def plan_reservoirs(prj: "Project") -> Plan:
    """Plan a reservoirs run for an already constructed ``Project``.

    SWOT is counted with one CMR search over the targets' envelope (before the
    prior-lake filter, so an upper bound). ICESat-2 and Sentinel are queried per
    reservoir during the download and are listed without a count.
    """
    gdf = prj.reservoirs.gdf.to_crs("EPSG:4326")
    bbox = [float(v) for v in gdf.total_bounds]
    plan = Plan("reservoirs", bbox=bbox, targets=len(gdf))

    for mission in prj.to_download:
        if mission == "swot":
            _login()
            count, size = _count_swot(
                SWOT_LAKE_SHORT_NAME, tuple(bbox), _as_list(prj.startdates[mission]),
                _as_list(prj.enddates[mission]),
            )
            plan.items.append(PlanItem(SWOT_LAKE_SHORT_NAME, count, size, "upper bound; prior lakes only are kept"))
            if count == 0:
                plan.warnings.append("No SWOT Lake granules found for this area and period.")
        else:
            plan.items.append(PlanItem(mission, None, None, "queried per reservoir during the download"))
    return _finish(plan)


def _as_list(date: Any) -> list[int]:
    """Config dates are ``[y, m, d]`` lists or ``datetime.date``; return a list."""
    return list(date) if isinstance(date, (list, tuple)) else [date.year, date.month, date.day]
