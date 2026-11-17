"""Unit tests for HydroEO.progress and HydroEO.planning (no network)."""

import datetime
from types import SimpleNamespace
from unittest.mock import patch

import geopandas as gpd
import pytest
from shapely.geometry import box

from HydroEO import planning, progress
from tests.unit.test_river_profile import _base_config, _make_chainage_shp, project_with_wse_tile  # noqa: F401

pytestmark = pytest.mark.unit


def test_report_without_callback_is_a_noop():
    progress.set_callback(None)
    progress.report("x", 0, 1)


def test_listen_restores_previous_callback():
    outer, inner = [], []
    with progress.listen(outer.append):
        with progress.listen(inner.append):
            progress.report("a", 1, 2)
        progress.report("b", 2, 2)
    progress.report("c", 0, 0)
    assert [e.stage for e in inner] == ["a"]
    assert [e.stage for e in outer] == ["b"]


def test_failing_callback_does_not_raise():
    def boom(_):
        raise RuntimeError("x")

    with progress.listen(boom):
        progress.report("a", 0, 1)


def test_river_profile_emits_download_and_process_events(project_with_wse_tile):  # noqa: F811
    from HydroEO.satellites.swot.river_profile import calculate_river_profile

    project_dir, chainage_path = project_with_wse_tile
    events = []
    with (
        patch("HydroEO.satellites.swot.river_profile.download_raster"),
        progress.listen(events.append),
    ):
        calculate_river_profile(_base_config(project_dir, chainage_path), str(project_dir))
    stages = [(e.stage, e.done, e.total) for e in events]
    assert stages[0] == ("download", 0, 1)
    assert ("download", 1, 1) in stages
    assert stages[-1] == ("process", 1, 1)


class _Granule:
    def size(self):
        return 10.0


def test_plan_river_profile_counts_granules(tmp_path):
    chainage = _make_chainage_shp(tmp_path)
    config = _base_config(tmp_path, chainage)
    with (
        patch("HydroEO.planning._login"),
        patch("HydroEO.planning._search", return_value=[_Granule()] * 3) as search,
    ):
        plan = planning.plan_river_profile(config)
    assert plan.total_granules == 3
    assert plan.total_size_mb == 30.0
    assert plan.targets == 25
    assert search.call_args.kwargs["short_name"] == config.get("product", "SWOT_L2_HR_Raster_D")
    assert not plan.warnings


def test_plan_river_profile_warns_when_nothing_found(tmp_path):
    config = _base_config(tmp_path, _make_chainage_shp(tmp_path))
    with patch("HydroEO.planning._login"), patch("HydroEO.planning._search", return_value=[]):
        plan = planning.plan_river_profile(config)
    assert plan.total_granules == 0
    assert "No SWOT raster granules" in plan.warnings[0]


def test_plan_reservoirs_counts_swot_and_lists_others():
    gdf = gpd.GeoDataFrame({"id": [1, 2]}, geometry=[box(0, 0, 1, 1), box(2, 2, 3, 3)], crs="EPSG:4326")
    prj = SimpleNamespace(
        reservoirs=SimpleNamespace(gdf=gdf),
        to_download=["swot", "icesat2"],
        startdates={"swot": [2024, 1, 1]},
        enddates={"swot": datetime.date(2024, 2, 1)},
    )
    with (
        patch("HydroEO.planning._login"),
        patch("HydroEO.planning._search", return_value=[_Granule()] * 2) as search,
    ):
        plan = planning.plan_reservoirs(prj)
    assert plan.targets == 2
    assert plan.bbox == [0.0, 0.0, 3.0, 3.0]
    assert [(i.label, i.granules) for i in plan.items] == [("SWOT_L2_HR_LakeSP_D", 2), ("icesat2", None)]
    assert search.call_args.kwargs["temporal"] == ("2024-01-01 00:00:00", "2024-02-01 23:59:59")


def test_large_plan_gets_warning():
    plan = planning.Plan("x", items=[planning.PlanItem("p", planning.LARGE_GRANULE_COUNT + 1)])
    assert planning._finish(plan).warnings


def test_granule_size_supports_attribute_and_method():
    class AsAttribute:
        size = 4.0

    assert planning._granule_size_mb([AsAttribute(), _Granule()]) == 14.0
    assert planning._granule_size_mb([object()]) is None
