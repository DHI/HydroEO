import datetime
from pathlib import Path

import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import Point

from HydroEO.waterbody import Rivers


@pytest.fixture
def empty_rivers_gdf():
    return gpd.GeoDataFrame({"geometry": []}, geometry="geometry", crs="EPSG:4326")


def make_rivers(tmp_path: Path, empty_rivers_gdf, id_key: str = "river_id"):
    rivers = Rivers(
        gdf=empty_rivers_gdf.copy(),
        id_key=id_key,
        dirs={
            "main": str(tmp_path),
            "swot": str(tmp_path / "swot"),
        },
    )
    rivers.mission_options = {
        "swot": {
            "hydrocron_fields": {
                "nodes": ["node_id", "node_q", "time_str", "wse"],
                "reaches": ["reach_id", "reach_q", "time_str", "wse"],
            },
            "quality_filters": {
                "nodes": {"max_q": 2},
                "reaches": {"max_q": 1},
            },
        }
    }
    return rivers


@pytest.mark.unit
def test_download_swot_hydrocron_single_node_configured_id_filters_quality(
    tmp_path, empty_rivers_gdf, monkeypatch
):
    """One node, configured_id fallback. Output: swot/rivers/<wb_id>/timeseries.csv."""
    rivers = make_rivers(tmp_path, empty_rivers_gdf)
    rivers.target_id_col = "node_id"
    rivers.target_ids = [101]
    rivers.configured_id = "gauja"
    captured = {}

    def fake_request(target_id, startdate, enddate, fields, feature):
        captured.update(
            {
                "target_id": target_id,
                "startdate": startdate,
                "enddate": enddate,
                "fields": fields,
                "feature": feature,
            }
        )
        return (
            200,
            {
                "results": {
                    "csv": "node_id,node_q,time_str,wse\n101,1,2024-01-01T00:00:00Z,10.0\n101,3,2024-01-02T00:00:00Z,11.0\n101,2,2024-01-03T00:00:00Z,12.0\n"
                }
            },
            {},
        )

    monkeypatch.setattr(rivers, "_request_hydrocron_timeseries", fake_request)

    summary = rivers.download_swot_hydrocron(
        startdate=[2024, 1, 1],
        enddate=[2024, 2, 1],
        update_existing=False,
    )

    assert summary == {
        "requested": 1,
        "successful": 1,
        "failed": 0,
        "empty_after_filter": 0,
    }
    assert captured == {
        "target_id": 101,
        "startdate": datetime.date(2024, 1, 1),
        "enddate": datetime.date(2024, 2, 1),
        "fields": ["node_id", "node_q", "time_str", "wse"],
        "feature": "Node",
    }

    output_path = tmp_path / "swot" / "rivers" / "gauja" / "timeseries.csv"
    saved = pd.read_csv(output_path)
    assert saved["node_q"].tolist() == [1, 2]
    assert saved["node_id"].tolist() == [101, 101]


@pytest.mark.unit
def test_download_swot_hydrocron_multiple_nodes_same_waterbody_produces_one_csv(
    tmp_path, empty_rivers_gdf, monkeypatch
):
    """Two nodes belonging to the same AOI feature must be merged into one CSV."""
    rivers = make_rivers(tmp_path, empty_rivers_gdf)
    rivers.target_id_col = "node_id"
    rivers.target_ids = [101, 102]
    rivers.target_features = gpd.GeoDataFrame(
        {
            "node_id": [101, 102],
            "river_id": ["gauja", "gauja"],
            "geometry": [Point(6.0, 46.0), Point(6.1, 46.1)],
        },
        geometry="geometry",
        crs="EPSG:4326",
    )

    def fake_request(target_id, startdate, enddate, fields, feature):
        _ = (startdate, enddate, fields, feature)
        return (
            200,
            {
                "results": {
                    "csv": f"node_id,node_q,time_str,wse\n{target_id},1,2024-01-01T00:00:00Z,10.0\n"
                }
            },
            {},
        )

    monkeypatch.setattr(rivers, "_request_hydrocron_timeseries", fake_request)

    summary = rivers.download_swot_hydrocron(
        startdate=[2024, 1, 1], enddate=[2024, 2, 1], update_existing=False
    )

    assert summary == {
        "requested": 1,
        "successful": 1,
        "failed": 0,
        "empty_after_filter": 0,
    }

    output_path = tmp_path / "swot" / "rivers" / "gauja" / "timeseries.csv"
    assert output_path.exists()
    df = pd.read_csv(output_path)
    assert sorted(df["node_id"].tolist()) == [101, 102]


@pytest.mark.unit
def test_download_swot_hydrocron_multiple_waterbodies_produce_separate_csvs(
    tmp_path, empty_rivers_gdf, monkeypatch
):
    """Nodes belonging to different AOI features each produce their own CSV."""
    rivers = make_rivers(tmp_path, empty_rivers_gdf)
    rivers.target_id_col = "node_id"
    rivers.target_ids = [101, 201]
    rivers.target_features = gpd.GeoDataFrame(
        {
            "node_id": [101, 201],
            "river_id": ["gauja", "daugava"],
            "geometry": [Point(6.0, 46.0), Point(6.5, 46.5)],
        },
        geometry="geometry",
        crs="EPSG:4326",
    )

    def fake_request(target_id, startdate, enddate, fields, feature):
        _ = (startdate, enddate, fields, feature)
        return (
            200,
            {
                "results": {
                    "csv": f"node_id,node_q,time_str,wse\n{target_id},1,2024-01-01T00:00:00Z,10.0\n"
                }
            },
            {},
        )

    monkeypatch.setattr(rivers, "_request_hydrocron_timeseries", fake_request)

    summary = rivers.download_swot_hydrocron(
        startdate=[2024, 1, 1], enddate=[2024, 2, 1], update_existing=False
    )

    assert summary == {
        "requested": 2,
        "successful": 2,
        "failed": 0,
        "empty_after_filter": 0,
    }
    gauja_path = tmp_path / "swot" / "rivers" / "gauja" / "timeseries.csv"
    daugava_path = tmp_path / "swot" / "rivers" / "daugava" / "timeseries.csv"
    assert gauja_path.exists()
    assert daugava_path.exists()
    assert pd.read_csv(gauja_path)["node_id"].iloc[0] == 101
    assert pd.read_csv(daugava_path)["node_id"].iloc[0] == 201


@pytest.mark.unit
def test_download_swot_hydrocron_reaches_uses_target_feature_mapping(
    tmp_path, empty_rivers_gdf, monkeypatch
):
    rivers = make_rivers(tmp_path, empty_rivers_gdf, id_key="river_name")
    rivers.target_id_col = "reach_id"
    rivers.target_ids = [2001]
    rivers.target_features = gpd.GeoDataFrame(
        {
            "reach_id": [2001],
            "river_name": ["mainstem"],
            "geometry": [Point(6.0, 46.0)],
        },
        geometry="geometry",
        crs="EPSG:4326",
    )

    def fake_request(target_id, startdate, enddate, fields, feature):
        _ = (target_id, startdate, enddate, fields, feature)
        return (
            200,
            {
                "results": {
                    "csv": "reach_id,reach_q,time_str,wse\n2001,1,2024-01-01T00:00:00Z,10.0\n2001,2,2024-01-02T00:00:00Z,11.0\n"
                }
            },
            {},
        )

    monkeypatch.setattr(rivers, "_request_hydrocron_timeseries", fake_request)

    summary = rivers.download_swot_hydrocron(
        startdate=[2024, 1, 1], enddate=[2024, 2, 1], update_existing=False
    )

    output_path = tmp_path / "swot" / "rivers" / "mainstem" / "timeseries.csv"
    assert summary == {
        "requested": 1,
        "successful": 1,
        "failed": 0,
        "empty_after_filter": 0,
    }
    saved = pd.read_csv(output_path)
    assert saved["reach_q"].tolist() == [1]
    assert saved["reach_id"].tolist() == [2001]


@pytest.mark.unit
def test_download_swot_hydrocron_failed_and_empty_waterbodies(
    tmp_path, empty_rivers_gdf, monkeypatch, caplog
):
    """Waterbody with all-failed requests counts as failed.
    Waterbody with data but all filtered out counts as empty_after_filter."""
    rivers = make_rivers(tmp_path, empty_rivers_gdf)
    rivers.target_id_col = "node_id"
    rivers.target_ids = [101, 102]
    rivers.target_features = gpd.GeoDataFrame(
        {
            "node_id": [101, 102],
            "river_id": ["alpha", "beta"],
            "geometry": [Point(6.0, 46.0), Point(6.1, 46.1)],
        },
        geometry="geometry",
        crs="EPSG:4326",
    )

    def fake_request(target_id, startdate, enddate, fields, feature):
        _ = (startdate, enddate, fields, feature)
        if target_id == 101:
            return 500, {}, {}
        return (
            200,
            {
                "results": {
                    "csv": "node_id,node_q,time_str,wse\n102,5,2024-01-02T00:00:00Z,11.0\n"
                }
            },
            {},
        )

    monkeypatch.setattr(rivers, "_request_hydrocron_timeseries", fake_request)

    with caplog.at_level("INFO"):
        summary = rivers.download_swot_hydrocron(
            startdate=[2024, 1, 1], enddate=[2024, 2, 1], update_existing=False
        )

    assert summary == {
        "requested": 2,
        "successful": 0,
        "failed": 1,
        "empty_after_filter": 1,
    }
    assert "no usable CSV payload" in caplog.text
    assert "empty after filtering" in caplog.text


@pytest.mark.unit
def test_download_swot_hydrocron_uses_existing_latest_obs_for_updates(
    tmp_path, empty_rivers_gdf, monkeypatch
):
    rivers = make_rivers(tmp_path, empty_rivers_gdf)
    rivers.target_id_col = "node_id"
    rivers.target_ids = [101]
    rivers.configured_id = "gauja"

    existing_path = tmp_path / "swot" / "rivers" / "gauja" / "timeseries.csv"
    existing_path.parent.mkdir(parents=True)
    existing_path.write_text(
        "node_id,node_q,time_str,wse\n101,1,2024-01-15T00:00:00Z,10.0\n",
        encoding="utf-8",
    )

    captured = {}

    def fake_request(target_id, startdate, enddate, fields, feature):
        captured["startdate"] = startdate
        _ = (target_id, enddate, fields, feature)
        return (
            200,
            {
                "results": {
                    "csv": "node_id,node_q,time_str,wse\n101,1,2024-01-16T00:00:00Z,10.5\n"
                }
            },
            {},
        )

    monkeypatch.setattr(rivers, "_request_hydrocron_timeseries", fake_request)

    summary = rivers.download_swot_hydrocron(
        startdate=[2024, 1, 1], enddate=[2024, 2, 1], update_existing=True
    )

    assert summary["successful"] == 1
    assert captured["startdate"] == datetime.date(2024, 1, 15)
