"""Tests for the river-profile interactive nodes map."""

import json
import logging
import re
from unittest.mock import patch

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import Point

from HydroEO.satellites.swot.river_profile import calculate_river_profile
from HydroEO.satellites.swot.river_profile_map import _label_to_iso, write_nodes_map
from tests.unit.test_river_profile import (  # noqa: F401
    NAME,
    _base_config,
    project_with_wse_tile,
)

pytestmark = pytest.mark.unit


def _embedded_data(html_path):
    page = html_path.read_text(encoding="utf-8")
    m = re.search(r'<script id="data" type="application/json">(.*?)</script>', page, re.DOTALL)
    return json.loads(m.group(1).replace("<\\/", "</"))


def _gdf(n=4):
    return gpd.GeoDataFrame(
        {
            "node_id": np.arange(n),
            "geometry": [Point(500_000.0, 2_600_000.0 + i * 200.0) for i in range(n)],
        },
        crs="EPSG:32645",
    )


@pytest.mark.parametrize(
    "label, expected",
    [
        ("20240115T000000", "2024-01-15T00:00:00"),
        ("20240115T013045_1", "2024-01-15T01:30:45"),
        ("2024-03-01_12-00-00", "2024-03-01T12:00:00"),
        ("some_granule_stem", None),
    ],
)
def test_label_to_iso(label, expected):
    assert _label_to_iso(label) == expected


def test_write_nodes_map_embeds_data_in_chronological_order(tmp_path):
    gdf = _gdf()
    x = np.array([0.0, 200.0, 400.0, 600.0])
    final_cols = {
        "20240301T000000": np.array([10.0, 10.1234567, np.nan, 10.3]),
        "20240115T000000": np.array([9.0, 9.1, 9.2, 9.3]),
    }
    out = tmp_path / "map.html"

    write_nodes_map(gdf, x, list(final_cols), final_cols, out, river_name="r</script>x")

    data = _embedded_data(out)
    assert data["node_id"] == [0, 1, 2, 3]
    assert data["km"] == [0.0, 0.2, 0.4, 0.6]
    assert data["labels"] == ["20240115T000000", "20240301T000000"]
    assert data["dates"] == ["2024-01-15T00:00:00", "2024-03-01T00:00:00"]
    assert data["wse"][1] == [10.0, 10.123, None, 10.3]
    assert data["name"] == "r</script>x"
    assert "Download data as CSV" in out.read_text(encoding="utf-8")
    assert all(-90 <= v <= 90 for v in data["lat"])


def test_write_nodes_map_keeps_unparseable_labels_last(tmp_path):
    final_cols = {"weird": np.ones(4), "20240115T000000": np.zeros(4)}
    out = tmp_path / "map.html"

    write_nodes_map(_gdf(), np.arange(4.0), list(final_cols), final_cols, out, river_name="r")

    data = _embedded_data(out)
    assert data["labels"] == ["20240115T000000", "weird"]
    assert data["dates"] == ["2024-01-15T00:00:00", None]


def test_calculate_river_profile_writes_nodes_map_and_node_ids(
    project_with_wse_tile,  # noqa: F811
    caplog,
):
    project_dir, chainage_path = project_with_wse_tile
    config = _base_config(project_dir, chainage_path)

    with (
        patch("HydroEO.satellites.swot.river_profile.download_raster"),
        caplog.at_level(logging.INFO, logger="HydroEO.satellites.swot.river_profile"),
    ):
        calculate_river_profile(config, project_dir=str(project_dir), global_crs="EPSG:4326")

    results_dir = project_dir / "results" / NAME
    map_path = results_dir / f"{NAME}_nodes_map.html"
    assert map_path.resolve().as_uri() in caplog.text
    df = pd.read_csv(results_dir / f"{NAME}_profiles_final.csv")
    assert list(df.columns[:2]) == ["node_id", "distance_along_river_m"]
    assert df["node_id"].tolist() == list(range(25))

    shp = next((results_dir / "profiles_final").glob("*.shp"))
    assert gpd.read_file(shp)["node_id"].tolist() == list(range(25))

    data = _embedded_data(results_dir / f"{NAME}_nodes_map.html")
    assert data["node_id"] == list(range(25))
    assert data["labels"] == [c for c in df.columns if c not in ("node_id", "distance_along_river_m")]
    assert np.allclose(
        [np.nan if v is None else v for v in data["wse"][0]],
        df[data["labels"][0]].round(3),
        equal_nan=True,
    )
