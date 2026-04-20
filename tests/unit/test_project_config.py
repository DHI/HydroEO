"""Stage 2 unit tests for Project config defaults and validation."""

from pathlib import Path

import geopandas as gpd
import pytest
import yaml
from shapely.geometry import box


@pytest.fixture
def _mock_reservoir_gdf():
    return gpd.GeoDataFrame(
        {"project": ["demo"], "geometry": [box(6.0, 46.2, 6.9, 46.6)]},
        crs="EPSG:4326",
    )


def _write_config(path: Path, data: dict):
    path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")


@pytest.mark.unit
def test_project_applies_stage2_defaults(tmp_path, monkeypatch, _mock_reservoir_gdf):
    """Optional Stage 2 mission keys should default when omitted."""
    from HydroEO.project import Project

    reservoir_path = tmp_path / "reservoirs.shp"
    reservoir_path.write_text("placeholder", encoding="utf-8")

    cfg_path = tmp_path / "config.yaml"
    _write_config(
        cfg_path,
        {
            "project": {"main_dir": str(tmp_path / "out")},
            "gis": {"global_crs": "EPSG:4326"},
            "reservoirs": {"path": str(reservoir_path), "id_key": "project"},
            "icesat2": {
                "download": False,
                "process": True,
                "startdate": [2024, 1, 1],
                "enddate": [2024, 2, 1],
            },
            "sentinel3": {
                "download": False,
                "process": False,
                "startdate": [2024, 1, 1],
                "enddate": [2024, 2, 1],
            },
        },
    )

    monkeypatch.setattr(
        "HydroEO.project.gpd.read_file",
        lambda *_args, **_kwargs: _mock_reservoir_gdf.copy(),
    )

    proj = Project(name="defaults", config=str(cfg_path))

    assert proj.config["icesat2"]["mad_threshold"] == 5.0
    assert proj.config["icesat2"]["track_keys"] == [
        "gt1l",
        "gt1r",
        "gt2l",
        "gt2r",
        "gt3l",
        "gt3r",
    ]
    assert proj.config["sentinel3"]["sigma0_max"] == 1e5
    assert "ht_ortho" in proj.config["icesat2"]["atl13_fields"]
    assert "delta_time" in proj.config["icesat2"]["atl13_fields"]


@pytest.mark.unit
def test_project_invalid_optional_values_raise_clear_error(
    tmp_path,
    monkeypatch,
    _mock_reservoir_gdf,
):
    """Invalid Stage 2 optional values should raise a descriptive validation error."""
    from HydroEO.project import Project

    reservoir_path = tmp_path / "reservoirs.shp"
    reservoir_path.write_text("placeholder", encoding="utf-8")

    cfg_path = tmp_path / "config.yaml"
    _write_config(
        cfg_path,
        {
            "project": {"main_dir": str(tmp_path / "out")},
            "reservoirs": {"path": str(reservoir_path), "id_key": "project"},
            "icesat2": {
                "download": False,
                "process": True,
                "startdate": [2024, 1, 1],
                "enddate": [2024, 2, 1],
                "atl13_fields": ["ht_ortho", "not_a_real_field"],
            },
        },
    )

    monkeypatch.setattr(
        "HydroEO.project.gpd.read_file",
        lambda *_args, **_kwargs: _mock_reservoir_gdf.copy(),
    )

    with pytest.raises(ValueError, match="Invalid ATL13 fields"):
        Project(name="bad-optional", config=str(cfg_path))


@pytest.mark.unit
def test_validate_config_reports_multiple_common_issues_at_once():
    """validate_config() should collect and report multiple errors in one exception."""
    from HydroEO.project import Project

    proj = Project.__new__(Project)
    proj.config = {
        "project": {"main_dir": ""},
        "reservoirs": {"path": "/path/does/not/exist.shp"},
        "sentinel3": {
            "download": "yes",
            "process": True,
            "startdate": [2024, 1],
            "enddate": [2024, 2, 1],
            "sigma0_max": -10,
        },
    }

    with pytest.raises(ValueError) as exc_info:
        proj.validate_config()

    msg = str(exc_info.value)
    assert "project.main_dir" in msg
    assert "reservoirs.id_key" in msg
    assert "sentinel3.download" in msg
    assert "sentinel3.startdate" in msg
    assert "sentinel3.sigma0_max" in msg
