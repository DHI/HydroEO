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


@pytest.fixture
def _mock_river_gdf():
    return gpd.GeoDataFrame(
        {
            "reach_id": [1001],
            "node_id": [1001],
            "geometry": [box(6.0, 46.2, 6.1, 46.3)],
        },
        crs="EPSG:4326",
    )


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
    # SlideRule returns core fields (height, lat/lon, date, rgt, cycle_number, beam)
    # by default — atl13_fields is empty unless ancillary fields are explicitly requested.
    assert proj.config["icesat2"]["atl13_fields"] == []
    # The atl13 sub-dict with SlideRule sub-parameters must be injected by defaults.
    assert "atl13" in proj.config["icesat2"]
    assert proj.config["icesat2"]["atl13"]["pass_invalid"] is False


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


@pytest.mark.unit
def test_validate_config_rejects_both_reservoirs_and_rivers_sections():
    from HydroEO.project import Project

    proj = Project.__new__(Project)
    proj.config = {
        "project": {"main_dir": "/tmp/hydroeo"},
        "reservoirs": {"path": "/tmp/reservoirs.shp", "id_key": "rid"},
        "rivers": {"feature_numbers": [1, 2, 3], "feature_type": "reaches", "id": "r"},
    }

    with pytest.raises(ValueError, match="mutually exclusive"):
        proj.validate_config()


@pytest.mark.unit
def test_validate_config_rejects_missing_waterbody_branch():
    from HydroEO.project import Project

    proj = Project.__new__(Project)
    proj.config = {"project": {"main_dir": "/tmp/hydroeo"}}

    with pytest.raises(
        ValueError,
        match="provide one of 'reservoirs', 'rivers', 'swot_raster', 'swot_pixc', or 'river_profile'",
    ):
        proj.validate_config()


@pytest.mark.unit
def test_project_accepts_rivers_aoi_branch(tmp_path, monkeypatch, _mock_river_gdf):
    from HydroEO.project import Project

    aoi_path = tmp_path / "river_aoi.gpkg"
    aoi_path.write_text("placeholder", encoding="utf-8")

    cfg_path = tmp_path / "config.yaml"
    _write_config(
        cfg_path,
        {
            "project": {"main_dir": str(tmp_path / "out")},
            "gis": {"global_crs": "EPSG:4326"},
            "rivers": {
                "aoi_path": str(aoi_path),
                "continent_key": "eu",
                "feature_type": "reaches",
                "id_key": "river_id",
                "buffer_meters": 500.0,
            },
        },
    )

    monkeypatch.setattr(
        "HydroEO.project.gpd.read_file",
        lambda *_args, **_kwargs: _mock_river_gdf.copy(),
    )

    def _mock_initialize_rivers(prj):
        prj.rivers.target_features = _mock_river_gdf.copy()
        prj.rivers.target_id_col = "reach_id"
        prj.rivers.target_ids = [1001]

    monkeypatch.setattr(
        "HydroEO.project.flows.initialize_rivers",
        _mock_initialize_rivers,
    )

    proj = Project(name="rivers-aoi", config=str(cfg_path))
    proj.initialize()
    assert hasattr(proj, "rivers")
    assert not hasattr(proj, "reservoirs")
    assert proj.rivers.target_ids == [1001]


@pytest.mark.unit
def test_project_accepts_rivers_node_number_branch(tmp_path):
    from HydroEO.project import Project

    cfg_path = tmp_path / "config.yaml"
    _write_config(
        cfg_path,
        {
            "project": {"main_dir": str(tmp_path / "out")},
            "gis": {"global_crs": "EPSG:4326"},
            "rivers": {
                "feature_numbers": [10, 11, 12],
                "feature_type": "nodes",
                "id": "demo-river",
            },
        },
    )

    proj = Project(name="rivers-feature-numbers-nodes", config=str(cfg_path))
    assert hasattr(proj, "rivers")
    assert proj.rivers.target_ids == [10, 11, 12]
    assert proj.rivers.target_id_col == "node_id"
    proj.initialize()
    assert proj.rivers.target_ids == [10, 11, 12]


@pytest.mark.unit
def test_initialize_skips_sword_preparation_for_non_aoi_rivers(tmp_path, monkeypatch):
    from HydroEO.project import Project

    cfg_path = tmp_path / "config.yaml"
    _write_config(
        cfg_path,
        {
            "project": {"main_dir": str(tmp_path / "out")},
            "gis": {"global_crs": "EPSG:4326"},
            "rivers": {
                "feature_numbers": [10, 11, 12],
                "feature_type": "nodes",
                "id": "demo-river",
            },
        },
    )

    def _unexpected_prepare(*_args, **_kwargs):
        raise AssertionError(
            "_prepare_rivers_from_sword should not be called for feature_numbers"
        )

    monkeypatch.setattr(
        "HydroEO.flows._river_init._prepare_rivers_from_sword",
        _unexpected_prepare,
    )

    proj = Project(name="rivers-feature-numbers", config=str(cfg_path))
    proj.initialize()
    assert proj.rivers.target_ids == [10, 11, 12]


@pytest.mark.unit
def test_validate_config_rejects_both_node_and_reach_numbers():
    from HydroEO.project import Project

    proj = Project.__new__(Project)
    proj.config = {
        "project": {"main_dir": "/tmp/hydroeo"},
        "rivers": {"node_numbers": [1], "reach_numbers": [2]},
    }

    with pytest.raises(ValueError, match="no longer supported"):
        proj.validate_config()


@pytest.mark.unit
def test_project_applies_swot_hydrocron_defaults_for_rivers(tmp_path):
    from HydroEO.project import Project

    cfg_path = tmp_path / "config.yaml"
    _write_config(
        cfg_path,
        {
            "project": {"main_dir": str(tmp_path / "out")},
            "gis": {"global_crs": "EPSG:4326"},
            "rivers": {
                "feature_numbers": [10],
                "feature_type": "nodes",
                "id": "demo-river",
            },
            "swot": {
                "download": True,
                "process": False,
                "startdate": [2024, 1, 1],
                "enddate": [2024, 2, 1],
            },
        },
    )

    proj = Project(name="rivers-swot-defaults", config=str(cfg_path))

    assert proj.config["swot"]["hydrocron_fields"]["nodes"][0] == "node_id"
    assert proj.config["swot"]["hydrocron_fields"]["reaches"][0] == "reach_id"
    assert proj.config["swot"]["quality_filters"]["nodes"]["max_q"] == 2
    assert proj.mission_options["swot"]["quality_filters"]["reaches"]["max_q"] == 2


@pytest.mark.unit
def test_validate_config_rejects_invalid_swot_hydrocron_shapes():
    from HydroEO.project import Project

    proj = Project.__new__(Project)
    proj.config = {
        "project": {"main_dir": "/tmp/hydroeo"},
        "rivers": {"feature_numbers": [1], "feature_type": "nodes", "id": "demo-river"},
        "swot": {
            "download": True,
            "process": False,
            "startdate": [2024, 1, 1],
            "enddate": [2024, 2, 1],
            "hydrocron_fields": {"nodes": "node_id", "invalid": []},
            "quality_filters": {"nodes": {"max_q": "2"}, "reaches": {}},
        },
    }

    with pytest.raises(ValueError) as exc_info:
        proj.validate_config()

    msg = str(exc_info.value)
    assert "swot.hydrocron_fields" in msg
    assert "swot.quality_filters.nodes.max_q" in msg
    assert "swot.quality_filters.reaches.max_q" in msg


@pytest.mark.unit
def test_validate_config_requires_rivers_id_for_number_inputs():
    from HydroEO.project import Project

    proj = Project.__new__(Project)
    proj.config = {
        "project": {"main_dir": "/tmp/hydroeo"},
        "rivers": {"feature_numbers": [2], "feature_type": "reaches"},
    }

    with pytest.raises(ValueError, match="rivers.id"):
        proj.validate_config()


@pytest.mark.unit
def test_validate_config_rejects_aoi_path_with_number_inputs():
    from HydroEO.project import Project

    proj = Project.__new__(Project)
    proj.config = {
        "project": {"main_dir": "/tmp/hydroeo"},
        "rivers": {
            "aoi_path": "/tmp/aoi.gpkg",
            "feature_numbers": [1],
            "feature_type": "nodes",
            "id": "r",
        },
    }

    with pytest.raises(ValueError, match="mutually exclusive"):
        proj.validate_config()


@pytest.mark.unit
def test_validate_config_requires_feature_keys_for_aoi_path():
    from HydroEO.project import Project

    proj = Project.__new__(Project)
    proj.config = {
        "project": {"main_dir": "/tmp/hydroeo"},
        "rivers": {"aoi_path": "/tmp/aoi.gpkg", "id_key": "river_id"},
    }

    with pytest.raises(ValueError, match="continent_key"):
        proj.validate_config()


@pytest.mark.unit
def test_validate_config_rejects_negative_river_buffer():
    from HydroEO.project import Project

    proj = Project.__new__(Project)
    proj.config = {
        "project": {"main_dir": "/tmp/hydroeo"},
        "rivers": {
            "aoi_path": "/tmp/aoi.gpkg",
            "continent_key": "eu",
            "feature_type": "nodes",
            "id_key": "river_id",
            "buffer_meters": -1,
        },
    }

    with pytest.raises(ValueError, match="buffer_meters"):
        proj.validate_config()


@pytest.mark.unit
def test_project_global_date_fallback(tmp_path, monkeypatch, _mock_reservoir_gdf):
    """Mission sections should inherit project-level startdate/enddate when omitted."""
    from HydroEO.project import Project

    reservoir_path = tmp_path / "reservoirs.shp"
    reservoir_path.write_text("placeholder", encoding="utf-8")

    cfg_path = tmp_path / "config.yaml"
    _write_config(
        cfg_path,
        {
            "project": {
                "main_dir": str(tmp_path / "out"),
                "startdate": [2024, 1, 1],
                "enddate": [2024, 12, 31],
            },
            "gis": {"global_crs": "EPSG:4326"},
            "reservoirs": {"path": str(reservoir_path), "id_key": "project"},
            "icesat2": {"download": False, "process": False},
            "sentinel3": {"download": False, "process": False},
        },
    )

    monkeypatch.setattr(
        "HydroEO.project.gpd.read_file",
        lambda *_args, **_kwargs: _mock_reservoir_gdf.copy(),
    )

    proj = Project(name="global-dates", config=str(cfg_path))
    assert proj.startdates["icesat2"] == [2024, 1, 1]
    assert proj.enddates["icesat2"] == [2024, 12, 31]
    assert proj.startdates["sentinel3"] == [2024, 1, 1]


@pytest.mark.unit
def test_project_no_warning_for_icesat2_with_rivers_configured(tmp_path):
    """ICESat-2 (and Sentinel-3/6) support rivers directly (see
    flows._download_rivers_icesat2/_download_rivers_sentinel), so no
    UserWarning should fire when a rivers section is present -- only when
    NEITHER reservoirs nor rivers is configured (see
    test_project_warns_incompatible_satellites_when_neither_mode_configured
    for that case). This replaces a stale test that expected a warning
    here from before ICESat-2/Sentinel gained river support."""
    from HydroEO.project import Project
    import warnings

    cfg_path = tmp_path / "config.yaml"
    _write_config(
        cfg_path,
        {
            "project": {"main_dir": str(tmp_path / "out")},
            "gis": {"global_crs": "EPSG:4326"},
            "rivers": {"feature_numbers": [10], "feature_type": "nodes", "id": "r"},
            "icesat2": {
                "download": True,
                "process": False,
                "startdate": [2024, 1, 1],
                "enddate": [2024, 2, 1],
            },
        },
    )

    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        Project(name="rivers-no-warn", config=str(cfg_path))


@pytest.mark.unit
def test_project_warns_incompatible_satellites_when_neither_mode_configured(tmp_path):
    """ICESat-2/Sentinel-3/6 configured for download/process without a
    reservoirs or rivers section genuinely has no effect (neither can
    spatially filter observations without one), so this is the one case
    that should still warn."""
    from HydroEO.project import Project

    cfg_path = tmp_path / "config.yaml"
    _write_config(
        cfg_path,
        {
            "project": {
                "main_dir": str(tmp_path / "out"),
                "startdate": [2024, 1, 1],
                "enddate": [2024, 2, 1],
            },
            "gis": {"global_crs": "EPSG:4326"},
            "swot_raster": {
                "aoi": {"name": "aoi", "type": "bbox", "bbox": [0, 0, 1, 1]},
                "product": "SWOT_L2_HR_Raster_D",
                "startdate": [2024, 1, 1],
                "enddate": [2024, 2, 1],
            },
            "icesat2": {
                "download": True,
                "process": False,
                "startdate": [2024, 1, 1],
                "enddate": [2024, 2, 1],
            },
        },
    )

    with pytest.warns(UserWarning, match="icesat2"):
        Project(name="neither-mode-warn", config=str(cfg_path))


@pytest.mark.unit
def test_project_enabled_false_skips_mode(tmp_path):
    """Setting enabled: false on a mode section should skip that mode."""
    from HydroEO.project import Project

    cfg_path = tmp_path / "config.yaml"
    _write_config(
        cfg_path,
        {
            "project": {"main_dir": str(tmp_path / "out")},
            "rivers": {
                "enabled": False,
                "feature_numbers": [10],
                "feature_type": "nodes",
                "id": "r",
            },
            "swot_raster": {
                "enabled": True,
                "aoi": {"name": "test", "type": "bbox", "bbox": [0.0, 0.0, 1.0, 1.0]},
                "product": "SWOT_L2_HR_Raster_D",
                "startdate": [2025, 1, 1],
                "enddate": [2025, 2, 1],
            },
        },
    )

    proj = Project(name="enabled-flag", config=str(cfg_path))
    assert not hasattr(proj, "rivers")
    assert hasattr(proj, "swot_raster_config")


@pytest.mark.unit
def test_swot_raster_falls_back_to_project_dates(tmp_path):
    """swot_raster may omit its own startdate/enddate and rely on
    project-level dates (as documented in configs/swot_raster.md); the
    resolved config passed to download_raster() must still carry concrete
    dates."""
    from HydroEO.project import Project

    cfg_path = tmp_path / "config.yaml"
    _write_config(
        cfg_path,
        {
            "project": {
                "main_dir": str(tmp_path / "out"),
                "startdate": [2024, 1, 1],
                "enddate": [2024, 12, 31],
            },
            "swot_raster": {
                "aoi": {"name": "test", "type": "bbox", "bbox": [0.0, 0.0, 1.0, 1.0]},
                "product": "SWOT_L2_HR_Raster_D",
            },
        },
    )

    proj = Project(name="date-fallback", config=str(cfg_path))
    assert proj.swot_raster_config["startdate"] == [2024, 1, 1]
    assert proj.swot_raster_config["enddate"] == [2024, 12, 31]


@pytest.mark.unit
def test_swot_pixc_falls_back_to_project_dates(tmp_path):
    """Same fallback as swot_raster, exercised for swot_pixc (documented in
    configs/swot_pixc.md) since the two sections are mutually exclusive and
    validated/backfilled by separate code paths."""
    from HydroEO.project import Project

    cfg_path = tmp_path / "config.yaml"
    _write_config(
        cfg_path,
        {
            "project": {
                "main_dir": str(tmp_path / "out"),
                "startdate": [2024, 1, 1],
                "enddate": [2024, 12, 31],
            },
            "swot_pixc": {
                "aoi": {"name": "test", "type": "bbox", "bbox": [0.0, 0.0, 1.0, 1.0]},
                "product": "SWOT_L2_HR_PIXC_D",
            },
        },
    )

    proj = Project(name="pixc-date-fallback", config=str(cfg_path))
    assert proj.swot_pixc_config["startdate"] == [2024, 1, 1]
    assert proj.swot_pixc_config["enddate"] == [2024, 12, 31]


@pytest.mark.unit
def test_validate_config_rejects_swot_raster_without_any_dates():
    from HydroEO.project import Project

    proj = Project.__new__(Project)
    proj.config = {
        "project": {"main_dir": "/tmp/hydroeo"},
        "swot_raster": {
            "aoi": {"name": "test", "type": "bbox", "bbox": [0.0, 0.0, 1.0, 1.0]},
            "product": "SWOT_L2_HR_Raster_D",
        },
    }

    with pytest.raises(ValueError, match="swot_raster.startdate"):
        proj.validate_config()


@pytest.fixture
def _mock_chainage_shp(tmp_path):
    from shapely.geometry import Point

    path = tmp_path / "chainage.shp"
    gdf = gpd.GeoDataFrame(
        {"cngmeters": [0.0, 100.0], "geometry": [Point(0, 0), Point(0, 100)]},
        crs="EPSG:32645",
    )
    gdf.to_file(path)
    return path


@pytest.mark.unit
def test_project_accepts_river_profile_branch(tmp_path, _mock_chainage_shp):
    """A valid, enabled river_profile section should populate
    river_profile_config and be picked up by download()."""
    from HydroEO.project import Project

    cfg_path = tmp_path / "config.yaml"
    _write_config(
        cfg_path,
        {
            "project": {"main_dir": str(tmp_path / "out")},
            "river_profile": {
                "name": "test_river",
                "chainage_path": str(_mock_chainage_shp),
                "startdate": [2024, 1, 1],
                "enddate": [2024, 2, 1],
            },
        },
    )

    proj = Project(name="river-profile", config=str(cfg_path))
    assert hasattr(proj, "river_profile_config")
    assert not hasattr(proj, "reservoirs")
    assert not hasattr(proj, "rivers")


@pytest.mark.unit
def test_river_profile_falls_back_to_project_dates(tmp_path, _mock_chainage_shp):
    """river_profile may omit its own startdate/enddate and rely on
    project-level dates, matching the swot_raster/swot_pixc fallback."""
    from HydroEO.project import Project

    cfg_path = tmp_path / "config.yaml"
    _write_config(
        cfg_path,
        {
            "project": {
                "main_dir": str(tmp_path / "out"),
                "startdate": [2024, 1, 1],
                "enddate": [2024, 12, 31],
            },
            "river_profile": {
                "name": "test_river",
                "chainage_path": str(_mock_chainage_shp),
            },
        },
    )

    proj = Project(name="river-profile-date-fallback", config=str(cfg_path))
    assert proj.river_profile_config["startdate"] == [2024, 1, 1]
    assert proj.river_profile_config["enddate"] == [2024, 12, 31]


@pytest.mark.unit
def test_validate_config_rejects_river_profile_invalid_calendar_date(_mock_chainage_shp):
    from HydroEO.project import Project

    proj = Project.__new__(Project)
    proj.config = {
        "project": {"main_dir": "/tmp/hydroeo"},
        "river_profile": {
            "chainage_path": str(_mock_chainage_shp),
            "startdate": [2024, 13, 40],  # not a real calendar date
            "enddate": [2024, 2, 1],
        },
    }

    with pytest.raises(ValueError, match="river_profile.startdate"):
        proj.validate_config()


@pytest.mark.unit
def test_validate_config_rejects_river_profile_unsupported_product(_mock_chainage_shp):
    from HydroEO.project import Project

    proj = Project.__new__(Project)
    proj.config = {
        "project": {"main_dir": "/tmp/hydroeo"},
        "river_profile": {
            "chainage_path": str(_mock_chainage_shp),
            "startdate": [2024, 1, 1],
            "enddate": [2024, 2, 1],
            "product": "SWOT_L2_LR_SSH_2.0",
        },
    }

    with pytest.raises(ValueError, match="river_profile.product"):
        proj.validate_config()


@pytest.mark.unit
def test_validate_config_rejects_river_profile_non_numeric_orbit_exclusion_bound(
    _mock_chainage_shp,
):
    from HydroEO.project import Project

    proj = Project.__new__(Project)
    proj.config = {
        "project": {"main_dir": "/tmp/hydroeo"},
        "river_profile": {
            "chainage_path": str(_mock_chainage_shp),
            "startdate": [2024, 1, 1],
            "enddate": [2024, 2, 1],
            "orbit_exclusions": [{"orbit": "467", "max_chainage_m": "50000"}],
        },
    }

    with pytest.raises(ValueError, match="orbit_exclusions.*max_chainage_m"):
        proj.validate_config()


@pytest.mark.unit
def test_validate_config_rejects_river_profile_reversed_date_range(_mock_chainage_shp):
    from HydroEO.project import Project

    proj = Project.__new__(Project)
    proj.config = {
        "project": {"main_dir": "/tmp/hydroeo"},
        "river_profile": {
            "chainage_path": str(_mock_chainage_shp),
            "startdate": [2024, 12, 31],
            "enddate": [2024, 1, 1],
        },
    }

    with pytest.raises(ValueError, match="river_profile.startdate.*cannot be after"):
        proj.validate_config()


@pytest.mark.unit
@pytest.mark.parametrize(
    "filters_override,expected_match",
    [
        ({"soft_clamp": {"bin_width_m": 0}}, "soft_clamp.bin_width_m"),
        ({"hampel_1": {"action": "delete"}}, "hampel_1.action"),
        ({"hampel_1": {"win_m": -5}}, "hampel_1.win_m"),
        ({"rolling_quantile": {"q": 1.5}}, "rolling_quantile.q"),
        ({"density_cull": {"low_pct": -1}}, "density_cull.low_pct"),
        ({"spline_fill": {"k": 0}}, "spline_fill.k"),
        ({"not_a_real_stage": {"enabled": True}}, "not_a_real_stage"),
    ],
)
def test_validate_config_rejects_bad_river_profile_filter_values(
    _mock_chainage_shp, filters_override, expected_match
):
    from HydroEO.project import Project

    proj = Project.__new__(Project)
    proj.config = {
        "project": {"main_dir": "/tmp/hydroeo"},
        "river_profile": {
            "chainage_path": str(_mock_chainage_shp),
            "startdate": [2024, 1, 1],
            "enddate": [2024, 2, 1],
            "filters": filters_override,
        },
    }

    with pytest.raises(ValueError, match=expected_match):
        proj.validate_config()


@pytest.mark.unit
def test_validate_config_rejects_river_profile_missing_chainage_path():
    from HydroEO.project import Project

    proj = Project.__new__(Project)
    proj.config = {
        "project": {"main_dir": "/tmp/hydroeo"},
        "river_profile": {"startdate": [2024, 1, 1], "enddate": [2024, 2, 1]},
    }

    with pytest.raises(ValueError, match="river_profile.chainage_path"):
        proj.validate_config()


@pytest.mark.unit
def test_validate_config_rejects_river_profile_and_reservoirs_together(_mock_chainage_shp):
    from HydroEO.project import Project

    proj = Project.__new__(Project)
    proj.config = {
        "project": {"main_dir": "/tmp/hydroeo"},
        "reservoirs": {"path": "/tmp/res.shp", "id_key": "id"},
        "river_profile": {
            "chainage_path": str(_mock_chainage_shp),
            "startdate": [2024, 1, 1],
            "enddate": [2024, 2, 1],
        },
    }

    with pytest.raises(ValueError, match="mutually exclusive"):
        proj.validate_config()
