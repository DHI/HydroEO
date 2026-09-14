"""Tests for the river_profile orchestration pipeline (mocked download)."""

from pathlib import Path
from unittest.mock import patch

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import Point

from HydroEO.satellites.swot.river_profile import (
    calculate_river_profile,
    _apply_orbit_exclusions,
    _load_chainage,
    _resolve_filters,
)

pytestmark = pytest.mark.unit

UTM_CRS = "EPSG:32645"
PRODUCT = "SWOT_L2_HR_Raster_D"
NAME = "test_river"


def _make_chainage_shp(tmp_path, n=25, spacing_m=200.0, reverse_field=False) -> Path:
    xs = np.full(n, 500_000.0)
    ys = 2_600_000.0 + np.arange(n) * spacing_m
    dist = np.arange(n) * spacing_m
    if reverse_field:
        dist = dist[::-1]
    gdf = gpd.GeoDataFrame(
        {"cngmeters": dist, "geometry": [Point(x, y) for x, y in zip(xs, ys)]},
        crs=UTM_CRS,
    )
    path = tmp_path / "chainage.shp"
    gdf.to_file(path)
    return path


def _make_wse_tif(path: Path, crs=UTM_CRS):
    path.parent.mkdir(parents=True, exist_ok=True)
    width, height = 40, 400
    transform = from_origin(499_900.0, 2_605_500.0, 10.0, 10.0)
    rows, cols = np.meshgrid(np.arange(height), np.arange(width), indexing="ij")
    data = (10.0 + 0.0001 * rows).astype(np.float32)
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=height,
        width=width,
        count=1,
        dtype="float32",
        crs=crs,
        transform=transform,
        nodata=np.nan,
    ) as dst:
        dst.write(data, 1)


def _base_config(tmp_path, chainage_path, **overrides):
    cfg = {
        "name": NAME,
        "chainage_path": str(chainage_path),
        "chainage_field": "cngmeters",
        "startdate": [2024, 1, 1],
        "enddate": [2024, 1, 31],
        "plot_enable": False,
    }
    cfg.update(overrides)
    return cfg


@pytest.fixture
def project_with_wse_tile(tmp_path):
    chainage_path = _make_chainage_shp(tmp_path)
    processed_dir = tmp_path / "processed" / "swot_raster" / NAME / PRODUCT
    wse_path = processed_dir / "SWOT_L2_HR_Raster_100m_UTM45N_x_x_x_001_001_20240115T000000_20240115T000020_PGC0_01_wse.tif"
    _make_wse_tif(wse_path)
    return tmp_path, chainage_path


def test_calculate_river_profile_writes_final_outputs_only_by_default(project_with_wse_tile):
    project_dir, chainage_path = project_with_wse_tile
    config = _base_config(project_dir, chainage_path)

    with patch("HydroEO.satellites.swot.river_profile.download_raster") as mock_download:
        calculate_river_profile(config, project_dir=str(project_dir), global_crs="EPSG:4326")

    assert mock_download.called
    call_kwargs = mock_download.call_args.kwargs
    assert call_kwargs["config"]["merge_tiles"] is False
    assert call_kwargs["config"]["aoi"]["type"] == "bbox"

    results_dir = project_dir / "results" / NAME
    final_shps = list((results_dir / "profiles_final").glob("*.shp"))
    assert len(final_shps) == 1

    csv_path = results_dir / f"{NAME}_profiles_final.csv"
    assert csv_path.exists()
    df = pd.read_csv(csv_path)
    assert "distance_along_river_m" in df.columns
    assert len(df) == 25

    assert (results_dir / "quality_report.csv").exists()

    # keep_intermediates defaults to False: no intermediate dirs written
    assert not (results_dir / "profiles_raw").exists()
    assert not (results_dir / "profiles_prefilter").exists()
    assert not (results_dir / "profiles_hampel1").exists()


def test_calculate_river_profile_keeps_intermediates_when_requested(project_with_wse_tile):
    project_dir, chainage_path = project_with_wse_tile
    config = _base_config(project_dir, chainage_path, keep_intermediates=True)

    with patch("HydroEO.satellites.swot.river_profile.download_raster"):
        calculate_river_profile(config, project_dir=str(project_dir), global_crs="EPSG:4326")

    results_dir = project_dir / "results" / NAME
    assert list((results_dir / "profiles_raw").glob("*.shp"))
    assert list((results_dir / "profiles_prefilter").glob("*.shp"))
    assert list((results_dir / "profiles_hampel1").glob("*.shp"))
    assert (results_dir / f"{NAME}_profiles_raw.csv").exists()


def test_calculate_river_profile_no_wse_tiles_is_a_noop(tmp_path):
    chainage_path = _make_chainage_shp(tmp_path)
    config = _base_config(tmp_path, chainage_path)

    with patch("HydroEO.satellites.swot.river_profile.download_raster"):
        calculate_river_profile(config, project_dir=str(tmp_path), global_crs="EPSG:4326")

    results_dir = tmp_path / "results" / NAME
    assert not list((results_dir / "profiles_final").glob("*.shp"))


def test_load_chainage_requires_configured_field(tmp_path):
    chainage_path = _make_chainage_shp(tmp_path)
    config = {"chainage_path": str(chainage_path), "chainage_field": "not_a_column"}

    with pytest.raises(ValueError, match="not_a_column"):
        _load_chainage(config)


def test_load_chainage_reverses_when_requested(tmp_path):
    chainage_path = _make_chainage_shp(tmp_path)
    config = {"chainage_path": str(chainage_path), "reverse_chainage": True}

    gdf, x = _load_chainage(config)

    assert x.min() == pytest.approx(0.0)
    # the point originally at distance 0 should now carry the max distance
    assert gdf.iloc[-1]["cngmeters"] == pytest.approx(0.0)


def test_apply_orbit_exclusions_masks_matching_range():
    x = np.array([0.0, 10_000.0, 60_000.0, 80_000.0])
    y = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float32)
    rules = [{"orbit": "467", "max_chainage_m": 50_000.0}]

    out, n_excluded = _apply_orbit_exclusions("SWOT_..._467_..._wse.tif", x, y, rules)

    assert n_excluded == 2
    assert np.isnan(out[0]) and np.isnan(out[1])
    assert out[2] == pytest.approx(3.0) and out[3] == pytest.approx(4.0)


def test_apply_orbit_exclusions_ignores_non_matching_orbit():
    x = np.array([0.0, 10_000.0])
    y = np.array([1.0, 2.0], dtype=np.float32)
    rules = [{"orbit": "467", "max_chainage_m": 50_000.0}]

    out, n_excluded = _apply_orbit_exclusions("SWOT_..._230_..._wse.tif", x, y, rules)

    assert n_excluded == 0
    np.testing.assert_array_equal(out, y)


def test_resolve_filters_merges_user_overrides_over_defaults():
    resolved = _resolve_filters({"hampel_1": {"enabled": False}, "unknown_stage": {"x": 1}})

    assert resolved["hampel_1"]["enabled"] is False
    assert resolved["hampel_1"]["win_m"] == 20_000.0  # untouched default preserved
    assert resolved["soft_clamp"]["enabled"] is True
    assert "unknown_stage" not in resolved
