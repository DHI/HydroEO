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
    _build_swot_raster_config,
    _load_chainage,
    _plot_combined,
    _process_one_date,
    _resolve_filters,
    _unique_label,
    _COMBINED_PLOT_MAX_LEGEND_ENTRIES,
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


def test_load_chainage_rejects_non_finite_values(tmp_path):
    xs = np.full(5, 500_000.0)
    ys = 2_600_000.0 + np.arange(5) * 200.0
    dist = np.array([0.0, 200.0, np.nan, 600.0, 800.0])
    gdf = gpd.GeoDataFrame(
        {"cngmeters": dist, "geometry": [Point(x, y) for x, y in zip(xs, ys)]},
        crs=UTM_CRS,
    )
    path = tmp_path / "chainage_with_nan.shp"
    gdf.to_file(path)

    with pytest.raises(ValueError, match="non-finite"):
        _load_chainage({"chainage_path": str(path), "chainage_field": "cngmeters"})


def test_load_chainage_reverses_when_requested(tmp_path):
    chainage_path = _make_chainage_shp(tmp_path)
    config = {"chainage_path": str(chainage_path), "reverse_chainage": True}

    gdf, x = _load_chainage(config)

    assert x.min() == pytest.approx(0.0)
    # the gdf's chainage column must be kept consistent with the reversed
    # working distance array, not left holding the pre-reversal values
    np.testing.assert_allclose(gdf["cngmeters"].to_numpy(), x)
    # the point originally at distance 0 (now sorted last) carries the max
    # reversed distance, not its original raw value
    assert gdf.iloc[-1]["cngmeters"] == pytest.approx(x.max())


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


def test_plot_combined_drops_legend_above_threshold(tmp_path):
    x = np.linspace(0, 1000, 10)
    n_dates = _COMBINED_PLOT_MAX_LEGEND_ENTRIES + 5
    series = {f"date_{i}": x for i in range(n_dates)}
    out_path = tmp_path / "combined.png"

    with patch("matplotlib.pyplot.legend") as mock_legend:
        _plot_combined(x, series, "title", out_path, ylim=None, dpi=72)

    mock_legend.assert_not_called()
    assert out_path.exists()


def test_plot_combined_keeps_legend_below_threshold(tmp_path):
    x = np.linspace(0, 1000, 10)
    n_dates = _COMBINED_PLOT_MAX_LEGEND_ENTRIES - 5
    series = {f"date_{i}": x for i in range(n_dates)}
    out_path = tmp_path / "combined.png"

    with patch("matplotlib.pyplot.legend") as mock_legend:
        _plot_combined(x, series, "title", out_path, ylim=None, dpi=72)

    mock_legend.assert_called_once()
    assert out_path.exists()


def test_calculate_river_profile_final_output_uses_y_final_not_y_spline(
    project_with_wse_tile,
):
    """Regression test: the emitted final CSV/shapefile must contain
    y_final (measured values, gaps spline-filled) not y_spline (the smooth
    curve evaluated everywhere, which silently overwrites real
    observations) - see docs/river_profile.md's 'never overwrites real
    data' contract for spline_fill."""
    project_dir, chainage_path = project_with_wse_tile
    config = _base_config(project_dir, chainage_path)
    label = "20240115T000000"
    n = 25
    fake_record = {
        "label": label,
        "y_raw": np.zeros(n, dtype=np.float32),
        "y_pref": np.zeros(n, dtype=np.float32),
        "y_h1": np.zeros(n, dtype=np.float32),
        "y_rq": np.zeros(n, dtype=np.float32),
        "y_after_density": np.zeros(n, dtype=np.float32),
        "y_h2": np.zeros(n, dtype=np.float32),
        "y_final": np.full(n, 1.0, dtype=np.float32),
        "y_spline": np.full(n, 999.0, dtype=np.float32),
        "quality": {
            "label": label,
            "file": "fake_wse.tif",
            "n_finite_raw": n,
            "n_orbit_excluded": 0,
            "n_soft_clamp_changed": 0,
            "n_hampel1_flagged": 0,
            "n_density_dropped": 0,
            "n_hampel2_flagged": 0,
            "n_finite_final": n,
        },
    }

    with (
        patch("HydroEO.satellites.swot.river_profile.download_raster"),
        patch(
            "HydroEO.satellites.swot.river_profile._process_one_date",
            return_value=fake_record,
        ),
    ):
        calculate_river_profile(config, project_dir=str(project_dir), global_crs="EPSG:4326")

    results_dir = project_dir / "results" / NAME
    df = pd.read_csv(results_dir / f"{NAME}_profiles_final.csv")
    values = df[label].to_numpy()
    assert np.allclose(values, 1.0)
    assert not np.allclose(values, 999.0)


def test_unique_label_disambiguates_collisions():
    seen: dict[str, int] = {}
    assert _unique_label("20240101T000000", seen) == "20240101T000000"
    assert _unique_label("20240101T000000", seen) == "20240101T000000_1"
    assert _unique_label("20240101T000000", seen) == "20240101T000000_2"
    assert _unique_label("other", seen) == "other"


def test_calculate_river_profile_disambiguates_same_timestamp_tiles(project_with_wse_tile):
    """Two granules that share the same extracted timestamp label (e.g. two
    tiles from the same pass split across a UTM zone boundary) must not
    silently overwrite each other's output."""
    project_dir, chainage_path = project_with_wse_tile
    processed_dir = project_dir / "processed" / "swot_raster" / NAME / PRODUCT
    second_path = (
        processed_dir
        / "SWOT_L2_HR_Raster_100m_UTM46N_x_x_x_001_002_20240115T000000_20240115T000020_PGD0_02_wse.tif"
    )
    _make_wse_tif(second_path)
    config = _base_config(project_dir, chainage_path)

    with patch("HydroEO.satellites.swot.river_profile.download_raster"):
        calculate_river_profile(config, project_dir=str(project_dir), global_crs="EPSG:4326")

    results_dir = project_dir / "results" / NAME
    final_shps = list((results_dir / "profiles_final").glob("*.shp"))
    assert len(final_shps) == 2

    df = pd.read_csv(results_dir / f"{NAME}_profiles_final.csv")
    assert "20240115T000000" in df.columns
    assert "20240115T000000_1" in df.columns

    # the quality report's 'label' column must match the disambiguated
    # labels used for the CSV columns/shapefiles above, not the raw
    # (colliding) timestamp extracted from the filename
    quality_df = pd.read_csv(results_dir / "quality_report.csv")
    assert sorted(quality_df["label"]) == ["20240115T000000", "20240115T000000_1"]


def test_calculate_river_profile_saves_geoid_regardless_of_keep_intermediates(
    project_with_wse_tile,
):
    """geoid is an explicitly requested output variable, not an
    intermediate filtering stage - it must be saved even with the default
    keep_intermediates: false."""
    project_dir, chainage_path = project_with_wse_tile
    processed_dir = project_dir / "processed" / "swot_raster" / NAME / PRODUCT
    geoid_path = processed_dir / (
        "SWOT_L2_HR_Raster_100m_UTM45N_x_x_x_001_001_"
        "20240115T000000_20240115T000020_PGC0_01_geoid.tif"
    )
    _make_wse_tif(geoid_path)
    config = _base_config(project_dir, chainage_path, variables=["wse", "geoid"])
    assert "keep_intermediates" not in config  # defaults to False

    with patch("HydroEO.satellites.swot.river_profile.download_raster"):
        calculate_river_profile(config, project_dir=str(project_dir), global_crs="EPSG:4326")

    results_dir = project_dir / "results" / NAME
    assert list((results_dir / "profiles_geoid").glob("*.shp"))


def test_build_swot_raster_config_handles_zero_buffer(tmp_path):
    """Regression test: buffering point geometries by 0 produces empty
    geometries with NaN bounds; a zero buffer must fall back to the raw
    extent instead of producing a NaN bbox."""
    chainage_path = _make_chainage_shp(tmp_path)
    gdf, _ = _load_chainage({"chainage_path": str(chainage_path), "chainage_field": "cngmeters"})
    config = {
        "startdate": [2024, 1, 1],
        "enddate": [2024, 1, 31],
        "aoi_buffer_meters": 0,
    }

    swot_cfg = _build_swot_raster_config(config, NAME, gdf, "EPSG:4326")

    assert not any(np.isnan(v) for v in swot_cfg["aoi"]["bbox"])


def test_calculate_river_profile_records_fully_orbit_excluded_date(project_with_wse_tile):
    """Regression test: if orbit_exclusions masks every sampled point for a
    date, that date must still get a quality_report row (recording the
    exclusion) instead of silently vanishing from the report."""
    project_dir, chainage_path = project_with_wse_tile
    config = _base_config(
        project_dir,
        chainage_path,
        orbit_exclusions=[{"orbit": "PGC0", "max_chainage_m": 100_000.0}],
    )

    with patch("HydroEO.satellites.swot.river_profile.download_raster"):
        calculate_river_profile(config, project_dir=str(project_dir), global_crs="EPSG:4326")

    results_dir = project_dir / "results" / NAME
    quality_df = pd.read_csv(results_dir / "quality_report.csv")
    assert len(quality_df) == 1
    row = quality_df.iloc[0]
    assert row["n_orbit_excluded"] == row["n_finite_raw"]
    assert row["n_finite_final"] == 0


def test_calculate_river_profile_clears_stale_output_between_runs(project_with_wse_tile):
    """results/<name> is a regenerated report - a rerun with
    keep_intermediates toggled off must not leave the previous run's
    intermediate directories behind."""
    project_dir, chainage_path = project_with_wse_tile
    results_dir = project_dir / "results" / NAME

    with patch("HydroEO.satellites.swot.river_profile.download_raster"):
        calculate_river_profile(
            _base_config(project_dir, chainage_path, keep_intermediates=True),
            project_dir=str(project_dir),
            global_crs="EPSG:4326",
        )
    assert (results_dir / "profiles_raw").exists()

    with patch("HydroEO.satellites.swot.river_profile.download_raster"):
        calculate_river_profile(
            _base_config(project_dir, chainage_path),
            project_dir=str(project_dir),
            global_crs="EPSG:4326",
        )

    assert not (results_dir / "profiles_raw").exists()
    assert list((results_dir / "profiles_final").glob("*.shp"))


def test_density_cull_does_not_drop_uniform_density_profile(project_with_wse_tile):
    """Regression test: when every point has identical local density (e.g. a
    small/short profile whose density_cull window is narrower than the point
    spacing), the percentile threshold equals every count, and a non-strict
    <= comparison used to mark the entire profile as low-density."""
    project_dir, chainage_path = project_with_wse_tile
    gdf, x = _load_chainage({"chainage_path": str(chainage_path), "chainage_field": "cngmeters"})
    wse_path = next(
        (project_dir / "processed" / "swot_raster" / NAME / PRODUCT).glob("*_wse.tif")
    )
    filters = _resolve_filters(
        {"density_cull": {"total_win_m": 50.0, "abs_min": 0}}  # narrower than 200m spacing
    )

    record = _process_one_date(
        wse_path=wse_path,
        gdf=gdf,
        x=x,
        filters=filters,
        orbit_exclusions=[],
        zero_is_nodata=True,
    )

    assert record is not None
    assert np.isfinite(record["y_after_density"]).any()
