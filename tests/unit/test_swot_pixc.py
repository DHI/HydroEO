"""Integration tests for SWOT Pixel Cloud download, preprocess, and rasterization workflows.

Tests use synthetic netCDF/GeoJSON/GeoTIFF files and mocked earthaccess to validate
the full pipeline without network calls.
"""

from unittest.mock import patch
import shutil

import geopandas as gpd
import numpy as np
import pytest
import rasterio
import xarray as xr
from shapely.geometry import box

from HydroEO.satellites.swot.pixc import (
    download_pixc,
    _download_granules,
    _preprocess_granules,
    _rasterize_granules,
)


# ============================================================================
# FIXTURES: Synthetic Data Generation
# ============================================================================


@pytest.fixture
def synthetic_pixc_netcdf(tmp_path):
    """Create a minimal valid SWOT PIXC netCDF file with pixel_cloud group."""
    # Use a realistic SWOT PIXC filename with a valid datetime
    nc_file = (
        tmp_path
        / "SWOT_L2_HR_PIXC_034_012_053L_20250609T035604_20250609T035612_PID0_01.nc"
    )

    # Create realistic PIXC data: WGS84 lat/lon coordinates
    n_points = 500

    # Simulated geographic extent (roughly over test AOI)
    lon = np.random.uniform(24.0, 26.0, n_points)
    lat = np.random.uniform(56.0, 58.0, n_points)

    # Height values
    height = np.random.uniform(0, 100, n_points)
    geoid = np.random.uniform(40, 50, n_points)
    solid_earth_tide = np.random.uniform(-0.5, 0.5, n_points)
    load_tide = np.random.uniform(-0.1, 0.1, n_points)
    pole_tide = np.random.uniform(-0.05, 0.05, n_points)

    # Classification: mostly open_water (4) and water_near_land (3)
    classification = np.random.choice(
        [3, 4, 4, 4], size=n_points
    )  # bias toward open_water

    # Optional fields
    water_frac = np.random.uniform(0, 1, n_points)
    phase_noise_std = np.random.uniform(0, 0.5, n_points)
    dheight_dphase = np.random.uniform(-0.1, 0.1, n_points)
    sig0 = np.random.uniform(0, 30, n_points)

    # Build pixel_cloud group
    pc_data = {
        "longitude": (["num_pixels"], lon),
        "latitude": (["num_pixels"], lat),
        "height": (["num_pixels"], height),
        "geoid": (["num_pixels"], geoid),
        "solid_earth_tide": (["num_pixels"], solid_earth_tide),
        "load_tide_fes": (["num_pixels"], load_tide),
        "pole_tide": (["num_pixels"], pole_tide),
        "classification": (["num_pixels"], classification),
        "water_frac": (["num_pixels"], water_frac),
        "phase_noise_std": (["num_pixels"], phase_noise_std),
        "dheight_dphase": (["num_pixels"], dheight_dphase),
        "sig0": (["num_pixels"], sig0),
    }

    ds = xr.Dataset(pc_data)

    # Write as netCDF with pixel_cloud group
    ds.to_netcdf(nc_file, group="pixel_cloud", engine="netcdf4")

    return nc_file


@pytest.fixture
def synthetic_pixc_netcdf_empty_aoi(tmp_path):
    """Create a PIXC netCDF that doesn't overlap the test AOI."""
    # Use a realistic SWOT PIXC filename with a valid datetime
    nc_file = (
        tmp_path
        / "SWOT_L2_HR_PIXC_034_012_053L_20250610T035604_20250610T035612_PID0_01.nc"
    )

    n_points = 100
    # Far away from test bbox
    lon = np.random.uniform(10.0, 12.0, n_points)
    lat = np.random.uniform(40.0, 42.0, n_points)

    height = np.random.uniform(0, 100, n_points)
    geoid = np.random.uniform(40, 50, n_points)
    solid_earth_tide = np.random.uniform(-0.5, 0.5, n_points)
    load_tide = np.random.uniform(-0.1, 0.1, n_points)
    pole_tide = np.random.uniform(-0.05, 0.05, n_points)
    classification = np.random.choice([3, 4], size=n_points)
    water_frac = np.random.uniform(0, 1, n_points)
    phase_noise_std = np.random.uniform(0, 0.5, n_points)
    dheight_dphase = np.random.uniform(-0.1, 0.1, n_points)
    sig0 = np.random.uniform(0, 30, n_points)

    pc_data = {
        "longitude": (["num_pixels"], lon),
        "latitude": (["num_pixels"], lat),
        "height": (["num_pixels"], height),
        "geoid": (["num_pixels"], geoid),
        "solid_earth_tide": (["num_pixels"], solid_earth_tide),
        "load_tide_fes": (["num_pixels"], load_tide),
        "pole_tide": (["num_pixels"], pole_tide),
        "classification": (["num_pixels"], classification),
        "water_frac": (["num_pixels"], water_frac),
        "phase_noise_std": (["num_pixels"], phase_noise_std),
        "dheight_dphase": (["num_pixels"], dheight_dphase),
        "sig0": (["num_pixels"], sig0),
    }

    ds = xr.Dataset(pc_data)
    ds.to_netcdf(nc_file, group="pixel_cloud", engine="netcdf4")

    return nc_file


@pytest.fixture
def trimmed_geojson(tmp_path):
    """Create a synthetic trimmed GeoJSON file."""
    # Real SWOT PIXC filename format: SWOT_L2_HR_PIXC_PPP_TTT_SSS_YYYYMMDDThhmmss_YYYYMMDDThhmmss_PIDX_VV_trimmed.geojson
    # Example: SWOT_L2_HR_PIXC_034_012_053L_20250609T035604_20250609T035612_PID0_01_trimmed.geojson
    gj_file = (
        tmp_path
        / "SWOT_L2_HR_PIXC_034_012_053L_20250601T120000_20250601T120030_PID0_01_trimmed.geojson"
    )

    # Create point data
    n_points = 100
    lon = np.random.uniform(24.5, 25.5, n_points)
    lat = np.random.uniform(56.5, 57.5, n_points)

    height = np.random.uniform(0, 100, n_points)
    heightEGM = height - np.random.uniform(40, 50, n_points)

    gdf = gpd.GeoDataFrame(
        {
            "height": height,
            "heightEGM": heightEGM,
        },
        geometry=gpd.points_from_xy(lon, lat),
        crs="EPSG:4326",
    )

    gdf.to_file(str(gj_file), driver="GeoJSON")
    return gj_file


@pytest.fixture
def pixc_bbox_config():
    """SWOT PIXC config with bbox AOI."""
    return {
        "aoi": {
            "name": "test_aoi",
            "type": "bbox",
            "bbox": [24.0, 56.0, 26.0, 58.0],  # Overlaps synthetic test data
        },
        "product": "SWOT_L2_HR_PIXC_D",
        "startdate": [2025, 6, 1],
        "enddate": [2025, 6, 10],
        "classes": ["open_water", "water_near_land"],
        "fields": ["heightEGM"],
        "grid_resolution": 1000,  # 1 km for testing
    }


@pytest.fixture
def pixc_shapefile_config(tmp_path):
    """SWOT PIXC config with shapefile AOI."""
    # Create a test shapefile overlapping the synthetic data
    gdf = gpd.GeoDataFrame(
        geometry=[box(24.5, 56.5, 25.5, 57.5)],
        crs="EPSG:4326",
    )
    shp_path = tmp_path / "aoi.shp"
    gdf.to_file(shp_path)

    return {
        "aoi": {
            "name": "test_aoi_shp",
            "type": "shapefile",
            "path": str(shp_path),
        },
        "product": "SWOT_L2_HR_PIXC_D",
        "startdate": [2025, 6, 1],
        "enddate": [2025, 6, 10],
        "classes": ["open_water"],
        "fields": ["heightEGM", "height"],
        "grid_resolution": 1000,
    }


@pytest.fixture
def _fake_pixc_result():
    """Mock earthaccess result object for PIXC granules."""

    class FakeEAResult(dict):
        def __init__(self):
            super().__init__()
            self["umm"] = {"GranuleUR": "G1234567890_SWOT_L2_HR_PIXC_D"}
            self.data_links = ["https://example.com/fake_pixc.nc"]

    return FakeEAResult()


# ============================================================================
# TESTS
# ============================================================================


@pytest.mark.unit
def test_download_granules_bbox_config(tmp_path, pixc_bbox_config, _fake_pixc_result):
    """Test download phase with bbox configuration."""
    raw_dir = tmp_path / "raw"
    raw_dir.mkdir(parents=True)

    with patch("HydroEO.satellites.swot.download.earthaccess") as mock_ea:
        mock_ea.search_data.return_value = [_fake_pixc_result]
        mock_ea.download.return_value = ["file1.nc"]
        mock_ea.login.return_value = None

        downloaded = _download_granules(
            pixc_bbox_config, raw_dir, set(), ("user", "pass")
        )

        assert len(downloaded) == 1
        mock_ea.search_data.assert_called_once()


@pytest.mark.unit
def test_download_granules_skip_processed(tmp_path, pixc_bbox_config):
    """Test that already-processed granules are skipped."""
    raw_dir = tmp_path / "raw"
    raw_dir.mkdir(parents=True)

    # Mark a granule as already processed
    processed = {"G1234567890_SWOT_L2_HR_PIXC_D"}

    with patch("HydroEO.satellites.swot.download.earthaccess") as mock_ea:
        mock_ea.search_data.return_value = []
        mock_ea.login.return_value = None

        downloaded = _download_granules(
            pixc_bbox_config, raw_dir, processed, ("user", "pass")
        )

        assert len(downloaded) == 0


@pytest.mark.unit
def test_preprocess_granules_extract_and_clip(
    tmp_path, pixc_bbox_config, synthetic_pixc_netcdf
):
    """Test extraction, filtering, and AOI clipping of PIXC data."""
    raw_dir = tmp_path / "raw"
    trimmed_dir = tmp_path / "trimmed"
    raw_dir.mkdir(parents=True)
    trimmed_dir.mkdir(parents=True)

    # Copy synthetic netCDF to raw dir
    shutil.copy(synthetic_pixc_netcdf, raw_dir / synthetic_pixc_netcdf.name)

    log_path = raw_dir / "downloaded.log"
    log_path.touch()

    _preprocess_granules(pixc_bbox_config, raw_dir, trimmed_dir, log_path)

    # Check that trimmed GeoJSON was created
    trimmed_files = list(trimmed_dir.glob("*_trimmed.geojson"))
    assert len(trimmed_files) > 0

    # Verify GeoJSON contents
    gdf = gpd.read_file(trimmed_files[0])
    assert len(gdf) > 0
    assert "heightEGM" in gdf.columns
    assert all(gdf.geometry.geom_type == "Point")


@pytest.mark.unit
def test_preprocess_granules_empty_after_clipping(
    tmp_path, pixc_bbox_config, synthetic_pixc_netcdf_empty_aoi
):
    """Test that GeoJSONs with no points after clipping are not written."""
    raw_dir = tmp_path / "raw"
    trimmed_dir = tmp_path / "trimmed"
    raw_dir.mkdir(parents=True)
    trimmed_dir.mkdir(parents=True)

    shutil.copy(
        synthetic_pixc_netcdf_empty_aoi, raw_dir / synthetic_pixc_netcdf_empty_aoi.name
    )

    log_path = raw_dir / "downloaded.log"
    log_path.touch()

    _preprocess_granules(pixc_bbox_config, raw_dir, trimmed_dir, log_path)

    # No trimmed files should exist (data is outside AOI)
    trimmed_files = list(trimmed_dir.glob("*_trimmed.geojson"))
    assert len(trimmed_files) == 0


@pytest.mark.unit
def test_rasterize_granules_heightegm_field(
    tmp_path, pixc_bbox_config, trimmed_geojson
):
    """Test rasterization of heightEGM field."""
    trimmed_dir = tmp_path / "trimmed"
    raster_dir = tmp_path / "raster"
    trimmed_dir.mkdir(parents=True)
    raster_dir.mkdir(parents=True)

    shutil.copy(trimmed_geojson, trimmed_dir / trimmed_geojson.name)

    _rasterize_granules(pixc_bbox_config, trimmed_dir, raster_dir)

    # Check that GeoTIFF was created
    tif_files = list(raster_dir.glob("*.tif"))
    assert len(tif_files) > 0

    # Verify raster properties
    with rasterio.open(tif_files[0]) as src:
        assert src.count == 1
        assert src.crs is not None
        data = src.read(1)
        assert data.shape[0] > 0
        assert data.shape[1] > 0


@pytest.mark.unit
def test_rasterize_granules_auto_detect_crs(
    tmp_path, pixc_bbox_config, trimmed_geojson
):
    """Test auto-detection of target CRS."""
    trimmed_dir = tmp_path / "trimmed"
    raster_dir = tmp_path / "raster"
    trimmed_dir.mkdir(parents=True)
    raster_dir.mkdir(parents=True)

    # Remove target_crs to force auto-detection
    config = pixc_bbox_config.copy()
    if "target_crs" in config:
        del config["target_crs"]

    shutil.copy(trimmed_geojson, trimmed_dir / trimmed_geojson.name)

    _rasterize_granules(config, trimmed_dir, raster_dir)

    tif_files = list(raster_dir.glob("*.tif"))
    assert len(tif_files) > 0


@pytest.mark.unit
def test_download_pixc_end_to_end(tmp_path, pixc_bbox_config, synthetic_pixc_netcdf):
    """Test full end-to-end PIXC download pipeline."""
    raw_dir = tmp_path / "swot_pixc" / "test_aoi" / "raw" / "SWOT_L2_HR_PIXC_D"
    raw_dir.mkdir(parents=True)

    # Pre-populate with synthetic NC
    shutil.copy(synthetic_pixc_netcdf, raw_dir / synthetic_pixc_netcdf.name)

    # Mock earthaccess at the download module level where _download_granules uses it
    with patch("HydroEO.satellites.swot.download.earthaccess") as mock_ea:
        mock_ea.login.return_value = None
        mock_ea.search_data.return_value = []
        mock_ea.download.return_value = []

        download_pixc(
            pixc_bbox_config,
            str(tmp_path),
            ("user", "pass"),
        )

    # Verify outputs
    raster_dir = tmp_path / "swot_pixc" / "test_aoi" / "raster"
    tif_files = list(raster_dir.glob("*.tif"))
    assert len(tif_files) > 0


@pytest.mark.unit
def test_preprocess_with_shapefile_aoi(
    tmp_path, pixc_shapefile_config, synthetic_pixc_netcdf
):
    """Test preprocessing with shapefile-based AOI."""
    raw_dir = tmp_path / "raw"
    trimmed_dir = tmp_path / "trimmed"
    raw_dir.mkdir(parents=True)
    trimmed_dir.mkdir(parents=True)

    shutil.copy(synthetic_pixc_netcdf, raw_dir / synthetic_pixc_netcdf.name)

    log_path = raw_dir / "downloaded.log"
    log_path.touch()

    _preprocess_granules(pixc_shapefile_config, raw_dir, trimmed_dir, log_path)

    trimmed_files = list(trimmed_dir.glob("*_trimmed.geojson"))
    assert len(trimmed_files) > 0

    gdf = gpd.read_file(trimmed_files[0])
    assert len(gdf) > 0


@pytest.mark.unit
def test_rasterize_multiple_fields(tmp_path, pixc_bbox_config, trimmed_geojson):
    """Test rasterization with multiple fields."""
    trimmed_dir = tmp_path / "trimmed"
    raster_dir = tmp_path / "raster"
    trimmed_dir.mkdir(parents=True)
    raster_dir.mkdir(parents=True)

    config = pixc_bbox_config.copy()
    config["fields"] = ["heightEGM", "height"]

    shutil.copy(trimmed_geojson, trimmed_dir / trimmed_geojson.name)

    _rasterize_granules(config, trimmed_dir, raster_dir)

    tif_files = sorted(raster_dir.glob("*.tif"))
    # Should have one date with two fields
    assert len(tif_files) >= 1


@pytest.mark.unit
def test_download_granules_no_results(tmp_path, pixc_bbox_config):
    """Test handling when no granules are found."""
    raw_dir = tmp_path / "raw"
    raw_dir.mkdir(parents=True)

    with patch("HydroEO.satellites.swot.download.earthaccess") as mock_ea:
        mock_ea.search_data.return_value = []
        mock_ea.login.return_value = None

        downloaded = _download_granules(
            pixc_bbox_config, raw_dir, set(), ("user", "pass")
        )

        assert len(downloaded) == 0
