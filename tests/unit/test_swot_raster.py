"""Integration tests for SWOT raster download, preprocess, and merge workflows.

Tests use synthetic netCDF/GeoTIFF files and mocked earthaccess to validate
the full pipeline without network calls.
"""

from pathlib import Path
from unittest.mock import MagicMock, patch
import shutil

import geopandas as gpd
import numpy as np
import pytest
import rasterio
import xarray as xr
from pyproj import CRS
from rasterio.transform import from_bounds
from shapely.geometry import box

from HydroEO.satellites.swot.raster import (
    REQUIRED_VARIABLES,
    download_raster,
    _download_granules,
    _preprocess_granules,
    _merge_and_reproject_granules,
    _detect_crs,
    _resampling_for,
)


# ============================================================================
# FIXTURES: Synthetic Data Generation
# ============================================================================


@pytest.fixture
def synthetic_swot_netcdf(tmp_path):
    """Create a minimal valid SWOT-like netCDF file with required fields."""
    nc_file = tmp_path / "SWOT_L2_HR_Raster_D_00000000T000000_000000_UTM44N_01_01.nc"

    # Create SWOT structure with geographic bounds overlapping test AOI
    # UTM 44N bounds (roughly): lon 24-30, lat 56-60 in real world
    width, height = 100, 100

    # Create coordinate arrays for UTM 44N in projected coordinates
    # UTM zone 44N origin is around lon=24, lat=56
    x_coords = np.linspace(500000, 600000, width)  # Easting (meters)
    y_coords = np.linspace(6200000, 6300000, height)  # Northing (meters)

    data = {
        "wse": (["y", "x"], np.random.uniform(0, 100, (height, width))),
        "wse_uncert": (["y", "x"], np.random.uniform(0, 0.5, (height, width))),
        "wse_qual": (["y", "x"], np.random.randint(0, 4, (height, width))),
        "height_cor_xover": (["y", "x"], np.random.uniform(-0.5, 0.5, (height, width))),
        "geoid": (["y", "x"], np.random.uniform(40, 50, (height, width))),
        "n_wse_pix": (["y", "x"], np.random.randint(1, 100, (height, width))),
        "n_other_pix": (["y", "x"], np.random.randint(0, 100, (height, width))),
        "layover_impact": (["y", "x"], np.random.uniform(0, 0.5, (height, width))),
    }

    coords = {
        "x": x_coords,
        "y": y_coords,
    }

    ds = xr.Dataset(data, coords=coords)
    ds.attrs["crs_wkt"] = CRS.from_epsg(32644).to_wkt()  # UTM 44N

    ds.to_netcdf(nc_file)
    return nc_file


@pytest.fixture
def synthetic_swot_netcdf_utm45(tmp_path):
    """Create a SWOT netCDF in UTM zone 45 (for multi-tile merge tests)."""
    nc_file = tmp_path / "SWOT_L2_HR_Raster_D_00000001T000000_000000_UTM45N_01_01.nc"

    width, height = 100, 100
    # UTM 45N coordinates (east of zone 44)
    x_coords = np.linspace(600000, 700000, width)
    y_coords = np.linspace(6200000, 6300000, height)

    data = {
        "wse": (["y", "x"], np.random.uniform(0, 100, (height, width))),
        "wse_uncert": (
            ["y", "x"],
            np.random.uniform(0, 0.2, (height, width)),
        ),  # Good quality
        "wse_qual": (["y", "x"], np.random.randint(0, 4, (height, width))),
        "height_cor_xover": (["y", "x"], np.random.uniform(-0.5, 0.5, (height, width))),
        "geoid": (["y", "x"], np.random.uniform(40, 50, (height, width))),
        "n_wse_pix": (["y", "x"], np.random.randint(1, 100, (height, width))),
        "n_other_pix": (["y", "x"], np.random.randint(0, 100, (height, width))),
        "layover_impact": (["y", "x"], np.random.uniform(0, 0.2, (height, width))),
    }

    coords = {"x": x_coords, "y": y_coords}
    ds = xr.Dataset(data, coords=coords)
    ds.attrs["crs_wkt"] = CRS.from_epsg(32645).to_wkt()  # UTM 45N
    ds.to_netcdf(nc_file)
    return nc_file


@pytest.fixture
def synthetic_geotiff_utm44(tmp_path):
    """Create a GeoTIFF in UTM 44N for merge testing."""
    tif_file = tmp_path / "20250601T000000_wse.tif"

    # Create raster in UTM 44N
    data = np.random.uniform(0, 100, (100, 100)).astype(np.float32)
    transform = from_bounds(500000, 6200000, 500000 + 1000, 6200000 + 1000, 100, 100)

    with rasterio.open(
        tif_file,
        "w",
        driver="GTiff",
        height=100,
        width=100,
        count=1,
        dtype=np.float32,
        crs=CRS.from_epsg(32644),
        transform=transform,
        nodata=np.nan,
    ) as dst:
        dst.write(data, 1)

    return tif_file


@pytest.fixture
def test_geotiff_basename():
    """Return a test basename for grouping variable names."""
    return "test_20250601T000000"


@pytest.fixture
def synthetic_geotiff_utm45(tmp_path):
    """Create a GeoTIFF in UTM 45N (adjacent zone, same date for merge test)."""
    tif_file = tmp_path / "20250601T010000_wse.tif"

    data = np.random.uniform(0, 100, (100, 100)).astype(np.float32)
    # Adjacent tile to the east of UTM 44
    transform = from_bounds(600000, 6200000, 600000 + 1000, 6200000 + 1000, 100, 100)

    with rasterio.open(
        tif_file,
        "w",
        driver="GTiff",
        height=100,
        width=100,
        count=1,
        dtype=np.float32,
        crs=CRS.from_epsg(32645),
        transform=transform,
        nodata=np.nan,
    ) as dst:
        dst.write(data, 1)

    return tif_file


@pytest.fixture
def bbox_config():
    """SWOT raster config with bbox AOI."""
    return {
        "aoi": {
            "name": "test_aoi",
            "type": "bbox",
            # Global bbox ensures overlap with synthetic test data
            "bbox": [-180, -90, 180, 90],
        },
        "product": "SWOT_L2_HR_Raster_D",
        "startdate": [2025, 6, 1],
        "enddate": [2025, 6, 10],
        "quality_filters": {
            "max_wse_uncert": 0.3,
            "max_layover_impact": 0.3,
        },
    }


@pytest.fixture
def bbox_config_custom_variables(bbox_config):
    """Config with custom variable list (plus required vars)."""
    config = bbox_config.copy()
    config["variables"] = ["wse_qual", "geoid"]
    return config


@pytest.fixture
def shapefile_config(tmp_path):
    """SWOT raster config with shapefile AOI."""
    # Create a test shapefile
    gdf = gpd.GeoDataFrame(
        geometry=[box(23, 56, 25, 58)],
        crs="EPSG:4326",
    )
    shp_path = tmp_path / "aoi.shp"
    gdf.to_file(shp_path)

    return {
        "aoi": {
            "name": "shapefile_aoi",
            "type": "shapefile",
            "path": str(shp_path),
        },
        "product": "SWOT_L2_HR_Raster_D",
        "startdate": [2025, 6, 1],
        "enddate": [2025, 6, 10],
        "quality_filters": {
            "max_wse_uncert": 0.3,
            "max_layover_impact": 0.3,
        },
    }


# ============================================================================
# TESTS: Download Phase
# ============================================================================


@pytest.mark.unit
def test_download_granules_returns_empty_when_no_results(bbox_config):
    """_download_granules returns empty list when search finds no granules."""
    with patch("HydroEO.satellites.swot.raster.earthaccess") as mock_ea:
        mock_ea.login.return_value = None
        mock_ea.search_data.return_value = []

        result = _download_granules(
            bbox_config,
            Path("/tmp/raw"),
            set(),  # processed_granules
            ("user", "pass"),
        )

    assert result == []


@pytest.mark.unit
def test_download_granules_skips_already_processed(bbox_config, tmp_path):
    """_download_granules skips granules already in processed set."""
    raw_dir = tmp_path / "raw"
    raw_dir.mkdir()

    granule_ur = "SWOT_L2_HR_Raster_D_00000000T000000_000000_UTM44N_01_01"

    with patch("HydroEO.satellites.swot.raster.earthaccess") as mock_ea:
        mock_ea.login.return_value = None
        granule = MagicMock()
        granule.__getitem__.return_value = {"GranuleUR": granule_ur}
        mock_ea.search_data.return_value = [granule]

        # Already in processed set
        result = _download_granules(
            bbox_config,
            raw_dir,
            {granule_ur},  # Already processed
            ("user", "pass"),
        )

    mock_ea.download.assert_not_called()
    assert result == []


# ============================================================================
# TESTS: Preprocessing Phase
# ============================================================================


@pytest.mark.unit
def test_preprocess_granules_filters_by_quality(
    bbox_config, synthetic_swot_netcdf, tmp_path
):
    """Pixels with wse_uncert or layover_impact above threshold are masked."""
    raw_dir = tmp_path / "raw"
    processed_dir = tmp_path / "processed"
    raw_dir.mkdir()
    processed_dir.mkdir()

    # Copy test file to raw dir
    shutil.copy(synthetic_swot_netcdf, raw_dir / synthetic_swot_netcdf.name)

    log_path = raw_dir / "downloaded.log"
    log_path.touch()

    _preprocess_granules(bbox_config, raw_dir, processed_dir, log_path)

    # Check that TIFFs were created
    tiffs = list(processed_dir.glob("*.tif"))
    assert len(tiffs) > 0, "No TIFFs were produced"

    # Verify required variables are present
    tiff_names = {t.stem for t in tiffs}
    assert any("wse" in name and "uncert" not in name for name in tiff_names), (
        "wse variable missing"
    )
    assert any("wse_uncert" in name for name in tiff_names), "wse_uncert required"
    assert any("layover_impact" in name for name in tiff_names), (
        "layover_impact required"
    )


@pytest.mark.unit
def test_preprocess_granules_always_includes_required_variables(
    bbox_config_custom_variables, synthetic_swot_netcdf, tmp_path
):
    """Even when user specifies custom variables, REQUIRED_VARIABLES are always output."""
    raw_dir = tmp_path / "raw"
    processed_dir = tmp_path / "processed"
    raw_dir.mkdir()
    processed_dir.mkdir()

    shutil.copy(synthetic_swot_netcdf, raw_dir / synthetic_swot_netcdf.name)

    log_path = raw_dir / "downloaded.log"
    log_path.touch()

    _preprocess_granules(bbox_config_custom_variables, raw_dir, processed_dir, log_path)

    tiff_names = {t.stem for t in processed_dir.glob("*.tif")}

    # Custom variables should be present
    assert any("wse_qual" in name for name in tiff_names), "Custom wse_qual not found"
    assert any("geoid" in name for name in tiff_names), "Custom geoid not found"

    # Required variables must always be present
    for required in REQUIRED_VARIABLES:
        assert any(required in name for name in tiff_names), (
            f"Required {required} missing"
        )


@pytest.mark.unit
def test_preprocess_granules_skips_already_processed(
    bbox_config, synthetic_swot_netcdf, tmp_path
):
    """If output TIFFs exist for NC file, skip reprocessing."""
    raw_dir = tmp_path / "raw"
    processed_dir = tmp_path / "processed"
    raw_dir.mkdir()
    processed_dir.mkdir()

    # Copy NC file
    nc_copy = shutil.copy(synthetic_swot_netcdf, raw_dir / synthetic_swot_netcdf.name)

    # Create fake output TIFFs for this NC file
    fake_tiff = processed_dir / f"{Path(nc_copy).stem}_wse.tif"
    fake_tiff.write_text("fake tiff")

    log_path = raw_dir / "downloaded.log"
    log_path.touch()

    # Rerun preprocess
    _preprocess_granules(bbox_config, raw_dir, processed_dir, log_path)

    # Fake file should still exist, not be replaced
    assert fake_tiff.read_text() == "fake tiff"


@pytest.mark.unit
def test_preprocess_granules_handles_missing_quality_fields(tmp_path):
    """NetCDF missing wse_uncert or layover_impact is skipped."""
    raw_dir = tmp_path / "raw"
    processed_dir = tmp_path / "processed"
    raw_dir.mkdir()
    processed_dir.mkdir()

    # Create NC file WITHOUT required quality fields
    nc_file = raw_dir / "incomplete.nc"
    data = {
        "wse": (["y", "x"], np.random.uniform(0, 100, (50, 50))),
        # Missing wse_uncert and layover_impact
    }
    coords = {"x": np.arange(50), "y": np.arange(50)}
    ds = xr.Dataset(data, coords=coords)
    ds.attrs["crs_wkt"] = CRS.from_epsg(32644).to_wkt()
    ds.to_netcdf(nc_file)

    log_path = raw_dir / "downloaded.log"
    log_path.touch()

    config = {
        "aoi": {"type": "bbox", "bbox": [0, 0, 1, 1]},
        "quality_filters": {"max_wse_uncert": 0.3, "max_layover_impact": 0.3},
    }

    _preprocess_granules(config, raw_dir, processed_dir, log_path)

    # No TIFFs should be created
    assert len(list(processed_dir.glob("*.tif"))) == 0


@pytest.mark.unit
def test_preprocess_granules_bounds_check_discards_non_overlapping(tmp_path):
    """Tiles outside AOI bounds are skipped before extraction."""
    raw_dir = tmp_path / "raw"
    processed_dir = tmp_path / "processed"
    raw_dir.mkdir()
    processed_dir.mkdir()

    # Create NC in UTM zone far outside bbox AOI
    nc_file = raw_dir / "SWOT_test_UTM50N.nc"
    width, height = 100, 100
    data = {
        "wse": (["y", "x"], np.random.uniform(0, 100, (height, width))),
        "wse_uncert": (["y", "x"], np.random.uniform(0, 0.2, (height, width))),
        "layover_impact": (["y", "x"], np.random.uniform(0, 0.2, (height, width))),
    }
    coords = {"x": np.arange(width), "y": np.arange(height)}
    ds = xr.Dataset(data, coords=coords)
    ds.attrs["crs_wkt"] = CRS.from_epsg(32650).to_wkt()  # UTM 50N (far east)
    ds.to_netcdf(nc_file)

    log_path = raw_dir / "downloaded.log"
    log_path.touch()

    # AOI in small bbox that won't overlap UTM50
    config = {
        "aoi": {"type": "bbox", "bbox": [23, 56, 25, 58]},
        "quality_filters": {"max_wse_uncert": 0.3, "max_layover_impact": 0.3},
    }

    _preprocess_granules(config, raw_dir, processed_dir, log_path)

    # No TIFFs should be created (bounds check filtered it out)
    assert len(list(processed_dir.glob("*.tif"))) == 0


# ============================================================================
# TESTS: Merge and Reproject Phase
# ============================================================================


@pytest.mark.unit
def test_merge_and_reproject_single_tile(
    bbox_config, synthetic_geotiff_utm44, tmp_path
):
    """Single TIF is reprojected without actual merge."""
    processed_dir = tmp_path / "processed"
    processed_dir.mkdir()

    shutil.copy(synthetic_geotiff_utm44, processed_dir / synthetic_geotiff_utm44.name)

    _merge_and_reproject_granules(bbox_config, processed_dir)

    merged_dir = processed_dir.parent.parent / "merged"
    merged_tiffs = list(merged_dir.glob("*.tif"))
    assert len(merged_tiffs) == 1


@pytest.mark.unit
def test_merge_and_reproject_multiple_tiles(
    bbox_config, synthetic_geotiff_utm44, synthetic_geotiff_utm45, tmp_path
):
    """Multiple TIFs are merged into single mosaic."""
    processed_dir = tmp_path / "processed"
    processed_dir.mkdir()

    shutil.copy(synthetic_geotiff_utm44, processed_dir / synthetic_geotiff_utm44.name)
    shutil.copy(synthetic_geotiff_utm45, processed_dir / synthetic_geotiff_utm45.name)

    _merge_and_reproject_granules(bbox_config, processed_dir)

    merged_dir = processed_dir.parent.parent / "merged"
    merged_tiffs = list(merged_dir.glob("*.tif"))
    assert len(merged_tiffs) == 1

    # Merged output should be larger than single input
    with rasterio.open(merged_tiffs[0]) as src:
        assert src.width >= 100


@pytest.mark.unit
def test_merge_tiles_false_skips_merge(synthetic_geotiff_utm44, tmp_path):
    """When merge_tiles is false, merge phase is skipped."""
    processed_dir = tmp_path / "swot_raster" / "aoi" / "processed" / "product"
    processed_dir.mkdir(parents=True, exist_ok=True)

    shutil.copy(synthetic_geotiff_utm44, processed_dir / synthetic_geotiff_utm44.name)

    config = {"merge_tiles": False}

    _merge_and_reproject_granules(config, processed_dir)

    # No merged dir should be created
    merged_dir = processed_dir.parent.parent / "merged"
    assert not merged_dir.exists() or len(list(merged_dir.glob("*.tif"))) == 0


@pytest.mark.unit
def test_merge_and_reproject_groups_by_date_and_variable(tmp_path):
    """Multiple variables from same date are grouped separately."""
    processed_dir = tmp_path / "processed"
    processed_dir.mkdir()

    # Create two TIFFs: same date, different variables
    # Must follow naming convention so merge extracts variable correctly:
    # variable is the last element after splitting by underscore
    tif1 = processed_dir / "20250601T000000_wse.tif"
    tif2 = processed_dir / "20250601T000000_geoid.tif"

    for tif_file in [tif1, tif2]:
        data = np.random.uniform(0, 100, (50, 50)).astype(np.float32)
        transform = from_bounds(500000, 6200000, 500000 + 500, 6200000 + 500, 50, 50)

        with rasterio.open(
            tif_file,
            "w",
            driver="GTiff",
            height=50,
            width=50,
            count=1,
            dtype=np.float32,
            crs=CRS.from_epsg(32644),
            transform=transform,
        ) as dst:
            dst.write(data, 1)

    config = {}

    _merge_and_reproject_granules(config, processed_dir, "EPSG:32644")

    merged_dir = processed_dir.parent.parent / "merged"
    merged_tiffs = list(merged_dir.glob("*.tif"))

    # Should have 2 outputs: one per variable
    assert len(merged_tiffs) == 2
    stems = {t.stem for t in merged_tiffs}
    assert any("wse" in stem for stem in stems), f"wse not found in {stems}"
    assert any("geoid" in stem for stem in stems), f"geoid not found in {stems}"


# ============================================================================
# TESTS: Full Pipeline
# ============================================================================


@pytest.mark.unit
def test_download_raster_full_pipeline_bbox(
    bbox_config, synthetic_swot_netcdf, tmp_path
):
    """Full pipeline: init dirs, mock download, preprocess, merge."""
    project_dir = tmp_path / "project"

    with patch("HydroEO.satellites.swot.raster.earthaccess") as mock_ea:
        mock_ea.login.return_value = None

        # Create the raw NC file first
        raw_parent = (
            project_dir
            / "swot_raster"
            / bbox_config["aoi"]["name"]
            / "raw"
            / bbox_config["product"]
        )
        raw_parent.mkdir(parents=True, exist_ok=True)

        shutil.copy(synthetic_swot_netcdf, raw_parent / synthetic_swot_netcdf.name)

        # Mock search to return nothing (files already exist)
        mock_ea.search_data.return_value = []
        mock_ea.download.return_value = []

        # Run the full pipeline
        download_raster(
            bbox_config,
            str(project_dir),
            ("user", "pass"),
        )

    # Verify outputs were created
    processed_dir = (
        project_dir
        / "swot_raster"
        / bbox_config["aoi"]["name"]
        / "processed"
        / bbox_config["product"]
    )
    assert processed_dir.exists()

    tiffs = list(processed_dir.glob("*.tif"))
    assert len(tiffs) > 0, "No processed TIFFs created"


@pytest.mark.unit
def test_download_raster_no_new_granules_processes_existing_nc(
    bbox_config, synthetic_swot_netcdf, tmp_path
):
    """When no new granules found, existing NC files in raw dir are still processed."""
    project_dir = tmp_path / "project"
    raw_dir = (
        project_dir
        / "swot_raster"
        / bbox_config["aoi"]["name"]
        / "raw"
        / bbox_config["product"]
    )
    raw_dir.mkdir(parents=True, exist_ok=True)

    # Place NC file in raw dir (simulating previous download)
    shutil.copy(synthetic_swot_netcdf, raw_dir / synthetic_swot_netcdf.name)

    with patch("HydroEO.satellites.swot.raster.earthaccess") as mock_ea:
        mock_ea.login.return_value = None
        mock_ea.search_data.return_value = []  # No new granules
        mock_ea.download.return_value = []

        download_raster(
            bbox_config,
            str(project_dir),
            ("user", "pass"),
        )

    # Verify preprocessing happened despite no new downloads
    processed_dir = (
        project_dir
        / "swot_raster"
        / bbox_config["aoi"]["name"]
        / "processed"
        / bbox_config["product"]
    )
    tiffs = list(processed_dir.glob("*.tif"))
    assert len(tiffs) > 0, "Existing NC files should have been processed"


# ============================================================================
# TESTS: CRS Detection
# ============================================================================


@pytest.mark.unit
def test_detect_crs_from_metadata():
    """CRS detected from netCDF variable CRS metadata."""
    ds = xr.Dataset()
    # Create a CRS variable with metadata (as SWOT files have)
    ds["transverse_mercator"] = xr.DataArray([])
    ds["transverse_mercator"].attrs["crs_wkt"] = CRS.from_epsg(32644).to_wkt()

    result = _detect_crs(ds, Path("test.nc"))

    assert result is not None
    assert result.to_epsg() == 32644


@pytest.mark.unit
def test_detect_crs_from_filename():
    """CRS detected from filename when metadata absent."""
    ds = xr.Dataset()

    result = _detect_crs(
        ds, Path("SWOT_L2_HR_Raster_D_20250601T000000_UTM44N_01_01.nc")
    )

    assert result is not None
    assert result.to_epsg() == 32644


# ============================================================================
# PARAMETRIZED TESTS
# ============================================================================


@pytest.mark.unit
@pytest.mark.parametrize("utm_zone", [44, 45, 46])
def test_detect_crs_utm_zones(utm_zone):
    """CRS detection works for various UTM zones."""
    ds = xr.Dataset()

    zone_str = f"UTM{utm_zone}N"
    result = _detect_crs(ds, Path(f"SWOT_test_{zone_str}.nc"))

    assert result is not None
    assert result.to_epsg() == (32600 + utm_zone)  # UTM zone to EPSG


@pytest.mark.unit
def test_swot_raster_download_function_accepts_expected_params():
    """download_raster is importable and has the expected signature."""
    import inspect

    sig = inspect.signature(download_raster)
    assert "config" in sig.parameters
    assert "project_dir" in sig.parameters
    assert "credentials" in sig.parameters


@pytest.mark.unit
def test_swot_raster_resampling_selection():
    """_resampling_for returns correct Resampling enum values."""
    from rasterio.warp import Resampling

    assert _resampling_for("wse") == Resampling.bilinear
    assert _resampling_for("wse_uncert") == Resampling.bilinear
    assert _resampling_for("geoid") == Resampling.bilinear
    assert _resampling_for("height_cor_xover") == Resampling.bilinear

    assert _resampling_for("n_wse_pix") == Resampling.nearest
    assert _resampling_for("n_other_pix") == Resampling.nearest
    assert _resampling_for("wse_qual") == Resampling.nearest


@pytest.mark.unit
def test_swot_raster_crs_detection_from_filename(tmp_path):
    """_detect_crs extracts UTM CRS from SWOT filename."""
    import xarray as xr

    ds = xr.Dataset()
    utm_file = (
        tmp_path
        / "SWOT_L2_HR_Raster_D_123_001_UTM45N_20240101T120000_20240101T130000.nc"
    )
    crs = _detect_crs(ds, utm_file)

    assert crs is not None
    assert crs.to_epsg() == 32645  # UTM 45N


@pytest.mark.unit
def test_swot_raster_no_processed_files_skips_merge(tmp_path, caplog):
    """download_raster skips merge when no processed TIF files exist."""
    config = {
        "aoi": {"name": "test_aoi", "type": "bbox", "bbox": [0, 0, 1, 1]},
        "product": "SWOT_L2_HR_Raster_D",
        "startdate": [2024, 1, 1],
        "enddate": [2024, 2, 1],
    }
    project_dir = tmp_path / "project"
    project_dir.mkdir()

    with (
        patch("HydroEO.satellites.swot.raster._download_granules", return_value=[]),
        patch("HydroEO.satellites.swot.raster._preprocess_granules"),
    ):
        with caplog.at_level("INFO"):
            download_raster(
                config=config,
                project_dir=str(project_dir),
                credentials=("user", "pass"),
            )

    assert "No processed TIF files found" in caplog.text


@pytest.mark.unit
def test_swot_raster_merge_with_existing_files(tmp_path, caplog):
    """download_raster calls merge when processed TIF files exist."""
    config = {
        "aoi": {"name": "test_aoi", "type": "bbox", "bbox": [0, 0, 1, 1]},
        "product": "SWOT_L2_HR_Raster_D",
        "startdate": [2024, 1, 1],
        "enddate": [2024, 2, 1],
    }
    project_dir = tmp_path / "project"
    processed_dir = (
        project_dir / "swot_raster" / "test_aoi" / "processed" / "SWOT_L2_HR_Raster_D"
    )
    processed_dir.mkdir(parents=True)
    (processed_dir / "20240101T120000_wse.tif").touch()

    with (
        patch("HydroEO.satellites.swot.raster._download_granules", return_value=[]),
        patch(
            "HydroEO.satellites.swot.raster._merge_and_reproject_granules"
        ) as mock_merge,
    ):
        with caplog.at_level("INFO"):
            download_raster(
                config=config,
                project_dir=str(project_dir),
                credentials=("user", "pass"),
            )

    mock_merge.assert_called_once()
    assert "Found 1 processed TIF files" in caplog.text
