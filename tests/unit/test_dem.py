"""Unit and integration tests for DEM download and processing (dem.py).

Tests cover:
- Bounding box geometry creation
- Tile clipping and merging workflows
- Full end-to-end download_cop_dem pipeline (with mocked network calls)
- Error handling for no tiles or no overlaps
"""

from __future__ import annotations

import zipfile
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest
import rasterio
from pyproj import CRS
from rasterio.transform import from_bounds
from shapely.geometry import box

from HydroEO.downloaders.dem import (
    download_cop_dem,
    clip_tif_to_aoi,
    merge_clipped_tifs,
    process,
)


# ============================================================================
# FIXTURES: Synthetic Data Generation
# ============================================================================


@pytest.fixture
def synthetic_dem_tif(tmp_path):
    """Create a minimal valid GeoTIFF with DEM-like data (EPSG:4326)."""
    tif_file = tmp_path / "test_dem.tif"

    # Simple DEM-like raster: 10x10 pixels covering bbox 10-11°E, 50-51°N
    width, height = 10, 10
    data = np.random.randint(100, 1000, (height, width), dtype=np.uint16)

    transform = from_bounds(10.0, 50.0, 11.0, 51.0, width, height)
    crs = CRS.from_epsg(4326)

    with rasterio.open(
        tif_file,
        "w",
        driver="GTiff",
        height=height,
        width=width,
        count=1,
        dtype=data.dtype,
        crs=crs,
        transform=transform,
        compress="lzw",
    ) as dst:
        dst.write(data, 1)

    return tif_file


@pytest.fixture
def synthetic_dem_zip(tmp_path, synthetic_dem_tif):
    """Create a ZIP archive containing a synthetic DEM tile (mimics CDSE download)."""
    zip_file = tmp_path / "COP_DEM_tile_001.zip"

    # Create a nested directory structure inside the ZIP (mimics CDSE format)
    with zipfile.ZipFile(zip_file, "w") as archive:
        # Add the DEM file with a nested path structure
        internal_path = "COP_DEM_GLO_30_srtm_srtm_rel/N50/E010/N50E010_0101_DEM.tif"
        archive.write(synthetic_dem_tif, arcname=internal_path)

    return zip_file


@pytest.fixture
def synthetic_dem_multiple_zips(tmp_path):
    """Create multiple synthetic DEM ZIP tiles covering a larger area."""
    zips = []

    # Create two overlapping DEM tiles for merge testing
    tile_specs = [
        {
            "zip_name": "COP_DEM_tile_001.zip",
            "bbox": (10.0, 50.0, 11.0, 51.0),  # E010 N050
            "dem_name": "N50E010_0101_DEM.tif",
        },
        {
            "zip_name": "COP_DEM_tile_002.zip",
            "bbox": (11.0, 50.0, 12.0, 51.0),  # E011 N050
            "dem_name": "N50E011_0101_DEM.tif",
        },
    ]

    for spec in tile_specs:
        zip_file = tmp_path / spec["zip_name"]

        # Create synthetic DEM for this tile
        width, height = 10, 10
        data = np.random.randint(100, 1000, (height, width), dtype=np.uint16)
        transform = from_bounds(*spec["bbox"], width, height)
        crs = CRS.from_epsg(4326)

        # Create temp DEM file
        temp_dem = tmp_path / f"temp_{spec['dem_name']}"
        with rasterio.open(
            temp_dem,
            "w",
            driver="GTiff",
            height=height,
            width=width,
            count=1,
            dtype=data.dtype,
            crs=crs,
            transform=transform,
            compress="lzw",
        ) as dst:
            dst.write(data, 1)

        # Create ZIP with nested path
        internal_path = f"COP_DEM_GLO_30_srtm_srtm_rel/N50/{spec['dem_name'].split('_')[1]}/{spec['dem_name']}"
        with zipfile.ZipFile(zip_file, "w") as archive:
            archive.write(temp_dem, arcname=internal_path)

        zips.append(zip_file)
        temp_dem.unlink()  # Clean up temp DEM

    return zips


# ============================================================================
# TESTS: Geometry and bbox handling
# ============================================================================


@pytest.mark.unit
def test_bbox_to_geojson_polygon():
    """Test that bbox coordinates are converted to a valid GeoJSON polygon."""
    from HydroEO.downloaders.dem import download_cop_dem

    minx, miny, maxx, maxy = 10.0, 50.0, 11.0, 51.0

    # We'll verify the bbox conversion by checking that the function creates
    # the correct geometry internally. We mock the download to avoid network calls.
    with patch("HydroEO.downloaders.dem.download_cop_dem_for_polygon") as mock_download:
        mock_download.return_value = []
        try:
            download_cop_dem(minx, miny, maxx, maxy, "/tmp/out")
        except RuntimeError:
            # Expected: no tiles found
            pass

        # Verify the aoi_geojson passed to download_cop_dem_for_polygon
        call_args = mock_download.call_args
        aoi = call_args.kwargs["aoi_geojson"]

        assert aoi["type"] == "Polygon"
        coords = aoi["coordinates"][0]
        # Check that polygon is closed (first and last points are the same)
        assert coords[0] == coords[-1]
        # Check corners
        assert coords[0] == [minx, miny]
        assert coords[1] == [maxx, miny]
        assert coords[2] == [maxx, maxy]
        assert coords[3] == [minx, maxy]


# ============================================================================
# TESTS: Clipping and merging
# ============================================================================


@pytest.mark.unit
def test_clip_tif_to_aoi(synthetic_dem_tif, tmp_path):
    """Test that a single GeoTIFF can be clipped to an AOI."""

    # Define AOI as a smaller box within the DEM tile
    aoi = box(10.2, 50.2, 10.8, 50.8)

    tmp_dir = tmp_path / "clipping"
    tmp_dir.mkdir()

    result = clip_tif_to_aoi(str(synthetic_dem_tif), aoi, str(tmp_dir))

    # Check that clipped file was created
    assert result is not None
    assert Path(result).exists()

    # Verify the clipped raster has smaller bounds
    with rasterio.open(result) as src_clipped:
        clipped_bounds = src_clipped.bounds
        # Clipped bounds should be smaller than original
        original_bounds = (10.0, 50.0, 11.0, 51.0)
        assert clipped_bounds.left > original_bounds[0]
        assert clipped_bounds.bottom > original_bounds[1]


@pytest.mark.unit
def test_clip_tif_no_overlap(synthetic_dem_tif, tmp_path):
    """Test that clip_tif_to_aoi returns None for non-overlapping AOI."""

    # Define AOI that doesn't overlap the DEM tile (15-16°E, 50-51°N)
    aoi = box(15.0, 50.0, 16.0, 51.0)

    tmp_dir = tmp_path / "clipping"
    tmp_dir.mkdir()

    result = clip_tif_to_aoi(str(synthetic_dem_tif), aoi, str(tmp_dir))

    # Should return None for non-overlapping tiles
    assert result is None


@pytest.mark.unit
def test_merge_clipped_tifs(tmp_path):
    """Test merging of multiple clipped GeoTIFFs."""

    # Create two synthetic clipped GeoTIFFs
    tif_paths = []
    for i, bbox in enumerate([(10.0, 50.0, 10.5, 50.5), (10.5, 50.0, 11.0, 50.5)]):
        tif_file = tmp_path / f"clipped_{i}.tif"
        width, height = 5, 5
        data = np.full((height, width), 100 + i * 10, dtype=np.uint16)
        transform = from_bounds(*bbox, width, height)
        crs = CRS.from_epsg(4326)

        with rasterio.open(
            tif_file,
            "w",
            driver="GTiff",
            height=height,
            width=width,
            count=1,
            dtype=data.dtype,
            crs=crs,
            transform=transform,
            compress="lzw",
        ) as dst:
            dst.write(data, 1)

        tif_paths.append(str(tif_file))

    # Merge the clipped TIFFs
    output_path = tmp_path / "merged_dem.tif"
    merge_clipped_tifs(tif_paths, str(output_path))

    # Check that merged file was created
    assert output_path.exists()

    # Verify merged raster properties
    with rasterio.open(output_path) as src_merged:
        assert src_merged.crs == CRS.from_epsg(4326)
        # Merged bounds should encompass both input tiles
        assert src_merged.bounds.left == 10.0
        assert src_merged.bounds.bottom == 50.0


# ============================================================================
# TESTS: End-to-end pipeline
# ============================================================================


@pytest.mark.unit
def test_process_with_valid_tifs(synthetic_dem_tif, tmp_path):
    """Test the full process() function with a single valid DEM tile."""

    aoi_geojson = {
        "type": "Polygon",
        "coordinates": [
            [
                [10.1, 50.1],
                [10.9, 50.1],
                [10.9, 50.9],
                [10.1, 50.9],
                [10.1, 50.1],
            ]
        ],
    }

    output_path = tmp_path / "merged_dem.tif"

    process(
        aoi_geojson=aoi_geojson,
        tif_paths=[str(synthetic_dem_tif)],
        output_path=output_path,
        target_crs="EPSG:4326",
    )

    # Check that output was created
    assert output_path.exists()

    # Verify output is a valid GeoTIFF
    with rasterio.open(output_path) as src:
        assert src.crs == CRS.from_epsg(4326)
        assert src.count == 1
        data = src.read(1)
        assert data.size > 0


@pytest.mark.unit
def test_process_no_overlapping_tiles(synthetic_dem_tif, tmp_path):
    """Test that process() raises error when no tiles overlap the AOI."""
    aoi_geojson = {
        "type": "Polygon",
        "coordinates": [
            [
                [15.0, 50.0],
                [16.0, 50.0],
                [16.0, 51.0],
                [15.0, 51.0],
                [15.0, 50.0],
            ]
        ],
    }

    output_path = tmp_path / "merged_dem.tif"

    with pytest.raises(RuntimeError, match="No input rasters overlapped the AOI"):
        process(
            aoi_geojson=aoi_geojson,
            tif_paths=[str(synthetic_dem_tif)],
            output_path=output_path,
            target_crs="EPSG:4326",
        )


@pytest.mark.unit
def test_download_cop_dem_mocked(tmp_path):
    """Test download_cop_dem with fully mocked CDSE operations."""
    output_dir = tmp_path / "dem_output"
    output_dir.mkdir()

    # Create a temporary DEM file to simulate a downloaded and unzipped tile
    test_dem = tmp_path / "test_dem_source.tif"
    width, height = 10, 10
    data = np.random.randint(100, 1000, (height, width), dtype=np.uint16)
    transform = from_bounds(10.0, 50.0, 11.0, 51.0, width, height)
    crs = CRS.from_epsg(4326)

    with rasterio.open(
        test_dem,
        "w",
        driver="GTiff",
        height=height,
        width=width,
        count=1,
        dtype=data.dtype,
        crs=crs,
        transform=transform,
        compress="lzw",
    ) as dst:
        dst.write(data, 1)

    # Mock the download and unzipping process
    with patch("HydroEO.downloaders.dem.download_cop_dem_for_polygon") as mock_download:
        # Create a temporary ZIP file
        mock_zip = tmp_path / "mock_tile.zip"
        with zipfile.ZipFile(mock_zip, "w") as zf:
            zf.write(test_dem, arcname="COP_DEM_GLO_30/N50/E010/N50E010_0101_DEM.tif")

        mock_download.return_value = [mock_zip]

        # Call download_cop_dem
        result = download_cop_dem(
            minx=10.0,
            miny=50.0,
            maxx=11.0,
            maxy=51.0,
            output_dir=str(output_dir),
            username="test_user",
            password="test_pass",
        )

        # Check that result file exists
        assert result.exists()
        assert result.name == "cop_dem_merged.tif"

        # Verify CDSE was called with correct parameters
        mock_download.assert_called_once()
        call_kwargs = mock_download.call_args.kwargs
        assert call_kwargs["username"] == "test_user"
        assert call_kwargs["password"] == "test_pass"


@pytest.mark.unit
def test_download_cop_dem_no_tiles(tmp_path):
    """Test download_cop_dem with no tiles found."""
    output_dir = tmp_path / "dem_output"
    output_dir.mkdir()

    # Mock the download to return no tiles
    with patch("HydroEO.downloaders.dem.download_cop_dem_for_polygon") as mock_download:
        mock_download.return_value = []

        with pytest.raises(RuntimeError, match="No COP-DEM tiles found"):
            download_cop_dem(
                minx=10.0,
                miny=50.0,
                maxx=11.0,
                maxy=51.0,
                output_dir=str(output_dir),
                username="test_user",
                password="test_pass",
            )


@pytest.mark.unit
def test_download_cop_dem_custom_dataset(tmp_path):
    """Test download_cop_dem with custom dataset parameter."""
    output_dir = tmp_path / "dem_output"
    output_dir.mkdir()

    test_dem = tmp_path / "test_dem_source.tif"
    width, height = 10, 10
    data = np.random.randint(100, 1000, (height, width), dtype=np.uint16)
    transform = from_bounds(10.0, 50.0, 11.0, 51.0, width, height)
    crs = CRS.from_epsg(4326)

    with rasterio.open(
        test_dem,
        "w",
        driver="GTiff",
        height=height,
        width=width,
        count=1,
        dtype=data.dtype,
        crs=crs,
        transform=transform,
        compress="lzw",
    ) as dst:
        dst.write(data, 1)

    with patch("HydroEO.downloaders.dem.download_cop_dem_for_polygon") as mock_download:
        mock_zip = tmp_path / "mock_tile.zip"
        with zipfile.ZipFile(mock_zip, "w") as zf:
            zf.write(test_dem, arcname="COP_DEM_GLO_90/N50/E010/N50E010_0101_DEM.tif")

        mock_download.return_value = [mock_zip]

        custom_dataset = "COP-DEM_GLO-90-DGED/2023_1"
        download_cop_dem(
            minx=10.0,
            miny=50.0,
            maxx=11.0,
            maxy=51.0,
            output_dir=str(output_dir),
            username="test_user",
            password="test_pass",
            dataset=custom_dataset,
        )

        # Check that the custom dataset was passed
        call_kwargs = mock_download.call_args.kwargs
        assert call_kwargs["dataset"] == custom_dataset


@pytest.mark.unit
def test_download_cop_dem_custom_output_filename(tmp_path):
    """Test download_cop_dem with custom output filename."""
    output_dir = tmp_path / "dem_output"
    output_dir.mkdir()

    test_dem = tmp_path / "test_dem_source.tif"
    width, height = 10, 10
    data = np.random.randint(100, 1000, (height, width), dtype=np.uint16)
    transform = from_bounds(10.0, 50.0, 11.0, 51.0, width, height)
    crs = CRS.from_epsg(4326)

    with rasterio.open(
        test_dem,
        "w",
        driver="GTiff",
        height=height,
        width=width,
        count=1,
        dtype=data.dtype,
        crs=crs,
        transform=transform,
        compress="lzw",
    ) as dst:
        dst.write(data, 1)

    with patch("HydroEO.downloaders.dem.download_cop_dem_for_polygon") as mock_download:
        mock_zip = tmp_path / "mock_tile.zip"
        with zipfile.ZipFile(mock_zip, "w") as zf:
            zf.write(test_dem, arcname="COP_DEM_GLO_30/N50/E010/N50E010_0101_DEM.tif")

        mock_download.return_value = [mock_zip]

        custom_filename = "custom_dem_output.tif"
        result = download_cop_dem(
            minx=10.0,
            miny=50.0,
            maxx=11.0,
            maxy=51.0,
            output_dir=str(output_dir),
            username="test_user",
            password="test_pass",
            output_filename=custom_filename,
        )

        assert result.name == custom_filename
        assert result.exists()


# ============================================================================
# TESTS: CLI integration
# ============================================================================


@pytest.mark.unit
def test_fetch_cop_dem_cli_help():
    """Test that the fetch cop-dem CLI command is available and has help."""
    from typer.testing import CliRunner
    from HydroEO.cli import app

    runner = CliRunner()
    result = runner.invoke(app, ["fetch", "cop-dem", "--help"])
    assert result.exit_code == 0
    assert "cop-dem" in result.output.lower()
    assert "--bbox" in result.output
    assert "--cdse-username" in result.output
    assert "--cdse-password" in result.output
