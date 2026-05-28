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
    _resolve_cdse_dataset,
    _validate_layers,
    clip_tif_to_aoi,
    download_cop_dem,
    merge_clipped_tifs,
    process,
    VALID_LAYERS,
    _CDSE_GLO30,
    _CDSE_GLO90,
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


def _make_tif(path: Path, bbox: tuple, dtype=np.uint16) -> None:
    """Write a minimal synthetic GeoTIFF to *path* covering *bbox*."""
    width, height = 10, 10
    info = np.iinfo(dtype) if np.dtype(dtype).kind in ("u", "i") else None
    high = min(1000, info.max) if info is not None else None
    if high is not None:
        data = np.random.randint(1, high, (height, width), dtype=dtype)
    else:
        data = np.random.rand(height, width).astype(dtype)
    transform = from_bounds(*bbox, width, height)
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=height,
        width=width,
        count=1,
        dtype=data.dtype,
        crs=CRS.from_epsg(4326),
        transform=transform,
        compress="lzw",
    ) as dst:
        dst.write(data, 1)


@pytest.fixture
def synthetic_auxfiles_zip(tmp_path):
    """
    Create a ZIP that mirrors the real CDSE COP-DEM archive structure:

        tile_name/
        ├── DEM/tile_name_DEM.tif
        └── AUXFILES/
            ├── tile_name_EDM.tif
            ├── tile_name_FLM.tif
            ├── tile_name_HEM.tif
            └── tile_name_WBM.tif
    """
    tile = "Copernicus_DSM_10_N50_00_E010_00"
    bbox = (10.0, 50.0, 11.0, 51.0)
    zip_file = tmp_path / "COP_DEM_auxfiles_tile.zip"

    layers = {
        f"{tile}/DEM/{tile}_DEM.tif": np.uint16,
        f"{tile}/AUXFILES/{tile}_EDM.tif": np.uint8,
        f"{tile}/AUXFILES/{tile}_FLM.tif": np.uint8,
        f"{tile}/AUXFILES/{tile}_HEM.tif": np.float32,
        f"{tile}/AUXFILES/{tile}_WBM.tif": np.uint8,
    }

    with zipfile.ZipFile(zip_file, "w") as archive:
        for internal_path, dtype in layers.items():
            tif_tmp = tmp_path / Path(internal_path).name
            _make_tif(tif_tmp, bbox, dtype=dtype)
            archive.write(tif_tmp, arcname=internal_path)
            tif_tmp.unlink()

    return zip_file


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

    test_dem = tmp_path / "test_dem_source.tif"
    _make_tif(test_dem, (10.0, 50.0, 11.0, 51.0))

    with patch("HydroEO.downloaders.dem.download_cop_dem_for_polygon") as mock_download:
        mock_zip = tmp_path / "mock_tile.zip"
        with zipfile.ZipFile(mock_zip, "w") as zf:
            zf.write(test_dem, arcname="COP_DEM_GLO_30/N50/E010/N50E010_0101_DEM.tif")

        mock_download.return_value = [mock_zip]

        result = download_cop_dem(
            minx=10.0,
            miny=50.0,
            maxx=11.0,
            maxy=51.0,
            output_dir=str(output_dir),
            username="test_user",
            password="test_pass",
            layers=["DEM30"],
        )

        # Returns a dict keyed by layer name
        assert isinstance(result, dict)
        assert "DEM30" in result
        assert result["DEM30"].exists()
        assert result["DEM30"].name == "cop_dem_merged_DEM30.tif"

        mock_download.assert_called_once()
        call_kwargs = mock_download.call_args.kwargs
        assert call_kwargs["username"] == "test_user"
        assert call_kwargs["password"] == "test_pass"
        # GLO-30 dataset used for DEM30
        assert call_kwargs["dataset"] == _CDSE_GLO30


@pytest.mark.unit
def test_download_cop_dem_no_tiles(tmp_path):
    """Test download_cop_dem with no tiles found."""
    output_dir = tmp_path / "dem_output"
    output_dir.mkdir()

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
                layers=["DEM30"],
            )


@pytest.mark.unit
def test_download_cop_dem_dem90(tmp_path):
    """Test download_cop_dem with DEM90 selects the GLO-90 CDSE dataset."""
    output_dir = tmp_path / "dem_output"
    output_dir.mkdir()

    test_dem = tmp_path / "test_dem_source.tif"
    _make_tif(test_dem, (10.0, 50.0, 11.0, 51.0))

    with patch("HydroEO.downloaders.dem.download_cop_dem_for_polygon") as mock_download:
        mock_zip = tmp_path / "mock_tile.zip"
        with zipfile.ZipFile(mock_zip, "w") as zf:
            zf.write(test_dem, arcname="COP_DEM_GLO_90/N50/E010/N50E010_0101_DEM.tif")

        mock_download.return_value = [mock_zip]

        result = download_cop_dem(
            minx=10.0,
            miny=50.0,
            maxx=11.0,
            maxy=51.0,
            output_dir=str(output_dir),
            username="test_user",
            password="test_pass",
            layers=["DEM90"],
        )

        # GLO-90 dataset must be selected when DEM90 is requested
        call_kwargs = mock_download.call_args.kwargs
        assert call_kwargs["dataset"] == _CDSE_GLO90

        assert "DEM90" in result
        assert result["DEM90"].name == "cop_dem_merged_DEM90.tif"


@pytest.mark.unit
def test_download_cop_dem_custom_output_basename(tmp_path):
    """Test download_cop_dem with a custom output_basename."""
    output_dir = tmp_path / "dem_output"
    output_dir.mkdir()

    test_dem = tmp_path / "test_dem_source.tif"
    _make_tif(test_dem, (10.0, 50.0, 11.0, 51.0))

    with patch("HydroEO.downloaders.dem.download_cop_dem_for_polygon") as mock_download:
        mock_zip = tmp_path / "mock_tile.zip"
        with zipfile.ZipFile(mock_zip, "w") as zf:
            zf.write(test_dem, arcname="COP_DEM_GLO_30/N50/E010/N50E010_0101_DEM.tif")

        mock_download.return_value = [mock_zip]

        result = download_cop_dem(
            minx=10.0,
            miny=50.0,
            maxx=11.0,
            maxy=51.0,
            output_dir=str(output_dir),
            username="test_user",
            password="test_pass",
            layers=["DEM30"],
            output_basename="my_dem",
        )

        # Basename + layer suffix
        assert "DEM30" in result
        assert result["DEM30"].name == "my_dem_DEM30.tif"
        assert result["DEM30"].exists()


@pytest.mark.unit
def test_download_cop_dem_custom_output_basename_with_extension(tmp_path):
    """Extension in output_basename should be stripped gracefully."""
    output_dir = tmp_path / "dem_output"
    output_dir.mkdir()

    test_dem = tmp_path / "test_dem_source.tif"
    _make_tif(test_dem, (10.0, 50.0, 11.0, 51.0))

    with patch("HydroEO.downloaders.dem.download_cop_dem_for_polygon") as mock_download:
        mock_zip = tmp_path / "mock_tile.zip"
        with zipfile.ZipFile(mock_zip, "w") as zf:
            zf.write(test_dem, arcname="COP_DEM_GLO_30/N50/E010/N50E010_0101_DEM.tif")

        mock_download.return_value = [mock_zip]

        result = download_cop_dem(
            minx=10.0,
            miny=50.0,
            maxx=11.0,
            maxy=51.0,
            output_dir=str(output_dir),
            username="test_user",
            password="test_pass",
            layers=["DEM30"],
            output_basename="my_dem.tif",  # extension supplied — should be stripped
        )

        assert result["DEM30"].name == "my_dem_DEM30.tif"


# ============================================================================
# TESTS: Layer validation helpers
# ============================================================================


@pytest.mark.unit
def test_validate_layers_valid():
    """No error raised for all individually valid layer names."""
    for layer in VALID_LAYERS:
        _validate_layers([layer])  # should not raise


@pytest.mark.unit
def test_validate_layers_unknown():
    """Unknown layer name raises ValueError."""
    with pytest.raises(ValueError, match="Unknown layer"):
        _validate_layers(["BOGUS"])


@pytest.mark.unit
def test_validate_layers_dem30_dem90_conflict():
    """Requesting DEM30 and DEM90 together raises ValueError."""
    with pytest.raises(ValueError, match="Cannot request DEM30 and DEM90 together"):
        _validate_layers(["DEM30", "DEM90"])


@pytest.mark.unit
def test_validate_layers_mixed_valid_invalid():
    """Mix of valid and invalid layer names still raises on the invalid one."""
    with pytest.raises(ValueError, match="Unknown layer"):
        _validate_layers(["DEM30", "NOPE"])


@pytest.mark.unit
def test_resolve_cdse_dataset_dem30():
    """DEM30 resolves to the GLO-30 CDSE dataset."""
    assert _resolve_cdse_dataset(["DEM30"]) == _CDSE_GLO30


@pytest.mark.unit
def test_resolve_cdse_dataset_dem90():
    """DEM90 resolves to the GLO-90 CDSE dataset."""
    assert _resolve_cdse_dataset(["DEM90"]) == _CDSE_GLO90


@pytest.mark.unit
def test_resolve_cdse_dataset_aux_only():
    """Aux-layer-only requests default to the GLO-30 CDSE dataset."""
    for layer in ("EDM", "FLM", "HEM", "WBM"):
        assert _resolve_cdse_dataset([layer]) == _CDSE_GLO30


@pytest.mark.unit
def test_resolve_cdse_dataset_dem90_with_aux():
    """DEM90 + aux layers still resolve to GLO-90."""
    assert _resolve_cdse_dataset(["DEM90", "WBM"]) == _CDSE_GLO90


# ============================================================================
# TESTS: Multi-layer pipeline
# ============================================================================


@pytest.mark.unit
def test_download_cop_dem_multi_layer(tmp_path, synthetic_auxfiles_zip):
    """Request DEM30 + WBM; both output files should be created from one ZIP."""
    output_dir = tmp_path / "dem_output"
    output_dir.mkdir()

    with patch("HydroEO.downloaders.dem.download_cop_dem_for_polygon") as mock_download:
        mock_download.return_value = [synthetic_auxfiles_zip]

        result = download_cop_dem(
            minx=10.0,
            miny=50.0,
            maxx=11.0,
            maxy=51.0,
            output_dir=str(output_dir),
            username="test_user",
            password="test_pass",
            layers=["DEM30", "WBM"],
        )

        assert set(result.keys()) == {"DEM30", "WBM"}
        assert result["DEM30"].name == "cop_dem_merged_DEM30.tif"
        assert result["WBM"].name == "cop_dem_merged_WBM.tif"
        assert result["DEM30"].exists()
        assert result["WBM"].exists()

        # ZIPs downloaded only once (single CDSE call)
        mock_download.assert_called_once()


@pytest.mark.unit
def test_download_cop_dem_all_aux_layers(tmp_path, synthetic_auxfiles_zip):
    """Request all four aux layers in one call."""
    output_dir = tmp_path / "dem_output"
    output_dir.mkdir()

    with patch("HydroEO.downloaders.dem.download_cop_dem_for_polygon") as mock_download:
        mock_download.return_value = [synthetic_auxfiles_zip]

        result = download_cop_dem(
            minx=10.0,
            miny=50.0,
            maxx=11.0,
            maxy=51.0,
            output_dir=str(output_dir),
            username="test_user",
            password="test_pass",
            layers=["EDM", "FLM", "HEM", "WBM"],
        )

        assert set(result.keys()) == {"EDM", "FLM", "HEM", "WBM"}
        for layer, path in result.items():
            assert path.exists(), f"Output for {layer} not found"
            assert path.name == f"cop_dem_merged_{layer}.tif"

        # Aux-only → GLO-30 dataset
        assert mock_download.call_args.kwargs["dataset"] == _CDSE_GLO30


@pytest.mark.unit
def test_download_cop_dem_aux_only_defaults_to_glo30(tmp_path, synthetic_auxfiles_zip):
    """Requesting only aux layers (no DEM30/DEM90) uses the GLO-30 dataset."""
    output_dir = tmp_path / "dem_output"
    output_dir.mkdir()

    with patch("HydroEO.downloaders.dem.download_cop_dem_for_polygon") as mock_download:
        mock_download.return_value = [synthetic_auxfiles_zip]

        download_cop_dem(
            minx=10.0,
            miny=50.0,
            maxx=11.0,
            maxy=51.0,
            output_dir=str(output_dir),
            username="u",
            password="p",
            layers=["WBM"],
        )

        assert mock_download.call_args.kwargs["dataset"] == _CDSE_GLO30


@pytest.mark.unit
def test_download_cop_dem_dem90_with_aux(tmp_path, synthetic_auxfiles_zip):
    """DEM90 + WBM uses the GLO-90 dataset for both layers."""
    output_dir = tmp_path / "dem_output"
    output_dir.mkdir()

    with patch("HydroEO.downloaders.dem.download_cop_dem_for_polygon") as mock_download:
        mock_download.return_value = [synthetic_auxfiles_zip]

        result = download_cop_dem(
            minx=10.0,
            miny=50.0,
            maxx=11.0,
            maxy=51.0,
            output_dir=str(output_dir),
            username="u",
            password="p",
            layers=["DEM90", "WBM"],
        )

        assert mock_download.call_args.kwargs["dataset"] == _CDSE_GLO90
        assert set(result.keys()) == {"DEM90", "WBM"}


@pytest.mark.unit
def test_download_cop_dem_default_layers(tmp_path):
    """Calling download_cop_dem without layers= defaults to DEM30."""
    output_dir = tmp_path / "dem_output"
    output_dir.mkdir()

    test_dem = tmp_path / "test_dem_source.tif"
    _make_tif(test_dem, (10.0, 50.0, 11.0, 51.0))

    with patch("HydroEO.downloaders.dem.download_cop_dem_for_polygon") as mock_download:
        mock_zip = tmp_path / "mock_tile.zip"
        with zipfile.ZipFile(mock_zip, "w") as zf:
            zf.write(test_dem, arcname="COP_DEM_GLO_30/N50/E010/N50E010_0101_DEM.tif")
        mock_download.return_value = [mock_zip]

        result = download_cop_dem(
            minx=10.0,
            miny=50.0,
            maxx=11.0,
            maxy=51.0,
            output_dir=str(output_dir),
            username="u",
            password="p",
        )

        assert list(result.keys()) == ["DEM30"]
        assert result["DEM30"].exists()


@pytest.mark.unit
def test_download_cop_dem_invalid_layer(tmp_path):
    """Invalid layer name raises ValueError before any download."""
    with pytest.raises(ValueError, match="Unknown layer"):
        download_cop_dem(
            minx=10.0,
            miny=50.0,
            maxx=11.0,
            maxy=51.0,
            output_dir=str(tmp_path),
            username="u",
            password="p",
            layers=["INVALID"],
        )


@pytest.mark.unit
def test_download_cop_dem_dem30_dem90_conflict(tmp_path):
    """Requesting DEM30 and DEM90 together raises ValueError."""
    with pytest.raises(ValueError, match="Cannot request DEM30 and DEM90 together"):
        download_cop_dem(
            minx=10.0,
            miny=50.0,
            maxx=11.0,
            maxy=51.0,
            output_dir=str(tmp_path),
            username="u",
            password="p",
            layers=["DEM30", "DEM90"],
        )
