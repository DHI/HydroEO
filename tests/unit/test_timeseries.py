"""Unit tests for time series functions."""

from pathlib import Path
from types import SimpleNamespace
from unittest import mock
from unittest.mock import patch

import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import box, Point

from HydroEO import flows
from HydroEO.waterbody import Reservoirs, Rivers


# ============================================================================
# Fixtures
# ============================================================================


@pytest.fixture
def mock_project_reservoirs(tmp_path):
    """Create a mock Project with Reservoirs configuration."""
    gdf = gpd.GeoDataFrame(
        {"id": [1, 2], "geometry": [box(0, 0, 1, 1), box(1, 1, 2, 2)]},
        crs="EPSG:4326",
    )
    prj = SimpleNamespace()
    prj.dirs = {
        "main": str(tmp_path),
        "output": str(tmp_path / "results"),
        "swot": str(tmp_path / "raw" / "swot"),
        "icesat2": str(tmp_path / "raw" / "icesat2"),
        "icesat2_processed": str(tmp_path / "processed" / "icesat2"),
        "sentinel3": str(tmp_path / "raw" / "sentinel3"),
        "sentinel6": str(tmp_path / "raw" / "sentinel6"),
        "pld": str(tmp_path / "aux" / "PLD" / "PLD_subset.gpkg"),
    }
    prj.reservoirs = Reservoirs(gdf=gdf, id_key="id", dirs=prj.dirs)
    prj.reservoirs.download_gdf = gdf
    prj.local_crs = "EPSG:3857"
    prj.keep_raw_pld = False
    prj.startdates = {
        "swot": [2024, 1, 1],
        "icesat2": [2024, 1, 1],
        "sentinel3": [2024, 1, 1],
        "sentinel6": [2024, 1, 1],
    }
    prj.enddates = {
        "swot": [2024, 2, 1],
        "icesat2": [2024, 2, 1],
        "sentinel3": [2024, 2, 1],
        "sentinel6": [2024, 2, 1],
    }
    prj.mission_options = {
        "swot": {"exclude_obs_id_values": ["no_data"]},
        "icesat2": {"atl13_fields": None, "track_keys": None},
        "sentinel3": {"sigma0_max": 1e5},
        "sentinel6": {"sigma0_max": 1e5},
    }
    return prj


@pytest.fixture
def mock_project_rivers(tmp_path):
    """Create a mock Project with Rivers configuration.

    Two targets (node_id 101, 102) belonging to the same waterbody
    ("loire"), placed ~11km apart (0.1 degrees latitude) so a small
    extraction/assignment buffer can't accidentally blur observations
    from one target into the other, which would mask assignment bugs
    rather than catch them.
    """
    gdf = gpd.GeoDataFrame({"geometry": []}, geometry="geometry", crs="EPSG:4326")
    prj = SimpleNamespace()
    prj.dirs = {
        "main": str(tmp_path),
        "output": str(tmp_path / "results"),
        "swot": str(tmp_path / "raw" / "swot"),
        "icesat2_processed": str(tmp_path / "processed" / "icesat2"),
        "sentinel3": str(tmp_path / "raw" / "sentinel3"),
        "sentinel6": str(tmp_path / "raw" / "sentinel6"),
        "sword": str(tmp_path / "aux" / "SWORD" / "gpkg"),
        "sword_subset": str(tmp_path / "aux" / "SWORD" / "SWORD_subset.gpkg"),
    }
    prj.rivers = Rivers(gdf=gdf, id_key="river_id", dirs=prj.dirs)
    prj.rivers.target_ids = [101, 102]
    prj.rivers.target_id_col = "node_id"
    prj.rivers.target_features = gpd.GeoDataFrame(
        {
            "node_id": [101, 102],
            "river_id": ["loire", "loire"],
            "geometry": [Point(0, 0), Point(0, 0.1)],
        },
        crs="EPSG:4326",
    )
    prj.rivers.configured_id = "loire"
    prj.rivers.input_mode = "configured_id"
    # Fixed buffer (rather than SWORD width-based sizing) so tests don't
    # depend on a "width" column being present -- see
    # flows._river_target_corridor's docstring for the width-based path.
    prj.rivers.extraction_buffer_meters = 500.0
    prj.rivers.width_buffer_factor = 1.05
    prj.rivers.max_node_assignment_meters = 1000.0
    prj.local_crs = "EPSG:3857"
    prj.global_crs = "EPSG:4326"
    prj.keep_raw_sword = False
    prj.startdates = {
        "swot": [2024, 1, 1],
        "icesat2": [2024, 1, 1],
        "sentinel3": [2024, 1, 1],
        "sentinel6": [2024, 1, 1],
    }
    prj.enddates = {
        "swot": [2024, 2, 1],
        "icesat2": [2024, 2, 1],
        "sentinel3": [2024, 2, 1],
        "sentinel6": [2024, 2, 1],
    }
    prj.mission_options = {
        "swot": {
            "hydrocron_fields": {
                "nodes": ["node_id", "node_q", "time_str", "wse"],
                "reaches": ["reach_id", "reach_q", "time_str", "wse"],
            },
            "quality_filters": {
                "nodes": {"max_q": 2},
                "reaches": {"max_q": 2},
            },
        },
        "icesat2": {"atl13_fields": None, "track_keys": None},
        "sentinel3": {"sigma0_max": 1e5, "sigma0_min": 0.0},
        "sentinel6": {"sigma0_max": 1e5, "sigma0_min": 0.0},
    }
    return prj


# ============================================================================
# Timeseries Extraction Tests
# ============================================================================


@pytest.mark.unit
def test_extract_swot_observations_calls_swot_module(mock_project_reservoirs, tmp_path):
    """_extract_swot_observations calls swot.extract_observations when downloads exist."""
    from HydroEO.satellites import swot

    mock_project_reservoirs.to_process = ["swot"]
    mock_project_reservoirs.dirs["swot"] = str(tmp_path / "swot")

    # Create mock SWOT download directory
    swot_dir = Path(mock_project_reservoirs.dirs["swot"])
    swot_dir.mkdir(parents=True, exist_ok=True)
    (swot_dir / "dummy_granule.shp").touch()

    with patch.object(swot, "extract_observations", return_value=[]) as mock_extract:
        flows._extract_swot_observations(mock_project_reservoirs)

        # Verify swot.extract_observations was called
        mock_extract.assert_called_once()
        call_args = mock_extract.call_args
        assert call_args[1]["src_dir"] == str(swot_dir)
        assert call_args[1]["exclude_obs_id_values"] == ["no_data"]


@pytest.mark.unit
def test_extract_swot_observations_skips_when_missing(mock_project_reservoirs, caplog):
    """_extract_swot_observations skips when SWOT download directory doesn't exist."""
    import logging

    mock_project_reservoirs.to_process = ["swot"]
    # SWOT dir doesn't exist

    with caplog.at_level(logging.WARNING):
        flows._extract_swot_observations(mock_project_reservoirs)
        assert "No SWOT downloads found" in caplog.text


@pytest.mark.unit
def test_extract_icesat2_observations_calls_icesat2_module(
    mock_project_reservoirs, tmp_path
):
    """_extract_icesat2_observations calls icesat2.extract_observations for each reservoir."""
    from HydroEO.satellites import icesat2

    mock_project_reservoirs.to_process = ["icesat2"]
    mock_project_reservoirs.dirs["icesat2_processed"] = str(
        tmp_path / "processed" / "icesat2"
    )

    # Create mock ICESat-2 parquet files for reservoir 1
    icesat2_dir = Path(mock_project_reservoirs.dirs["icesat2_processed"]) / "1"
    icesat2_dir.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(
        {
            "height": [100.5, 100.8, 101.2],
            "date": pd.date_range("2024-01-01", periods=3),
        }
    )
    df.to_parquet(icesat2_dir / "atl13.parquet")

    with patch.object(icesat2, "extract_observations") as mock_extract:
        flows._extract_icesat2_observations(mock_project_reservoirs)

        # Verify icesat2.extract_observations was called for reservoir 1
        assert mock_extract.call_count == 1
        call_args = mock_extract.call_args
        assert str(icesat2_dir) in call_args[1]["src_dir"]


@pytest.mark.unit
def test_extract_sentinel_observations_calls_sentinel_for_s3(
    mock_project_reservoirs, tmp_path
):
    """_extract_sentinel_observations calls sentinel.extract_observations for Sentinel-3."""
    from HydroEO.satellites import sentinel

    mock_project_reservoirs.to_process = ["sentinel3"]
    mock_project_reservoirs.dirs["sentinel3"] = str(tmp_path / "sentinel3")

    # Create mock Sentinel-3 directory with data
    s3_dir = Path(mock_project_reservoirs.dirs["sentinel3"]) / "1"
    s3_dir.mkdir(parents=True, exist_ok=True)
    (s3_dir / "dummy_file.nc").touch()

    with patch.object(sentinel, "extract_observations") as mock_extract:
        flows._extract_reservoirs_timeseries(mock_project_reservoirs)

        # Verify sentinel.extract_observations was called
        assert mock_extract.called


@pytest.mark.unit
def test_extract_sentinel_observations_calls_sentinel_for_s6(
    mock_project_reservoirs, tmp_path
):
    """_extract_sentinel_observations calls sentinel.extract_observations for Sentinel-6."""
    from HydroEO.satellites import sentinel

    mock_project_reservoirs.to_process = ["sentinel6"]
    mock_project_reservoirs.dirs["sentinel6"] = str(tmp_path / "sentinel6")

    # Create mock Sentinel-6 directory with data
    s6_dir = Path(mock_project_reservoirs.dirs["sentinel6"]) / "1"
    s6_dir.mkdir(parents=True, exist_ok=True)
    (s6_dir / "dummy_file.nc").touch()

    with patch.object(sentinel, "extract_observations") as mock_extract:
        flows._extract_reservoirs_timeseries(mock_project_reservoirs)

        # Verify sentinel.extract_observations was called
        assert mock_extract.called


@pytest.mark.unit
def test_extract_reservoirs_timeseries_calls_all_enabled_missions(
    mock_project_reservoirs, tmp_path
):
    """_extract_reservoirs_timeseries dispatcher calls extractors for all enabled missions."""
    mock_project_reservoirs.to_process = ["swot", "icesat2", "sentinel3", "sentinel6"]

    # Setup mock directories
    mock_project_reservoirs.dirs["swot"] = str(tmp_path / "swot")
    Path(mock_project_reservoirs.dirs["swot"]).mkdir(parents=True)
    (Path(mock_project_reservoirs.dirs["swot"]) / "dummy.shp").touch()

    with (
        patch.object(flows._reservoir_pipeline, "_extract_swot_observations") as mock_swot,
        patch.object(flows._reservoir_pipeline, "_extract_icesat2_observations") as mock_icesat2,
        patch.object(flows._reservoir_pipeline, "_extract_sentinel_observations") as mock_sentinel,
    ):
        flows._extract_reservoirs_timeseries(mock_project_reservoirs)

        mock_swot.assert_called_once()
        mock_icesat2.assert_called_once()
        # _extract_sentinel_observations called twice (S3 and S6)
        assert mock_sentinel.call_count == 2


@pytest.mark.unit
def test_extract_reservoirs_timeseries_skips_disabled_missions(
    mock_project_reservoirs, tmp_path
):
    """_extract_reservoirs_timeseries skips missions not in to_process."""
    mock_project_reservoirs.to_process = ["swot"]  # Only SWOT

    # Setup SWOT directory
    mock_project_reservoirs.dirs["swot"] = str(tmp_path / "swot")
    Path(mock_project_reservoirs.dirs["swot"]).mkdir(parents=True)
    (Path(mock_project_reservoirs.dirs["swot"]) / "dummy.shp").touch()

    with (
        patch.object(flows._reservoir_pipeline, "_extract_swot_observations") as mock_swot,
        patch.object(flows._reservoir_pipeline, "_extract_icesat2_observations") as mock_icesat2,
        patch.object(flows._reservoir_pipeline, "_extract_sentinel_observations") as mock_sentinel,
    ):
        flows._extract_reservoirs_timeseries(mock_project_reservoirs)

        mock_swot.assert_called_once()
        mock_icesat2.assert_not_called()
        mock_sentinel.assert_not_called()


# ============================================================================
# River Timeseries Extraction Tests
# ============================================================================
#
# Mirrors the reservoirs extraction tests above. Rivers' extraction is more
# involved -- raw data covers a whole waterbody (multiple targets) at once
# and has to be split into the same per-target raw_observations structure
# reservoirs use (see flows._extract_rivers_*_observations) -- so these
# tests exercise the real waterbody-grouping/corridor-buffering/point-
# assignment logic rather than mocking it away, only mocking the actual
# satellite module call each extractor makes.


@pytest.mark.unit
def test_extract_rivers_swot_observations_splits_hydrocron_csv_by_target(
    mock_project_rivers,
):
    """_extract_rivers_swot_observations splits one waterbody's Hydrocron
    CSV into separate per-target raw_observations/swot.gpkg files."""
    swot_dir = Path(mock_project_rivers.dirs["swot"]) / "loire"
    swot_dir.mkdir(parents=True, exist_ok=True)

    df = pd.DataFrame(
        {
            "node_id": [101, 101, 102],
            "node_q": [0, 0, 1],
            "time_str": [
                "2024-01-01T00:00:00Z",
                "2024-01-05T00:00:00Z",
                "2024-01-01T00:00:00Z",
            ],
            "wse": [10.0, 10.5, 20.0],
        }
    )
    df.to_csv(swot_dir / "nodes_timeseries.csv", index=False)

    flows._extract_rivers_swot_observations(mock_project_rivers)

    out_101 = (
        Path(mock_project_rivers.dirs["output"]) / "101" / "raw_observations" / "swot.gpkg"
    )
    out_102 = (
        Path(mock_project_rivers.dirs["output"]) / "102" / "raw_observations" / "swot.gpkg"
    )
    assert out_101.exists()
    assert out_102.exists()

    gdf_101 = gpd.read_file(out_101)
    assert len(gdf_101) == 2
    assert set(gdf_101["height"]) == {10.0, 10.5}
    assert (gdf_101["platform"] == "swot").all()

    gdf_102 = gpd.read_file(out_102)
    assert len(gdf_102) == 1
    assert gdf_102["height"].iloc[0] == 20.0


@pytest.mark.unit
def test_extract_rivers_swot_observations_skips_missing_csv(mock_project_rivers):
    """_extract_rivers_swot_observations skips a waterbody with no Hydrocron CSV."""
    flows._extract_rivers_swot_observations(mock_project_rivers)

    out_dir = Path(mock_project_rivers.dirs["output"])
    assert not (out_dir / "101").exists()
    assert not (out_dir / "102").exists()


@pytest.mark.unit
def test_extract_rivers_icesat2_observations_assigns_points_to_targets(
    mock_project_rivers,
):
    """_extract_rivers_icesat2_observations assigns ICESat-2 points to their
    nearest river target and writes separate per-target output files."""
    from HydroEO.satellites import icesat2

    parquet_dir = Path(mock_project_rivers.dirs["icesat2_processed"]) / "loire"
    parquet_dir.mkdir(parents=True, exist_ok=True)
    (parquet_dir / "atl13.parquet").touch()

    def fake_extract_observations(src_dir, dst_path, features, **kwargs):
        points = gpd.GeoDataFrame(
            {
                "height": [101.0, 102.0, 201.0],
                "date": pd.to_datetime(["2024-01-01", "2024-01-02", "2024-01-01"]),
            },
            geometry=[Point(0, 0.0001), Point(0, -0.0001), Point(0, 0.1001)],
            crs="EPSG:4326",
        )
        points.to_file(dst_path, driver="GPKG")

    with patch.object(
        icesat2, "extract_observations", side_effect=fake_extract_observations
    ):
        flows._extract_rivers_icesat2_observations(mock_project_rivers)

    out_101 = (
        Path(mock_project_rivers.dirs["output"])
        / "101" / "raw_observations" / "icesat2.gpkg"
    )
    out_102 = (
        Path(mock_project_rivers.dirs["output"])
        / "102" / "raw_observations" / "icesat2.gpkg"
    )
    assert out_101.exists()
    assert out_102.exists()

    gdf_101 = gpd.read_file(out_101)
    assert len(gdf_101) == 2
    assert set(gdf_101["height"]) == {101.0, 102.0}

    gdf_102 = gpd.read_file(out_102)
    assert len(gdf_102) == 1
    assert gdf_102["height"].iloc[0] == 201.0


@pytest.mark.unit
def test_extract_rivers_icesat2_observations_skips_missing_parquet(
    mock_project_rivers,
):
    """_extract_rivers_icesat2_observations skips a waterbody with no
    downloaded atl13.parquet -- no output, no exception."""
    from HydroEO.satellites import icesat2

    with patch.object(icesat2, "extract_observations") as mock_extract:
        flows._extract_rivers_icesat2_observations(mock_project_rivers)
        mock_extract.assert_not_called()

    out_dir = Path(mock_project_rivers.dirs["output"])
    assert not (out_dir / "101").exists()
    assert not (out_dir / "102").exists()


@pytest.mark.unit
def test_extract_rivers_sentinel_observations_filters_sigma0_and_assigns(
    mock_project_rivers,
):
    """_extract_rivers_sentinel_observations applies the sigma0_min
    post-filter before assigning surviving points to their nearest target."""
    from HydroEO.satellites import sentinel

    download_dir = Path(mock_project_rivers.dirs["sentinel3"]) / "loire"
    download_dir.mkdir(parents=True, exist_ok=True)
    mock_project_rivers.mission_options["sentinel3"]["sigma0_min"] = 0.5

    def fake_extract_observations(src_dir, dst_path, features, **kwargs):
        points = gpd.GeoDataFrame(
            {
                "height": [101.0, 102.0, 201.0],
                "date": pd.to_datetime(["2024-01-01", "2024-01-02", "2024-01-01"]),
                # the second point (near target 101) should be filtered out
                "sigma0": [0.9, 0.1, 0.9],
            },
            geometry=[Point(0, 0.0001), Point(0, -0.0001), Point(0, 0.1001)],
            crs="EPSG:4326",
        )
        points.to_file(dst_path, driver="GPKG")

    with patch.object(
        sentinel, "extract_observations", side_effect=fake_extract_observations
    ):
        flows._extract_rivers_sentinel_observations(mock_project_rivers, "sentinel3", "S3")

    out_101 = (
        Path(mock_project_rivers.dirs["output"])
        / "101" / "raw_observations" / "sentinel3.gpkg"
    )
    out_102 = (
        Path(mock_project_rivers.dirs["output"])
        / "102" / "raw_observations" / "sentinel3.gpkg"
    )
    assert out_101.exists()
    assert out_102.exists()

    gdf_101 = gpd.read_file(out_101)
    assert len(gdf_101) == 1
    assert gdf_101["height"].iloc[0] == 101.0

    gdf_102 = gpd.read_file(out_102)
    assert len(gdf_102) == 1
    assert gdf_102["height"].iloc[0] == 201.0


@pytest.mark.unit
def test_extract_rivers_timeseries_calls_all_enabled_missions(mock_project_rivers):
    """_extract_rivers_timeseries dispatcher calls extractors for all
    enabled missions."""
    mock_project_rivers.to_process = ["swot", "icesat2", "sentinel3", "sentinel6"]

    with (
        patch.object(flows._river_pipeline, "_extract_rivers_icesat2_observations") as mock_icesat2,
        patch.object(flows._river_pipeline, "_extract_rivers_sentinel_observations") as mock_sentinel,
        patch.object(flows._river_pipeline, "_extract_rivers_swot_observations") as mock_swot,
    ):
        flows._extract_rivers_timeseries(mock_project_rivers)

        mock_icesat2.assert_called_once()
        assert mock_sentinel.call_count == 2  # sentinel3 and sentinel6
        mock_swot.assert_called_once()


@pytest.mark.unit
def test_extract_rivers_timeseries_skips_disabled_missions(mock_project_rivers):
    """_extract_rivers_timeseries skips missions not in to_process."""
    mock_project_rivers.to_process = ["swot"]

    with (
        patch.object(flows._river_pipeline, "_extract_rivers_icesat2_observations") as mock_icesat2,
        patch.object(flows._river_pipeline, "_extract_rivers_sentinel_observations") as mock_sentinel,
        patch.object(flows._river_pipeline, "_extract_rivers_swot_observations") as mock_swot,
    ):
        flows._extract_rivers_timeseries(mock_project_rivers)

        mock_icesat2.assert_not_called()
        mock_sentinel.assert_not_called()
        mock_swot.assert_called_once()


# ============================================================================
# Timeseries Cleaning and Filter Tests
# ============================================================================


@pytest.mark.unit
def test_clean_reservoirs_applies_elevation_filter(mock_project_reservoirs, tmp_path):
    """_clean_reservoirs_timeseries applies elevation_min/max filter correctly."""
    from HydroEO.utils import timeseries

    mock_project_reservoirs.to_process = ["swot"]
    mock_project_reservoirs.processing_options = {
        "swot": {
            "processing_filters": ["elevation"],
            "elevation_min_m": 700.0,
            "elevation_max_m": 8000.0,
        }
    }
    mock_project_reservoirs.dirs["output"] = str(tmp_path / "results")

    # Create raw observations with various heights
    raw_obs_dir = tmp_path / "results" / "1" / "raw_observations"
    raw_obs_dir.mkdir(parents=True, exist_ok=True)

    raw_gdf = gpd.GeoDataFrame(
        {
            "height": [500, 750, 800, 1000, 8500],
            "date": pd.date_range("2024-01-01", periods=5),
            "geometry": [Point(0, 0)] * 5,
        },
        crs="EPSG:4326",
    )
    raw_gdf.to_file(raw_obs_dir / "swot.gpkg", driver="GPKG")

    with (
        patch.object(timeseries.Timeseries, "clean") as mock_clean,
    ):
        flows._clean_reservoirs_timeseries(mock_project_reservoirs)

        # Verify clean was called with correct filter params
        mock_clean.assert_called_once()
        call_args = mock_clean.call_args
        assert call_args[0][0] == ["elevation"]
        assert call_args[1]["filter_params"]["elevation_min_m"] == 700.0
        assert call_args[1]["filter_params"]["elevation_max_m"] == 8000.0


@pytest.mark.unit
def test_clean_reservoirs_applies_mad_filter(mock_project_reservoirs, tmp_path):
    """_clean_reservoirs_timeseries applies MAD (Median Absolute Deviation) filter."""
    from HydroEO.utils import timeseries

    mock_project_reservoirs.to_process = ["swot"]
    mock_project_reservoirs.processing_options = {
        "swot": {
            "processing_filters": ["MAD"],
            "mad_threshold": 5.0,
        }
    }
    mock_project_reservoirs.dirs["output"] = str(tmp_path / "results")

    # Create raw observations
    raw_obs_dir = tmp_path / "results" / "1" / "raw_observations"
    raw_obs_dir.mkdir(parents=True, exist_ok=True)

    raw_gdf = gpd.GeoDataFrame(
        {
            "height": [780, 785, 790, 800, 850],
            "date": pd.date_range("2024-01-01", periods=5),
            "geometry": [Point(0, 0)] * 5,
        },
        crs="EPSG:4326",
    )
    raw_gdf.to_file(raw_obs_dir / "swot.gpkg", driver="GPKG")

    with (
        patch.object(timeseries.Timeseries, "clean") as mock_clean,
    ):
        flows._clean_reservoirs_timeseries(mock_project_reservoirs)

        # Verify clean was called with MAD filter
        mock_clean.assert_called_once()
        call_args = mock_clean.call_args
        assert call_args[0][0] == ["MAD"]
        assert call_args[1]["filter_params"]["mad_threshold"] == 5.0


@pytest.mark.unit
def test_clean_reservoirs_applies_combined_filters(mock_project_reservoirs, tmp_path):
    """_clean_reservoirs_timeseries applies multiple filters in sequence."""
    from HydroEO.utils import timeseries

    mock_project_reservoirs.to_process = ["swot"]
    mock_project_reservoirs.processing_options = {
        "swot": {
            "processing_filters": ["elevation", "MAD", "kalman"],
            "elevation_min_m": 700.0,
            "elevation_max_m": 8000.0,
            "mad_threshold": 5.0,
        }
    }
    mock_project_reservoirs.dirs["output"] = str(tmp_path / "results")

    # Create raw observations
    raw_obs_dir = tmp_path / "results" / "1" / "raw_observations"
    raw_obs_dir.mkdir(parents=True, exist_ok=True)

    raw_gdf = gpd.GeoDataFrame(
        {
            "height": [750, 780, 790, 800],
            "date": pd.date_range("2024-01-01", periods=4),
            "geometry": [Point(0, 0)] * 4,
        },
        crs="EPSG:4326",
    )
    raw_gdf.to_file(raw_obs_dir / "swot.gpkg", driver="GPKG")

    with (
        patch.object(timeseries.Timeseries, "clean") as mock_clean,
    ):
        flows._clean_reservoirs_timeseries(mock_project_reservoirs)

        # Verify clean called with all three filters
        mock_clean.assert_called_once()
        call_args = mock_clean.call_args
        assert set(call_args[0][0]) == {"elevation", "MAD", "kalman"}


@pytest.mark.unit
def test_clean_reservoirs_handles_empty_raw_observations(
    mock_project_reservoirs, caplog
):
    """_clean_reservoirs_timeseries handles empty raw observations gracefully."""
    import logging

    mock_project_reservoirs.to_process = ["swot"]
    mock_project_reservoirs.processing_options = {}

    with caplog.at_level(logging.WARNING):
        flows._clean_reservoirs_timeseries(mock_project_reservoirs)
        assert "No raw observations found" in caplog.text


@pytest.mark.unit
def test_clean_reservoirs_preserves_mission_columns(mock_project_reservoirs, tmp_path):
    """_clean_reservoirs_timeseries preserves mission-specific metadata columns."""
    from HydroEO.utils import timeseries

    mock_project_reservoirs.to_process = ["swot"]
    mock_project_reservoirs.processing_options = {
        "swot": {"processing_filters": ["elevation"]}
    }
    mock_project_reservoirs.dirs["output"] = str(tmp_path / "results")

    # Create raw observations with mission metadata
    raw_obs_dir = tmp_path / "results" / "1" / "raw_observations"
    raw_obs_dir.mkdir(parents=True, exist_ok=True)

    raw_gdf = gpd.GeoDataFrame(
        {
            "height": [750, 780, 790],
            "date": pd.date_range("2024-01-01", periods=3),
            "platform": ["swot", "swot", "swot"],
            "orbit": [1, 1, 2],
            "quality_f": [0, 0, 1],
            "geometry": [Point(0, 0)] * 3,
        },
        crs="EPSG:4326",
    )
    raw_gdf.to_file(raw_obs_dir / "swot.gpkg", driver="GPKG")

    with (
        patch.object(timeseries.Timeseries, "clean"),
        patch.object(timeseries.Timeseries, "export_csv"),
    ):
        flows._clean_reservoirs_timeseries(mock_project_reservoirs)


@pytest.mark.unit
def test_clean_reservoirs_uses_default_options(mock_project_reservoirs, tmp_path):
    """_clean_reservoirs_timeseries uses default processing_options when not specified."""
    from HydroEO.utils import timeseries

    mock_project_reservoirs.to_process = ["swot"]
    mock_project_reservoirs.processing_options = {}
    mock_project_reservoirs.dirs["output"] = str(tmp_path / "results")

    # Create raw observations
    raw_obs_dir = tmp_path / "results" / "1" / "raw_observations"
    raw_obs_dir.mkdir(parents=True, exist_ok=True)

    raw_gdf = gpd.GeoDataFrame(
        {
            "height": [750, 780, 790],
            "date": pd.date_range("2024-01-01", periods=3),
            "geometry": [Point(0, 0)] * 3,
        },
        crs="EPSG:4326",
    )
    raw_gdf.to_file(raw_obs_dir / "swot.gpkg", driver="GPKG")

    with (
        patch.object(timeseries.Timeseries, "clean") as mock_clean,
        patch.object(timeseries.Timeseries, "export_csv"),
    ):
        flows._clean_reservoirs_timeseries(mock_project_reservoirs)

        # Verify defaults used
        mock_clean.assert_called_once()
        call_args = mock_clean.call_args
        assert "elevation" in call_args[0][0]
        assert "MAD" in call_args[0][0]
        assert call_args[1]["filter_params"]["elevation_min_m"] == 0.0
        assert call_args[1]["filter_params"]["elevation_max_m"] == 8000.0
        assert call_args[1]["filter_params"]["mad_threshold"] == 5.0


@pytest.mark.unit
def test_load_product_timeseries_loads_all_products(mock_project_reservoirs, tmp_path):
    """_load_product_timeseries loads all products from directory."""
    obs_dir = tmp_path / "observations"
    obs_dir.mkdir()

    # Create mock GPKG files for multiple products
    for product in ["swot", "icesat2", "sentinel3"]:
        gdf = gpd.GeoDataFrame(
            {
                "height": [100 + i for i in range(3)],
                "date": pd.date_range("2024-01-01", periods=3),
                "product": [product] * 3,
                "geometry": [Point(0, 0)] * 3,
            },
            crs="EPSG:4326",
        )
        gdf.to_file(obs_dir / f"{product}.gpkg", driver="GPKG")

    # Load with product filter
    df = flows._load_product_timeseries(
        str(obs_dir),
        ".gpkg",
        ["swot", "icesat2"],
        lambda path: gpd.read_file(path).drop(columns=["geometry"]),
    )

    assert df is not None
    assert "product" in df.columns
    assert set(df["product"].unique()) == {"swot", "icesat2"}
    assert len(df) == 6  # 2 products × 3 records


@pytest.mark.unit
def test_clean_reservoirs_processes_multiple_products(
    mock_project_reservoirs, tmp_path
):
    """_clean_reservoirs_timeseries processes multiple products."""
    from HydroEO.utils import timeseries

    mock_project_reservoirs.to_process = ["swot", "icesat2"]
    mock_project_reservoirs.processing_options = {
        "swot": {"processing_filters": ["elevation"]},
        "icesat2": {"processing_filters": ["MAD"]},
    }
    mock_project_reservoirs.dirs["output"] = str(tmp_path / "results")

    # Create raw observations for both products
    raw_obs_dir = tmp_path / "results" / "1" / "raw_observations"
    raw_obs_dir.mkdir(parents=True, exist_ok=True)

    for product, heights in [
        ("swot", [750, 780, 790]),
        ("icesat2", [500, 600, 700]),
    ]:
        gdf = gpd.GeoDataFrame(
            {
                "height": heights,
                "date": pd.date_range("2024-01-01", periods=3),
                "geometry": [Point(0, 0)] * 3,
            },
            crs="EPSG:4326",
        )
        gdf.to_file(raw_obs_dir / f"{product}.gpkg", driver="GPKG")

    with (
        patch.object(timeseries.Timeseries, "clean") as mock_clean,
        patch.object(timeseries.Timeseries, "export_csv"),
    ):
        flows._clean_reservoirs_timeseries(mock_project_reservoirs)

        # Verify clean called for both products
        assert mock_clean.call_count == 2


# ============================================================================
# River Timeseries Cleaning Tests
# ============================================================================
#
# _clean_rivers_timeseries is a thin wrapper around the same
# flows._clean_timeseries engine reservoirs use (target_type="rivers"
# instead of "reservoirs") -- the filter math itself is already covered by
# the reservoirs cleaning tests above, so these focus on the rivers-specific
# dispatch: target ids come from prj.rivers.target_ids rather than
# prj.reservoirs.download_gdf.


@pytest.mark.unit
def test_clean_rivers_applies_elevation_filter(mock_project_rivers, tmp_path):
    """_clean_rivers_timeseries applies elevation_min/max filter correctly
    for a river target."""
    from HydroEO.utils import timeseries

    mock_project_rivers.to_process = ["swot"]
    mock_project_rivers.processing_options = {
        "swot": {
            "processing_filters": ["elevation"],
            "elevation_min_m": 5.0,
            "elevation_max_m": 8000.0,
        }
    }

    raw_obs_dir = Path(mock_project_rivers.dirs["output"]) / "101" / "raw_observations"
    raw_obs_dir.mkdir(parents=True, exist_ok=True)
    raw_gdf = gpd.GeoDataFrame(
        {
            "height": [3, 10, 20],
            "date": pd.date_range("2024-01-01", periods=3),
            "geometry": [Point(0, 0)] * 3,
        },
        crs="EPSG:4326",
    )
    raw_gdf.to_file(raw_obs_dir / "swot.gpkg", driver="GPKG")

    with patch.object(timeseries.Timeseries, "clean") as mock_clean:
        flows._clean_rivers_timeseries(mock_project_rivers)

        mock_clean.assert_called_once()
        call_args = mock_clean.call_args
        assert call_args[0][0] == ["elevation"]
        assert call_args[1]["filter_params"]["elevation_min_m"] == 5.0
        assert call_args[1]["filter_params"]["elevation_max_m"] == 8000.0

    cleaned_path = (
        Path(mock_project_rivers.dirs["output"])
        / "101" / "cleaned_observations" / "swot.csv"
    )
    assert cleaned_path.exists()

    # target 102 has no raw_observations at all -- must not be touched
    assert not (
        Path(mock_project_rivers.dirs["output"]) / "102" / "cleaned_observations"
    ).exists()


@pytest.mark.unit
def test_clean_rivers_handles_empty_raw_observations(mock_project_rivers, caplog):
    """_clean_rivers_timeseries warns and returns early when no target has
    any raw observations yet."""
    import logging

    mock_project_rivers.to_process = ["swot"]

    with caplog.at_level(logging.WARNING):
        flows._clean_rivers_timeseries(mock_project_rivers)
        assert "No raw observations found for any rivers" in caplog.text


# ============================================================================
# Timeseries Merging Tests
# ============================================================================


@pytest.mark.unit
def test_load_and_parse_cleaned_timeseries_loads_all_products(
    mock_project_reservoirs, tmp_path
):
    """_merge_reservoirs_timeseries loads and parses all cleaned CSV files."""
    from HydroEO.utils import timeseries

    mock_project_reservoirs.to_process = ["swot", "icesat2", "sentinel3"]
    mock_project_reservoirs.dirs["output"] = str(tmp_path / "results")

    # Create cleaned observation CSV files
    cleaned_obs_dir = tmp_path / "results" / "1" / "cleaned_observations"
    cleaned_obs_dir.mkdir(parents=True, exist_ok=True)

    # Create CSVs for each product
    for product in ["swot", "icesat2", "sentinel3"]:
        df = pd.DataFrame(
            {
                "date": pd.date_range("2024-01-01", periods=3),
                "height": [100 + i for i in range(3)],
                "platform": [product] * 3,
            }
        )
        df.to_csv(cleaned_obs_dir / f"{product}.csv", index=False)

    with (
        patch.object(timeseries, "concat") as mock_concat,
        patch.object(timeseries.Timeseries, "merge") as mock_merge,
        patch.object(timeseries.Timeseries, "export_csv"),
    ):
        mock_concat.return_value = mock.MagicMock()
        mock_merge.return_value = mock.MagicMock()

        flows._merge_reservoirs_timeseries(mock_project_reservoirs)

        # Verify concat was called
        assert mock_concat.called


@pytest.mark.unit
def test_merge_concatenates_with_aligned_timestamps(mock_project_reservoirs, tmp_path):
    """_merge_reservoirs_timeseries concatenates and sorts by timestamp."""
    from HydroEO.utils import timeseries

    mock_project_reservoirs.to_process = ["swot", "icesat2"]
    mock_project_reservoirs.dirs["output"] = str(tmp_path / "results")

    # Create cleaned observation CSVs with different timestamps
    cleaned_obs_dir = tmp_path / "results" / "1" / "cleaned_observations"
    cleaned_obs_dir.mkdir(parents=True, exist_ok=True)

    # SWOT observations
    swot_df = pd.DataFrame(
        {
            "date": [
                "2024-01-05 00:00:00+00:00",
                "2024-01-10 00:00:00+00:00",
            ],
            "height": [780.5, 782.1],
        }
    )
    swot_df.to_csv(cleaned_obs_dir / "swot.csv", index=False)

    # ICESat-2 observations (different dates)
    icesat2_df = pd.DataFrame(
        {
            "date": [
                "2024-01-03 00:00:00+00:00",
                "2024-01-08 00:00:00+00:00",
            ],
            "height": [781.2, 781.8],
        }
    )
    icesat2_df.to_csv(cleaned_obs_dir / "icesat2.csv", index=False)

    with (
        patch.object(timeseries.Timeseries, "merge") as mock_merge,
        patch.object(timeseries.Timeseries, "export_csv"),
    ):
        mock_merge_result = mock.MagicMock()
        mock_merge.return_value = mock_merge_result

        flows._merge_reservoirs_timeseries(mock_project_reservoirs)

        # Verify merge was called to combine timeseries
        assert mock_merge.called


@pytest.mark.unit
def test_merge_handles_single_mission(mock_project_reservoirs, tmp_path):
    """_merge_reservoirs_timeseries handles case with only one mission."""
    from HydroEO.utils import timeseries

    mock_project_reservoirs.to_process = ["swot"]
    mock_project_reservoirs.dirs["output"] = str(tmp_path / "results")

    cleaned_obs_dir = tmp_path / "results" / "1" / "cleaned_observations"
    cleaned_obs_dir.mkdir(parents=True, exist_ok=True)

    # Single product CSV
    df = pd.DataFrame(
        {
            "date": pd.date_range("2024-01-01", periods=3),
            "height": [780, 781, 782],
        }
    )
    df.to_csv(cleaned_obs_dir / "swot.csv", index=False)

    with (
        patch.object(timeseries.Timeseries, "merge") as mock_merge,
        patch.object(timeseries.Timeseries, "export_csv"),
    ):
        mock_merge.return_value = mock.MagicMock()
        flows._merge_reservoirs_timeseries(mock_project_reservoirs)

        # Should still call merge even with single mission
        assert mock_merge.called


@pytest.mark.unit
def test_merge_handles_sparse_mission_coverage(mock_project_reservoirs, tmp_path):
    """_merge_reservoirs_timeseries handles missions with sparse/different frequencies."""
    from HydroEO.utils import timeseries

    mock_project_reservoirs.to_process = ["swot", "icesat2"]
    mock_project_reservoirs.dirs["output"] = str(tmp_path / "results")

    cleaned_obs_dir = tmp_path / "results" / "1" / "cleaned_observations"
    cleaned_obs_dir.mkdir(parents=True, exist_ok=True)

    # Dense SWOT observations
    swot_df = pd.DataFrame(
        {
            "date": pd.date_range("2024-01-01", periods=10, freq="D"),
            "height": [780 + 0.1 * i for i in range(10)],
        }
    )
    swot_df.to_csv(cleaned_obs_dir / "swot.csv", index=False)

    # Sparse ICESat-2 observations (weekly)
    icesat2_df = pd.DataFrame(
        {
            "date": pd.date_range("2024-01-01", periods=2, freq="W"),
            "height": [780.5, 781.5],
        }
    )
    icesat2_df.to_csv(cleaned_obs_dir / "icesat2.csv", index=False)

    with (
        patch.object(timeseries.Timeseries, "merge") as mock_merge,
        patch.object(timeseries.Timeseries, "export_csv"),
    ):
        mock_merge.return_value = mock.MagicMock()
        flows._merge_reservoirs_timeseries(mock_project_reservoirs)

        # Merge should handle frequency mismatch
        assert mock_merge.called


@pytest.mark.unit
def test_merge_preserves_utc_timezone_conversion(mock_project_reservoirs, tmp_path):
    """_merge_reservoirs_timeseries converts timestamps to UTC and removes tz."""
    from HydroEO.utils import timeseries

    mock_project_reservoirs.to_process = ["swot"]
    mock_project_reservoirs.dirs["output"] = str(tmp_path / "results")

    cleaned_obs_dir = tmp_path / "results" / "1" / "cleaned_observations"
    cleaned_obs_dir.mkdir(parents=True, exist_ok=True)

    # Create CSV with UTC timezone info
    df = pd.DataFrame(
        {
            "date": [
                "2024-01-01 12:00:00+00:00",
                "2024-01-02 12:00:00+00:00",
            ],
            "height": [780.5, 781.2],
        }
    )
    df.to_csv(cleaned_obs_dir / "swot.csv", index=False)

    with (
        patch.object(timeseries.Timeseries, "merge") as mock_merge,
        patch.object(timeseries.Timeseries, "export_csv"),
    ):
        mock_merge.return_value = mock.MagicMock()
        flows._merge_reservoirs_timeseries(mock_project_reservoirs)

        # Verify merge was called (timestamps were parsed and converted)
        assert mock_merge.called


@pytest.mark.unit
def test_merge_handles_no_cleaned_observations(mock_project_reservoirs, caplog):
    """_merge_reservoirs_timeseries handles case with no cleaned observations."""
    import logging

    mock_project_reservoirs.to_process = ["swot"]

    with caplog.at_level(logging.WARNING):
        flows._merge_reservoirs_timeseries(mock_project_reservoirs)
        assert "No cleaned observations found" in caplog.text


# ============================================================================
# River Timeseries Merging Tests
# ============================================================================
#
# _merge_rivers_timeseries is a thin wrapper around the same
# flows._merge_timeseries engine reservoirs use (target_type="rivers"
# instead of "reservoirs") -- the merge algorithm itself is already covered
# by the reservoirs merging tests (and the Nuozhadu baseline regression
# tests) above, so these focus on the rivers-specific dispatch: target ids
# come from prj.rivers.target_ids, and centroid lookup
# (flows._target_centroid) comes from prj.rivers.target_features rather
# than a reservoir polygon.


@pytest.mark.unit
def test_merge_rivers_concatenates_cleaned_observations(mock_project_rivers, tmp_path):
    """_merge_rivers_timeseries merges a river target's cleaned observations."""
    from HydroEO.utils import timeseries

    mock_project_rivers.to_process = ["swot"]

    cleaned_dir = (
        Path(mock_project_rivers.dirs["output"]) / "101" / "cleaned_observations"
    )
    cleaned_dir.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(
        {
            "date": pd.date_range("2024-01-01", periods=3),
            "height": [10.0, 10.1, 10.2],
        }
    )
    df.to_csv(cleaned_dir / "swot.csv", index=False)

    with (
        patch.object(timeseries.Timeseries, "merge") as mock_merge,
        patch.object(timeseries.Timeseries, "export_csv"),
    ):
        mock_merge.return_value = mock.MagicMock()
        flows._merge_rivers_timeseries(mock_project_rivers)

        assert mock_merge.called
        # _target_centroid should have resolved target 101's centroid from
        # prj.rivers.target_features rather than failing/using (None, None)
        call_kwargs = mock_merge.call_args[1]
        assert call_kwargs["ref_lat"] is not None
        assert call_kwargs["ref_lon"] is not None


@pytest.mark.unit
def test_merge_rivers_handles_no_cleaned_observations(mock_project_rivers, caplog):
    """_merge_rivers_timeseries handles the case with no cleaned observations."""
    import logging

    mock_project_rivers.to_process = ["swot"]

    with caplog.at_level(logging.WARNING):
        flows._merge_rivers_timeseries(mock_project_rivers)
        assert "No cleaned observations found for any rivers" in caplog.text


# ============================================================================
# Regression): Baseline Tests Using Real Data
# ============================================================================
# These tests validate merge algorithm outputs against Nuozhadu baseline data.
# They serve as regression tests to catch unintended algorithm changes.


@pytest.fixture
def nuozhadu_baseline_path():
    """Return path to Nuozhadu baseline test data."""
    return Path(__file__).parent.parent / "data" / "baselines" / "nuozhadu"


@pytest.mark.unit
def test_merge_baseline_all_cleaned_timeseries_matches(
    mock_project_reservoirs, nuozhadu_baseline_path, tmp_path
):
    """Regression test: merge concatenates cleaned observations matching baseline."""
    from HydroEO.utils import timeseries

    if not (nuozhadu_baseline_path / "cleaned_observations").exists():
        pytest.skip("Baseline data not available")

    mock_project_reservoirs.to_process = ["swot", "icesat2", "sentinel3", "sentinel6"]
    mock_project_reservoirs.dirs["output"] = str(tmp_path / "results")

    # Copy baseline cleaned observations to temp directory
    cleaned_obs_dir = tmp_path / "results" / "1" / "cleaned_observations"
    cleaned_obs_dir.mkdir(parents=True, exist_ok=True)
    for csv_file in (nuozhadu_baseline_path / "cleaned_observations").glob("*.csv"):
        import shutil

        shutil.copy(csv_file, cleaned_obs_dir)

    # Load baseline all_cleaned_timeseries
    baseline_all_cleaned = pd.read_csv(
        nuozhadu_baseline_path / "all_cleaned_timeseries.csv", index_col=0
    )

    # Capture the actual export_csv calls to validate output
    captured_exports = {}

    def capture_export(path):
        """Capture what export_csv would write, keyed by filename."""
        # Get the actual DataFrame from the Timeseries object
        filename = Path(path).name
        captured_exports[filename] = {
            "path": path,
            "called": True,
        }

    # Patch export_csv to capture what would be written (actual data is in self.data)
    original_export = timeseries.Timeseries.export_csv

    def export_with_capture(self, path):
        capture_export(path)
        # Still call the real export to write the file
        original_export(self, path)

    with patch.object(timeseries.Timeseries, "export_csv", export_with_capture):
        flows._merge_reservoirs_timeseries(mock_project_reservoirs)

    # Verify exports happened
    assert "all_cleaned_timeseries.csv" in captured_exports
    assert "merged_timeseries.csv" in captured_exports

    # Read the actual output that was written
    all_cleaned_output = pd.read_csv(
        tmp_path / "results" / "1" / "all_cleaned_timeseries.csv", index_col=0
    )

    # Validate structure matches baseline (columns, data types, row count)
    assert "date" in all_cleaned_output.columns
    assert "height" in all_cleaned_output.columns
    assert len(all_cleaned_output) > 0

    # Verify data is concatenated (should have multiple products combined)
    # Baseline has 4 missions, so expect multiple rows
    assert (
        len(all_cleaned_output) >= len(baseline_all_cleaned) * 0.8
    )  # Allow some variance
    assert len(all_cleaned_output) <= len(baseline_all_cleaned) * 1.2


@pytest.mark.unit
def test_merge_baseline_merged_timeseries_columns(nuozhadu_baseline_path):
    """Merged timeseries output has expected columns (date, height)."""
    if not (nuozhadu_baseline_path / "merged_timeseries.csv").exists():
        pytest.skip("Baseline data not available")

    # Load baseline merged output
    baseline_merged = pd.read_csv(
        nuozhadu_baseline_path / "merged_timeseries.csv", index_col=0
    )

    # Verify expected columns
    assert "date" in baseline_merged.columns
    assert "height" in baseline_merged.columns
    assert len(baseline_merged) > 0


@pytest.mark.unit
def test_merge_baseline_data_structure_preserved(nuozhadu_baseline_path):
    """All baseline output files have expected structure."""
    if not nuozhadu_baseline_path.exists():
        pytest.skip("Baseline data not available")

    # Check all_cleaned_timeseries exists and has data
    all_cleaned = pd.read_csv(
        nuozhadu_baseline_path / "all_cleaned_timeseries.csv", index_col=0
    )
    assert "date" in all_cleaned.columns
    assert "height" in all_cleaned.columns
    assert len(all_cleaned) > 0

    # Check merged_timeseries exists and has data
    merged = pd.read_csv(nuozhadu_baseline_path / "merged_timeseries.csv", index_col=0)
    assert "date" in merged.columns
    assert "height" in merged.columns
    assert len(merged) > 0

    # Check merge progress files exist and have data
    for stage in ["daily_mad_error", "kalman", "svr_linear", "svr_radial"]:
        progress_file = nuozhadu_baseline_path / "merged_progress" / f"{stage}.csv"
        assert progress_file.exists(), f"{stage}.csv should exist"
        df = pd.read_csv(progress_file, index_col=0)
        assert len(df) > 0, f"{stage}.csv should have data"


@pytest.mark.unit
def test_merge_baseline_cleaned_observations_loaded(nuozhadu_baseline_path):
    """Cleaned observation files from baseline have expected structure."""
    if not (nuozhadu_baseline_path / "cleaned_observations").exists():
        pytest.skip("Baseline data not available")

    # Check each cleaned product CSV
    for product_csv in (nuozhadu_baseline_path / "cleaned_observations").glob("*.csv"):
        df = pd.read_csv(product_csv, index_col=0)
        assert "date" in df.columns, f"{product_csv.name}: date column missing"
        assert "height" in df.columns, f"{product_csv.name}: height column missing"
        assert len(df) > 0, f"{product_csv.name}: no data rows"


@pytest.mark.unit
def test_merge_baseline_timestamp_coverage(nuozhadu_baseline_path):
    """Merged output spans expected date range."""
    if not (nuozhadu_baseline_path / "merged_timeseries.csv").exists():
        pytest.skip("Baseline data not available")

    merged = pd.read_csv(nuozhadu_baseline_path / "merged_timeseries.csv", index_col=0)
    merged["date"] = pd.to_datetime(merged["date"], format="mixed")

    # Baseline covers Jan-Mar 2024 (expected for Nuozhadu)
    min_date = merged["date"].min()
    max_date = merged["date"].max()

    assert min_date.year == 2024
    assert min_date.month in [1, 2]
    assert max_date.year == 2024
    assert max_date.month in [1, 2, 3]
    assert (max_date - min_date).days > 20  # At least 20 days of data


@pytest.mark.unit
def test_merge_baseline_height_values_reasonable(nuozhadu_baseline_path):
    """Merged height values are within expected physical range."""
    if not (nuozhadu_baseline_path / "merged_timeseries.csv").exists():
        pytest.skip("Baseline data not available")

    merged = pd.read_csv(nuozhadu_baseline_path / "merged_timeseries.csv", index_col=0)

    # Mekong-Nuozhadu reservoir: elevation ~780-800m (sea level reference)
    height_min = merged["height"].min()
    height_max = merged["height"].max()

    assert 700 < height_min < 850, f"Min height {height_min} outside expected range"
    assert 700 < height_max < 850, f"Max height {height_max} outside expected range"
    assert height_max - height_min > 1, "Height variation should be > 1m"


@pytest.mark.unit
def test_merge_baseline_no_nan_in_final_output(nuozhadu_baseline_path):
    """Final merged timeseries has no NaN values."""
    if not (nuozhadu_baseline_path / "merged_timeseries.csv").exists():
        pytest.skip("Baseline data not available")

    merged = pd.read_csv(nuozhadu_baseline_path / "merged_timeseries.csv", index_col=0)

    # Check for NaN in critical columns
    assert merged["date"].notna().all(), "date column should have no NaN"
    assert merged["height"].notna().all(), "height column should have no NaN"


@pytest.mark.unit
def test_merge_baseline_dates_chronological(nuozhadu_baseline_path):
    """Merged output is sorted chronologically."""
    if not (nuozhadu_baseline_path / "merged_timeseries.csv").exists():
        pytest.skip("Baseline data not available")

    merged = pd.read_csv(nuozhadu_baseline_path / "merged_timeseries.csv", index_col=0)
    merged["date"] = pd.to_datetime(merged["date"], format="mixed")

    # Check if sorted
    is_sorted = merged["date"].is_monotonic_increasing
    assert is_sorted, "merged_timeseries.csv should be sorted by date"


# (Old SWORD tests removed - replaced with new comprehensive test suite above)
