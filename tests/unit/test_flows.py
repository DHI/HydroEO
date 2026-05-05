"""Unit tests for flows.py orchestration functions."""

import datetime
import os
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
        "output": str(tmp_path / "reservoirs"),
        "swot": str(tmp_path / "swot"),
        "icesat2": str(tmp_path / "icesat2"),
        "sentinel3": str(tmp_path / "sentinel3"),
        "sentinel6": str(tmp_path / "sentinel6"),
        "pld": str(tmp_path / "pld.shp"),
    }
    prj.reservoirs = Reservoirs(gdf=gdf, id_key="id", dirs=prj.dirs)
    prj.reservoirs.download_gdf = gdf
    prj.local_crs = "EPSG:3857"
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
    """Create a mock Project with Rivers configuration."""
    gdf = gpd.GeoDataFrame({"geometry": []}, geometry="geometry", crs="EPSG:4326")
    prj = SimpleNamespace()
    prj.dirs = {
        "main": str(tmp_path),
        "output": str(tmp_path / "rivers"),
        "swot": str(tmp_path / "swot"),
        "sword": str(tmp_path / "SWORD_v17b_gpkg" / "gpkg"),
    }
    prj.rivers = Rivers(gdf=gdf, id_key="river_id", dirs=prj.dirs)
    prj.rivers.target_ids = [101, 102]
    prj.rivers.target_id_col = "node_id"
    prj.rivers.target_features = None
    prj.rivers.configured_id = "loire"
    prj.rivers.input_mode = "configured_id"
    prj.local_crs = "EPSG:3857"
    prj.startdates = {"swot": [2024, 1, 1]}
    prj.enddates = {"swot": [2024, 2, 1]}
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
        }
    }
    return prj


# ============================================================================
# Initialization Tests
# ============================================================================


@pytest.mark.unit
def test_initialize_reservoirs_with_pld_download(mock_project_reservoirs, monkeypatch):
    """initialize_reservoirs calls _download_pld when enabled."""
    mock_project_reservoirs.to_download = ["swot"]
    mock_project_reservoirs.to_process = []

    with (
        patch.object(flows, "_download_pld") as mock_download,
        patch.object(flows, "_assign_pld_id") as mock_assign,
        patch.object(flows, "_flag_missing_priors") as mock_flag,
    ):
        flows.initialize_reservoirs(mock_project_reservoirs)

        # Verify all helper functions are called
        mock_download.assert_called_once_with(mock_project_reservoirs)
        mock_assign.assert_called_once_with(mock_project_reservoirs)
        mock_flag.assert_called_once_with(mock_project_reservoirs)


@pytest.mark.unit
def test_initialize_rivers_aoi_branch(mock_project_rivers):
    """initialize_rivers calls _prepare_rivers_from_sword for aoi_path mode."""
    mock_project_rivers.rivers.input_mode = "aoi_path"
    mock_project_rivers.rivers.aoi_gdf = gpd.GeoDataFrame(
        {"river_id": ["a"], "geometry": [box(0, 0, 1, 1)]},
        crs="EPSG:4326",
    )
    mock_project_rivers.rivers.continent_key = "eu"
    mock_project_rivers.rivers.feature_type = "nodes"
    mock_project_rivers.rivers.buffer_meters = 500

    with patch.object(flows, "_prepare_rivers_from_sword") as mock_prepare:
        flows.initialize_rivers(mock_project_rivers)
        mock_prepare.assert_called_once_with(mock_project_rivers)


@pytest.mark.unit
def test_initialize_rivers_configured_id_branch(mock_project_rivers):
    """initialize_rivers skips SWORD for configured_id mode."""
    mock_project_rivers.rivers.input_mode = "configured_id"

    with patch.object(flows, "_prepare_rivers_from_sword") as mock_prepare:
        flows.initialize_rivers(mock_project_rivers)
        mock_prepare.assert_not_called()


# ============================================================================
# Download Orchestration Tests
# ============================================================================


@pytest.mark.unit
def test_download_reservoirs_dispatches_swot(mock_project_reservoirs):
    """download_reservoirs calls _download_reservoirs_swot when enabled."""
    mock_project_reservoirs.to_download = ["swot"]

    with (
        patch.object(flows, "_download_reservoirs_swot") as mock_swot,
        patch.object(flows, "_download_reservoirs_icesat2") as mock_ice,
        patch.object(flows, "_download_reservoirs_sentinel") as mock_sent,
    ):
        flows.download_reservoirs(mock_project_reservoirs)

        mock_swot.assert_called_once_with(mock_project_reservoirs)
        mock_ice.assert_not_called()
        mock_sent.assert_not_called()


@pytest.mark.unit
def test_download_reservoirs_dispatches_all_missions(mock_project_reservoirs):
    """download_reservoirs dispatches all enabled satellite missions."""
    mock_project_reservoirs.to_download = [
        "swot",
        "icesat2",
        "sentinel3",
        "sentinel6",
    ]

    with (
        patch.object(flows, "_download_reservoirs_swot") as mock_swot,
        patch.object(flows, "_download_reservoirs_icesat2") as mock_ice,
        patch.object(flows, "_download_reservoirs_sentinel") as mock_sent,
    ):
        flows.download_reservoirs(mock_project_reservoirs)

        mock_swot.assert_called_once()
        mock_ice.assert_called_once()
        assert mock_sent.call_count == 2  # sentinel3 and sentinel6


@pytest.mark.unit
def test_download_reservoirs_skips_disabled_missions(mock_project_reservoirs):
    """download_reservoirs skips missions not in to_download."""
    mock_project_reservoirs.to_download = ["swot"]

    with (
        patch.object(flows, "_download_reservoirs_swot") as mock_swot,
        patch.object(flows, "_download_reservoirs_icesat2") as mock_ice,
    ):
        flows.download_reservoirs(mock_project_reservoirs)

        mock_swot.assert_called_once()
        mock_ice.assert_not_called()


@pytest.mark.unit
def test_download_rivers_calls_hydrocron(mock_project_rivers):
    """download_rivers calls _download_swot_hydrocron_timeseries."""
    mock_project_rivers.to_download = ["swot"]

    with patch.object(flows, "_download_swot_hydrocron_timeseries") as mock_hydrocron:
        flows.download_rivers(mock_project_rivers)

        mock_hydrocron.assert_called_once()
        call_args = mock_hydrocron.call_args
        assert call_args[0][0] == mock_project_rivers
        assert isinstance(call_args[0][1], datetime.date)
        assert isinstance(call_args[0][2], datetime.date)


@pytest.mark.unit
def test_download_rivers_skips_when_swot_not_in_to_download(mock_project_rivers):
    """download_rivers returns early if swot not in to_download."""
    mock_project_rivers.to_download = ["icesat2"]

    with patch.object(flows, "_download_swot_hydrocron_timeseries") as mock_hydrocron:
        flows.download_rivers(mock_project_rivers)
        mock_hydrocron.assert_not_called()


# ============================================================================
# Timeseries Processing Tests
# ============================================================================


@pytest.mark.unit
def test_create_reservoirs_timeseries_orchestrates_steps(mock_project_reservoirs):
    """create_reservoirs_timeseries calls extract, clean, and merge in order."""
    mock_project_reservoirs.to_process = ["swot", "icesat2"]
    mock_project_reservoirs.processing_options = {}

    call_order = []

    def track_extract(prj):
        call_order.append("extract")

    def track_clean(prj):
        call_order.append("clean")

    def track_merge(prj):
        call_order.append("merge")

    with (
        patch.object(
            flows, "_extract_reservoirs_timeseries", side_effect=track_extract
        ),
        patch.object(flows, "_clean_reservoirs_timeseries", side_effect=track_clean),
        patch.object(flows, "_merge_reservoirs_timeseries", side_effect=track_merge),
    ):
        flows.create_reservoirs_timeseries(mock_project_reservoirs)

        # Verify all called and in order
        assert call_order == ["extract", "clean", "merge"]


@pytest.mark.unit
def test_create_reservoirs_timeseries_calls_export_when_toggled(
    mock_project_reservoirs,
):
    """create_reservoirs_timeseries calls _export_cleaned_to_dfs0 when enabled."""
    mock_project_reservoirs.to_process = ["swot"]
    mock_project_reservoirs.processing_options = {}
    mock_project_reservoirs.reservoirs.export_to_dfs0 = True

    call_order = []

    def track_extract(prj):
        call_order.append("extract")

    def track_clean(prj):
        call_order.append("clean")

    def track_export(prj):
        call_order.append("export")

    def track_merge(prj):
        call_order.append("merge")

    with (
        patch.object(
            flows, "_extract_reservoirs_timeseries", side_effect=track_extract
        ),
        patch.object(flows, "_clean_reservoirs_timeseries", side_effect=track_clean),
        patch.object(flows, "_export_cleaned_to_dfs0", side_effect=track_export),
        patch.object(flows, "_merge_reservoirs_timeseries", side_effect=track_merge),
    ):
        flows.create_reservoirs_timeseries(mock_project_reservoirs)

        # Verify export is called in correct position
        assert call_order == ["extract", "clean", "export", "merge"]


@pytest.mark.unit
def test_create_reservoirs_timeseries_skips_export_when_disabled(
    mock_project_reservoirs,
):
    """create_reservoirs_timeseries skips export when toggle is false."""
    mock_project_reservoirs.to_process = ["swot"]
    mock_project_reservoirs.processing_options = {}
    mock_project_reservoirs.reservoirs.export_to_dfs0 = False

    with (
        patch.object(flows, "_extract_reservoirs_timeseries"),
        patch.object(flows, "_clean_reservoirs_timeseries"),
        patch.object(flows, "_export_cleaned_to_dfs0") as mock_export,
        patch.object(flows, "_merge_reservoirs_timeseries"),
    ):
        flows.create_reservoirs_timeseries(mock_project_reservoirs)

        # Verify export is not called
        mock_export.assert_not_called()


@pytest.mark.unit
def test_export_cleaned_to_dfs0_writes_files(
    mock_project_reservoirs, tmp_path, monkeypatch
):
    """_export_cleaned_to_dfs0 writes dfs0 files for each product."""
    # Create mock cleaned_observations directory structure
    mock_project_reservoirs.to_process = ["swot", "icesat2"]
    output_dir = tmp_path / "reservoirs" / "1" / "cleaned_observations"
    output_dir.mkdir(parents=True, exist_ok=True)
    mock_project_reservoirs.dirs["output"] = str(tmp_path / "reservoirs")

    # Create sample CSVs
    dates = pd.date_range("2024-01-01", periods=3, freq="D")
    for product in ["swot", "icesat2"]:
        csv_path = output_dir / f"{product}.csv"
        df = pd.DataFrame(
            {
                "date": dates,
                "height": [100.5, 101.2, 100.8],
            }
        )
        df.to_csv(csv_path, index=False)

    # Mock mikeio module functions
    mock_ds = mock.MagicMock()

    with (
        patch.object(flows.mikeio, "from_pandas", return_value=mock_ds),
        patch.object(flows.mikeio, "ItemInfo", return_value=mock.MagicMock()),
        patch.object(
            flows.mikeio, "EUMType", mock.MagicMock(Water_Level="Water_Level")
        ),
    ):
        flows._export_cleaned_to_dfs0(mock_project_reservoirs)

        # Verify from_pandas was called for each product
        assert flows.mikeio.from_pandas.call_count >= 2
        # Verify to_dfs was called for each
        assert mock_ds.to_dfs.call_count >= 2


@pytest.mark.unit
def test_export_cleaned_to_dfs0_skips_missing_csv(mock_project_reservoirs, tmp_path):
    """_export_cleaned_to_dfs0 skips products with no CSV file."""
    mock_project_reservoirs.to_process = ["swot", "icesat2"]
    output_dir = tmp_path / "reservoirs" / "1" / "cleaned_observations"
    output_dir.mkdir(parents=True, exist_ok=True)
    mock_project_reservoirs.dirs["output"] = str(tmp_path / "reservoirs")

    # Create only one CSV (swot missing)
    dates = pd.date_range("2024-01-01", periods=3, freq="D")
    df = pd.DataFrame(
        {
            "date": dates,
            "height": [100.5, 101.2, 100.8],
        }
    )
    (output_dir / "icesat2.csv").write_text(df.to_csv(index=False))

    mock_ds = mock.MagicMock()

    with (
        patch.object(flows.mikeio, "from_pandas", return_value=mock_ds),
        patch.object(flows.mikeio, "ItemInfo", return_value=mock.MagicMock()),
        patch.object(
            flows.mikeio, "EUMType", mock.MagicMock(Water_Level="Water_Level")
        ),
    ):
        flows._export_cleaned_to_dfs0(mock_project_reservoirs)

        # Verify from_pandas was called only for icesat2
        assert flows.mikeio.from_pandas.call_count == 1


@pytest.mark.unit
def test_generate_reservoirs_summaries_iterates_ids(mock_project_reservoirs):
    """generate_reservoirs_summaries loops over download_gdf IDs."""
    mock_project_reservoirs.reservoirs.download_gdf = pd.DataFrame(
        {"id": [1, 2], "geometry": [None, None]}
    )

    with (
        patch.object(flows, "_load_product_timeseries"),
        patch("HydroEO.flows.plotting.plot_crossings") as mock_plot,
        patch("HydroEO.flows.plotting.plot_cleaning"),
        patch("HydroEO.flows.plotting.plot_merging"),
    ):
        flows.generate_reservoirs_summaries(
            mock_project_reservoirs, show=False, save=False
        )

        # Verify plotting functions were called
        assert mock_plot.call_count >= 2


# ============================================================================
# Helper Function Tests
# ============================================================================


@pytest.mark.unit
def test_group_river_targets_by_waterbody_from_features(mock_project_rivers):
    """_group_river_targets_by_waterbody uses target_features when available."""
    mock_project_rivers.rivers.target_features = gpd.GeoDataFrame(
        {
            "node_id": [101, 102, 201],
            "river_id": ["Loire", "Loire", "Rhine"],
            "geometry": [None, None, None],
        }
    )
    mock_project_rivers.rivers.target_id_col = "node_id"
    mock_project_rivers.rivers.id_key = "river_id"

    result = flows._group_river_targets_by_waterbody(mock_project_rivers)

    assert result == {"Loire": [101, 102], "Rhine": [201]}


@pytest.mark.unit
def test_group_river_targets_by_waterbody_from_configured_id(mock_project_rivers):
    """_group_river_targets_by_waterbody falls back to configured_id."""
    mock_project_rivers.rivers.target_features = None
    mock_project_rivers.rivers.configured_id = "loire"
    mock_project_rivers.rivers.target_ids = [101, 102]

    result = flows._group_river_targets_by_waterbody(mock_project_rivers)

    assert result == {"loire": [101, 102]}


@pytest.mark.unit
def test_group_river_targets_by_waterbody_raises_when_unconfigured(mock_project_rivers):
    """_group_river_targets_by_waterbody raises ValueError when no config."""
    mock_project_rivers.rivers.target_features = None
    mock_project_rivers.rivers.configured_id = None

    with pytest.raises(ValueError, match="Unable to group river targets"):
        flows._group_river_targets_by_waterbody(mock_project_rivers)


@pytest.mark.unit
def test_get_latest_hydrocron_obs_date_returns_none_when_missing(tmp_path):
    """_get_latest_hydrocron_obs_date returns None when file doesn't exist."""
    missing_path = tmp_path / "nonexistent.csv"

    result = flows._get_latest_hydrocron_obs_date(str(missing_path))

    assert result is None


@pytest.mark.unit
def test_get_latest_hydrocron_obs_date_parses_existing_csv(tmp_path):
    """_get_latest_hydrocron_obs_date extracts latest timestamp from CSV."""
    csv_path = tmp_path / "nodes_timeseries.csv"
    df = pd.DataFrame(
        {
            "time_str": [
                "2024-01-01T00:00:00Z",
                "2024-01-05T00:00:00Z",
                "2024-01-03T00:00:00Z",
            ]
        }
    )
    df.to_csv(csv_path, index=False)

    result = flows._get_latest_hydrocron_obs_date(str(csv_path))

    assert result == datetime.date(2024, 1, 5)


@pytest.mark.unit
def test_get_latest_hydrocron_obs_date_handles_invalid_csv(tmp_path):
    """_get_latest_hydrocron_obs_date handles malformed CSV gracefully."""
    csv_path = tmp_path / "broken.csv"
    csv_path.write_text("this is not valid csv, broken here\n")

    result = flows._get_latest_hydrocron_obs_date(str(csv_path))

    assert result is None


@pytest.mark.unit
def test_get_latest_hydrocron_obs_date_handles_missing_time_str_column(tmp_path):
    """_get_latest_hydrocron_obs_date returns None if time_str column missing."""
    csv_path = tmp_path / "no_time_str.csv"
    df = pd.DataFrame({"node_id": [101, 102], "wse": [10.5, 11.0]})
    df.to_csv(csv_path, index=False)

    result = flows._get_latest_hydrocron_obs_date(str(csv_path))

    assert result is None


# ============================================================================
# PLD Workflow Tests
# ============================================================================


@pytest.mark.unit
def test_download_pld_skips_when_exists(mock_project_reservoirs, tmp_path, caplog):
    """_download_pld skips download when PLD file exists."""
    import logging

    pld_path = mock_project_reservoirs.dirs["pld"]
    Path(pld_path).parent.mkdir(parents=True, exist_ok=True)
    Path(pld_path).touch()  # Create the file

    with (
        caplog.at_level(logging.INFO),
        patch("HydroEO.downloaders.hydroweb.download_PLD") as mock_dl,
    ):
        flows._download_pld(mock_project_reservoirs)
        mock_dl.assert_not_called()
        assert "PLD located" in caplog.text


@pytest.mark.unit
def test_download_pld_downloads_when_missing(mock_project_reservoirs, caplog):
    """_download_pld downloads PLD when file doesn't exist."""
    import logging

    with (
        caplog.at_level(logging.INFO),
        patch("HydroEO.downloaders.hydroweb.download_PLD") as mock_dl,
    ):
        flows._download_pld(mock_project_reservoirs)
        mock_dl.assert_called_once()
        assert "Downloading PLD" in caplog.text


@pytest.mark.unit
def test_assign_pld_id_updates_gdf(mock_project_reservoirs):
    """_assign_pld_id joins PLD data and updates reservoirs gdf."""
    # Mock PLD GeoDataFrame
    pld_gdf = gpd.GeoDataFrame(
        {
            "lake_id": [1001, 1002],
            "res_id": [501, 502],
            "geometry": [Point(0.5, 0.5), Point(1.5, 1.5)],
        },
        crs="EPSG:4326",
    )

    with (
        patch("geopandas.read_file", return_value=pld_gdf),
        patch("geopandas.sjoin_nearest") as mock_sjoin,
    ):
        # Mock the sjoin result
        joined_gdf = mock_project_reservoirs.reservoirs.gdf.copy()
        joined_gdf["prior_lake_id"] = [1001, 1002]
        joined_gdf["prior_res_id"] = [501, 502]
        mock_sjoin.return_value = joined_gdf

        flows._assign_pld_id(mock_project_reservoirs)

        mock_sjoin.assert_called_once()


# ============================================================================
# SWORD/Rivers Workflow Tests
# ============================================================================


@pytest.mark.unit
def test_ensure_sword_database_returns_existing_path(tmp_path):
    """_ensure_sword_database returns path if SWORD database exists."""
    mock_prj = SimpleNamespace()
    sword_dir = tmp_path / "SWORD_v17b_gpkg" / "gpkg"
    sword_dir.mkdir(parents=True)
    mock_prj.dirs = {"main": str(tmp_path)}

    result = flows._ensure_sword_database(mock_prj)

    assert result == str(sword_dir)
    assert os.path.isdir(result)


@pytest.mark.unit
def test_ensure_sword_database_downloads_when_missing(tmp_path, caplog):
    """_ensure_sword_database downloads SWORD when database doesn't exist."""
    import logging

    mock_prj = SimpleNamespace()
    mock_prj.dirs = {"main": str(tmp_path)}

    with (
        caplog.at_level(logging.WARNING),
        patch("urllib.request.urlretrieve") as mock_retrieve,
        patch("zipfile.ZipFile"),
        patch("os.remove"),
    ):
        flows._ensure_sword_database(mock_prj)

        mock_retrieve.assert_called_once()
        assert "Downloading" in caplog.text or "SWORD_v17b" in caplog.text


@pytest.mark.unit
def test_prepare_rivers_from_sword_buffers_aoi(mock_project_rivers, tmp_path):
    """_prepare_rivers_from_sword buffers AOI in local CRS."""
    sword_gpkg = gpd.GeoDataFrame(
        {
            "node_id": [101, 102],
            "geometry": [Point(0.5, 0.5), Point(1.5, 1.5)],
        },
        crs="EPSG:4326",
    )

    mock_project_rivers.rivers.aoi_gdf = gpd.GeoDataFrame(
        {
            "river_id": ["Loire"],
            "geometry": [box(0, 0, 1, 1)],
        },
        crs="EPSG:4326",
    )
    mock_project_rivers.rivers.continent_key = "eu"
    mock_project_rivers.rivers.feature_type = "nodes"
    mock_project_rivers.rivers.buffer_meters = 1000

    sword_dir = tmp_path / "SWORD_v17b_gpkg" / "gpkg"
    sword_dir.mkdir(parents=True)
    gpkg_path = sword_dir / "eu_sword_nodes_v17b.gpkg"
    gpkg_path.touch()  # Create fake file

    with (
        patch.object(flows, "_ensure_sword_database", return_value=str(sword_dir)),
        patch("geopandas.read_file", return_value=sword_gpkg),
        patch("geopandas.sjoin", return_value=sword_gpkg) as mock_sjoin,
    ):
        flows._prepare_rivers_from_sword(mock_project_rivers)

        mock_sjoin.assert_called_once()
        assert mock_project_rivers.rivers.target_id_col == "node_id"
