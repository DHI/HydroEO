"""Unit tests for internal flow orchestration modules."""

import pandas as pd
import pytest
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock

from HydroEO.flows import (
    ReservoirDownloadFlow,
    RiverDownloadFlow,
    SWOTRasterDownloadFlow,
    PlottingFlow,
    PreprocessFlow,
)


class DummyReservoirs:
    def __init__(self):
        self.calls = []
        self.id_key = "rid"
        self.download_gdf = pd.DataFrame({"rid": ["A", "B"]})

    def download_altimetry(self, **kwargs):
        self.calls.append(("download_altimetry", kwargs))

    def download_swot_hydrocron(self, **kwargs):
        self.calls.append(("download_swot_hydrocron", kwargs))
        return {
            "requested": 1,
            "successful": 1,
            "failed": 0,
            "empty_after_filter": 0,
        }

    def extract_product_timeseries(self, products):
        self.calls.append(("extract_product_timeseries", products))

    def clean_product_timeseries(self, products, filter_options_by_product):
        self.calls.append(
            (
                "clean_product_timeseries",
                {
                    "products": products,
                    "filter_options_by_product": filter_options_by_product,
                },
            )
        )

    def merge_product_timeseries(self, products):
        self.calls.append(("merge_product_timeseries", products))

    def summarize_crossings_by_id(self, rid, show, save):
        self.calls.append(("summarize_crossings_by_id", rid, show, save))

    def summarize_cleaning_by_id(self, rid, show, save):
        self.calls.append(("summarize_cleaning_by_id", rid, show, save))

    def summarize_merging_by_id(self, rid, show, save):
        self.calls.append(("summarize_merging_by_id", rid, show, save))


@pytest.mark.unit
def test_download_flow_dispatches_products_and_credentials():
    reservoirs = DummyReservoirs()

    flow = ReservoirDownloadFlow(
        reservoirs=reservoirs,
        to_download=["swot", "icesat2", "sentinel3", "sentinel6"],
        startdates={
            "swot": [2024, 1, 1],
            "icesat2": [2024, 1, 2],
            "sentinel3": [2024, 1, 3],
            "sentinel6": [2024, 1, 4],
        },
        enddates={
            "swot": [2024, 2, 1],
            "icesat2": [2024, 2, 2],
            "sentinel3": [2024, 2, 3],
            "sentinel6": [2024, 2, 4],
        },
        earthdata_credentials=("edl-user", "edl-pass"),
        creodias_credentials_provider=lambda: ("creo-user", "creo-pass"),
    )

    flow.run(update_existing=True, enddate_overrides={"swot": [2024, 3, 1]})

    calls = [c for c in reservoirs.calls if c[0] == "download_altimetry"]
    assert len(calls) == 4

    swot = calls[0][1]
    assert swot["product"] == "SWOT_LAKE"
    assert swot["enddate"] == [2024, 3, 1]
    assert swot["update_existing"] is True

    icesat2 = calls[1][1]
    assert icesat2["product"] == "ATL13"
    assert icesat2["credentials"] == ("edl-user", "edl-pass")

    s3 = calls[2][1]
    assert s3["product"] == "S3"
    assert s3["credentials"] == ("creo-user", "creo-pass")

    s6 = calls[3][1]
    assert s6["product"] == "S6"
    assert s6["credentials"] == ("creo-user", "creo-pass")


@pytest.mark.unit
def test_preprocess_flow_runs_extract_clean_merge_in_order():
    reservoirs = DummyReservoirs()

    flow = PreprocessFlow(
        reservoirs=reservoirs,
        to_process=["swot", "icesat2"],
        processing_options={"swot": {"processing_filters": ["elevation"]}},
    )

    flow.run()

    assert reservoirs.calls[0] == ("extract_product_timeseries", ["swot", "icesat2"])
    assert reservoirs.calls[1][0] == "clean_product_timeseries"
    assert reservoirs.calls[2] == ("merge_product_timeseries", ["swot", "icesat2"])


@pytest.mark.unit
def test_plotting_flow_runs_all_summary_steps_for_each_reservoir():
    reservoirs = DummyReservoirs()

    flow = PlottingFlow(reservoirs=reservoirs)
    flow.run(show=False, save=True)

    summarize_calls = [c for c in reservoirs.calls if c[0].startswith("summarize_")]
    assert len(summarize_calls) == 6

    expected = [
        ("summarize_crossings_by_id", "A", False, True),
        ("summarize_cleaning_by_id", "A", False, True),
        ("summarize_merging_by_id", "A", False, True),
        ("summarize_crossings_by_id", "B", False, True),
        ("summarize_cleaning_by_id", "B", False, True),
        ("summarize_merging_by_id", "B", False, True),
    ]
    assert summarize_calls == expected


@pytest.mark.unit
def test_river_download_flow_dispatches_swot_and_skips_other_missions(caplog):
    rivers = DummyReservoirs()

    flow = RiverDownloadFlow(
        rivers=rivers,
        to_download=["swot", "icesat2"],
        startdates={"swot": [2024, 1, 1], "icesat2": [2024, 1, 1]},
        enddates={"swot": [2024, 2, 1], "icesat2": [2024, 2, 1]},
    )

    with caplog.at_level("WARNING"):
        flow.run(update_existing=True, enddate_overrides={"swot": [2024, 3, 1]})

    assert rivers.calls == [
        (
            "download_swot_hydrocron",
            {
                "startdate": [2024, 1, 1],
                "enddate": [2024, 3, 1],
                "update_existing": True,
            },
        )
    ]
    assert "Skipping unsupported river mission" in caplog.text


@pytest.mark.unit
def test_swot_raster_download_flow_initialization():
    """Test SWOTRasterDownloadFlow instantiation with valid config."""
    config = {
        "aoi": {"name": "test_aoi", "type": "bbox", "bbox": [0, 0, 1, 1]},
        "product": "SWOT_L2_HR_Raster_D",
        "startdate": [2024, 1, 1],
        "enddate": [2024, 2, 1],
        "granule_filter": None,
        "target_crs": "EPSG:4326",
    }

    flow = SWOTRasterDownloadFlow(
        swot_raster_config=config,
        project_dir="/tmp/test_project",
        earthdata_credentials=("user", "pass"),
    )

    assert flow.swot_raster_config == config
    assert flow.project_dir == "/tmp/test_project"
    assert flow.earthdata_credentials == ("user", "pass")


@pytest.mark.unit
def test_swot_raster_download_flow_resampling_selection():
    """Test correct resampling method selection for different variable types."""
    config = {
        "aoi": {"name": "test_aoi", "type": "bbox", "bbox": [0, 0, 1, 1]},
        "product": "SWOT_L2_HR_Raster_D",
        "startdate": [2024, 1, 1],
        "enddate": [2024, 2, 1],
    }

    flow = SWOTRasterDownloadFlow(
        swot_raster_config=config,
        project_dir="/tmp/test_project",
        earthdata_credentials=("user", "pass"),
    )

    from rasterio.warp import Resampling

    # Test continuous variables get bilinear resampling
    assert flow._resampling_for("wse") == Resampling.bilinear
    assert flow._resampling_for("wse_uncert") == Resampling.bilinear
    assert flow._resampling_for("geoid") == Resampling.bilinear
    assert flow._resampling_for("height_cor_xover") == Resampling.bilinear

    # Test discrete variables get nearest resampling
    assert flow._resampling_for("n_wse_pix") == Resampling.nearest
    assert flow._resampling_for("n_other_pix") == Resampling.nearest
    assert flow._resampling_for("wse_qual") == Resampling.nearest


@pytest.mark.unit
def test_swot_raster_download_flow_crs_detection_from_filename(tmp_path):
    """Test CRS detection from SWOT filename pattern."""
    import xarray as xr

    config = {
        "aoi": {"name": "test_aoi", "type": "bbox", "bbox": [0, 0, 1, 1]},
        "product": "SWOT_L2_HR_Raster_D",
        "startdate": [2024, 1, 1],
        "enddate": [2024, 2, 1],
    }

    flow = SWOTRasterDownloadFlow(
        swot_raster_config=config,
        project_dir=str(tmp_path),
        earthdata_credentials=("user", "pass"),
    )

    # Create a dummy dataset
    ds = xr.Dataset()

    # Test UTM zone extraction from filename
    from pyproj import CRS

    utm_file = (
        tmp_path
        / "SWOT_L2_HR_Raster_D_123_001_UTM45N_20240101T120000_20240101T130000.nc"
    )
    crs = flow._detect_crs(ds, utm_file)

    # Should detect UTM 45 North
    assert crs is not None
    assert crs.to_epsg() == 32645  # EPSG code for UTM 45N


@pytest.mark.unit
def test_swot_raster_download_flow_no_processed_files_skips_merge(tmp_path, caplog):
    """Test that merge phase is skipped when no processed files exist."""
    config = {
        "aoi": {"name": "test_aoi", "type": "bbox", "bbox": [0, 0, 1, 1]},
        "product": "SWOT_L2_HR_Raster_D",
        "startdate": [2024, 1, 1],
        "enddate": [2024, 2, 1],
        "target_crs": "EPSG:4326",
    }

    project_dir = tmp_path / "project"
    project_dir.mkdir()

    with patch.object(SWOTRasterDownloadFlow, "_download_granules", return_value=[]):
        flow = SWOTRasterDownloadFlow(
            swot_raster_config=config,
            project_dir=str(project_dir),
            earthdata_credentials=("user", "pass"),
        )

        with caplog.at_level("INFO"):
            flow.run()

        # Should log that no processed TIFs were found
        assert "No processed TIF files found" in caplog.text


@pytest.mark.unit
def test_swot_raster_download_flow_merge_with_existing_files(tmp_path, caplog):
    """Test that merge phase runs when processed files exist."""
    config = {
        "aoi": {"name": "test_aoi", "type": "bbox", "bbox": [0, 0, 1, 1]},
        "product": "SWOT_L2_HR_Raster_D",
        "startdate": [2024, 1, 1],
        "enddate": [2024, 2, 1],
        "target_crs": "EPSG:4326",
    }

    project_dir = (
        tmp_path
        / "project"
        / "swot_raster"
        / "test_aoi"
        / "processed"
        / "SWOT_L2_HR_Raster_D"
    )
    project_dir.mkdir(parents=True)

    # Create a dummy TIF file
    dummy_tif = project_dir / "20240101T120000_wse.tif"
    dummy_tif.touch()

    with patch.object(SWOTRasterDownloadFlow, "_download_granules", return_value=[]):
        with patch.object(
            SWOTRasterDownloadFlow, "_merge_and_reproject_granules"
        ) as mock_merge:
            flow = SWOTRasterDownloadFlow(
                swot_raster_config=config,
                project_dir=str(tmp_path / "project"),
                earthdata_credentials=("user", "pass"),
            )

            with caplog.at_level("INFO"):
                flow.run()

            # Merge should have been called
            mock_merge.assert_called_once()
            assert "Found 1 processed TIF files" in caplog.text
