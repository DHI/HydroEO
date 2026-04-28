"""Unit tests for waterbody orchestration methods and SWOT raster pipeline."""

import pandas as pd
import pytest
from unittest.mock import patch, MagicMock, call

from HydroEO.waterbody import Reservoirs, Rivers
from HydroEO.satellites.swot.raster import _resampling_for, _detect_crs, download_raster


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_STARTDATES = {
    "swot": [2024, 1, 1],
    "icesat2": [2024, 1, 2],
    "sentinel3": [2024, 1, 3],
    "sentinel6": [2024, 1, 4],
}

_ENDDATES = {
    "swot": [2024, 2, 1],
    "icesat2": [2024, 2, 2],
    "sentinel3": [2024, 2, 3],
    "sentinel6": [2024, 2, 4],
}


# ---------------------------------------------------------------------------
# Reservoirs.download()
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_reservoirs_download_dispatches_products_and_credentials():
    mock_res = MagicMock(spec=Reservoirs)

    Reservoirs.download(
        mock_res,
        to_download=["swot", "icesat2", "sentinel3", "sentinel6"],
        startdates=_STARTDATES,
        enddates=_ENDDATES,
        earthdata_credentials=("edl-user", "edl-pass"),
        creodias_credentials_provider=lambda: ("creo-user", "creo-pass"),
        update_existing=True,
        enddate_overrides={"swot": [2024, 3, 1]},
    )

    calls = mock_res.download_altimetry.call_args_list
    assert len(calls) == 4

    swot_kw = calls[0].kwargs
    assert swot_kw["product"] == "SWOT_LAKE"
    assert swot_kw["enddate"] == [2024, 3, 1]
    assert swot_kw["update_existing"] is True

    icesat2_kw = calls[1].kwargs
    assert icesat2_kw["product"] == "ATL13"
    assert icesat2_kw["credentials"] == ("edl-user", "edl-pass")

    s3_kw = calls[2].kwargs
    assert s3_kw["product"] == "S3"
    assert s3_kw["credentials"] == ("creo-user", "creo-pass")

    s6_kw = calls[3].kwargs
    assert s6_kw["product"] == "S6"
    assert s6_kw["credentials"] == ("creo-user", "creo-pass")


# ---------------------------------------------------------------------------
# Reservoirs.process()
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_reservoirs_process_runs_extract_clean_merge_in_order():
    mock_res = MagicMock(spec=Reservoirs)
    call_order = []

    mock_res.extract_product_timeseries.side_effect = lambda p: call_order.append(
        "extract"
    )
    mock_res.clean_product_timeseries.side_effect = (
        lambda products, filter_options_by_product: call_order.append("clean")
    )
    mock_res.merge_product_timeseries.side_effect = lambda products: call_order.append(
        "merge"
    )

    Reservoirs.process(
        mock_res,
        to_process=["swot", "icesat2"],
        processing_options={"swot": {"processing_filters": ["elevation"]}},
    )

    assert call_order == ["extract", "clean", "merge"]
    mock_res.extract_product_timeseries.assert_called_once_with(["swot", "icesat2"])
    mock_res.merge_product_timeseries.assert_called_once_with(
        products=["swot", "icesat2"]
    )


# ---------------------------------------------------------------------------
# Reservoirs.summarize()
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_reservoirs_summarize_runs_all_steps_for_each_reservoir():
    mock_res = MagicMock(spec=Reservoirs)
    mock_res.download_gdf = pd.DataFrame({"rid": ["A", "B"]})
    mock_res.id_key = "rid"

    Reservoirs.summarize(mock_res, show=False, save=True)

    mock_res.summarize_crossings_by_id.assert_has_calls(
        [call("A", show=False, save=True), call("B", show=False, save=True)]
    )
    mock_res.summarize_cleaning_by_id.assert_has_calls(
        [call("A", show=False, save=True), call("B", show=False, save=True)]
    )
    mock_res.summarize_merging_by_id.assert_has_calls(
        [call("A", show=False, save=True), call("B", show=False, save=True)]
    )
    assert mock_res.summarize_crossings_by_id.call_count == 2
    assert mock_res.summarize_cleaning_by_id.call_count == 2
    assert mock_res.summarize_merging_by_id.call_count == 2


# ---------------------------------------------------------------------------
# Rivers.download()
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_rivers_download_dispatches_swot_and_skips_other_missions(caplog):
    mock_rivers = MagicMock(spec=Rivers)
    mock_rivers.download_swot_hydrocron.return_value = {
        "requested": 1,
        "successful": 1,
        "failed": 0,
        "empty_after_filter": 0,
    }

    with caplog.at_level("WARNING"):
        Rivers.download(
            mock_rivers,
            to_download=["swot", "icesat2"],
            startdates={"swot": [2024, 1, 1], "icesat2": [2024, 1, 1]},
            enddates={"swot": [2024, 2, 1], "icesat2": [2024, 2, 1]},
            update_existing=True,
            enddate_overrides={"swot": [2024, 3, 1]},
        )

    mock_rivers.download_swot_hydrocron.assert_called_once_with(
        startdate=[2024, 1, 1],
        enddate=[2024, 3, 1],
        update_existing=True,
    )
    assert "Skipping unsupported river mission" in caplog.text


# ---------------------------------------------------------------------------
# SWOT raster free functions
# ---------------------------------------------------------------------------


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
        "target_crs": "EPSG:4326",
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
        "target_crs": "EPSG:4326",
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
