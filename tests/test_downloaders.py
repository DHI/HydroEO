"""Unit tests for download routines.

Covers:
  (a) SWOT sample download path (mocked earthaccess)
  (b) Sentinel-3 sample download path (mocked requests / CDSE)
  (c) Invalid / unrecognised product name handling
  (d) Authentication failure handling (CreoDIAS token endpoint)

All tests in this file are mocked (no live network calls).
Live tests live in test_integration.py.
"""

import datetime
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from tests.conftest import (
    CDSE_EMPTY_RESPONSE,
    CDSE_ONE_RESULT_RESPONSE,
    SWOT_AOI,
    TEST_END,
    TEST_START,
)


# ---------------------------------------------------------------------------
# (a) SWOT download path
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_swot_query_uses_baseline_d_short_name():
    """swot.query() must pass SWOT_L2_HR_LakeSP_D to earthaccess.search_data."""
    from HydroEO.satellites.swot import SWOT_LAKE_SHORT_NAME, query

    captured = {}

    def fake_search(**kwargs):
        captured.update(kwargs)
        return []

    with patch("HydroEO.satellites.swot.earthaccess") as mock_ea:
        mock_ea.search_data.side_effect = fake_search
        mock_ea.login.return_value = None

        query(
            aoi=SWOT_AOI,
            startdate=TEST_START,
            enddate=TEST_END,
            earthdata_credentials=("user", "pass"),
        )

    assert captured.get("short_name") == SWOT_LAKE_SHORT_NAME


@pytest.mark.unit
def test_swot_download_skips_already_logged(tmp_path, fake_ea_prior_result):
    """swot.download() must skip files already listed in the log."""
    from HydroEO.satellites.swot import download

    log = tmp_path / "downloaded.log"
    file_stem = fake_ea_prior_result.data_links()[0].split("/")[-1].split(".")[0]
    log.write_text(file_stem + "\n")

    with patch("HydroEO.satellites.swot.earthaccess") as mock_ea:
        mock_ea.download.return_value = []
        result = download([fake_ea_prior_result], str(tmp_path))

    mock_ea.download.assert_not_called()
    assert result == []


@pytest.mark.unit
def test_swot_download_queues_new_file(tmp_path, fake_ea_prior_result):
    """swot.download() must request a download for files not in the log."""
    from HydroEO.satellites.swot import download

    with patch("HydroEO.satellites.swot.earthaccess") as mock_ea:
        fake_path = str(tmp_path / "SWOT_file.zip")
        mock_ea.download.return_value = [fake_path]
        result = download([fake_ea_prior_result], str(tmp_path))

    mock_ea.download.assert_called_once()
    assert result == [fake_path]


@pytest.mark.unit
def test_swot_download_handles_path_objects(tmp_path, fake_ea_prior_result):
    """swot.download() must accept pathlib.Path objects from earthaccess.download."""
    from HydroEO.satellites.swot import download

    with patch("HydroEO.satellites.swot.earthaccess") as mock_ea:
        fake_path = tmp_path / "SWOT_file.zip"
        mock_ea.download.return_value = [fake_path]
        result = download([fake_ea_prior_result], str(tmp_path))

    assert result == [fake_path]
    log_path = tmp_path / "downloaded.log"
    assert "SWOT_file" in log_path.read_text()


# ---------------------------------------------------------------------------
# SWOT granule filter — baseline D naming conventions
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_swot_prior_granule_selected(fake_ea_prior_result):
    """Granules whose URL contains '_prior_' must be selected by the filter."""
    link = fake_ea_prior_result.data_links()[0]
    filename = link.split("/")[-1].lower()
    # Replicate the system.py filter logic
    selected = "_prior_" in filename or "prior" in filename.split("_")
    assert selected


@pytest.mark.unit
def test_swot_obs_granule_not_selected(fake_ea_obs_result):
    """Observed-lake granules must NOT be selected by the prior filter."""
    link = fake_ea_obs_result.data_links()[0]
    filename = link.split("/")[-1].lower()
    selected = "_prior_" in filename or "prior" in filename.split("_")
    assert not selected


# ---------------------------------------------------------------------------
# (b) Sentinel-3 download path
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_sentinel3_query_returns_ids():
    """sentinel.query() for S3 product must return a list of IDs from CDSE."""
    from HydroEO.satellites.sentinel import query

    with patch("HydroEO.satellites.sentinel.creodias") as mock_cdse:
        mock_cdse.query.return_value = {"abc-123": {"id": "abc-123"}}
        ids = query(aoi=SWOT_AOI, startdate=TEST_START, enddate=TEST_END, product="S3")

    assert ids == ["abc-123"]
    call_kwargs = mock_cdse.query.call_args.kwargs
    assert call_kwargs["collection"] == "Sentinel3"
    assert call_kwargs["productType"] == "SR_2_LAN_HY"


@pytest.mark.unit
def test_sentinel6_query_returns_ids():
    """sentinel.query() for S6 product must return a list of IDs from CDSE."""
    from HydroEO.satellites.sentinel import query

    with patch("HydroEO.satellites.sentinel.creodias") as mock_cdse:
        mock_cdse.query.return_value = {"def-456": {"id": "def-456"}}
        ids = query(aoi=SWOT_AOI, startdate=TEST_START, enddate=TEST_END, product="S6")

    assert ids == ["def-456"]
    call_kwargs = mock_cdse.query.call_args.kwargs
    assert call_kwargs["collection"] == "Sentinel6"


@pytest.mark.unit
def test_sentinel3_download_skips_already_logged(tmp_path):
    """sentinel.download() must skip IDs already present in the log file."""
    from HydroEO.satellites.sentinel import download

    uid = "abc-123"
    log = tmp_path / "downloaded.log"
    log.write_text(uid + "\n")

    with patch("HydroEO.satellites.sentinel.creodias") as mock_cdse:
        mock_cdse.download_list.return_value = (None, None)
        download(
            ids=[uid],
            download_directory=str(tmp_path),
            creodias_credentials=("user", "pass"),
        )

    # download_list must be called with an empty list of UIDs to download
    call_args = mock_cdse.download_list.call_args
    uids_arg = call_args.args[0] if call_args.args else call_args.kwargs.get("uids", call_args.kwargs.get("ids", []))
    assert uid not in uids_arg


# ---------------------------------------------------------------------------
# (c) Invalid / unrecognised product name
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_invalid_product_raises_in_sentinel_query():
    """sentinel.query() must raise ValueError for unrecognised product labels."""
    from HydroEO.satellites.sentinel import query

    with pytest.raises(ValueError, match="unrecognized"):
        query(
            aoi=SWOT_AOI,
            startdate=TEST_START,
            enddate=TEST_END,
            product="MODIS",
        )


@pytest.mark.unit
def test_invalid_product_raises_in_system_download_altimetry():
    """System.download_altimetry() must raise ValueError for unknown product strings."""
    from HydroEO.system import System
    import geopandas as gpd
    from shapely.geometry import box

    gdf = gpd.GeoDataFrame(
        {"id": [1], "geometry": [box(6.0, 46.2, 6.9, 46.6)]},
        crs="EPSG:4326",
    )
    sys = System.__new__(System)
    sys.gdf = gdf
    sys.id_key = "id"
    sys.dirs = {"main": "/tmp"}
    sys.download_gdf = gdf
    sys.type = "reservoirs"

    with pytest.raises(ValueError, match="not accepted"):
        sys.download_altimetry(
            product="LANDSAT",
            startdate=(2024, 1, 1),
            enddate=(2024, 3, 30),
        )


# ---------------------------------------------------------------------------
# (d) Authentication failure
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_creodias_get_token_raises_on_bad_credentials():
    """_get_token() must raise RuntimeError when the token endpoint returns no access_token."""
    from HydroEO.downloaders.creodias import _get_token

    bad_response = MagicMock()
    bad_response.json.return_value = {"error": "invalid_grant", "error_description": "Invalid credentials"}

    with patch("HydroEO.downloaders.creodias.requests.post", return_value=bad_response):
        with pytest.raises(RuntimeError, match="Unable to get token"):
            _get_token("wrong_user", "wrong_pass")


@pytest.mark.unit
def test_creodias_query_raises_on_http_error():
    """creodias.query() must surface HTTP errors raised by requests."""
    import requests
    from HydroEO.downloaders.creodias import query

    mock_response = MagicMock()
    mock_response.raise_for_status.side_effect = requests.HTTPError("401 Unauthorized")

    with patch("HydroEO.downloaders.creodias.requests.get", return_value=mock_response):
        with pytest.raises(requests.HTTPError):
            query(
                collection="Sentinel3",
                start_date=TEST_START,
                end_date=TEST_END,
            )
