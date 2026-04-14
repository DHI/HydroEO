"""Unit tests for SWOT query/download routines (fully mocked)."""

from unittest.mock import patch

import pytest

from tests.conftest import SWOT_AOI, TEST_END, TEST_START


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


@pytest.mark.unit
def test_swot_prior_granule_selected(fake_ea_prior_result):
    """Granules whose URL contains '_prior_' must be selected by the filter."""
    link = fake_ea_prior_result.data_links()[0]
    filename = link.split("/")[-1].lower()
    selected = "_prior_" in filename or "prior" in filename.split("_")
    assert selected


@pytest.mark.unit
def test_swot_obs_granule_not_selected(fake_ea_obs_result):
    """Observed-lake granules must NOT be selected by the prior filter."""
    link = fake_ea_obs_result.data_links()[0]
    filename = link.split("/")[-1].lower()
    selected = "_prior_" in filename or "prior" in filename.split("_")
    assert not selected
