"""Unit tests for CREODIAS client/auth helpers (fully mocked)."""

from unittest.mock import MagicMock, patch

import pytest

from tests.conftest import TEST_END, TEST_START


@pytest.mark.unit
def test_creodias_get_token_raises_on_bad_credentials():
    """_get_token() must raise RuntimeError when the token endpoint returns no access_token."""
    from HydroEO.downloaders.creodias import _get_token

    bad_response = MagicMock()
    bad_response.json.return_value = {
        "error": "invalid_grant",
        "error_description": "Invalid credentials",
    }

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
