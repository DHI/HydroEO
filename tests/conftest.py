"""Shared pytest fixtures for HydroEO unit and integration tests."""

import datetime
import pytest

# ---------------------------------------------------------------------------
# Sample identifiers (fixed for Stage 1 — do NOT change without updating CI)
# ---------------------------------------------------------------------------
SWOT_LAKE_ID = "6340048832"       # Lake Geneva in PLD
SWOT_AOI = [                       # small bounding box around Lake Geneva
    (6.0, 46.2),
    (6.9, 46.2),
    (6.9, 46.6),
    (6.0, 46.6),
    (6.0, 46.2),
]
TEST_START = datetime.date(2024, 1, 1)
TEST_END = datetime.date(2024, 3, 30)

S3_LAKE_NAME = "TITICACA"


# ---------------------------------------------------------------------------
# Lightweight earthaccess mock
# ---------------------------------------------------------------------------

class _FakeEAResult:
    """Minimal stand-in for an earthaccess DataGranule result."""

    def __init__(self, filename: str, url: str | None = None):
        self._filename = filename
        self._url = url or f"https://example.com/swot/{filename}.zip"

    def data_links(self):
        return [self._url]


@pytest.fixture
def fake_ea_prior_result():
    """One baseline-D prior-lake granule result."""
    return _FakeEAResult(
        "SWOT_L2_HR_LakeSP_D_2_20240101T000000_20240101T235959_PIC0_01_prior_0001",
        "https://example.com/swot/SWOT_L2_HR_LakeSP_D_2_20240101T000000_20240101T235959_PIC0_01_prior_0001.zip",
    )


@pytest.fixture
def fake_ea_obs_result():
    """One baseline-D observed-lake granule result (no 'prior' in name)."""
    return _FakeEAResult(
        "SWOT_L2_HR_LakeSP_D_2_20240101T000000_20240101T235959_PIC0_01_obs_0001",
        "https://example.com/swot/SWOT_L2_HR_LakeSP_D_2_20240101T000000_20240101T235959_PIC0_01_obs_0001.zip",
    )


# ---------------------------------------------------------------------------
# CDSE mock helpers
# ---------------------------------------------------------------------------

CDSE_EMPTY_RESPONSE = {
    "features": [],
    "properties": {
        "id": "test",
        "totalResults": 0,
        "links": [],
    },
}

CDSE_ONE_RESULT_RESPONSE = {
    "features": [
        {
            "id": "abc-123",
            "properties": {"startDate": "2024-01-01T00:00:00Z"},
        }
    ],
    "properties": {
        "id": "test",
        "totalResults": 1,
        "links": [],   # no 'next' page
    },
}
