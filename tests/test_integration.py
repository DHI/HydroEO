"""Live integration tests — require valid API credentials.

These tests make real HTTP calls and are gated by environment variables / GitHub
Actions secrets. A test is automatically *skipped* if the required secret is
absent, so CI remains green on forks or when credentials are not configured.

Marks: all tests in this file carry both `integration` and `unit` markers.
To run locally with credentials:

    export EDL_USERNAME=...
    export EDL_PASSWORD=...
    export CREODIAS_USERNAME=...
    export CREODIAS_PASSWORD=...
    pytest tests/test_integration.py -m integration -v
"""

import datetime
import os

import pytest

# ---------------------------------------------------------------------------
# Shared coordinates and date window
# ---------------------------------------------------------------------------

SWOT_LAKE_ID = "6340048832"       # Lake Geneva in PLD
SWOT_AOI = [
    (6.0, 46.2),
    (6.9, 46.2),
    (6.9, 46.6),
    (6.0, 46.6),
    (6.0, 46.2),
]
TEST_START = datetime.date(2024, 1, 1)
TEST_END = datetime.date(2024, 3, 30)

TITICACA_AOI = [
    (-70.5, -16.5),
    (-68.5, -16.5),
    (-68.5, -15.0),
    (-70.5, -15.0),
    (-70.5, -16.5),
]

# ---------------------------------------------------------------------------
# Skip helpers
# ---------------------------------------------------------------------------

_has_edl = pytest.mark.skipif(
    not (os.environ.get("EDL_USERNAME") and os.environ.get("EDL_PASSWORD")),
    reason="EDL_USERNAME / EDL_PASSWORD not set — skipping live SWOT test",
)

_has_creodias = pytest.mark.skipif(
    not (os.environ.get("CREODIAS_USERNAME") and os.environ.get("CREODIAS_PASSWORD")),
    reason="CREODIAS_USERNAME / CREODIAS_PASSWORD not set — skipping live Sentinel test",
)


# ---------------------------------------------------------------------------
# (a) SWOT live sample search
# ---------------------------------------------------------------------------


@pytest.mark.integration
@_has_edl
def test_swot_live_query_returns_results():
    """Live SWOT query for Lake Geneva AOI must return ≥ 1 baseline-D granule."""
    from HydroEO.satellites.swot import SWOT_LAKE_SHORT_NAME, query

    results = query(
        aoi=SWOT_AOI,
        startdate=TEST_START,
        enddate=TEST_END,
        earthdata_credentials=(
            os.environ["EDL_USERNAME"],
            os.environ["EDL_PASSWORD"],
        ),
    )

    assert len(results) > 0, (
        f"Expected ≥ 1 {SWOT_LAKE_SHORT_NAME} granule for Lake Geneva between "
        f"{TEST_START} and {TEST_END}, got 0. "
        "Check that baseline D data is available for this period."
    )


@pytest.mark.integration
@_has_edl
def test_swot_live_prior_granules_present():
    """At least one prior-lake granule must be present in the baseline-D results."""
    from HydroEO.satellites.swot import query

    results = query(
        aoi=SWOT_AOI,
        startdate=TEST_START,
        enddate=TEST_END,
        earthdata_credentials=(
            os.environ["EDL_USERNAME"],
            os.environ["EDL_PASSWORD"],
        ),
    )

    prior_results = [
        r for r in results
        if "_prior_" in r.data_links()[0].split("/")[-1].lower()
    ]

    assert len(prior_results) > 0, (
        "No prior-lake granules found. Check that the granule filter matches "
        "baseline-D file naming conventions."
    )


# ---------------------------------------------------------------------------
# (b) Sentinel-3 live sample search
# ---------------------------------------------------------------------------


@pytest.mark.integration
@_has_creodias
def test_sentinel3_live_query_returns_results():
    """Live Sentinel-3 query for Lake Titicaca must return ≥ 1 granule ID."""
    from HydroEO.satellites.sentinel import query

    ids = query(
        aoi=TITICACA_AOI,
        startdate=TEST_START,
        enddate=TEST_END,
        product="S3",
    )

    assert len(ids) > 0, (
        f"Expected ≥ 1 Sentinel-3 SR_2_LAN_HY granule for Lake Titicaca between "
        f"{TEST_START} and {TEST_END}, got 0. "
        "Verify the productType and collection name against the CDSE catalogue."
    )


# ---------------------------------------------------------------------------
# (d) Auth failure — live token endpoint with bad credentials
# ---------------------------------------------------------------------------


@pytest.mark.integration
def test_creodias_bad_credentials_raises():
    """Attempting to get a CDSE token with wrong credentials must raise RuntimeError."""
    from HydroEO.downloaders.creodias import _get_token

    with pytest.raises(RuntimeError, match="Unable to get token"):
        _get_token("definitely_not_a_real_user", "definitely_wrong_password")
