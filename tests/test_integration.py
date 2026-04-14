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
import re

import pytest

# ---------------------------------------------------------------------------
# Shared coordinates and date window
# ---------------------------------------------------------------------------

SWOT_LAKE_ID = "6340048832"  # Lake Geneva in PLD
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
    reason="EDL_USERNAME / EDL_PASSWORD not set — skipping EDL-backed live tests",
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

    os.environ["EARTHDATA_USERNAME"] = os.environ["EDL_USERNAME"]
    os.environ["EARTHDATA_PASSWORD"] = os.environ["EDL_PASSWORD"]

    results = query(
        aoi=SWOT_AOI,
        startdate=TEST_START,
        enddate=TEST_END,
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

    os.environ["EARTHDATA_USERNAME"] = os.environ["EDL_USERNAME"]
    os.environ["EARTHDATA_PASSWORD"] = os.environ["EDL_PASSWORD"]

    results = query(
        aoi=SWOT_AOI,
        startdate=TEST_START,
        enddate=TEST_END,
    )

    prior_results = [
        r for r in results if "_prior_" in r.data_links()[0].split("/")[-1].lower()
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


# ---------------------------------------------------------------------------
# (c) SWOT & Sentinel schema validation
# ---------------------------------------------------------------------------


@pytest.mark.integration
@_has_edl
def test_swot_granule_has_data_links():
    """Live: Each SWOT granule result must have data_links() returning URLs."""
    from HydroEO.satellites.swot import query

    os.environ["EARTHDATA_USERNAME"] = os.environ["EDL_USERNAME"]
    os.environ["EARTHDATA_PASSWORD"] = os.environ["EDL_PASSWORD"]

    results = query(
        aoi=SWOT_AOI,
        startdate=TEST_START,
        enddate=TEST_END,
    )

    assert len(results) > 0, "No SWOT results to validate"

    for granule in results[:1]:  # Test first granule only
        links = granule.data_links()
        assert isinstance(links, list), "data_links() must return a list"
        assert len(links) > 0, "granule.data_links() returned empty list"
        assert any(url.endswith((".zip", ".nc")) for url in links), (
            f"Expected .zip or .nc URLs, got {links}"
        )


@pytest.mark.integration
@_has_edl
def test_swot_prior_granule_naming_convention():
    """Live: Verify at least one granule uses '_prior_' naming (guards filename filter)."""
    from HydroEO.satellites.swot import query

    os.environ["EARTHDATA_USERNAME"] = os.environ["EDL_USERNAME"]
    os.environ["EARTHDATA_PASSWORD"] = os.environ["EDL_PASSWORD"]

    results = query(
        aoi=SWOT_AOI,
        startdate=TEST_START,
        enddate=TEST_END,
    )

    prior_granules = [
        r for r in results if "_prior_" in r.data_links()[0].split("/")[-1].lower()
    ]

    assert len(prior_granules) > 0, (
        "Expected ≥ 1 granule with '_prior_' in filename. "
        "This guards the granule filtering logic — if missing, the filter may need updating."
    )


@pytest.mark.integration
@_has_creodias
def test_sentinel3_response_shape():
    """Live: Sentinel-3 query response must have expected OData schema."""
    from HydroEO.downloaders.creodias import query

    results = query(
        collection="SENTINEL-3",
        start_date=TEST_START,
        end_date=TEST_END,
        geometry=TITICACA_AOI,
        productType="SR_2_LAN_HY",
    )

    assert len(results) > 0, "No Sentinel-3 results — cannot validate schema"

    for result_id, result_data in list(results.items())[:1]:
        assert "Id" in result_data, f"Missing 'Id' key in response: {result_data}"
        assert "Name" in result_data, f"Missing 'Name' key in response: {result_data}"
        assert isinstance(result_data.get("Id"), str), "Id should be a string"
        assert isinstance(result_data.get("Name"), str), "Name should be a string"


@pytest.mark.integration
@_has_creodias
def test_sentinel6_live_query_returns_results():
    """Live: Sentinel-6 query must return ≥ 1 granule ID for Lake Titicaca."""
    from HydroEO.satellites.sentinel import query

    ids = query(
        aoi=TITICACA_AOI,
        startdate=TEST_START,
        enddate=TEST_END,
        product="S6",
    )

    assert len(ids) > 0, (
        f"Expected ≥ 1 Sentinel-6 granule for Lake Titicaca between "
        f"{TEST_START} and {TEST_END}, got 0. "
        "Verify S6 productType 'P4_2__LR_____' is available in CDSE."
    )


@pytest.mark.integration
@_has_creodias
def test_sentinel6_response_shape():
    """Live: Sentinel-6 query response must match OData schema (same as S3)."""
    from HydroEO.downloaders.creodias import query

    results = query(
        collection="SENTINEL-6",
        start_date=TEST_START,
        end_date=TEST_END,
        geometry=TITICACA_AOI,
        productType="P4_2__LR_____",
    )

    assert len(results) > 0, "No Sentinel-6 results — cannot validate schema"

    for result_id, result_data in list(results.items())[:1]:
        assert "Id" in result_data, f"Missing 'Id' key in S6 response: {result_data}"
        assert "Name" in result_data, (
            f"Missing 'Name' key in S6 response: {result_data}"
        )


@pytest.mark.integration
@_has_edl
def test_icesat2_live_harmony_query_downloads_atl13(tmp_path):
    """Live: ATL13 request should run through Harmony and download at least one file."""
    from HydroEO.satellites import icesat2

    os.environ["EARTHDATA_USERNAME"] = os.environ["EDL_USERNAME"]
    os.environ["EARTHDATA_PASSWORD"] = os.environ["EDL_PASSWORD"]

    download_dir = tmp_path / "icesat2_harmony"
    download_dir.mkdir(parents=True, exist_ok=True)

    try:
        _ = icesat2.query(
            aoi=TITICACA_AOI,
            startdate=TEST_START,
            enddate=TEST_END,
            earthdata_credentials=(
                os.environ["EDL_USERNAME"],
                os.environ["EDL_PASSWORD"],
            ),
            download_directory=str(download_dir),
            product="ATL13",
        )
    except Exception as exc:
        message = str(exc)

        # Live endpoints can intermittently return 5xx; treat this as infra flakiness.
        if re.search(r"HTTP\s+5\d\d", message) or any(
            code in message for code in [" 500", " 502", " 503", " 504"]
        ):
            pytest.skip(f"Skipping due to transient Harmony/EDL service issue: {message}")

        # Keep credential/auth problems as failures because they indicate setup issues.
        if "incorrect or missing credentials" in message or "HTTP 401" in message:
            pytest.fail(f"EDL credentials rejected during Harmony auth: {message}")

        raise

    downloaded_files = [p for p in download_dir.rglob("*") if p.is_file()]
    assert len(downloaded_files) > 0, (
        "Harmony ATL13 query completed but no files were downloaded. "
        "Check Harmony availability, ATL13 coverage, or auth configuration."
    )
