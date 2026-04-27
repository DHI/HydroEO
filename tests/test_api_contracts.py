"""API contract tests — credential-free endpoint & schema validation.

These tests check:
  1. Endpoints are reachable (no DNS/firewall issues)
  2. Response structures match expected schemas (detects API versioning changes)
  3. Auth flows use expected protocols (OIDC, etc.)

No credentials required for most. Marked `api_contract` for CI.
"""

import datetime
import os

import pytest

# Lake Geneva (SWOT test area)
SWOT_AOI = [
    (6.0, 46.2),
    (6.9, 46.2),
    (6.9, 46.6),
    (6.0, 46.6),
    (6.0, 46.2),
]
TEST_START = datetime.date(2024, 1, 1)
TEST_END = datetime.date(2024, 3, 30)

_has_edl = pytest.mark.skipif(
    not (os.environ.get("EDL_USERNAME") and os.environ.get("EDL_PASSWORD")),
    reason="EDL_USERNAME / EDL_PASSWORD not set",
)

_has_creodias = pytest.mark.skipif(
    not (os.environ.get("CREODIAS_USERNAME") and os.environ.get("CREODIAS_PASSWORD")),
    reason="CREODIAS_USERNAME / CREODIAS_PASSWORD not set",
)

_has_hydroweb = pytest.mark.skipif(
    not os.environ.get("HYDROWEB_API_KEY"),
    reason="HYDROWEB_API_KEY not set",
)


# ---------------------------------------------------------------------------
# CDSE OData endpoint checks (credential-free)
# ---------------------------------------------------------------------------


@pytest.mark.api_contract
def test_cdse_odata_endpoint_reachable():
    """CDSE OData Products endpoint must be reachable and return valid collection."""
    import requests

    # Query a specific collection endpoint (more reliable than root)
    url = "https://catalogue.dataspace.copernicus.eu/odata/v1/Products?$filter=Collection/Name%20eq%20%27SENTINEL-3%27&$top=1"

    try:
        resp = requests.get(url, timeout=10)
        # Should be reachable; may be 401 without auth, 200 with auth, 400 for bad query
        assert resp.status_code in (200, 401, 400), (
            f"CDSE OData endpoint returned {resp.status_code}. "
            "Endpoint may be down or URL changed."
        )
        if resp.status_code == 200:
            data = resp.json()
            assert "value" in data or "features" in data, (
                f"CDSE response missing expected 'value' or 'features' key: {list(data.keys())}"
            )
    except requests.exceptions.ConnectionError as e:
        pytest.fail(f"CDSE OData endpoint unreachable: {e}")
    except requests.exceptions.Timeout:
        pytest.fail("CDSE OData endpoint timed out (>10s)")


@pytest.mark.api_contract
def test_cdse_token_endpoint_bad_credentials_returns_4xx():
    """CDSE token endpoint must reject bad credentials with 400-range status."""
    import requests

    url = "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token"

    payload = {
        "grant_type": "password",
        "username": "invalid@test.local",
        "password": "invalid_password",
        "client_id": "cdse-public",
    }

    try:
        resp = requests.post(url, data=payload, timeout=10)
        assert 400 <= resp.status_code < 500, (
            f"Expected 4xx on bad creds, got {resp.status_code}. "
            "Token endpoint protocol may have changed."
        )
        data = resp.json()
        assert "error" in data, (
            f"CDSE token error response missing 'error' field: {list(data.keys())}"
        )
    except requests.exceptions.ConnectionError as e:
        pytest.fail(f"CDSE token endpoint unreachable: {e}")
    except requests.exceptions.Timeout:
        pytest.fail("CDSE token endpoint timed out (>10s)")


# ---------------------------------------------------------------------------
# earthaccess product existence checks (requires EDL creds)
# ---------------------------------------------------------------------------


@pytest.mark.api_contract
@_has_edl
def test_earthaccess_swot_product_exists():
    """earthaccess must find SWOT_L2_HR_LakeSP_D product."""
    import earthaccess

    os.environ["EARTHDATA_USERNAME"] = os.environ["EDL_USERNAME"]
    os.environ["EARTHDATA_PASSWORD"] = os.environ["EDL_PASSWORD"]

    try:
        results = earthaccess.search_data(
            short_name="SWOT_L2_HR_LakeSP_D",
            count=1,
        )
        assert len(results) > 0, (
            "No SWOT_L2_HR_LakeSP_D product found in earthaccess. "
            "Product may have been retired or short_name changed."
        )
    except Exception as e:
        pytest.fail(f"earthaccess SWOT product search failed: {e}")


@pytest.mark.api_contract
@_has_edl
def test_earthaccess_icesat2_product_exists():
    """earthaccess must find ATL13 product."""
    import earthaccess

    os.environ["EARTHDATA_USERNAME"] = os.environ["EDL_USERNAME"]
    os.environ["EARTHDATA_PASSWORD"] = os.environ["EDL_PASSWORD"]

    try:
        results = earthaccess.search_data(
            short_name="ATL13",
            count=1,
        )
        assert len(results) > 0, (
            "No ATL13 product found in earthaccess. "
            "Product may have been renamed or data granules removed."
        )
    except Exception as e:
        pytest.fail(f"earthaccess ATL13 product search failed: {e}")


# ---------------------------------------------------------------------------
# HydroWeb API key validation
# ---------------------------------------------------------------------------


@pytest.mark.api_contract
@_has_hydroweb
def test_hydroweb_api_key_valid():
    """HydroWeb API key must be accepted (no 401/403 on basic client init)."""
    try:
        from py_hydroweb import Client

        client = Client(api_key=os.environ["HYDROWEB_API_KEY"])
        # If initialization succeeds without raising an auth error, the key is valid
        assert client is not None
    except ImportError:
        pytest.skip("py_hydroweb not installed")
    except Exception as e:
        pytest.fail(f"HydroWeb client initialization failed: {e}")


# ---------------------------------------------------------------------------
# Hydrocron (SWOT river timeseries) endpoint checks — credential-free
# ---------------------------------------------------------------------------

# Loire River node 23227000010171 — confirmed working in Hydrocron with data from 2023+.
_HYDROCRON_TEST_NODE_ID = "23227000010171"
_HYDROCRON_URL = "https://soto.podaac.earthdatacloud.nasa.gov/hydrocron/v1/timeseries"


@pytest.mark.api_contract
def test_hydrocron_endpoint_is_reachable():
    """Hydrocron timeseries endpoint must be reachable and return a JSON body."""
    import requests

    params = {
        "feature": "Node",
        "feature_id": _HYDROCRON_TEST_NODE_ID,
        "start_time": "2024-01-01T00:00:00Z",
        "end_time": "2024-01-15T00:00:00Z",
        "output": "csv",
        "fields": "node_id,node_q,time_str,wse",
    }

    try:
        resp = requests.get(_HYDROCRON_URL, params=params, timeout=30)
    except requests.exceptions.ConnectionError as exc:
        pytest.fail(f"Hydrocron endpoint unreachable: {exc}")
    except requests.exceptions.Timeout:
        pytest.fail("Hydrocron endpoint timed out (>30 s)")

    assert resp.status_code == 200, (
        f"Expected 200 from Hydrocron, got {resp.status_code}. "
        "Endpoint URL or API version may have changed."
    )

    payload = resp.json()
    assert isinstance(payload, dict), (
        f"Hydrocron response must be a JSON object, got {type(payload)}"
    )


@pytest.mark.api_contract
def test_hydrocron_response_schema_has_results_csv_key():
    """Hydrocron timeseries response must include 'results.csv' for csv output mode."""
    import requests

    params = {
        "feature": "Node",
        "feature_id": _HYDROCRON_TEST_NODE_ID,
        "start_time": "2024-01-01T00:00:00Z",
        "end_time": "2024-02-01T00:00:00Z",
        "output": "csv",
        "fields": "node_id,node_q,time_str,wse",
    }

    try:
        resp = requests.get(_HYDROCRON_URL, params=params, timeout=30)
    except (requests.exceptions.ConnectionError, requests.exceptions.Timeout) as exc:
        pytest.skip(f"Hydrocron unreachable — skipping schema validation: {exc}")

    if resp.status_code != 200:
        pytest.skip(f"Hydrocron returned {resp.status_code} — cannot validate schema")

    payload = resp.json()
    assert "results" in payload, (
        f"Top-level 'results' key missing from Hydrocron response. Keys: {list(payload.keys())}"
    )
    assert "csv" in payload["results"], (
        f"'results.csv' key missing. Keys under 'results': {list(payload['results'].keys())}"
    )
    assert isinstance(payload["results"]["csv"], str), (
        "'results.csv' must be a string (the raw CSV text)"
    )
