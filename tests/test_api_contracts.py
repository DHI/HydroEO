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
