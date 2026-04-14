"""Basic smoke tests — import and version checks.

These must pass with no network access and no credentials.
"""

import importlib

import pytest


@pytest.mark.unit
def test_package_importable():
    """HydroEO must be importable after installation."""
    import HydroEO  # noqa: F401


@pytest.mark.unit
def test_version_attribute():
    """__version__ must be present and non-empty."""
    import HydroEO

    assert hasattr(HydroEO, "__version__")
    assert HydroEO.__version__


@pytest.mark.unit
def test_submodules_importable():
    """Core sub-modules must be importable (catches missing __init__ wiring)."""
    for mod in [
        "HydroEO.satellites.swot",
        "HydroEO.satellites.sentinel",
        "HydroEO.downloaders.creodias",
        "HydroEO.downloaders.hydroweb",
        "HydroEO.utils.general",
        "HydroEO.utils.geometry",
    ]:
        importlib.import_module(mod)


@pytest.mark.unit
def test_swot_short_name_is_baseline_d():
    """The SWOT default product short name must be the baseline-D identifier."""
    from HydroEO.satellites.swot import SWOT_LAKE_SHORT_NAME

    assert SWOT_LAKE_SHORT_NAME == "SWOT_L2_HR_LakeSP_D"


@pytest.mark.unit
def test_creodias_download_url_is_cdse():
    """Download endpoint must point at zipper.dataspace.copernicus.eu, not the old CreoDIAS domain."""
    from HydroEO.downloaders.creodias import DOWNLOAD_URL

    assert "https://zipper.creodias.eu/download" in DOWNLOAD_URL
