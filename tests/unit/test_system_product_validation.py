"""Unit tests for system-level product validation paths."""

import geopandas as gpd
import pytest
from shapely.geometry import box


@pytest.mark.unit
def test_invalid_product_raises_in_system_download_altimetry():
    """System.download_altimetry() must raise ValueError for unknown product strings."""
    from HydroEO.system import System

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
