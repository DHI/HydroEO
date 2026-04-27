"""Unit tests for system-level product validation paths."""

import geopandas as gpd
import pytest
from shapely.geometry import box


@pytest.mark.unit
def test_invalid_product_raises_in_system_download_altimetry():
    """Reservoirs.download_altimetry() must raise ValueError for unknown product strings."""
    from HydroEO.waterbody import Reservoirs

    gdf = gpd.GeoDataFrame(
        {"id": [1], "geometry": [box(6.0, 46.2, 6.9, 46.6)]},
        crs="EPSG:4326",
    )
    res = Reservoirs.__new__(Reservoirs)
    res.gdf = gdf
    res.id_key = "id"
    res.dirs = {"main": "/tmp"}
    res.download_gdf = gdf
    res.type = "reservoirs"

    with pytest.raises(ValueError, match="not accepted"):
        res.download_altimetry(
            product="LANDSAT",
            startdate=(2024, 1, 1),
            enddate=(2024, 3, 30),
        )
