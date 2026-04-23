"""Unit tests for ICESat-2 ATL13 SlideRule-based download routine."""

import datetime

import geopandas as gpd
import pandas as pd
import pytest
import shapely
from shapely.geometry import Point, box


@pytest.mark.unit
def test_query_parms_construction(monkeypatch):
    """sliderule.run must receive correct atl13.coord, poly, t0, t1, atl13_fields."""
    from HydroEO.satellites import icesat2
    from HydroEO.waterbody import HydroEODownloadError  # noqa: F401 — verify importable

    captured = {}

    def fake_run(endpoint, parms):
        captured["endpoint"] = endpoint
        captured["parms"] = parms
        # Return a non-empty GDF so query() doesn't raise HydroEODownloadError
        gdf = gpd.GeoDataFrame(
            {
                "ht_ortho": [1885.0],
                "cycle": [5],
                "gt": ["gt1l"],
                "rgt": [30],
                "ht_water_surf": [1885.0],
                "stdev_water_surf": [0.05],
                "water_depth": [3.0],
                "spot": [1],
                "srcid": [0],
                "segment_id_beg": [100],
            },
            geometry=[Point(36.35, -0.79)],
            crs="EPSG:4326",
            index=pd.DatetimeIndex(["2019-03-01T12:00:00Z"]),
        )
        return gdf

    monkeypatch.setattr(icesat2.sliderule, "run", fake_run)

    aoi = [
        (36.28, -0.87),
        (36.43, -0.87),
        (36.43, -0.72),
        (36.28, -0.72),
        (36.28, -0.87),
    ]
    startdate = datetime.date(2019, 1, 1)
    enddate = datetime.date(2019, 12, 31)
    fields = ["segment_quality"]

    icesat2.query(
        aoi=aoi,
        startdate=startdate,
        enddate=enddate,
        download_directory=None,
        atl13_fields=fields,
    )

    poly = shapely.Polygon(aoi)
    expected_centroid = poly.centroid

    parms = captured["parms"]
    assert captured["endpoint"] == "atl13x"
    assert "atl13" in parms
    assert "coord" in parms["atl13"]
    assert abs(parms["atl13"]["coord"]["lon"] - expected_centroid.x) < 1e-9
    assert abs(parms["atl13"]["coord"]["lat"] - expected_centroid.y) < 1e-9
    assert "poly" not in parms, (
        "poly must not be passed (causes empty results in SR v5)"
    )
    assert parms["t0"] == "2019-01-01T00:00:00Z"
    assert parms["t1"] == "2019-12-31T00:00:00Z"
    assert parms.get("atl13_fields") == fields


@pytest.mark.unit
def test_column_mapping():
    """SlideRule column names must be renamed to HydroEO schema names."""
    from HydroEO.satellites.icesat2 import (
        SR_ANCILLARY_COLUMN_MAP,
        SR_DEFAULT_COLUMN_MAP,
    )

    raw = gpd.GeoDataFrame(
        {
            "ht_ortho": [1885.0],
            "cycle": [5],
            "gt": ["gt1l"],
            "rgt": [30],
            "segment_quality": [0],
        },
        geometry=[Point(36.35, -0.79)],
        crs="EPSG:4326",
    )

    rename_map = {
        **{k: v for k, v in SR_DEFAULT_COLUMN_MAP.items() if k in raw.columns},
        **{k: v for k, v in SR_ANCILLARY_COLUMN_MAP.items() if k in raw.columns},
    }
    renamed = raw.rename(columns=rename_map)

    assert "height" in renamed.columns
    assert "cycle_number" in renamed.columns
    assert "beam" in renamed.columns
    assert "rgt" in renamed.columns
    assert "quality_seg" in renamed.columns
    assert "ht_ortho" not in renamed.columns
    assert "cycle" not in renamed.columns
    assert "gt" not in renamed.columns


@pytest.mark.unit
def test_empty_result_raises(monkeypatch):
    """Empty GeoDataFrame from sliderule.run must raise HydroEODownloadError mentioning HydroLAKES."""
    from HydroEO.satellites import icesat2
    from HydroEO.waterbody import HydroEODownloadError

    monkeypatch.setattr(icesat2.sliderule, "run", lambda *a, **kw: gpd.GeoDataFrame())

    aoi = [
        (36.28, -0.87),
        (36.43, -0.87),
        (36.43, -0.72),
        (36.28, -0.72),
        (36.28, -0.87),
    ]
    with pytest.raises(HydroEODownloadError, match="HydroLAKES"):
        icesat2.query(
            aoi=aoi,
            startdate=datetime.date(2019, 1, 1),
            enddate=datetime.date(2019, 12, 31),
            download_directory=None,
        )


@pytest.mark.unit
def test_connection_error_raises(monkeypatch):
    """Exception from sliderule.run must propagate as HydroEODownloadError."""
    from HydroEO.satellites import icesat2
    from HydroEO.waterbody import HydroEODownloadError

    def bad_run(*a, **kw):
        raise ConnectionError("Network unreachable")

    monkeypatch.setattr(icesat2.sliderule, "run", bad_run)

    aoi = [
        (36.28, -0.87),
        (36.43, -0.87),
        (36.43, -0.72),
        (36.28, -0.72),
        (36.28, -0.87),
    ]
    with pytest.raises(HydroEODownloadError, match="SlideRule atl13x request failed"):
        icesat2.query(
            aoi=aoi,
            startdate=datetime.date(2019, 1, 1),
            enddate=datetime.date(2019, 12, 31),
            download_directory=None,
        )


@pytest.mark.unit
def test_ams_400_raises_descriptive_error(monkeypatch):
    """AMS HTTP 400 (water body not registered) must raise HydroEODownloadError mentioning AMS."""
    from HydroEO.satellites import icesat2
    from HydroEO.waterbody import HydroEODownloadError

    def ams_failure(*a, **kw):
        # SlideRule raises FileNotFoundError on a /tmp/ path when the AMS returns 400.
        raise FileNotFoundError("/tmp/tmpXXXXXX")

    monkeypatch.setattr(icesat2.sliderule, "run", ams_failure)

    aoi = [
        (36.28, -0.87),
        (36.43, -0.87),
        (36.43, -0.72),
        (36.28, -0.72),
        (36.28, -0.87),
    ]
    with pytest.raises(HydroEODownloadError, match="SlideRule lookup failed"):
        icesat2.query(
            aoi=aoi,
            startdate=datetime.date(2019, 1, 1),
            enddate=datetime.date(2019, 12, 31),
            download_directory=None,
        )


@pytest.mark.unit
def test_invalid_atl13_field_name_raises_descriptive_error(tmp_path):
    """Invalid ATL13 fields must be rejected at config validation time."""
    import yaml

    from HydroEO.project import Project

    reservoir_path = tmp_path / "reservoirs.shp"
    gpd.GeoDataFrame(
        {"project": ["r1"], "geometry": [box(36.0, -1.0, 36.5, -0.5)]},
        crs="EPSG:4326",
    ).to_file(str(reservoir_path))

    cfg = {
        "project": {"main_dir": str(tmp_path / "out")},
        "reservoirs": {"path": str(reservoir_path), "id_key": "project"},
        "icesat2": {
            "download": False,
            "process": False,
            "startdate": [2019, 1, 1],
            "enddate": [2019, 12, 31],
            "atl13_fields": [
                "segment_quality",
                "photon_rate",
            ],  # photon_rate is invalid
        },
    }
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text(yaml.dump(cfg))

    with pytest.raises(ValueError, match="Invalid ATL13 fields"):
        Project(name="test", config=str(cfg_path))
