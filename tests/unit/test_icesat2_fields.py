"""Unit tests for configurable ICESat-2 ATL13 field selection."""

import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import box


@pytest.mark.unit
def test_extract_observations_passes_custom_atl13_fields(monkeypatch, tmp_path):
    """Custom atl13_fields should be passed through to ATL13 reader instances."""
    from HydroEO.satellites import icesat2

    src_dir = tmp_path / "src"
    src_dir.mkdir(parents=True, exist_ok=True)
    (src_dir / "sample.h5").write_text("placeholder", encoding="utf-8")

    features = gpd.GeoDataFrame(
        {"id": [1], "geometry": [box(6.0, 46.2, 6.9, 46.6)]},
        crs="EPSG:4326",
    )

    captured = []

    class FakeATL13:
        def __init__(self, infile, track_key, fields=None):
            captured.append(list(fields or []))
            self._track_key = track_key

        def check_height_data(self):
            return True

        def read(self):
            return pd.DataFrame(
                {
                    "height": [10.0],
                    "lat": [46.4],
                    "lon": [6.4],
                    "date": [pd.Timestamp("2024-01-01")],
                    "beam": [self._track_key],
                }
            )

    monkeypatch.setattr(icesat2, "ATL13", FakeATL13)
    monkeypatch.setattr(gpd.GeoDataFrame, "to_file", lambda self, _dst: None)

    custom_fields = [
        "ht_ortho",
        "segment_lat",
        "segment_lon",
        "delta_time",
        "segment_quality",
    ]
    track_keys = ["gt1l"]

    icesat2.extract_observations(
        src_dir=str(src_dir),
        dst_path=str(tmp_path / "out.shp"),
        features=features,
        atl13_fields=custom_fields,
        track_keys=track_keys,
    )

    assert captured, "Expected ATL13 reader to be called at least once"
    assert captured[0] == custom_fields


@pytest.mark.unit
def test_invalid_atl13_field_name_raises_descriptive_error(tmp_path):
    """Invalid ATL13 fields should be rejected at config validation time."""
    import yaml
    from HydroEO.project import Project

    reservoir_path = tmp_path / "reservoirs.shp"
    reservoir_path.write_text("placeholder", encoding="utf-8")

    cfg = {
        "project": {"main_dir": str(tmp_path / "out")},
        "reservoirs": {"path": str(reservoir_path), "id_key": "project"},
        "icesat2": {
            "download": False,
            "process": True,
            "startdate": [2024, 1, 1],
            "enddate": [2024, 2, 1],
            "atl13_fields": ["ht_ortho", "photon_rate"],
        },
    }

    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")

    with pytest.raises(ValueError, match="Invalid ATL13 fields"):
        Project(name="invalid-atl13", config=str(cfg_path))
