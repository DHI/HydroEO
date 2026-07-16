"""Unit tests for HydroEO.downloaders.hydroweb.download_PLD."""

from unittest.mock import patch

import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import Polygon

from HydroEO.downloaders import hydroweb


@pytest.mark.unit
def test_download_pld_keeps_lakes_that_extend_past_bbox_edge(tmp_path):
    """Regression test: download_PLD's bounds filter previously used
    .within(), which requires the ENTIRE lake polygon to be contained in the
    bounding box -- silently dropping any lake that extends even slightly
    past its edge. This is exactly the failure mode for a single-reservoir
    project (bounds == that one reservoir's own tight bbox) where the real
    PLD lake polygon doesn't align neatly with a buffered input polygon's
    bounding box. Fixed to use .intersects() instead, which correctly keeps
    any lake that overlaps the bbox at all.
    """
    # A pre-extracted PLD directory (raw_pld_path as a directory containing
    # .sqlite files directly) skips the download/zip-extraction path
    # entirely and goes straight to the buggy bounds-filtering logic.
    extracted_dir = tmp_path / "pld_extracted"
    extracted_dir.mkdir()
    (extracted_dir / "lakes.sqlite").touch()  # content irrelevant; gpd.read_file is mocked

    bounds = [0, 0, 1, 1]

    # lake_ok: fully inside the bbox. lake_edge: pokes 0.05 past the top edge
    # -- exactly the scenario a buffered reservoir polygon's bbox can miss.
    lake_ok = Polygon([(0.2, 0.2), (0.8, 0.2), (0.8, 0.8), (0.2, 0.8)])
    lake_edge = Polygon([(0.2, 0.2), (0.8, 0.2), (0.8, 1.05), (0.2, 1.05)])

    fake_sqlite_gdf = gpd.GeoDataFrame(
        {"res_id": [101, 102]},
        geometry=[lake_ok, lake_edge],
        crs="EPSG:4326",
    )
    # download_PLD sets gdf.index.name = "lake_id" after reading with
    # fid_as_index=True -- give the mock a named index to match.
    fake_sqlite_gdf.index = pd.Index([1, 2], name="lake_id")

    download_dir = tmp_path / "output"

    with patch.object(gpd, "read_file", return_value=fake_sqlite_gdf):
        hydroweb.download_PLD(
            download_dir=str(download_dir),
            bounds=bounds,
            raw_pld_path=str(extracted_dir),
            keep_raw=True,
        )

    result = gpd.read_file(download_dir / "PLD_subset.gpkg")
    assert set(result["res_id"]) == {101, 102}, (
        "lake_edge (which pokes past the bbox edge) was incorrectly dropped -- "
        "the within()-vs-intersects() bug regressed"
    )
