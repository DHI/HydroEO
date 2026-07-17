"""Unit tests for HydroEO.downloaders.hydroweb.download_PLD."""

from unittest.mock import patch

import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import Point, Polygon

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


@pytest.mark.unit
def test_download_pld_backfills_res_id_from_tile_files(tmp_path, caplog):
    """When a '_light' file (lake_id + geometry only) is present alongside
    full-schema continent-tile files, coverage should come from '_light'
    (guaranteed complete) while res_id gets backfilled per-lake_id from
    whichever tiles are present -- lakes whose continent tile isn't
    available should end up with res_id left NA, plus a warning naming how
    many, rather than silently missing that information."""
    import logging

    extracted_dir = tmp_path / "pld_extracted"
    extracted_dir.mkdir()
    (extracted_dir / "SWOT_LakeDatabase_light.gpkg").touch()
    (extracted_dir / "SWOT_LakeDatabase_AS.gpkg").touch()

    bounds = [0, 0, 3, 3]

    # Light file: covers everything (lake_id 1 and 2), no res_id column.
    light_gdf = gpd.GeoDataFrame(
        {"lake_id": [1, 2]},
        geometry=[Point(1, 1), Point(2, 2)],
        crs="EPSG:4326",
    )
    # AS tile: only has lake 1 (lake 2's continent tile isn't in this
    # download), but does carry res_id.
    as_tile_gdf = gpd.GeoDataFrame(
        {"lake_id": [1], "res_id": [9001]},
        geometry=[Point(1, 1)],
        crs="EPSG:4326",
    )

    def fake_read_file(filepath, *args, **kwargs):
        if "light" in str(filepath).lower():
            return light_gdf
        return as_tile_gdf

    with (
        patch.object(gpd, "read_file", side_effect=fake_read_file),
        caplog.at_level(logging.WARNING),
    ):
        hydroweb.download_PLD(
            download_dir=str(tmp_path / "output"),
            bounds=bounds,
            raw_pld_path=str(extracted_dir),
            keep_raw=True,
        )

    result = gpd.read_file(tmp_path / "output" / "PLD_subset.gpkg")
    result = result.set_index("lake_id")

    assert result.loc[1, "res_id"] == 9001.0, "lake 1's res_id should be backfilled from the AS tile"
    assert pd.isna(result.loc[2, "res_id"]), "lake 2's res_id should stay NA (no matching tile)"
    assert "1 of 2 PLD lake(s)" in caplog.text


@pytest.mark.unit
def test_download_pld_continent_codes_filters_backfill_tiles(tmp_path, caplog):
    """hydroweb.continent_codes should restrict which tile files are used
    for res_id backfill, and warn if a requested code has no matching file
    among what was actually downloaded."""
    import logging

    extracted_dir = tmp_path / "pld_extracted"
    extracted_dir.mkdir()
    (extracted_dir / "SWOT_LakeDatabase_light.gpkg").touch()
    (extracted_dir / "SWOT_LakeDatabase_AS.gpkg").touch()

    bounds = [0, 0, 3, 3]
    light_gdf = gpd.GeoDataFrame(
        {"lake_id": [1]}, geometry=[Point(1, 1)], crs="EPSG:4326"
    )
    as_tile_gdf = gpd.GeoDataFrame(
        {"lake_id": [1], "res_id": [9001]}, geometry=[Point(1, 1)], crs="EPSG:4326"
    )

    def fake_read_file(filepath, *args, **kwargs):
        return light_gdf if "light" in str(filepath).lower() else as_tile_gdf

    with (
        patch.object(gpd, "read_file", side_effect=fake_read_file),
        caplog.at_level(logging.WARNING),
    ):
        hydroweb.download_PLD(
            download_dir=str(tmp_path / "output"),
            bounds=bounds,
            raw_pld_path=str(extracted_dir),
            keep_raw=True,
            continent_codes=["AU"],  # doesn't match the AS tile that's present
        )

    result = gpd.read_file(tmp_path / "output" / "PLD_subset.gpkg").set_index("lake_id")
    assert pd.isna(result.loc[1, "res_id"]), (
        "res_id should NOT be backfilled from the AS tile since continent_codes "
        "only requested AU"
    )
    assert "['AU']" in caplog.text or "AU" in caplog.text
    assert "no matching tile file" in caplog.text
 