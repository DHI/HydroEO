"""Unit tests for the HydroEO CLI (hydroeo command group)."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from typer.testing import CliRunner

from HydroEO.cli import app

runner = CliRunner()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


@pytest.fixture
def config_file(tmp_path: Path) -> Path:
    """Write a minimal valid YAML config so the CLI can resolve the path."""
    cfg = tmp_path / "test.yaml"
    cfg.write_text("project:\n  main_dir: /tmp/out\n", encoding="utf-8")
    return cfg


# ---------------------------------------------------------------------------
# Global help / structure
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_help_exits_zero():
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    assert "hydroeo" in result.output.lower()


@pytest.mark.unit
def test_fetch_help_exits_zero():
    result = runner.invoke(app, ["fetch", "--help"])
    assert result.exit_code == 0
    assert "swot-raster" in result.output
    assert "swot-lake" in result.output
    assert "icesat2" in result.output
    assert "sentinel" in result.output


# ---------------------------------------------------------------------------
# Pipeline commands
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.parametrize(
    "cmd,method",
    [
        (["run"], ["initialize", "download"]),
        (["initialize"], ["initialize"]),
        (["download"], ["download"]),
        (["update"], ["update"]),
        (["timeseries"], ["create_timeseries"]),
        (["summaries"], ["generate_summaries"]),
    ],
)
def test_pipeline_command_calls_project(config_file, cmd, method):
    mock_prj = MagicMock()
    with patch("HydroEO.cli._make_project", return_value=mock_prj):
        result = runner.invoke(app, cmd + [str(config_file)])
    assert result.exit_code == 0, result.output
    for m in method:
        getattr(mock_prj, m).assert_called_once()


@pytest.mark.unit
def test_pipeline_uses_config_stem_as_default_name(config_file):
    """Project name defaults to the config file stem when --name is omitted."""
    mock_prj = MagicMock()
    with patch("HydroEO.cli._make_project", return_value=mock_prj) as mock_make:
        result = runner.invoke(app, ["initialize", str(config_file)])
    assert result.exit_code == 0
    mock_make.assert_called_once_with(config_file, None)


@pytest.mark.unit
def test_pipeline_custom_name(config_file):
    mock_prj = MagicMock()
    with patch("HydroEO.cli._make_project", return_value=mock_prj) as mock_make:
        result = runner.invoke(
            app, ["initialize", str(config_file), "--name", "my_project"]
        )
    assert result.exit_code == 0
    mock_make.assert_called_once_with(config_file, "my_project")


@pytest.mark.unit
def test_pipeline_missing_config_exits_nonzero():
    result = runner.invoke(app, ["run", "/nonexistent/config.yaml"])
    assert result.exit_code != 0


# ---------------------------------------------------------------------------
# fetch swot-raster
# ---------------------------------------------------------------------------

_COMMON_FETCH_ARGS = [
    "--bbox",
    "-10 40 10 60",
    "--start",
    "2023-01-01",
    "--end",
    "2023-06-01",
]


@pytest.mark.unit
def test_fetch_swot_raster_calls_download_raster(tmp_path):
    with patch("HydroEO.satellites.swot.raster.download_raster") as mock_dl:
        result = runner.invoke(
            app,
            [
                "fetch",
                "swot-raster",
                *_COMMON_FETCH_ARGS,
                "--output",
                str(tmp_path),
            ],
        )
    assert result.exit_code == 0, result.output
    mock_dl.assert_called_once()
    call_kwargs = mock_dl.call_args.kwargs
    assert call_kwargs["project_dir"] == str(tmp_path)
    cfg = call_kwargs["config"]
    assert cfg["aoi"]["type"] == "bbox"
    assert cfg["aoi"]["bbox"] == [-10.0, 40.0, 10.0, 60.0]


@pytest.mark.unit
def test_fetch_swot_raster_credentials_from_env(tmp_path, monkeypatch):
    monkeypatch.setenv("EARTHACCESS_USERNAME", "envuser")
    monkeypatch.setenv("EARTHACCESS_PASSWORD", "envpass")
    with patch("HydroEO.satellites.swot.raster.download_raster") as mock_dl:
        result = runner.invoke(
            app,
            [
                "fetch",
                "swot-raster",
                *_COMMON_FETCH_ARGS,
                "--output",
                str(tmp_path),
            ],
        )
    assert result.exit_code == 0, result.output
    creds = mock_dl.call_args.kwargs["credentials"]
    assert creds == ("envuser", "envpass")


@pytest.mark.unit
def test_fetch_swot_raster_credentials_from_args(tmp_path):
    with patch("HydroEO.satellites.swot.raster.download_raster") as mock_dl:
        result = runner.invoke(
            app,
            [
                "fetch",
                "swot-raster",
                *_COMMON_FETCH_ARGS,
                "--output",
                str(tmp_path),
                "--username",
                "arguser",
                "--password",
                "argpass",
            ],
        )
    assert result.exit_code == 0, result.output
    creds = mock_dl.call_args.kwargs["credentials"]
    assert creds == ("arguser", "argpass")


@pytest.mark.unit
@pytest.mark.parametrize("bad_date_flag", ["--start", "--end"])
def test_fetch_swot_raster_bad_date_exits_nonzero(tmp_path, bad_date_flag):
    args = [
        "fetch",
        "swot-raster",
        "--bbox",
        "-10 40 10 60",
        "--start",
        "2023-01-01",
        "--end",
        "2023-06-01",
        "--output",
        str(tmp_path),
    ]
    # Replace the date value for the targeted flag
    flag_idx = args.index(bad_date_flag)
    args[flag_idx + 1] = "not-a-date"
    with patch("HydroEO.satellites.swot.raster.download_raster"):
        result = runner.invoke(app, args)
    assert result.exit_code != 0


@pytest.mark.unit
@pytest.mark.parametrize("bad_bbox", ["1 2 3", "1 2 3 4 5", "a b c d", ""])
def test_fetch_swot_raster_bad_bbox_exits_nonzero(tmp_path, bad_bbox):
    with patch("HydroEO.satellites.swot.raster.download_raster"):
        result = runner.invoke(
            app,
            [
                "fetch",
                "swot-raster",
                "--bbox",
                bad_bbox,
                "--start",
                "2023-01-01",
                "--end",
                "2023-06-01",
                "--output",
                str(tmp_path),
            ],
        )
    assert result.exit_code != 0


# ---------------------------------------------------------------------------
# fetch swot-lake
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_fetch_swot_lake_calls_query_and_download(tmp_path):
    mock_results = MagicMock()
    with (
        patch("HydroEO.satellites.swot.query", return_value=mock_results) as mock_q,
        patch("HydroEO.satellites.swot.download") as mock_dl,
    ):
        result = runner.invoke(
            app,
            [
                "fetch",
                "swot-lake",
                *_COMMON_FETCH_ARGS,
                "--output",
                str(tmp_path),
            ],
        )
    assert result.exit_code == 0, result.output
    mock_q.assert_called_once()
    mock_dl.assert_called_once()
    assert mock_dl.call_args.kwargs["download_directory"] == str(tmp_path)


# ---------------------------------------------------------------------------
# fetch icesat2
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_fetch_icesat2_calls_query(tmp_path):
    import geopandas as gpd
    from shapely.geometry import Point

    fake_gdf = gpd.GeoDataFrame(
        {"height": [100.0]},
        geometry=[Point(0, 0)],
        crs="EPSG:4326",
    )
    with patch("HydroEO.satellites.icesat2.query", return_value=fake_gdf) as mock_q:
        result = runner.invoke(
            app,
            [
                "fetch",
                "icesat2",
                *_COMMON_FETCH_ARGS,
                "--output",
                str(tmp_path),
            ],
        )
    assert result.exit_code == 0, result.output
    mock_q.assert_called_once()
    call_kwargs = mock_q.call_args.kwargs
    assert call_kwargs["download_directory"] == str(tmp_path)
    assert "1 records" in result.output


# ---------------------------------------------------------------------------
# fetch sentinel
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_fetch_sentinel_calls_query_and_download(tmp_path, monkeypatch):
    monkeypatch.setenv("CREODIAS_USERNAME", "cuser")
    monkeypatch.setenv("CREODIAS_PASSWORD", "cpass")
    with (
        patch(
            "HydroEO.satellites.sentinel.query", return_value=["id1", "id2"]
        ) as mock_q,
        patch("HydroEO.satellites.sentinel.download") as mock_dl,
    ):
        result = runner.invoke(
            app,
            [
                "fetch",
                "sentinel",
                *_COMMON_FETCH_ARGS,
                "--output",
                str(tmp_path),
            ],
        )
    assert result.exit_code == 0, result.output
    mock_q.assert_called_once()
    mock_dl.assert_called_once()
    assert mock_dl.call_args.kwargs["creodias_credentials"] == ("cuser", "cpass")


@pytest.mark.unit
def test_fetch_sentinel_missing_creodias_creds_exits_nonzero(tmp_path, monkeypatch):
    monkeypatch.delenv("CREODIAS_USERNAME", raising=False)
    monkeypatch.delenv("CREODIAS_PASSWORD", raising=False)
    result = runner.invoke(
        app,
        [
            "fetch",
            "sentinel",
            *_COMMON_FETCH_ARGS,
            "--output",
            str(tmp_path),
        ],
    )
    assert result.exit_code != 0


@pytest.mark.unit
def test_fetch_sentinel_invalid_product_exits_nonzero(tmp_path, monkeypatch):
    monkeypatch.setenv("CREODIAS_USERNAME", "u")
    monkeypatch.setenv("CREODIAS_PASSWORD", "p")
    result = runner.invoke(
        app,
        [
            "fetch",
            "sentinel",
            *_COMMON_FETCH_ARGS,
            "--output",
            str(tmp_path),
            "--product",
            "S9",
        ],
    )
    assert result.exit_code != 0


@pytest.mark.unit
def test_fetch_sentinel_s6_product(tmp_path, monkeypatch):
    monkeypatch.setenv("CREODIAS_USERNAME", "u")
    monkeypatch.setenv("CREODIAS_PASSWORD", "p")
    with (
        patch("HydroEO.satellites.sentinel.query", return_value=[]) as mock_q,
    ):
        result = runner.invoke(
            app,
            [
                "fetch",
                "sentinel",
                *_COMMON_FETCH_ARGS,
                "--output",
                str(tmp_path),
                "--product",
                "S6",
            ],
        )
    assert result.exit_code == 0, result.output
    assert mock_q.call_args.kwargs["product"] == "S6"
