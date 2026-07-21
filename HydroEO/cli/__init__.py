"""HydroEO command-line interface.

Two command groups:
  hydroeo <pipeline-cmd> CONFIG   — config-driven pipeline (wraps Project)
  hydroeo fetch <satellite> ...   — direct per-satellite downloads (no config file)
"""

from __future__ import annotations

import datetime
import logging
import os
from pathlib import Path
from typing import List, Optional

import typer

from HydroEO.logging_config import setup_logging

app = typer.Typer(
    name="hydroeo",
    help="HydroEO: satellite altimetry for water resource applications.",
    no_args_is_help=True,
    add_completion=False,
)

fetch_app = typer.Typer(
    name="fetch",
    help="Direct per-satellite downloads (no config file required).",
    no_args_is_help=True,
    add_completion=False,
)
app.add_typer(fetch_app, name="fetch")


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


def _configure_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    setup_logging(level=level)


def _parse_date(value: str, label: str) -> datetime.date:
    try:
        return datetime.date.fromisoformat(value)
    except ValueError:
        typer.echo(
            f"Error: {label} must be in YYYY-MM-DD format, got '{value}'.", err=True
        )
        raise typer.Exit(code=1)


def _resolve_earthaccess_creds(
    username: Optional[str],
    password: Optional[str],
) -> tuple:
    u = username or os.environ.get("EARTHACCESS_USERNAME")
    p = password or os.environ.get("EARTHACCESS_PASSWORD")
    return u, p


def _resolve_creodias_creds(
    username: Optional[str],
    password: Optional[str],
) -> tuple:
    u = username or os.environ.get("CREODIAS_USERNAME")
    p = password or os.environ.get("CREODIAS_PASSWORD")
    return u, p


def _parse_bbox(value: str) -> List[float]:
    """Parse a quoted bbox string 'MINLON MINLAT MAXLON MAXLAT' into floats."""
    try:
        parts = value.split()
        if len(parts) != 4:
            raise ValueError
        return [float(p) for p in parts]
    except ValueError:
        typer.echo(
            f"Error: --bbox must be 'MINLON MINLAT MAXLON MAXLAT', got '{value}'.",
            err=True,
        )
        raise typer.Exit(code=1)


def _bbox_to_coord_list(bbox: List[float]) -> list:
    minlon, minlat, maxlon, maxlat = bbox
    return [
        (minlon, minlat),
        (maxlon, minlat),
        (maxlon, maxlat),
        (minlon, maxlat),
        (minlon, minlat),
    ]


def _make_project(config: Path, name: Optional[str]):
    from HydroEO.project import Project

    project_name = name or config.stem
    return Project(name=project_name, config=str(config))


# ---------------------------------------------------------------------------
# Pipeline commands (config-driven)
# ---------------------------------------------------------------------------


@app.command("run")
def cmd_run(
    config: Path = typer.Argument(
        ...,
        exists=True,
        readable=True,
        resolve_path=True,
        help="Path to a HydroEO YAML config file.",
    ),
    name: Optional[str] = typer.Option(
        None, "--name", "-n", help="Project name (defaults to config file stem)."
    ),
    verbose: bool = typer.Option(
        False, "--verbose", "-v", help="Enable debug logging."
    ),
):
    """Initialize and download (mirrors run.py)."""
    _configure_logging(verbose)
    prj = _make_project(config, name)
    prj.initialize()
    prj.download()


@app.command("initialize")
def cmd_initialize(
    config: Path = typer.Argument(
        ...,
        exists=True,
        readable=True,
        resolve_path=True,
        help="Path to a HydroEO YAML config file.",
    ),
    name: Optional[str] = typer.Option(
        None, "--name", "-n", help="Project name (defaults to config file stem)."
    ),
    verbose: bool = typer.Option(
        False, "--verbose", "-v", help="Enable debug logging."
    ),
):
    """Validate config, download PLD, assign waterbody IDs."""
    _configure_logging(verbose)
    prj = _make_project(config, name)
    prj.initialize()


@app.command("download")
def cmd_download(
    config: Path = typer.Argument(
        ...,
        exists=True,
        readable=True,
        resolve_path=True,
        help="Path to a HydroEO YAML config file.",
    ),
    name: Optional[str] = typer.Option(
        None, "--name", "-n", help="Project name (defaults to config file stem)."
    ),
    verbose: bool = typer.Option(
        False, "--verbose", "-v", help="Enable debug logging."
    ),
):
    """Download satellite data for the configured waterbodies."""
    _configure_logging(verbose)
    prj = _make_project(config, name)
    prj.download()


@app.command("update")
def cmd_update(
    config: Path = typer.Argument(
        ...,
        exists=True,
        readable=True,
        resolve_path=True,
        help="Path to a HydroEO YAML config file.",
    ),
    name: Optional[str] = typer.Option(
        None, "--name", "-n", help="Project name (defaults to config file stem)."
    ),
    verbose: bool = typer.Option(
        False, "--verbose", "-v", help="Enable debug logging."
    ),
):
    """Re-download data with end-date set to today."""
    _configure_logging(verbose)
    prj = _make_project(config, name)
    prj.update()


@app.command("timeseries")
def cmd_timeseries(
    config: Path = typer.Argument(
        ...,
        exists=True,
        readable=True,
        resolve_path=True,
        help="Path to a HydroEO YAML config file.",
    ),
    name: Optional[str] = typer.Option(
        None, "--name", "-n", help="Project name (defaults to config file stem)."
    ),
    verbose: bool = typer.Option(
        False, "--verbose", "-v", help="Enable debug logging."
    ),
):
    """Process downloaded data into time-series (reservoirs only)."""
    _configure_logging(verbose)
    prj = _make_project(config, name)
    prj.create_timeseries()


@app.command("summaries")
def cmd_summaries(
    config: Path = typer.Argument(
        ...,
        exists=True,
        readable=True,
        resolve_path=True,
        help="Path to a HydroEO YAML config file.",
    ),
    name: Optional[str] = typer.Option(
        None, "--name", "-n", help="Project name (defaults to config file stem)."
    ),
    verbose: bool = typer.Option(
        False, "--verbose", "-v", help="Enable debug logging."
    ),
):
    """Generate plots and summaries (reservoirs and rivers)."""
    _configure_logging(verbose)
    prj = _make_project(config, name)
    prj.generate_summaries()


# ---------------------------------------------------------------------------
# fetch swot-raster
# ---------------------------------------------------------------------------


@fetch_app.command("swot-raster")
def fetch_swot_raster(
    bbox: str = typer.Option(
        ...,
        "--bbox",
        help="Bounding box as a quoted string: 'MINLON MINLAT MAXLON MAXLAT'.",
    ),
    start: str = typer.Option(..., "--start", help="Start date (YYYY-MM-DD)."),
    end: str = typer.Option(..., "--end", help="End date (YYYY-MM-DD)."),
    aoi_name: str = typer.Option("aoi", "--aoi-name", help="Label for the AOI."),
    product: str = typer.Option(
        "SWOT_L2_HR_PIXC_2.0", "--product", help="SWOT raster product short name."
    ),
    output: str = typer.Option(
        "hydroeo_output", "--output", "-o", help="Output directory (created if absent)."
    ),
    username: Optional[str] = typer.Option(
        None, "--username", help="EarthAccess username (or set EARTHACCESS_USERNAME)."
    ),
    password: Optional[str] = typer.Option(
        None,
        "--password",
        help="EarthAccess password (or set EARTHACCESS_PASSWORD).",
        hide_input=True,
    ),
    verbose: bool = typer.Option(
        False, "--verbose", "-v", help="Enable debug logging."
    ),
):
    """Download SWOT raster data for a bounding box and time range."""
    _configure_logging(verbose)
    parsed_bbox = _parse_bbox(bbox)
    startdate = _parse_date(start, "--start")
    enddate = _parse_date(end, "--end")
    u, p = _resolve_earthaccess_creds(username, password)
    if u:
        os.environ["EARTHDATA_USERNAME"] = u
    if p:
        os.environ["EARTHDATA_PASSWORD"] = p

    config = {
        "aoi": {
            "name": aoi_name,
            "type": "bbox",
            "bbox": parsed_bbox,
        },
        "product": product,
        "startdate": [startdate.year, startdate.month, startdate.day],
        "enddate": [enddate.year, enddate.month, enddate.day],
    }

    from HydroEO.satellites.swot.raster import download_raster

    download_raster(config=config, project_dir=output)
    typer.echo(f"SWOT raster download complete. Output: {output}")


# ---------------------------------------------------------------------------
# fetch swot-pixc
# ---------------------------------------------------------------------------


@fetch_app.command("swot-pixc")
def fetch_swot_pixc(
    bbox: str = typer.Option(
        ...,
        "--bbox",
        help="Bounding box as a quoted string: 'MINLON MINLAT MAXLON MAXLAT'.",
    ),
    start: str = typer.Option(..., "--start", help="Start date (YYYY-MM-DD)."),
    end: str = typer.Option(..., "--end", help="End date (YYYY-MM-DD)."),
    aoi_name: str = typer.Option("aoi", "--aoi-name", help="Label for the AOI."),
    product: str = typer.Option(
        "SWOT_L2_HR_PIXC_D", "--product", help="SWOT pixel cloud product short name."
    ),
    output: str = typer.Option(
        "hydroeo_output", "--output", "-o", help="Output directory (created if absent)."
    ),
    classes: Optional[str] = typer.Option(
        None,
        "--classes",
        help="Comma-separated water classes to include (e.g., 'open_water,water_near_land'). "
        "Default: open_water,water_near_land.",
    ),
    fields: Optional[str] = typer.Option(
        None,
        "--fields",
        help="Comma-separated fields to grid (e.g., 'heightEGM,height'). Default: heightEGM.",
    ),
    grid_resolution: int = typer.Option(
        100, "--grid-resolution", help="Grid resolution in meters. Default: 100."
    ),
    stat_method: str = typer.Option(
        "median",
        "--stat-method",
        help="Binning statistic method (median, mean, min, max, count). Default: median.",
    ),
    username: Optional[str] = typer.Option(
        None, "--username", help="EarthAccess username (or set EARTHACCESS_USERNAME)."
    ),
    password: Optional[str] = typer.Option(
        None,
        "--password",
        help="EarthAccess password (or set EARTHACCESS_PASSWORD).",
        hide_input=True,
    ),
    verbose: bool = typer.Option(
        False, "--verbose", "-v", help="Enable debug logging."
    ),
):
    """Download SWOT pixel cloud data, preprocess, and grid to rasters.

    Point cloud data is filtered by water class, extracted with heightEGM computation,
    and gridded to regular GeoTIFF rasters using specified statistics.
    """
    _configure_logging(verbose)
    parsed_bbox = _parse_bbox(bbox)
    startdate = _parse_date(start, "--start")
    enddate = _parse_date(end, "--end")
    u, p = _resolve_earthaccess_creds(username, password)
    if u:
        os.environ["EARTHDATA_USERNAME"] = u
    if p:
        os.environ["EARTHDATA_PASSWORD"] = p

    # Parse optional classes and fields
    classes_list = (
        [c.strip() for c in classes.split(",")]
        if classes
        else ["open_water", "water_near_land"]
    )
    fields_list = [f.strip() for f in fields.split(",")] if fields else ["heightEGM"]

    config = {
        "aoi": {
            "name": aoi_name,
            "type": "bbox",
            "bbox": parsed_bbox,
        },
        "product": product,
        "startdate": [startdate.year, startdate.month, startdate.day],
        "enddate": [enddate.year, enddate.month, enddate.day],
        "classes": classes_list,
        "fields": fields_list,
        "grid_resolution": grid_resolution,
        "stat_method": stat_method,
    }

    from HydroEO.satellites.swot.pixc import download_pixc

    download_pixc(config=config, project_dir=output)
    typer.echo(f"SWOT pixel cloud download complete. Output: {output}")


# ---------------------------------------------------------------------------
# fetch swot-lake
# ---------------------------------------------------------------------------


@fetch_app.command("swot-lake")
def fetch_swot_lake(
    bbox: str = typer.Option(
        ...,
        "--bbox",
        help="Bounding box as a quoted string: 'MINLON MINLAT MAXLON MAXLAT'.",
    ),
    start: str = typer.Option(..., "--start", help="Start date (YYYY-MM-DD)."),
    end: str = typer.Option(..., "--end", help="End date (YYYY-MM-DD)."),
    output: str = typer.Option(
        "hydroeo_output", "--output", "-o", help="Output directory (created if absent)."
    ),
    username: Optional[str] = typer.Option(
        None, "--username", help="EarthAccess username (or set EARTHACCESS_USERNAME)."
    ),
    password: Optional[str] = typer.Option(
        None,
        "--password",
        help="EarthAccess password (or set EARTHACCESS_PASSWORD).",
        hide_input=True,
    ),
    verbose: bool = typer.Option(
        False, "--verbose", "-v", help="Enable debug logging."
    ),
):
    """Download SWOT lake (LakeSP) data for a bounding box and time range."""
    _configure_logging(verbose)
    parsed_bbox = _parse_bbox(bbox)
    startdate = _parse_date(start, "--start")
    enddate = _parse_date(end, "--end")

    u, p = _resolve_earthaccess_creds(username, password)
    if u:
        os.environ["EARTHDATA_USERNAME"] = u
    if p:
        os.environ["EARTHDATA_PASSWORD"] = p

    aoi = _bbox_to_coord_list(parsed_bbox)
    os.makedirs(output, exist_ok=True)

    from HydroEO.satellites import swot

    results = swot.query(aoi=aoi, startdate=startdate, enddate=enddate)
    swot.download(results=results, download_directory=output)
    typer.echo(f"SWOT lake download complete. Output: {output}")


# ---------------------------------------------------------------------------
# fetch icesat2
# ---------------------------------------------------------------------------


@fetch_app.command("icesat2")
def fetch_icesat2(
    bbox: str = typer.Option(
        ...,
        "--bbox",
        help="Bounding box as a quoted string: 'MINLON MINLAT MAXLON MAXLAT'.",
    ),
    start: str = typer.Option(..., "--start", help="Start date (YYYY-MM-DD)."),
    end: str = typer.Option(..., "--end", help="End date (YYYY-MM-DD)."),
    output: str = typer.Option(
        "hydroeo_output", "--output", "-o", help="Output directory (created if absent)."
    ),
    verbose: bool = typer.Option(
        False, "--verbose", "-v", help="Enable debug logging."
    ),
):
    """Download ICESat-2 ATL13 water surface elevations via SlideRule.

    The bounding box centroid is used for the ATL13x AMS water-body lookup.
    Output is written as atl13.parquet inside OUTPUT.
    """
    _configure_logging(verbose)
    parsed_bbox = _parse_bbox(bbox)
    startdate = _parse_date(start, "--start")
    enddate = _parse_date(end, "--end")

    aoi = _bbox_to_coord_list(parsed_bbox)
    os.makedirs(output, exist_ok=True)

    from HydroEO.satellites import icesat2

    gdf = icesat2.query(
        aoi=aoi,
        startdate=startdate,
        enddate=enddate,
        download_directory=output,
    )
    typer.echo(
        f"ICESat-2 download complete. {len(gdf)} records written to "
        f"{os.path.join(output, 'atl13.parquet')}"
    )


# ---------------------------------------------------------------------------
# fetch sentinel
# ---------------------------------------------------------------------------


@fetch_app.command("sentinel")
def fetch_sentinel(
    bbox: str = typer.Option(
        ...,
        "--bbox",
        help="Bounding box as a quoted string: 'MINLON MINLAT MAXLON MAXLAT'.",
    ),
    start: str = typer.Option(..., "--start", help="Start date (YYYY-MM-DD)."),
    end: str = typer.Option(..., "--end", help="End date (YYYY-MM-DD)."),
    output: str = typer.Option(
        "hydroeo_output", "--output", "-o", help="Output directory (created if absent)."
    ),
    product: str = typer.Option(
        "S3", "--product", help="Sentinel product: S3 (Sentinel-3) or S6 (Sentinel-6)."
    ),
    creodias_username: Optional[str] = typer.Option(
        None,
        "--creodias-username",
        help="CREODIAS username (or set CREODIAS_USERNAME).",
    ),
    creodias_password: Optional[str] = typer.Option(
        None,
        "--creodias-password",
        help="CREODIAS password (or set CREODIAS_PASSWORD).",
        hide_input=True,
    ),
    verbose: bool = typer.Option(
        False, "--verbose", "-v", help="Enable debug logging."
    ),
):
    """Download Sentinel-3 or Sentinel-6 data for a bounding box and time range."""
    _configure_logging(verbose)
    parsed_bbox = _parse_bbox(bbox)

    product = product.upper()
    if product not in ("S3", "S6"):
        typer.echo("Error: --product must be 'S3' or 'S6'.", err=True)
        raise typer.Exit(code=1)

    startdate = _parse_date(start, "--start")
    enddate = _parse_date(end, "--end")

    u, p = _resolve_creodias_creds(creodias_username, creodias_password)
    if not u or not p:
        typer.echo(
            "Error: CREODIAS credentials required. Pass --creodias-username/--creodias-password "
            "or set CREODIAS_USERNAME / CREODIAS_PASSWORD environment variables.",
            err=True,
        )
        raise typer.Exit(code=1)

    aoi = _bbox_to_coord_list(parsed_bbox)
    os.makedirs(output, exist_ok=True)

    from HydroEO.satellites import sentinel

    ids = sentinel.query(
        aoi=aoi,
        startdate=startdate,
        enddate=enddate,
        product=product,
        creodias_credentials=(u, p),
    )
    typer.echo(f"Found {len(ids)} Sentinel-{product[1]} products.")

    sentinel.download(
        ids=ids,
        download_directory=output,
        creodias_credentials=(u, p),
    )
    typer.echo(f"Sentinel download complete. Output: {output}")


# ---------------------------------------------------------------------------
# fetch cop-dem
# ---------------------------------------------------------------------------


@fetch_app.command("cop-dem")
def fetch_cop_dem(
    bbox: str = typer.Option(
        ...,
        "--bbox",
        help="Bounding box as a quoted string: 'MINLON MINLAT MAXLON MAXLAT'.",
    ),
    output: str = typer.Option(
        "hydroeo_output", "--output", "-o", help="Output directory (created if absent)."
    ),
    cdse_username: Optional[str] = typer.Option(
        None,
        "--cdse-username",
        help="CDSE username (or set CDSE_USERNAME).",
    ),
    cdse_password: Optional[str] = typer.Option(
        None,
        "--cdse-password",
        help="CDSE password (or set CDSE_PASSWORD).",
        hide_input=True,
    ),
    dataset: str = typer.Option(
        "DEM30",
        "--dataset",
        help=(
            "Comma-separated list of product layers to download. "
            "Valid values: DEM30 (30 m elevation, default), DEM90 (90 m elevation), "
            "EDM (Editing Mask), FLM (Filling Mask), HEM (Height Error Mask), "
            "WBM (Water Body Mask). "
            "DEM30 and DEM90 cannot be combined. "
            "Example: --dataset DEM30,WBM,EDM"
        ),
    ),
    output_filename: str = typer.Option(
        "cop_dem_merged",
        "--output-filename",
        help=(
            "Base name for the output GeoTIFFs (without extension). "
            "Each layer produces {output_filename}_{LAYER}.tif. "
            "Default: cop_dem_merged  →  cop_dem_merged_DEM30.tif, cop_dem_merged_WBM.tif, …"
        ),
    ),
    verbose: bool = typer.Option(
        False, "--verbose", "-v", help="Enable debug logging."
    ),
):
    """Download, merge, and clip COP-DEM layers for a bounding box.

    Downloads Copernicus DEM tiles from CDSE, merges them, and clips to the
    specified bounding box.  One GeoTIFF is produced per requested layer.

    Available layers: DEM30, DEM90, EDM, FLM, HEM, WBM.
    ZIP archives are downloaded once and all layers are extracted in a single
    pass, so requesting multiple layers incurs no extra download cost.

    Requires CDSE credentials (register at https://dataspace.copernicus.eu/
    and request CCM data access).
    """
    _configure_logging(verbose)
    parsed_bbox = _parse_bbox(bbox)

    u = cdse_username or os.environ.get("CDSE_USERNAME")
    p = cdse_password or os.environ.get("CDSE_PASSWORD")

    if not u or not p:
        typer.echo(
            "Error: CDSE credentials required. Pass --cdse-username/--cdse-password "
            "or set CDSE_USERNAME / CDSE_PASSWORD environment variables.",
            err=True,
        )
        raise typer.Exit(code=1)

    from HydroEO.downloaders.dem import VALID_LAYERS, download_cop_dem

    layer_list = [s.strip().upper() for s in dataset.split(",") if s.strip()]

    # Validate early so the user gets a clear error before any download starts
    unknown = set(layer_list) - VALID_LAYERS
    if unknown:
        typer.echo(
            f"Error: Unknown layer(s): {sorted(unknown)}. "
            f"Valid options: {sorted(VALID_LAYERS)}",
            err=True,
        )
        raise typer.Exit(code=1)
    if "DEM30" in layer_list and "DEM90" in layer_list:
        typer.echo(
            "Error: DEM30 and DEM90 cannot be requested together — "
            "they target different CDSE datasets (GLO-30 vs GLO-90).",
            err=True,
        )
        raise typer.Exit(code=1)

    minx, miny, maxx, maxy = parsed_bbox
    os.makedirs(output, exist_ok=True)

    output_paths = download_cop_dem(
        minx=minx,
        miny=miny,
        maxx=maxx,
        maxy=maxy,
        output_dir=output,
        username=u,
        password=p,
        layers=layer_list,
        output_basename=output_filename,
    )
    for layer, path in output_paths.items():
        typer.echo(f"COP-DEM [{layer}] complete. Output: {path}")


if __name__ == "__main__":
    app()
