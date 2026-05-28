"""
Download COP-DEM GLO-30 (or GLO-90) tiles from CDSE for a given polygon.

Prerequisites:
    pip install requests shapely

You need a CDSE account with CCM data access:
    https://dataspace.copernicus.eu/ → register, then request CCM access.
"""

import json
import logging
import os
import shutil
import tempfile
import zipfile
from pathlib import Path

import rasterio
import rioxarray  # noqa: F401 — activates the .rio accessor on xarray objects
import requests
import xarray as xr
from rasterio.crs import CRS
from rasterio.mask import mask as rasterio_mask
from rasterio.merge import merge as rasterio_merge
from shapely.geometry import box, mapping, shape
from shapely.ops import unary_union
from shapely.wkt import dumps as wkt_dumps
from tqdm import tqdm

logger = logging.getLogger(__name__)

# ── Configuration ─────────────────────────────────────────────────────────────

CDSE_USERNAME = os.environ.get("CDSE_USERNAME", "your_username")
CDSE_PASSWORD = os.environ.get("CDSE_PASSWORD", "your_password")

# CDSE dataset identifiers by resolution
_CDSE_GLO30 = "COP-DEM_GLO-30-DGED/2023_1"
_CDSE_GLO90 = "COP-DEM_GLO-90-DGED/2023_1"

# Product layers available in each COP-DEM ZIP archive.
# Each layer maps to the filename suffix used inside the archive:
#   DEM30/DEM90 → DEM/ subfolder   (elevation raster)
#   EDM/FLM/HEM/WBM → AUXFILES/ subfolder (quality masks)
VALID_LAYERS: set[str] = {"DEM30", "DEM90", "EDM", "FLM", "HEM", "WBM"}

LAYER_SUFFIXES: dict[str, str] = {
    "DEM30": "_DEM.tif",
    "DEM90": "_DEM.tif",
    "EDM": "_EDM.tif",  # Editing Mask
    "FLM": "_FLM.tif",  # Filling Mask
    "HEM": "_HEM.tif",  # Height Error Mask
    "WBM": "_WBM.tif",  # Water Body Mask
}

DEFAULT_LAYERS: list[str] = ["DEM30"]

OUTPUT_DIR = Path("/home/kaza/eodata/hydroeo/stage1-test/cop-dem")  # TODO

# Your area of interest as a GeoJSON-style dict (WGS84 / EPSG:4326)
AOI_GEOJSON = {
    "type": "Polygon",
    "coordinates": [
        [
            [12.0, 47.0],
            [13.5, 47.0],
            [13.5, 48.0],
            [12.0, 48.0],
            [12.0, 47.0],
        ]
    ],
}

# ── API endpoints ──────────────────────────────────────────────────────────────

TOKEN_URL = "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token"
CATALOGUE_URL = "https://catalogue.dataspace.copernicus.eu/odata/v1/Products"
DOWNLOAD_URL = "https://download.dataspace.copernicus.eu/odata/v1/Products"


# ── Step 1: Authenticate ───────────────────────────────────────────────────────


def get_access_token(username: str, password: str) -> str:
    """Obtain a short-lived access token from CDSE Keycloak."""
    response = requests.post(
        TOKEN_URL,
        data={
            "client_id": "cdse-public",
            "grant_type": "password",
            "username": username,
            "password": password,
        },
    )
    response.raise_for_status()
    token = response.json()["access_token"]
    logger.info("✓ Authentication successful")
    return token


# ── Step 2: Search catalogue ───────────────────────────────────────────────────


def search_cop_dem(aoi_wkt: str, dataset: str) -> list[dict]:
    """
    Query the CDSE OData catalogue for COP-DEM tiles intersecting the AOI.

    The collection name is 'COP-DEM'; the specific resolution/format is
    selected via the 'dataset' attribute filter.
    """
    # OData filter: collection + dataset attribute + spatial intersection
    odata_filter = (
        f"Collection/Name eq 'CCM' "
        f"and Attributes/OData.CSC.StringAttribute/any("
        f"att:att/Name eq 'dataset' and att/OData.CSC.StringAttribute/Value eq '{dataset}') "
        f"and OData.CSC.Intersects(area=geography'SRID=4326;{aoi_wkt}')"
    )

    params = {
        "$filter": odata_filter,
        "$expand": "Attributes",
        "$orderby": "Name asc",
        "$top": 1000,
    }

    response = requests.get(CATALOGUE_URL, params=params)
    response.raise_for_status()
    products = response.json().get("value", [])
    logger.info(f"✓ Found {len(products)} COP-DEM tile(s) intersecting the AOI")
    return products


# ── Step 3: Download tiles ─────────────────────────────────────────────────────


def download_tile(product: dict, token: str, output_dir: Path) -> Path:
    """Download a single COP-DEM product (tile) as a ZIP archive."""
    product_id = product["Id"]
    product_name = product["Name"]
    out_path = output_dir / f"{product_name}.zip"

    if out_path.exists():
        return out_path

    url = f"{DOWNLOAD_URL}({product_id})/$value"
    headers = {"Authorization": f"Bearer {token}"}

    with requests.get(url, headers=headers, stream=True, allow_redirects=True) as r:
        r.raise_for_status()
        downloaded = 0
        with open(out_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):  # 1 MB chunks
                f.write(chunk)
                downloaded += len(chunk)
    return out_path


# ── Main ───────────────────────────────────────────────────────────────────────


def download_cop_dem_for_polygon(
    aoi_geojson: dict,
    dataset: str = _CDSE_GLO30,
    username: str = CDSE_USERNAME,
    password: str = CDSE_PASSWORD,
    output_dir: Path = OUTPUT_DIR,
) -> list[Path]:
    """
    End-to-end: authenticate → search → download all COP-DEM tiles
    that intersect the given GeoJSON polygon.

    Returns a list of downloaded ZIP file paths.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Convert GeoJSON polygon to WKT (required by OData spatial filter)
    geom = shape(aoi_geojson)
    aoi_wkt = wkt_dumps(geom, rounding_precision=6)

    token = get_access_token(username, password)
    products = search_cop_dem(aoi_wkt, dataset)

    if not products:
        logger.info("No tiles found — check your AOI or dataset name.")
        return []

    downloaded_paths = []
    for product in tqdm(products, desc="Downloading COP DEM tiles"):
        path = download_tile(product, token, output_dir)
        downloaded_paths.append(path)

    return downloaded_paths


# ---------------------------------------------------------------------------
#
# MERGING AND CLIPPING
#
# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def load_aoi_geometry(geojson_path: str) -> dict:
    """
    Load an AOI from a GeoJSON file and return a single (merged) Shapely geometry
    plus its GeoJSON-serialisable dict representation.

    Supports FeatureCollection, Feature, or raw geometry objects.
    """
    with open(geojson_path) as f:
        data = json.load(f)

    geom_type = data.get("type")

    if geom_type == "FeatureCollection":
        shapes = [shape(feat["geometry"]) for feat in data["features"]]
    elif geom_type == "Feature":
        shapes = [shape(data["geometry"])]
    else:
        # Assume it is already a raw geometry (Polygon, MultiPolygon, etc.)
        shapes = [shape(data)]

    # Dissolve to a single geometry so the mask works cleanly
    aoi = unary_union(shapes)
    return aoi


def reproject_aoi_to_raster_crs(aoi_shapely, raster_crs: CRS):
    """
    Reproject a Shapely geometry from EPSG:4326 (assumed) into the raster's CRS
    using pyproj / shapely.

    If the raster is already in EPSG:4326 this is a no-op.
    """
    from pyproj import Transformer
    from shapely.ops import transform as shp_transform

    target_crs = CRS(raster_crs)

    # We assume the AOI GeoJSON is in WGS-84 (EPSG:4326) — the GeoJSON spec default.
    src_epsg = "EPSG:4326"
    dst_epsg = target_crs.to_epsg()

    if dst_epsg is None:
        # Fall back to WKT-based authority string if EPSG code is unavailable
        dst_authority = target_crs.to_wkt()
    else:
        dst_authority = f"EPSG:{dst_epsg}"

    if src_epsg == dst_authority:
        return aoi_shapely  # Already in the right CRS

    transformer = Transformer.from_crs(src_epsg, dst_authority, always_xy=True)
    return shp_transform(transformer.transform, aoi_shapely)


def reproject_raster_to_crs(tif_path: str, target_crs: str, tmp_dir: str) -> str:
    """
    Reproject a GeoTIFF to a common target CRS (e.g. "EPSG:4326") using rioxarray.

    Both the raster and the AOI are brought into the same CRS before clipping,
    so the AOI geometry never needs to be reprojected — it stays in target_crs
    throughout the whole pipeline.

    Parameters
    ----------
    tif_path   : path to the source GeoTIFF
    target_crs : EPSG string for the common CRS, e.g. "EPSG:4326"
    tmp_dir    : directory to write the reprojected temporary file

    Returns
    -------
    Path to the reprojected GeoTIFF written inside tmp_dir.
    If the raster is already in target_crs the original path is returned as-is
    (no unnecessary copy is made).
    """
    out_path = Path(tmp_dir) / f"reprojected_{Path(tif_path).name}"

    da = xr.open_dataset(tif_path, engine="rasterio")

    src_crs = CRS.from_user_input(da.rio.crs)
    dst_crs = CRS.from_string(target_crs)

    if src_crs == dst_crs:
        da.close()
        return tif_path

    da_reprojected = da.rio.reproject(target_crs)
    da_reprojected.rio.to_raster(str(out_path), driver="GTiff", compress="lzw")
    da.close()

    return str(out_path)


def clip_tif_to_aoi(
    tif_path: str,
    aoi_shapely,
    tmp_dir: str,
    target_crs: str = "EPSG:4326",
) -> str:
    """
    Clip a single GeoTIFF to the AOI geometry.

    Parameters
    ----------
    tif_path   : path to the source GeoTIFF
    aoi_shapely: Shapely geometry in EPSG:4326
    tmp_dir    : directory to write the clipped temporary file

    Returns
    -------
    Path to the clipped GeoTIFF written inside tmp_dir.
    """
    out_path = Path(tmp_dir) / f"clipped_{Path(tif_path).name}"
    if out_path.exists():
        return str(out_path)

    reprojected_path = reproject_raster_to_crs(tif_path, target_crs, tmp_dir)

    with rasterio.open(reprojected_path) as src:
        # 2. Bounding-box check — AOI and raster are now in the same CRS
        raster_bounds = box(*src.bounds)
        if not aoi_shapely.intersects(raster_bounds):
            return None

        # 3. Clip to the AOI
        clipped, transform = rasterio_mask(
            src,
            [mapping(aoi_shapely)],
            crop=True,  # Crop to AOI bounding box
            all_touched=True,  # Include pixels touching the boundary
            filled=True,  # Fill outside pixels with nodata
        )

        profile = src.profile.copy()
        profile.update(
            height=clipped.shape[1],
            width=clipped.shape[2],
            transform=transform,
            driver="GTiff",
            compress="lzw",
        )

        with rasterio.open(out_path, "w", **profile) as dst:
            dst.write(clipped)

    return str(out_path)


def merge_clipped_tifs(clipped_paths: list[str], output_path: str):
    """
    Merge a list of clipped GeoTIFFs into a single output raster using
    rasterio.merge (first-on-top strategy by default).

    The result is also validated by opening it with xarray / rioxarray.
    """
    datasets = [rasterio.open(p) for p in clipped_paths]

    try:
        mosaic, out_transform = rasterio_merge(datasets)
    finally:
        for ds in datasets:
            ds.close()

    # Re-open one source to copy the profile
    with rasterio.open(clipped_paths[0]) as ref:
        profile = ref.profile.copy()

    profile.update(
        height=mosaic.shape[1],
        width=mosaic.shape[2],
        transform=out_transform,
        driver="GTiff",
        compress="lzw",
    )

    with rasterio.open(output_path, "w", **profile) as dst:
        dst.write(mosaic)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def process(
    aoi_geojson: dict,
    tif_paths: list[str],
    output_path: Path,
    target_crs: str = "EPSG:4326",
):
    """
    Full pipeline:
      1. Load AOI geometry from GeoJSON
      2. Clip every GeoTIFF to the AOI (skipping non-overlapping tiles)
      3. Merge clipped tiles into a single output GeoTIFF
    """
    aoi = shape(aoi_geojson)

    with tempfile.TemporaryDirectory() as tmp_dir:
        clipped = []
        for tif in tqdm(tif_paths, desc="Clipping COM DEM tiles"):
            result = clip_tif_to_aoi(tif, aoi, tmp_dir, target_crs)
            if result:
                clipped.append(result)

        if not clipped:
            raise RuntimeError(
                "No input rasters overlapped the AOI. "
                "Check that your GeoTIFs and AOI share geographic coverage."
            )

        merge_clipped_tifs(clipped, output_path)


# ---------------------------------------------------------------------------
# Layer validation helpers
# ---------------------------------------------------------------------------


def _validate_layers(layers: list[str]) -> None:
    """
    Validate that all requested layer names are recognised and that DEM30 and
    DEM90 are not requested simultaneously (they target different CDSE datasets).

    Parameters
    ----------
    layers : list[str]
        Layer names, already normalised to upper-case.

    Raises
    ------
    ValueError
        On unknown layer names or a DEM30 + DEM90 conflict.
    """
    unknown = set(layers) - VALID_LAYERS
    if unknown:
        raise ValueError(
            f"Unknown layer(s): {sorted(unknown)}. "
            f"Valid options are: {sorted(VALID_LAYERS)}"
        )
    if "DEM30" in layers and "DEM90" in layers:
        raise ValueError(
            "Cannot request DEM30 and DEM90 together — they target different "
            "CDSE datasets (GLO-30 vs GLO-90). Request one DEM resolution at a time."
        )


def _resolve_cdse_dataset(layers: list[str]) -> str:
    """
    Return the CDSE dataset identifier string for the requested layers.

    If DEM90 is present, the GLO-90 dataset is used (aux layers, if also
    requested, are extracted from the same GLO-90 ZIP archives).
    Otherwise, the GLO-30 dataset is used.

    Parameters
    ----------
    layers : list[str]
        Validated, upper-case layer names.

    Returns
    -------
    str
        CDSE OData dataset identifier.
    """
    return _CDSE_GLO90 if "DEM90" in layers else _CDSE_GLO30


def download_cop_dem(
    minx: float,
    miny: float,
    maxx: float,
    maxy: float,
    output_dir: str,
    username: str = CDSE_USERNAME,
    password: str = CDSE_PASSWORD,
    layers: list[str] | None = None,
    output_basename: str = "cop_dem_merged",
) -> dict[str, Path]:
    """
    Download, merge, and clip one or more COP-DEM product layers for a bounding
    box (bbox).

    The same set of ZIP tiles is downloaded once from CDSE and all requested
    layers are extracted in a single pass, then each layer is independently
    clipped to the AOI and merged into a separate output GeoTIFF.

    Parameters
    ----------
    minx : float
        Minimum longitude (western edge).
    miny : float
        Minimum latitude (southern edge).
    maxx : float
        Maximum longitude (eastern edge).
    maxy : float
        Maximum latitude (northern edge).
    output_dir : str
        Directory where the merged output GeoTIFFs will be written.
    username : str
        CDSE username (or set CDSE_USERNAME env var).
    password : str
        CDSE password (or set CDSE_PASSWORD env var).
    layers : list[str] or None
        Product layers to extract.  Each entry must be one of:
        ``DEM30`` (30 m elevation), ``DEM90`` (90 m elevation),
        ``EDM`` (Editing Mask), ``FLM`` (Filling Mask),
        ``HEM`` (Height Error Mask), ``WBM`` (Water Body Mask).
        ``DEM30`` and ``DEM90`` cannot be combined (they target different
        CDSE datasets).  When only aux layers are requested (no DEM30/DEM90),
        the GLO-30 dataset is used by default.
        Defaults to ``["DEM30"]``.
    output_basename : str
        Base name (without extension) for the output files.  Each layer
        produces ``{output_basename}_{LAYER}.tif``.
        Defaults to ``"cop_dem_merged"``.

    Returns
    -------
    dict[str, Path]
        Mapping of layer name → path to the merged and clipped GeoTIFF.
        E.g. ``{"DEM30": Path(".../cop_dem_merged_DEM30.tif"), "WBM": ...}``.

    Raises
    ------
    ValueError
        On invalid layer names or a DEM30 + DEM90 conflict.
    RuntimeError
        If no tiles are found or if no rasters overlap the AOI for a layer.
    """
    if layers is None:
        layers = list(DEFAULT_LAYERS)
    layers = [lyr.upper() for lyr in layers]
    _validate_layers(layers)

    cdse_dataset = _resolve_cdse_dataset(layers)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Strip extension from basename so each layer gets its own suffix
    stem = Path(output_basename).stem
    output_paths: dict[str, Path] = {
        layer: output_dir / f"{stem}_{layer}.tif" for layer in layers
    }

    # Idempotency: return immediately if every output already exists
    if all(p.exists() for p in output_paths.values()):
        logger.info("All requested COP-DEM outputs already exist.")
        return output_paths

    # Create GeoJSON polygon from bbox
    aoi_geojson = {
        "type": "Polygon",
        "coordinates": [
            [
                [minx, miny],
                [maxx, miny],
                [maxx, maxy],
                [minx, maxy],
                [minx, miny],
            ]
        ],
    }

    # Step 1: Download ZIP tiles (single download for all layers)
    logger.info(
        f"Downloading COP-DEM tiles for bbox: ({minx}, {miny}, {maxx}, {maxy}) "
        f"[dataset: {cdse_dataset}, layers: {layers}]"
    )
    downloaded_zips = download_cop_dem_for_polygon(
        aoi_geojson=aoi_geojson,
        dataset=cdse_dataset,
        username=username,
        password=password,
        output_dir=output_dir / "zips",
    )

    if not downloaded_zips:
        raise RuntimeError(
            f"No COP-DEM tiles found for bbox ({minx}, {miny}, {maxx}, {maxy}). "
            "Check that your bbox is valid and within coverage."
        )

    # Step 2: Unzip tiles — extract each requested layer to its own subdirectory
    tif_dir = output_dir / "tifs"
    tif_dir.mkdir(parents=True, exist_ok=True)

    layer_tif_dirs: dict[str, Path] = {}
    for layer in layers:
        layer_dir = tif_dir / layer
        layer_dir.mkdir(parents=True, exist_ok=True)
        layer_tif_dirs[layer] = layer_dir

    layer_tifs: dict[str, list[Path]] = {layer: [] for layer in layers}

    for in_file in tqdm(downloaded_zips, desc="Unzipping COP-DEM tiles"):
        with zipfile.ZipFile(in_file, "r") as archive:
            for file in archive.namelist():
                for layer in layers:
                    suffix = LAYER_SUFFIXES[layer]
                    if not file.endswith(suffix):
                        continue

                    dest_dir = layer_tif_dirs[layer]
                    filename = Path(file).name
                    out_path = dest_dir / filename

                    if out_path.exists():
                        layer_tifs[layer].append(out_path)
                    else:
                        archive.extract(file, dest_dir)
                        extracted_path = dest_dir / file
                        if extracted_path.exists() and extracted_path != out_path:
                            shutil.move(str(extracted_path), str(out_path))

                        # Clean up any empty subdirectories left by extraction
                        for subdir in dest_dir.glob("*/"):
                            if subdir.is_dir() and not list(subdir.iterdir()):
                                shutil.rmtree(subdir, ignore_errors=True)

                        layer_tifs[layer].append(out_path)

    # Validate that at least one TIF was found per layer
    for layer, tifs in layer_tifs.items():
        if not tifs:
            raise RuntimeError(
                f"No {LAYER_SUFFIXES[layer]} files found for layer '{layer}' "
                "in the downloaded ZIP archives. "
                "The archive format may have changed or the layer is unavailable "
                f"for the selected dataset ({cdse_dataset})."
            )

    # Remove ZIPs — all layers have been extracted
    for f in downloaded_zips:
        f.unlink()
    shutil.rmtree(output_dir / "zips", ignore_errors=True)

    # Step 3: Clip and merge each layer independently
    for layer, tif_paths in layer_tifs.items():
        out_path = output_paths[layer]
        if out_path.exists():
            logger.info(f"Layer '{layer}' output already exists, skipping: {out_path}")
            continue
        logger.info(f"Processing layer '{layer}' → {out_path}")
        process(
            aoi_geojson=aoi_geojson,
            tif_paths=[str(p) for p in tif_paths],
            output_path=out_path,
            target_crs="EPSG:4326",
        )

    # Clean up extracted TIFs
    for layer in layers:
        for f in layer_tifs[layer]:
            f.unlink(missing_ok=True)
    shutil.rmtree(tif_dir, ignore_errors=True)

    logger.info(
        "COP-DEM processing complete. Outputs: "
        + ", ".join(f"{lyr}: {p}" for lyr, p in output_paths.items())
    )
    return output_paths
