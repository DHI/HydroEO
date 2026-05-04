"""
Download COP-DEM GLO-30 (or GLO-90) tiles from CDSE for a given polygon.

Prerequisites:
    pip install requests shapely

You need a CDSE account with CCM data access:
    https://dataspace.copernicus.eu/ → register, then request CCM access.
"""

import os
import requests
import logging
from pathlib import Path
from shapely.geometry import shape
from shapely.wkt import dumps as wkt_dumps

logger = logging.getLogger(__name__)

# ── Configuration ─────────────────────────────────────────────────────────────

CDSE_USERNAME = os.environ.get("CDSE_USERNAME", "your_username")
CDSE_PASSWORD = os.environ.get("CDSE_PASSWORD", "your_password")

# Dataset options:
#   "COP-DEM_GLO-30-DGED/2023_1"   ← 30 m, full-resolution GeoTIFF
#   "COP-DEM_GLO-90-DGED/2023_1"   ← 90 m, GeoTIFF
DATASET = "COP-DEM_GLO-30-DGED/2023_1"

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
    print("✓ Authentication successful")
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
    print(f"✓ Found {len(products)} COP-DEM tile(s) intersecting the AOI")
    return products


# ── Step 3: Download tiles ─────────────────────────────────────────────────────


def download_tile(product: dict, token: str, output_dir: Path) -> Path:
    """Download a single COP-DEM product (tile) as a ZIP archive."""
    product_id = product["Id"]
    product_name = product["Name"]
    out_path = output_dir / f"{product_name}.zip"

    if out_path.exists():
        print(f"  — Skipping (already downloaded): {product_name}")
        return out_path

    url = f"{DOWNLOAD_URL}({product_id})/$value"
    headers = {"Authorization": f"Bearer {token}"}

    with requests.get(url, headers=headers, stream=True, allow_redirects=True) as r:
        r.raise_for_status()
        total = int(r.headers.get("Content-Length", 0))
        downloaded = 0
        with open(out_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):  # 1 MB chunks
                f.write(chunk)
                downloaded += len(chunk)
                if total:
                    pct = 100 * downloaded / total
                    print(f"\r  Downloading {product_name}: {pct:.1f}%", end="")
    print(f"\n  ✓ Saved: {out_path}")
    return out_path


# ── Main ───────────────────────────────────────────────────────────────────────


def download_cop_dem_for_polygon(
    aoi_geojson: dict,
    dataset: str = DATASET,
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
        print("No tiles found — check your AOI or dataset name.")
        return []

    print("\nTiles to download:")
    for p in products:
        size_mb = p.get("ContentLength", 0) / 1e6
        print(f"  • {p['Name']}  ({size_mb:.1f} MB)")
        print(f"  • {p.get('S3Path')}")

    downloaded_paths = []
    for product in products:
        path = download_tile(product, token, output_dir)
        downloaded_paths.append(path)

    print(f"\n✓ Done. {len(downloaded_paths)} tile(s) saved to '{output_dir}'")
    return downloaded_paths


if __name__ == "__main__":
    download_cop_dem_for_polygon(
        aoi_geojson=AOI_GEOJSON,
        dataset=DATASET,
        output_dir=OUTPUT_DIR,
    )
