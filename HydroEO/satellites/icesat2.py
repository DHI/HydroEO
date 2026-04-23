"""ICESat-2 ATL13 download via SlideRule's ``atl13x`` endpoint.

Architecture
------------
SlideRule's ``atl13x`` endpoint is **water-body-centric**, not polygon-centric.
It uses an internal ATL13 Metadata Service (AMS) to identify which granules
contain a given water body.  The ``atl13.coord`` sub-key is **required** in
every ``atl13x`` parameter dict — without it the server returns a 400 error.
Passing only a ``poly`` and time range (as the former Harmony approach did) does
not work.

For each reservoir in HydroEO's workflow:

1. The **centroid** of the reservoir polygon is used as the AMS lookup
   coordinate (``atl13.coord``).  SlideRule resolves this coordinate to the
   registered water body; it is *not* a spatial filter.
2. The polygon itself is converted to ``poly`` format via
   ``sliderule.toregion()`` and used for **spatial subsetting** of the returned
   segments (useful for large lakes).
3. SlideRule streams back a ``GeoDataFrame`` directly — no local HDF5 files are
   written.  If ``download_directory`` is provided the GeoDataFrame is cached
   as ``atl13.parquet`` for later use by ``extract_observations()``.

Column notes
------------
* ``ht_ortho`` is already **orthometric (EGM2008 geoid-corrected)** height,
  applied by SlideRule server-side.  Do not apply a second geoid correction.
* ``orbit_number`` is not available from ``atl13x``; use ``rgt + cycle_number``
  as the equivalent orbit identifier.
* ``file_name`` is not returned (no local files); the ``srcid`` -> granule name
  mapping is available in ``gdf.attrs["meta"]["srctbl"]`` if needed.
* Latitude and longitude are encoded in the ``geometry`` column.  Explicit
  ``lat`` and ``lon`` columns are extracted from geometry for downstream compat.
"""

from __future__ import annotations

import datetime
import os

import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon

import sliderule
from HydroEO.utils import general

# Initialise SlideRule client once at module load time.
# sliderule.init() only sets the server URL — it does not open a network
# connection — so it is safe to call at import time.
sliderule.init("slideruleearth.io")

# ── Column mapping: SlideRule atl13x default output -> HydroEO schema ────────
# Default columns always present in the atl13x GeoDataFrame result:
#   cycle, geometry, gt, ht_ortho, ht_water_surf, rgt, segment_id_beg,
#   spot, srcid, stdev_water_surf, water_depth
#
# New columns returned by SlideRule not present in the old Harmony schema
# (kept, they are useful): ht_water_surf, stdev_water_surf, water_depth,
# spot, srcid, segment_id_beg.
SR_DEFAULT_COLUMN_MAP: dict[str, str] = {
    "ht_ortho": "height",  # orthometric water surface height (EGM2008, server-side)
    "cycle": "cycle_number",
    "gt": "beam",
    # rgt is kept as-is (no rename needed)
}

# ── Ancillary beam-group fields fetchable via parms["atl13_fields"] ───────────
# These additional HDF5 beam-group fields are fetched server-side by SlideRule
# when specified in parms["atl13_fields"].  They become extra columns in the
# returned GeoDataFrame.
SR_ANCILLARY_COLUMN_MAP: dict[str, str] = {
    "inland_water_body_id": "wb_id",
    "inland_water_body_size": "wb_size",
    "inland_water_body_type": "wb_type",
    "segment_slope_trk_bdy": "wb_slope",
    "err_ht_water_surf": "height_err",
    "segment_quality": "quality_seg",
    "segment_geoid": "geoid_track",
    "segment_geoid_free2mean": "geoid_corr_track",
    "segment_dem_ht": "dem",
    "segment_near_sat_fract": "sat_frac_track",
}

# Valid ancillary field names for config validation (keys of SR_ANCILLARY_COLUMN_MAP).
SR_ATL13_VALID_ANCILLARY_FIELDS: set[str] = set(SR_ANCILLARY_COLUMN_MAP)

# Kept for backward API compatibility — no ancillary fields are fetched by default.
ATL13_DEFAULT_FIELDS: list[str] = []

ATL13_CORE_FIELDS = ["height", "lat", "lon", "date"]


def query(
    aoi: list,
    startdate: datetime.date,
    enddate: datetime.date,
    download_directory: str | None,
    atl13_options: dict | None = None,
    atl13_fields: list[str] | None = None,
) -> gpd.GeoDataFrame:
    """Download ATL13 water surface elevations via SlideRule's ``atl13x`` endpoint.

    SlideRule's ``atl13x`` endpoint is water-body-centric.  The **centroid** of
    the reservoir polygon is used to identify the water body in the ATL13
    Metadata Service (AMS) via ``atl13.coord`` — no polygon spatial filter is
    sent to the server (doing so causes empty results in SlideRule v5).
    Precise spatial filtering is applied client-side in ``extract_observations()``.

    Parameters
    ----------
    aoi:
        List of ``(lon, lat)`` coordinate tuples defining the reservoir polygon.
    startdate, enddate:
        Download date range.
    download_directory:
        If provided, the returned GeoDataFrame is cached as ``atl13.parquet``
        inside this directory so that ``extract_observations()`` can read it
        later.  If ``None``, no file is written.
    atl13_options:
        Optional dict of ``atl13`` sub-key overrides, e.g.
        ``{"pass_invalid": True, "beams": ["gt1l", "gt1r"]}``.
        The ``coord`` key must **not** be set here — it is always derived from
        the polygon centroid at runtime.
    atl13_fields:
        Optional list of ancillary ATL13 beam-group HDF5 field names to fetch
        server-side (e.g. ``["segment_quality", "segment_geoid"]``).  Valid
        values are the keys of ``SR_ANCILLARY_COLUMN_MAP``.

    Returns
    -------
    gpd.GeoDataFrame
        Observations with HydroEO-schema column names and an explicit ``date``
        column derived from the SlideRule datetime index.

    Raises
    ------
    HydroEODownloadError
        On SlideRule connection / HTTP errors, or when the centroid coordinate
        does not fall inside any water body registered in the AMS database
        (GRWL / HydroLAKES).
    """
    from HydroEO.system import HydroEODownloadError

    poly = Polygon(aoi)
    centroid = poly.representative_point()

    atl13_sub: dict = {"coord": {"lon": centroid.x, "lat": centroid.y}}
    # atl13_sub: dict = {"coord": {"lon": 100.32201928690762, "lat": 22.858474628948642}}   # TODO
    atl13_sub.update(atl13_options or {})

    parms: dict = {
        "atl13": atl13_sub,  # required — water-body AMS lookup
        "t0": startdate.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "t1": enddate.strftime("%Y-%m-%dT%H:%M:%SZ"),
    }

    if atl13_fields:
        parms["atl13_fields"] = atl13_fields

    try:
        print("This will take time..")
        gdf = sliderule.run("atl13x", parms)
    except FileNotFoundError as exc:
        raise HydroEODownloadError(
            f"SlideRule lookup failed for centroid "
            f"(lon={centroid.x:.4f}, lat={centroid.y:.4f})."
        ) from exc
    except Exception as exc:
        raise HydroEODownloadError(f"SlideRule atl13x request failed: {exc}") from exc

    if gdf is None or gdf.empty:
        raise HydroEODownloadError(
            f"SlideRule atl13x returned no data for centroid "
            f"(lon={centroid.x:.4f}, lat={centroid.y:.4f}). "
            "This usually means the reservoir is not registered in the AMS "
            "database (GRWL / HydroLAKES). Verify that the water body appears "
            "in HydroLAKES and that the polygon centroid falls inside it."
        )

    # ── Apply HydroEO column schema ───────────────────────────────────────────
    rename_map = {k: v for k, v in SR_DEFAULT_COLUMN_MAP.items() if k in gdf.columns}
    if atl13_fields:
        rename_map.update(
            {k: v for k, v in SR_ANCILLARY_COLUMN_MAP.items() if k in gdf.columns}
        )
    gdf = gdf.rename(columns=rename_map)

    # Extract explicit lat/lon columns from geometry for downstream compatibility.
    gdf["lat"] = gdf.geometry.y
    gdf["lon"] = gdf.geometry.x

    # ht_ortho is already geoid-corrected (EGM2008) by SlideRule server-side.
    # Do not apply a second geoid correction here.

    # orbit_number is not available from atl13x.
    # rgt + cycle_number is the equivalent orbit identifier.
    # file_name is not returned (no local files); use
    # gdf.attrs["meta"]["srctbl"][srcid] to map srcid -> granule name if needed.

    gdf["date"] = gdf.index
    gdf = gdf.reset_index(drop=True)

    # ── Date filter ───────────────────────────────────────────────────────────
    # SlideRule t0/t1 may not be enforced server-side; filter client-side so
    # only observations within [startdate, enddate] are cached and returned.
    start_dt = pd.Timestamp(startdate).tz_localize("UTC")
    end_dt = pd.Timestamp(enddate).tz_localize("UTC")
    date_col = pd.to_datetime(gdf["date"], utc=True)
    gdf = gdf.loc[(date_col >= start_dt) & (date_col <= end_dt)].reset_index(drop=True)

    # ── Optional Parquet cache ────────────────────────────────────────────────
    if download_directory is not None:
        general.ifnotmakedirs(download_directory)
        gdf.to_parquet(os.path.join(download_directory, "atl13.parquet"))

    return gdf


def extract_observations(
    src_dir,
    dst_path,
    features,
    atl13_fields=None,  # kept for API compat; no-op in SlideRule path
    track_keys=None,  # kept for API compat; no-op in SlideRule path
):
    """Extract ICESat-2 ATL13 observations from a cached Parquet file.

    Reads the ``atl13.parquet`` file written by :func:`query` from *src_dir*,
    filters observations to those that fall within the *features* geometry, and
    saves the result as a shapefile at *dst_path*.

    ``atl13_fields`` and ``track_keys`` are accepted for API compatibility but
    are no-ops in the SlideRule path: field selection is handled at download
    time via the ``atl13_fields`` config key, and beam filtering via
    ``atl13.beams``.
    """
    parquet_path = os.path.join(src_dir, "atl13.parquet")
    if not os.path.exists(parquet_path):
        return

    gdf = gpd.read_parquet(parquet_path)
    if gdf.empty:
        return

    # Filter observations to those inside the reservoir geometry.
    gdf = gdf.loc[gdf.within(features.unary_union)].reset_index(drop=True)
    if len(gdf) == 0:
        return

    gdf["platform"] = "icesat2"
    gdf["product"] = "ATL13"

    # Ensure CRS matches features; convert if necessary.
    if gdf.crs is None:
        gdf = gdf.set_crs(features.crs)
    elif gdf.crs != features.crs:
        gdf = gdf.to_crs(features.crs)

    # Strip timezone info so shapefiles can represent the date column.
    if "date" in gdf.columns:
        gdf["date"] = pd.to_datetime(gdf["date"]).dt.tz_localize(None)

    gdf = _to_shapefile_safe_columns(gdf)
    gdf.to_file(dst_path)


def _to_shapefile_safe_columns(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Truncate column names to 10 chars (ESRI Shapefile field name limit)."""
    rename_map: dict[str, str] = {}
    used: set[str] = set()

    for col in gdf.columns:
        if col == "geometry":
            continue

        base = str(col)[:10]
        candidate = base
        suffix = 1
        while candidate in used:
            suffix_str = str(suffix)
            candidate = f"{base[: 10 - len(suffix_str)]}{suffix_str}"
            suffix += 1

        used.add(candidate)
        if candidate != col:
            rename_map[col] = candidate

    return gdf.rename(columns=rename_map)


def get_latest_obs_date(data_dir):
    shp_path = os.path.join(data_dir, "icesat2.shp")
    if os.path.exists(shp_path):
        gdf = gpd.read_file(shp_path)
        last_obs_date = max(gdf.date.values).astype(datetime.date)
        return last_obs_date
