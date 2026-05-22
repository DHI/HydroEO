# SWOT Pixel Cloud — Configuration Reference

Download SWOT L2 PIXC data and grid point clouds to regular rasters. Points are filtered by water classification, clipped to the AOI, and aggregated via binned statistics (median, mean, etc.). Requires NASA Earthdata credentials.

## Products

| Product | Description |
| --- | --- |
| `SWOT_L2_HR_PIXC_D` | High-resolution pixel cloud (default) |
| `SWOT_L2_HR_PIXC_2.0` | High-resolution pixel cloud v2.0 format |

Credentials: Earthdata account via `EARTHDATA_USERNAME` / `EARTHDATA_PASSWORD` (or config `earthaccess` section).

## Config reference

Start from [`configs/swot_pixc.yaml`](swot_pixc.yaml).

### `swot_pixc.aoi`

| Key | Type | Required | Description |
| --- | --- | --- | --- |
| `name` | string | ✅ | Identifier used in output directory names |
| `type` | string | ✅ | `bbox` \| `shapefile` \| `geopackage` |
| `bbox` | list | ✅ for `bbox` | `[lon_min, lat_min, lon_max, lat_max]` |
| `path` | string | ✅ for `shapefile`/`geopackage` | Path to geometry file |

### `swot_pixc`

| Key | Default | Description |
| --- | --- | --- |
| `product` | `SWOT_L2_HR_PIXC_D` | Product short name (see table above) |
| `startdate` / `enddate` | project dates | Per-section date override |
| `classes` | `[open_water, water_near_land]` | Water classification filter (see below) |
| `fields` | `[heightEGM]` | Point fields to extract and grid (see below) |
| `grid_resolution` | `100` | Output pixel size in metres |
| `stat_method` | `median` | Binning statistic: `median` \| `mean` \| `max` \| `min` |
| `target_crs` | `gis.global_crs` | Override output CRS (e.g. `"EPSG:32645"`) |

### Water classes

| Class | Description |
| --- | --- |
| `land` | |
| `land_near_water` | |
| `water_near_land` | ✅ default |
| `open_water` | ✅ default |
| `dark_water` | |
| `low_coh_water_near_land` | |
| `open_low_coh_water` | |

### Available fields

`heightEGM` (geoid-corrected height, default), `height`, `geoid`, `water_frac`, `phase_noise_std`, `dheight_dphase`, `sig0`, and others from the PIXC product.

## Processing pipeline

1. **Download** — queries Earthdata for PIXC granules within bbox and date range; tracks already-downloaded granules in `raw/<product>/downloaded.log`.
2. **Preprocess** — extracts points matching configured water classes, computes derived heights (e.g. `heightEGM` from geoid corrections), clips to AOI bounds, saves to GeoJSON. Raw netCDF files are retained for archival.
3. **Rasterize** — bins point data into regular grids using the specified statistic and resolution; outputs GeoTIFFs by date and variable.

## Output structure

```
{main_dir}/
  raw/
    swot_pixc/<aoi_name>/<product>/   # raw netCDF granule files
  processed/
    swot_pixc/<aoi_name>/             # trimmed GeoJSON point data (filtered by water class)
  results/
    <aoi_name>/                       # gridded rasters by date/field
```

## Direct download (no config file)

```sh
hydroeo fetch swot-pixc \
  --bbox "-10 40 10 60" \
  --start 2023-01-01 --end 2023-06-01 \
  --output ./output \
  --username <earthaccess-user> \
  --password <earthaccess-pass>

# Or set credentials as environment variables:
export EARTHDATA_USERNAME=...
export EARTHDATA_PASSWORD=...
hydroeo fetch swot-pixc --bbox "-10 40 10 60" --start 2023-01-01 --end 2023-06-01 --output ./output
```
