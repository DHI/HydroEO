# SWOT Raster — Configuration Reference

Download and process SWOT L2 HR/LR raster products for an arbitrary area of interest. Tiles are extracted, quality-filtered, clipped to the AOI, and merged by date/variable. Requires NASA Earthdata credentials.

## Products

| Product | Description |
| --- | --- |
| `SWOT_L2_HR_Raster_D` | High-resolution raster (default) |
| `SWOT_L2_LR_SSH_2.0` | Low-resolution sea surface height |
| `SWOT_L2_HR_RIVERSP_2.0` | High-resolution river surface product |

Credentials: Earthdata account via `EARTHDATA_USERNAME` / `EARTHDATA_PASSWORD` (or config `earthaccess` section).

## Config reference

Start from [`configs/swot_raster.yaml`](swot_raster.yaml).

### `swot_raster.aoi`

| Key | Type | Required | Description |
| --- | --- | --- | --- |
| `name` | string | ✅ | Identifier used in output directory names |
| `type` | string | ✅ | `bbox` \| `shapefile` \| `geopackage` |
| `bbox` | list | ✅ for `bbox` | `[lon_min, lat_min, lon_max, lat_max]` |
| `path` | string | ✅ for `shapefile`/`geopackage` | Path to geometry file |

### `swot_raster`

| Key | Default | Description |
| --- | --- | --- |
| `product` | `SWOT_L2_HR_Raster_D` | Product short name (see table above) |
| `startdate` / `enddate` | project dates | Per-section date override |
| `granule_filter` | — | Optional glob pattern to filter granules (e.g. `"*100m*"` for 100 m tiles) |
| `merge_tiles` | `true` | Merge tiles by date/variable and reproject; `false` = keep individual clipped tiles |
| `target_crs` | `gis.global_crs` | Output CRS for the merge phase (e.g. `"EPSG:32645"` for UTM 45N) |
| `variables` | all 8 | List of variables to extract; omit to extract all available |
| `quality_filters.max_wse_uncert` | `0.3` | Mask pixels with `wse_uncert ≥` this value (metres) |
| `quality_filters.max_layover_impact` | `0.3` | Mask pixels with `layover_impact ≥` this value (metres) |

**Available variables:** `wse`, `wse_uncert`, `wse_qual`, `height_cor_xover`, `geoid`, `n_wse_pix`, `n_other_pix`, `layover_impact`.

> `wse`, `wse_uncert`, and `layover_impact` are always extracted regardless of the `variables` list — they are required for quality masking.

## Processing pipeline

1. **Download** — queries Earthdata for granules within bbox and date range; tracks already-downloaded granules in `raw/<product>/downloaded.log` to avoid re-fetching.
2. **Preprocess** — extracts configured variables, applies quality filters, clips each tile to the AOI; outputs GeoTIFFs. Out-of-area tiles are discarded. Raw netCDF files are deleted after preprocessing.
3. **Merge** — groups clipped tiles by date and variable, merges multiple passes, reprojects to target CRS; outputs mosaics.

## Output structure

```
{main_dir}/
  raw/
    swot_raster/<aoi_name>/<product>/   # raw netCDF granule files (deleted after preprocessing)
  processed/
    swot_raster/<aoi_name>/<product>/   # extracted/clipped GeoTIFFs per granule per variable
  results/
    <aoi_name>/                         # merged mosaics by date/variable
```

## Direct download (no config file)

```sh
hydroeo fetch swot-raster \
  --bbox "-10 40 10 60" \
  --start 2023-01-01 --end 2023-06-01 \
  --output ./output \
  --username <earthaccess-user> \
  --password <earthaccess-pass>

# Or set credentials as environment variables:
export EARTHDATA_USERNAME=...
export EARTHDATA_PASSWORD=...
hydroeo fetch swot-raster --bbox "-10 40 10 60" --start 2023-01-01 --end 2023-06-01 --output ./output
```

Optional flags: `--aoi-name <label>` (default `aoi`), `--product <short-name>` (default `SWOT_L2_HR_Raster_D`).
