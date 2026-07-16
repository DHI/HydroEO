# Reservoirs & Lakes — Configuration Reference

Multi-satellite water surface elevation timeseries for reservoirs and lakes defined by polygon input. Supports SWOT Lake SP, ICESat-2 ATL13, Sentinel-3, and Sentinel-6.

## Products & credentials

| Product | Satellite | Credentials required |
| --- | --- | --- |
| `SWOT_L2_HR_LakeSP_D` | SWOT | Earthdata account + HydroWeb API key |
| `ATL13` | ICESat-2 | None (SlideRule public API) |
| — | Sentinel-3A/B | CREODIAS account |
| — | Sentinel-6 | CREODIAS account (default, Low Rate product); or Earthdata account for the High Rate product — see `source` under [`sentinel3` / `sentinel6`](#sentinel3--sentinel6) |

Credentials can be set in the config file or as environment variables:

| Service | Environment variables |
| --- | --- |
| Earthdata | `EARTHDATA_USERNAME`, `EARTHDATA_PASSWORD` |
| HydroWeb | `EODAG__HYDROWEB_NEXT__AUTH__CREDENTIALS__APIKEY` |
| CREODIAS | `CREODIAS_USERNAME`, `CREODIAS_PASSWORD` |

> On Windows, credentials may need to be set via `os.environ` before calling `project.download()` if config-file loading fails.

## Config reference

Start from [`configs/reservoirs.yaml`](reservoirs.yaml).

### `reservoirs`

| Key | Type | Required | Description |
| --- | --- | --- | --- |
| `enabled` | bool | ✅ | Must be `true` |
| `path` | string | ✅ | Path to `.shp` or `.gpkg` with one feature per reservoir |
| `id_key` | string | ✅ | Column name used as unique reservoir identifier |
| `export_to_dfs0` | bool | — | Export cleaned observations to dfs0 format (default `false`; requires `mikeio`) |

### `swot`

| Key | Default | Description |
| --- | --- | --- |
| `download` | `false`* | Download SWOT Lake SP granules |
| `process` | `false`* | Extract observations for each reservoir |
| `startdate` / `enddate` | project dates | Per-satellite date override |
| `pld_match_max_distance_m` | `100.0` | Max nearest-neighbour distance to match reservoirs to PLD lakes |
| `exclude_obs_id_values` | `["no_data"]` | SWOT `obs_id` values to drop during extraction |
| `processing_filters` | `["elevation", "MAD"]` | Ordered list of cleaning filters (see [Cleaning filters](#cleaning-filters)) |
| `elevation_min_m` | `0.0` | Lower bound for elevation filter (metres) |
| `elevation_max_m` | `8000.0` | Upper bound for elevation filter (metres) |
| `mad_threshold` | `5.0` | Outlier multiplier for MAD filter |
| `download_dir` | `{main_dir}/raw/swot` | Override download location |

*`false` is the code-level default if omitted entirely — you must explicitly set `download: true`/`process: true` to enable SWOT. The shipped [`configs/reservoirs.yaml`](reservoirs.yaml) template sets both to `true` out of the box.

### `icesat2`

| Key | Default | Description |
| --- | --- | --- |
| `download` / `process` | `false` | Enable/disable download and processing |
| `atl13_fields` | `[]` | Extra SlideRule ancillary fields to request (empty by design -- SlideRule's `atl13x` endpoint already returns the core fields listed below without any extra request) |
| `atl13.pass_invalid` | `false` | Include segments flagged as invalid |
| `atl13.beams` | `[]` | Subset by beam name (e.g. `["gt1l", "gt2r"]`); empty = all |
| `atl13.spots` | `[]` | Subset by spot number (e.g. `[1, 3]`); empty = all |
| `track_keys` | all 6 tracks | Validated against `["gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"]` if set, but currently a **no-op** in the SlideRule extraction path -- kept for API compatibility, doesn't yet filter anything. Use `atl13.beams` above for actual beam filtering. |
| `processing_filters` | `["elevation", "MAD"]` | Same options as SWOT |
| `download_dir` | `{main_dir}/raw/icesat2` | Override download location |

Valid `atl13_fields` values: `inland_water_body_id`, `inland_water_body_size`, `inland_water_body_type`, `segment_slope_trk_bdy`, `err_ht_water_surf`, `segment_quality`, `segment_geoid`, `segment_geoid_free2mean`, `segment_dem_ht`, `segment_near_sat_fract`.

Common output columns: `height` (EGM2008-corrected), `cycle_number`, `beam`, `rgt`, `ht_water_surf`, `stdev_water_surf`, `water_depth`, `spot`, `srcid`.

### `sentinel3` / `sentinel6`

| Key | Default | Description |
| --- | --- | --- |
| `download` / `process` | `false` | Enable/disable |
| `sigma0_max` | `100000.0` | Maximum accepted sigma0 during extraction |
| `download_threads` | `1` | Parallel CREODIAS download threads |
| `subset_file_id` | `"enhanced_measurement.nc"` | Filename identifying the file to extract from each downloaded product zip |
| `source` (sentinel6 only) | `"creodias"` | CREODIAS only distributes the Sentinel-6 Low Rate product. Set to `"earthdata"` for the High Rate (20Hz Ku-band) product via PO.DAAC/EarthData instead — needs `EARTHDATA_USERNAME`/`PASSWORD` rather than CREODIAS credentials. |
| `processing_filters` | `["elevation", "MAD"]` | Same filter options as SWOT |
| `download_dir` | `{main_dir}/raw/sentinel3` or `sentinel6` | Override download location |

### `hydroweb` (Prior Lake Database)

The PLD matches input reservoirs to SWOT lake identifiers. Required for SWOT extraction.

| Key | Default | Description |
| --- | --- | --- |
| `api_key` | — | HydroWeb API key (or env var) |
| `raw_pld_path` | — | Optional path to existing PLD zip or folder; skips download |
| `keep_raw_pld` | `false` | Retain raw PLD zip and temp folder after subset creation |
| `continent_codes` | — | Optional list of continent-tile suffixes (e.g. `["AS", "AU"]`) identifying which full-schema PLD files to use for backfilling `res_id`. Inspect your downloaded PLD folder (`aux/PLD/PLD_temp/...`) to see which ones are actually present. If omitted, all full-schema tiles found are used. |

PLD downloads include a global `_light` file (lake matching only, no `res_id`) alongside full-schema per-continent tiles (which additionally carry `res_id` — a cross-reference to external reservoir databases like GeoDAR, present only for PLD lakes classified as man-made reservoirs). HydroEO uses the `_light` file for coverage (guaranteed complete, avoids silently missing a continent your PLD download didn't include) and backfills `res_id` from whichever tile files are present, matched by `lake_id`. 
If no `_light` file is present at all, HydroEO falls back to using the tile files directly for coverage too, with a warning that this risks missing continents the download didn't include.

The PLD subset is written to `{main_dir}/aux/PLD/PLD_subset.gpkg` with QA/QC files `present_in_pld.gpkg` (matched) and `missing_in_pld.gpkg` (unmatched).

## Cleaning filters

Applied during `create_timeseries()`. Configured per mission under `processing_filters`:

| Filter | What it removes |
| --- | --- |
| `elevation` | Observations outside `[elevation_min_m, elevation_max_m]` |
| `MAD` | Outliers beyond `mad_threshold × MAD` from the median |
| `daily_mean` | Reduces multiple same-day observations to their mean |
| `hampel` | Spike detection using a sliding median window |
| `rolling_median` | Smoothing pass using a rolling median |

## Merging

Applied during `create_timeseries()`, after cleaning — combines all enabled missions' cleaned observations into one `merged_timeseries.csv` per reservoir using bias correction, SVR, and a Kalman filter. Override any key in `flows.DEFAULT_RESERVOIR_MERGING_OPTIONS` via a `merging_options` dict under `reservoirs:` in your config:

```yaml
reservoirs:
  merging_options:
    svr_radial_err: 1.0
    svr_radial_gamma: 0.00438
```

**If you use this pipeline, please also cite Schwatke et al. (2015)** — see the root [README](../README.md)'s Citation section. (DAHITI covers the windowed ADM-based outlier/uncertainty weighting, SVR, and Kalman filter specifically; the bias correction and multi-mission-combination logic are HydroEO's own contributions, not sourced from that paper.)

## Output structure

```
{main_dir}/
  aux/PLD/
    PLD_subset.gpkg            # SWOT Prior Lake Database subset
    present_in_pld.gpkg        # reservoirs matched to PLD lakes
    missing_in_pld.gpkg        # reservoirs NOT matched to PLD lakes
  raw/
    swot/                      # raw SWOT Lake SP granule downloads
    icesat2/                   # ICESat-2 download directory (empty; SlideRule writes to processed/)
    sentinel3/<id>/            # raw Sentinel-3 subsets per reservoir
    sentinel6/<id>/            # raw Sentinel-6 subsets per reservoir
  processed/
    icesat2/<id>/
      atl13.parquet            # ICESat-2 SlideRule output per reservoir
  results/
    <reservoir_id>/
      raw_observations/        # extracted, spatially filtered shapefiles per mission
      cleaned_observations/    # filter-cleaned CSVs per mission (+ dfs0 if enabled)
      merged_progress/         # intermediate merged CSVs
      merged_timeseries.csv    # all missions merged (wide format)
      all_cleaned_timeseries.csv
      crossing_summary.png
      cleaning_summary.png
      merging_summary.png
```

## Direct download (no config file)

Download individual satellite data using a bounding box without a config file:

```sh
# SWOT Lake SP
hydroeo fetch swot-lake \
  --bbox "-10 40 10 60" --start 2023-01-01 --end 2023-06-01 \
  --output ./output --username <user> --password <pass>

# ICESat-2 (no credentials needed; outputs atl13.parquet)
hydroeo fetch icesat2 \
  --bbox "-10 40 10 60" --start 2023-01-01 --end 2023-06-01 \
  --output ./output

# Sentinel-3 or Sentinel-6
hydroeo fetch sentinel \
  --bbox "-10 40 10 60" --start 2023-01-01 --end 2023-06-01 \
  --product S3 --output ./output \
  --creodias-username <user> --creodias-password <pass>
# Use --product S6 for Sentinel-6

# Copernicus DEM (free CDSE account required)
hydroeo fetch cop-dem \
  --bbox "-10 40 10 60" --output ./output \
  --cdse-username <user> --cdse-password <pass>
# Optional: --dataset "DEM90" (90 m; default is DEM30 = 30 m)
# Add quality layers: --dataset "DEM30,WBM,EDM"   (EDM, FLM, HEM, WBM also available)
```

Credentials can also be set as environment variables. Note the standalone `fetch` CLI commands use different variable names than the config-driven `Project` workflow above for Earthdata specifically:

| Command | Environment variables |
| --- | --- |
| `fetch swot-lake` / `swot-raster` / `swot-pixc` | `EARTHACCESS_USERNAME`, `EARTHACCESS_PASSWORD` |
| `fetch sentinel` | `CREODIAS_USERNAME`, `CREODIAS_PASSWORD` |
| `fetch cop-dem` | `CDSE_USERNAME`, `CDSE_PASSWORD` |

> **Note:** ICESat-2, Sentinel-3, and Sentinel-6 require reservoir or river targets for spatial filtering. Configuring them with `download: true` while neither a `reservoirs` nor a `rivers` section is present emits a `UserWarning` — they no longer require reservoirs specifically, since ICESat-2/Sentinel-3/6 also support rivers directly (see [rivers.md](rivers.md)).