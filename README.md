![](images/HydroEO.png)
# Easy access Earth observation data for water resource applications

Repo to allow users with little EO (Earth Observation) knowledge to access and download altimetry over reservoirs and lakes for integration into larger water resource projects.

> [!CAUTION]
> HydroEO has had limited testing and further developments are likely to come. Please report any bugs or issues here: https://github.com/DHI/HydroEO/issues

## Installation
> Note: Ensure that `uv` and `git` are installed on your system before running the following command.

**System dependency:** If you plan to download SWOT data, you must have **spatialite** installed at the OS level for spatial SQLite query support:

```sh
# macOS (Homebrew)
brew install spatialite

# Ubuntu/Debian
sudo apt-get install libspatialite-dev libspatialite7

# Fedora/RHEL
sudo dnf install spatialite-libs spatialite-devel
```

Then install HydroEO:

```sh
git clone https://github.com/DHI/HydroEO.git
cd HydroEO
uv sync
```
#### Python versions
HydroEO currently runs on Python 3.9 - 3.12.

## Quick start
Usage examples are available in [notebooks](./notebooks), with sample configurations for [reservoirs](./notebooks/example_config_reservoirs.yaml), [rivers](./notebooks/example_config_rivers.yaml), and [SWOT rasters](./notebooks/example_config_swot_raster.yaml).

HydroEO projects are configured with:
- project output location and CRS
- reservoir polygons and unique reservoir ID key
- provider credentials (Earthdata, CREODIAS, HydroWeb API key)
- per-product download/process flags and date ranges

Minimal workflow:

```python
from HydroEO.project import Project

project = Project(name="my_altimetry_project", config="config.yaml")
project.initialize()
project.download()
project.create_timeseries()
project.generate_summaries()
```

## Project capabilities

A HydroEO project is built around a single **water body branch** — reservoirs/lakes, rivers, or SWOT rasters for arbitrary areas. All three branches are mutually exclusive within one project config; a project targets one type only.

### Water body branches

| Branch | Status | Products | Description |
| --- | --- | --- | --- |
| `reservoirs` | ✅ available | SWOT Lake, ICESat-2, Sentinel-3, Sentinel-6 | Closed water bodies: lakes, reservoirs. Defined by polygon shapefile + unique ID key. Supports multi-mission downloads and full timeseries processing. |
| `rivers` | 🧪 partial | SWOT Hydrocron | River projects support initialization plus SWOT Hydrocron download (public API, no credentials needed). With `rivers.aoi_path`, `initialize()` downloads SWORD v17b if needed, loads the continent-specific nodes/reaches gpkg, optionally buffers the AOI, subsets SWORD to the AOI, and maps each node/reach to the AOI `id_key` for output naming. Explicit `node_numbers` / `reach_numbers` inputs are also supported for SWOT downloads when `rivers.id` is provided as the output folder key. River preprocessing and plotting remain unimplemented. |
| `swot_raster` | ✅ available | SWOT Rasters (L2 HR/LR products) | Arbitrary areas of interest defined by bounding box or shapefile/geopackage geometries. Downloads SWOT raster products (SWOT_L2_HR_Raster_D, SWOT_L2_LR_SSH_2.0, SWOT_L2_HR_RIVERSP_2.0) via Earthdata. Extracts selected variables, applies quality filters, clips to AOI, and merges/reprojects tiles by date. Requires Earthdata credentials. |

### Project lifecycle

Each project exposes six top-level operations:

```
report() → initialize() → download() → update() → create_timeseries() → generate_summaries()
```

| Method | What it does |
| --- | --- |
| `project.report()` | Logs the project name, logs the number of loaded water bodies in the active branch, and returns a preview of the branch GeoDataFrame (`gdf.head()`). Useful as a quick inspection step after loading a project. |
| `project.initialize()` | Loads and validates config, resolves credentials, prepares output directories, and reads the water body input. For rivers with `aoi_path`, it ensures the SWORD v17b database exists under `project.main_dir`, then subsets the selected continent/layer to the AOI, using `rivers.buffer_meters` before the spatial filter when configured. For SWOT rasters, validates AOI geometry and product selection. Reports all config issues in one message before doing any I/O. |
| `project.download()` | Downloads raw satellite data within the configured date range. **Reservoirs:** downloads from all enabled missions (SWOT, ICESat-2, Sentinel-3, Sentinel-6) keeping their existing raw-download layout. **Rivers:** SWOT Hydrocron API calls write filtered CSV outputs to `swot/rivers/<per_id>/<node_or_reach_id>/timeseries.csv`, where `<per_id>` comes from `rivers.id_key` value or fallback `rivers.id`. **SWOT Rasters:** downloads SWOT raster granules matching the AOI bounds and temporal range using Earthdata credentials, outputs to `swot_raster/<aoi_name>/raw/<product>/`. |
| `project.update()` | Extends existing downloads from the latest observation up to today. Safe to run repeatedly — already-downloaded data is not re-fetched. |
| `project.create_timeseries()` | Extracts observations spatially matched to each water body, applies the configured cleaning filters, and writes per-body timeseries shapefiles. |
| `project.generate_summaries(show, save)` | Produces diagnostic plots for each water body and mission: raw crossings, filter effect, and merged multi-mission timeseries. |

### Supported missions

#### Reservoirs / Lakes
A reservoir project can enable any combination of the four missions below. Each mission is independently toggled (`download: true/false`, `process: true/false`) and has its own date range.

| Mission key | Product key | Satellite | Archive start | Required credentials |
| --- | --- | --- | --- | --- |
| `swot` | `SWOT_LAKE` | SWOT Lake SP | 2023 | Earthdata account + HydroWeb API key (PLD) |
| `icesat2` | `ATL13` | ICESat-2 ATL13 | 2018 | None (SlideRule — no auth needed) |
| `sentinel3` | `S3` | Sentinel-3A/B | 2016 | CREODIAS account |
| `sentinel6` | `S6` | Sentinel-6 | 2020 | CREODIAS account |

#### Rivers
River projects support SWOT data only, via the Hydrocron public API (no credentials required):

| Product | Source | Archive start | Credentials |
| --- | --- | --- | --- |
| SWOT Hydrocron (nodes/reaches) | SWOT | 2023 | None (public API) |

#### SWOT Rasters
SWOT raster projects support three high-resolution and low-resolution SWOT Level-2 products for arbitrary areas:

| Product key | Description | Archive start | Credentials |
| --- | --- | --- | --- |
| `SWOT_L2_HR_Raster_D` | High-resolution raster data (~100m pixels) | 2024 | Earthdata account |
| `SWOT_L2_LR_SSH_2.0` | Low-resolution sea-surface-height raster (~250m pixels) | 2023 | Earthdata account |
| `SWOT_L2_HR_RIVERSP_2.0` | High-resolution river surface product (~100m pixels) | 2024 | Earthdata account |

### Output structure

All outputs are written under `project.main_dir`. Structure varies by branch:

#### Reservoirs
```
main_dir/
  <water_body_id>/
    raw_observations/      # extracted, spatially filtered shapefiles per mission
    cleaned_timeseries/    # filter-cleaned shapefiles per mission
    merged_timeseries/     # single merged CSV across all enabled missions
    plots/                 # PNG diagnostics (crossings, cleaning, merging)
  swot/                    # raw SWOT Lake SP downloads
  icesat2/                 # raw ICESat-2 downloads
  sentinel3/               # raw Sentinel-3 downloads
  sentinel6/               # raw Sentinel-6 downloads
```

#### Rivers
```
main_dir/
  swot/rivers/
    <per_id>/<node_or_reach_id>/
      timeseries.csv       # SWOT Hydrocron data per node/reach
```

#### SWOT Rasters
```
main_dir/
  swot_raster/<aoi_name>/
    raw/<product>/                # raw netCDF granule files
    processed/<product>/          # extracted/clipped GeoTIFFs per granule (temporary, deleted after merge)
    merged/                       # merged mosaics by date/variable
      20231004_wse_merged.tif
      20231004_wse_uncert_merged.tif
      20231005_wse_merged.tif
      ...
```

### Cleaning filters

Applied during `create_timeseries()`. Configured per mission under `processing_filters`:

| Filter | What it removes |
| --- | --- |
| `elevation` | Observations outside `[elevation_min_m, elevation_max_m]` |
| `MAD` | Outliers beyond `mad_threshold × MAD` from the median |
| `daily_mean` | Reduces multiple same-day observations to their mean |
| `hampel` | Spike detection using a sliding median window |
| `rolling_median` | Smoothing pass using a rolling median |

Configuration validation happens during `project.initialize()`, which checks required sections, input paths, date ranges, and optional parameter values before any downloads start.

## Configuration controls
HydroEO exposes mission-level optional parameters with defaults that preserve prior behaviour. Existing configs continue to run unchanged.

### Common per-mission processing options
Available under each mission section (`swot`, `icesat2`, `sentinel3`, `sentinel6`):
- `processing_filters`: list of filters. Allowed values: `elevation`, `MAD`, `daily_mean`, `hampel`, `rolling_median`.
- `elevation_min_m`: lower bound used by elevation filter.
- `elevation_max_m`: upper bound used by elevation filter.
- `mad_threshold`: MAD outlier multiplier.

### SWOT
- `pld_match_max_distance_m`: max nearest-neighbour distance for PLD matching.
- `exclude_obs_id_values`: list of SWOT `obs_id` values to exclude during extraction.
- `hydrocron_fields.nodes` and `hydrocron_fields.reaches`: optional field lists for river Hydrocron requests.
- `quality_filters.nodes.max_q` and `quality_filters.reaches.max_q`: optional Hydrocron river quality thresholds. The default `2` keeps records where `node_q <= 2` or `reach_q <= 2`.
- River SWOT downloads currently stop at `project.download()` and write CSV outputs under `swot/rivers/<per_id>/<node_or_reach_id>/timeseries.csv`. Each node/reach gets its own subfolder within the AOI-feature folder, so multiple nodes from the same AOI feature coexist safely.

### SWOT Rasters
SWOT raster workflow configuration:
- `aoi.name`: identifier for this AOI (used in output directory paths and log files).
- `aoi.type`: one of `"bbox"`, `"shapefile"`, or `"geopackage"`. For non-bbox types, provide `aoi.path`.
- `aoi.bbox`: for type `"bbox"`, array of `[lon_min, lat_min, lon_max, lat_max]`.
- `aoi.path`: for type `"shapefile"` or `"geopackage"`, path to the geometry file.
- `product`: one of `"SWOT_L2_HR_Raster_D"`, `"SWOT_L2_LR_SSH_2.0"`, or `"SWOT_L2_HR_RIVERSP_2.0"`.
- `startdate` and `enddate`: temporal range as `[year, month, day]`.
- `granule_filter`: optional filename pattern to filter granules (e.g., `"*100m*"` for 100m products).
- `merge_tiles`: optional boolean (default `true`). When `true`, processed tiles are merged by date/variable and reprojected. Set to `false` to keep individual clipped tiles without merging.
- `target_crs`: advanced optional override — output EPSG code for the merge/reproject phase (e.g., `"EPSG:32645"` for UTM 45N). By default the merge phase uses `gis.global_crs`. Only set this when you need a different output CRS than the project-level CRS.
- `variables`: optional list of variables to extract. Defaults to all 8 variables: `wse`, `wse_uncert`, `wse_qual`, `height_cor_xover`, `geoid`, `n_wse_pix`, `n_other_pix`, `layover_impact`. Regardless of the user-supplied list, `wse`, `wse_uncert`, and `layover_impact` are always extracted as they are required for quality masking.
- `quality_filters.max_wse_uncert`: pixels with `wse_uncert` at or above this value are masked. Default `0.3` (metres).
- `quality_filters.max_layover_impact`: pixels with `layover_impact` at or above this value are masked. Default `0.3` (metres).

The download phase tracks already-processed granules in `raw/<product>/downloaded.log` to avoid re-downloading on config reruns. Preprocessing extracts the configured variables, applies quality filters, clips each tile to the AOI, and outputs GeoTIFFs. The merge phase groups clipped tiles by date and variable, merges multiple passes, and reprojects to the target CRS. Only tiles that spatially overlap the AOI are processed — out-of-area tiles (e.g., distant UTM zones fetched by the search API) are discarded before extraction. Raw netCDF files are deleted after preprocessing to save space.

### Sentinel-3 and Sentinel-6
- `subset_file_id`: expected source netCDF filename inside each package.
- `sigma0_max`: maximum accepted `sigma0` during extraction.
- `download_threads`: thread count used for CREODIAS download batching.

### ICESat-2 ATL13 (SlideRule)

ICESat-2 ATL13 data is downloaded via [SlideRule](https://slideruleearth.io) using its
`atl13x` endpoint.  SlideRule streams results directly as a GeoDataFrame — no local HDF5
files are written.

**Important:** `atl13x` is **water-body-centric**, not polygon-centric.  SlideRule uses
an internal ATL13 Metadata Service (AMS) to identify which granules contain a given water
body.  The polygon centroid is passed as `atl13.coord` for the AMS lookup — this is the
only spatial parameter sent to the server.  The polygon is **not** passed as a spatial
filter (doing so causes empty results in SlideRule v5 because track segments that cross
outside a tight bounding box are excluded).  Precise spatial filtering is applied
client-side in `extract_observations()` using `gdf.within()` after download.
The reservoir must be registered in HydroLAKES or GRWL for the download to return data.

SlideRule does **not** require Earthdata authentication — no credentials are needed for
the ICESat-2 download path.

If `download_dir` is provided, the returned GeoDataFrame is cached as `atl13.parquet`
inside that directory.  If `download_dir` is omitted, the directory defaults to
`{main_dir}/icesat2/` and no extra cache is written between download and process steps.

#### Column changes from prior Harmony-based pipeline

| Column | Status | Notes |
| --- | --- | --- |
| `height` | ✅ present | From `ht_ortho` (EGM2008 geoid-corrected server-side — do **not** apply a second correction) |
| `cycle_number` | ✅ present | From SlideRule `cycle` column |
| `beam` | ✅ present | From SlideRule `gt` column |
| `rgt` | ✅ present | Unchanged |
| `ht_water_surf` | 🆕 new | Instantaneous water surface height |
| `stdev_water_surf` | 🆕 new | Standard deviation of water surface height |
| `water_depth` | 🆕 new | Estimated water depth |
| `spot` | 🆕 new | ICESat-2 spot number |
| `srcid` | 🆕 new | Source granule ID (map to name via `gdf.attrs["meta"]["srctbl"]`) |
| `orbit_number` | ❌ removed | Not available from `atl13x`; use `rgt + cycle_number` as equivalent orbit identifier |
| `file_name` | ❌ removed | No local files; use `gdf.attrs["meta"]["srctbl"][srcid]` if the granule name is needed |

#### Supported `atl13_fields` values

Use `icesat2.atl13_fields` to request additional ancillary beam-group fields from
SlideRule.  These are fetched server-side and returned as extra columns.
Field names are validated at config load time.

Supported ancillary field names:
- `inland_water_body_id`
- `inland_water_body_size`
- `inland_water_body_type`
- `segment_slope_trk_bdy`
- `err_ht_water_surf`
- `segment_quality`
- `segment_geoid`
- `segment_geoid_free2mean`
- `segment_dem_ht`
- `segment_near_sat_fract`

Field mapping from config key to output column:

| Config field name | Output column | Description |
| --- | --- | --- |
| `inland_water_body_id` | `wb_id` | Unique inland water-body identifier. |
| `inland_water_body_size` | `wb_size` | Estimated water-body size class. |
| `inland_water_body_type` | `wb_type` | Inland water-body type classification. |
| `segment_slope_trk_bdy` | `wb_slope` | Along-track boundary slope estimate. |
| `err_ht_water_surf` | `height_err` | Estimated uncertainty of water-surface height. |
| `segment_quality` | `quality_seg` | Segment-level quality flag. |
| `segment_geoid` | `geoid_track` | Geoid height along the segment track. |
| `segment_geoid_free2mean` | `geoid_corr_track` | Free-to-mean geoid correction term. |
| `segment_dem_ht` | `dem` | Reference DEM elevation sampled at segment location. |
| `segment_near_sat_fract` | `sat_frac_track` | Fraction of pulses flagged as near-saturated. |

Example:

```yaml
icesat2:
  atl13_fields:
    - inland_water_body_id
    - segment_quality
  atl13:
    pass_invalid: false
    beams: ["gt1l", "gt2l", "gt3l"]  # left beams only
    spots: []
```

## Available products

### Lake/Reservoir workflows
HydroEO supports 4 products for reservoir workflows. These are the product identifiers you will see in extracted outputs, merged timeseries, and lower-level processing methods:

| Product key | Source mission | Typical coverage | Key retrievable parameters |
| --- | --- | --- | --- |
| `SWOT_LAKE` | SWOT Lake SP | 2023-present | water surface elevation, water area, storage-related metrics, quality flags |
| `ATL13` | ICESat-2 | 2018-present | inland water surface elevation, along-track geometry, quality metadata |
| `S3` | Sentinel-3A/B | 2016-present | latitude, longitude, elevation, waveform, backscatter (`sig0`), time, altitude, tracker/range variables, geoid |
| `S6` | Sentinel-6 | 2020-present | latitude, longitude, altitude, `range_ocog`, wet/dry tropospheric corrections, backscatter, time, geoid |

Common across reservoir products: geolocation, observation time, water level/elevation measurements, and quality metadata.

### River workflows
River projects use SWOT Hydrocron, a public API delivering SWOT observations for predefined river nodes and reaches. Returns timeseries CSV with water surface elevation, discharge estimates, and quality flags.

### SWOT Raster workflows
SWOT raster projects download three L2 raster products via Earthdata. Each product is distributed as geotiff granules, one per satellite pass:

| Product | Resolution | Coverage | Key parameters |
| --- | --- | --- | --- |
| `SWOT_L2_HR_Raster_D` | ~100m | Land/water classification, elevation | digital elevation model (DEM), water surface elevation, water area fraction |
| `SWOT_L2_LR_SSH_2.0` | ~250m | Sea surface height, open water | sea surface height, geoid height, wind speed, wave height |
| `SWOT_L2_HR_RIVERSP_2.0` | ~100m | River surface | river channel elevation, width, slope, quality flags |

## Running tests
The development environment with all testing tools is set up via:

```sh
uv sync --all-extras
```

Once set up, run test suites with:

```sh
# Run all tests
make test

# Unit tests (fast, mocked)
uv run pytest -m unit

# API contract tests (endpoint/schema checks)
uv run pytest -m api_contract

# Integration tests (live API calls)
uv run pytest -m integration -v

# Coverage report
make coverage
```

For just the core package (without dev/test/notebook dependencies), use `uv sync` as shown in Installation.

### Pytest markers in this repo
Markers are defined in [pyproject.toml](./pyproject.toml):

- `unit`: fast mocked tests with no external network dependency.
- `integration`: live tests that call external services and require valid credentials.
- `api_contract`: tests focused on endpoint reachability and response schema expectations.

### Credentials for tests
Both Earthdata variable naming schemes are needed in practice:
- some tests and wrappers use `EDL_USERNAME` / `EDL_PASSWORD`
- earthaccess login paths use `EARTHDATA_USERNAME` / `EARTHDATA_PASSWORD`

Set both pairs to the same values when running the full suite:

```sh
export EDL_USERNAME="your_earthdata_username"
export EDL_PASSWORD="your_earthdata_password"
export EARTHDATA_USERNAME="$EDL_USERNAME"
export EARTHDATA_PASSWORD="$EDL_PASSWORD"
export CREODIAS_USERNAME="your_creodias_username"
export CREODIAS_PASSWORD="your_creodias_password"
export HYDROWEB_API_KEY="your_hydroweb_api_key"
```

## Project structure
Key modules:
- `HydroEO/project.py`: high-level project workflow orchestration.
- `HydroEO/waterbody.py`: shared water body base class and reservoir branch implementation.
- `HydroEO/flows/`: orchestration split into download, preprocessing, and plotting flows.
- `HydroEO/plotting.py`: shared plotting helpers for crossings, cleaning, and merged outputs.
- `HydroEO/satellites/swot.py`: public SWOT facade delegating to mission-specific download and preprocess modules.
- `HydroEO/satellites/icesat2.py`: public ICESat-2 facade delegating to mission-specific download and preprocess modules.
- `HydroEO/satellites/sentinel.py`: public Sentinel facade delegating to Sentinel-3 and Sentinel-6 download and preprocess logic.
- `HydroEO/downloaders/creodias.py`: CREODIAS/CDSE queries and download helpers.
- `HydroEO/downloaders/hydroweb.py`: HydroWeb integration and PLD support.

## CI/CD and GitHub setup

**Automated checks** run on every push and pull request via GitHub Actions (see [.github/workflows/ci.yml](.github/workflows/ci.yml)):

1. **Lint** (ruff) — checks code style and imports
2. **Unit tests** — runs on Ubuntu and Windows, Python 3.9 and 3.12; publishes coverage to Codecov
3. **Integration tests** — runs live API calls (only on official repo when enabled)

**For repo maintainers:** To enable live integration tests, set these GitHub secrets and variables:

*Secrets:*
- `EDL_USERNAME` — Earthdata Login username for NASA SWOT/ICESat-2 access
- `EDL_PASSWORD` — Earthdata Login password
- `CREODIAS_USERNAME` — Copernicus CDSE username for Sentinel-3/6 access
- `CREODIAS_PASSWORD` — Copernicus CDSE password

*Repository variables:*
- `RUN_INTEGRATION_TESTS` = `true` — enables integration tests on push (forks and PRs will skip them automatically)

Integration tests require live credentials and are gated to prevent exposure on untrusted branches.