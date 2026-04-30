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
A single unified configuration template covers all four use cases: [notebooks/example_config.yaml](./notebooks/example_config.yaml).

HydroEO projects are configured with:
- project output location, CRS, and shared date range
- the active water body branch (`reservoirs`, `rivers`, `swot_raster`, or `swot_pixc`) with `enabled: true`
- provider credentials (Earthdata, CREODIAS, HydroWeb API key) — or equivalent environment variables
- per-satellite download/process flags; dates default to `project.startdate`/`project.enddate` unless overridden

Minimal workflow:

```python
from HydroEO.project import Project

project = Project(name="my_altimetry_project", config="config.yaml")
project.initialize()
project.download()
project.create_timeseries()
project.generate_summaries()
```

## CLI

After installation, a `hydroeo` command is available on the PATH.

### Config-driven pipeline

Run any project lifecycle step directly from the terminal against a YAML config file:

```sh
# Full pipeline (initialize + download)
hydroeo run config.yaml

# Individual steps
hydroeo initialize config.yaml
hydroeo download    config.yaml
hydroeo update      config.yaml   # extend to today
hydroeo timeseries  config.yaml   # reservoirs only
hydroeo summaries   config.yaml   # reservoirs only
```

Use `--name` to override the project name (defaults to the config file stem) and `--verbose` / `-v` for debug logging:

```sh
hydroeo run config.yaml --name "Lake Geneva" --verbose
```

### Direct satellite downloads (no config file)

Download data for any satellite using a bounding box and date range, without writing a config file. Credentials can be passed as flags or via environment variables.

#### SWOT raster

```sh
hydroeo fetch swot-raster \
  --bbox "-10 40 10 60" \
  --start 2023-01-01 \
  --end   2023-06-01 \
  --output ./output \
  --username <earthaccess-user> \
  --password <earthaccess-pass>

# Or set credentials as environment variables:
export EARTHACCESS_USERNAME=...
export EARTHACCESS_PASSWORD=...
hydroeo fetch swot-raster --bbox "-10 40 10 60" --start 2023-01-01 --end 2023-06-01
```

Optional flags: `--aoi-name <label>` (default `aoi`), `--product <short-name>` (default `SWOT_L2_HR_Raster_D`).

#### SWOT pixel cloud

```sh
hydroeo fetch swot-pixc \
  --bbox "-10 40 10 60" \
  --start 2023-01-01 \
  --end   2023-06-01 \
  --output ./output \
  --username <earthaccess-user> \
  --password <earthaccess-pass>

# Or set credentials as environment variables:
export EARTHACCESS_USERNAME=...
export EARTHACCESS_PASSWORD=...
hydroeo fetch swot-pixc --bbox "-10 40 10 60" --start 2023-01-01 --end 2023-06-01
```

#### SWOT lake

```sh
hydroeo fetch swot-lake \
  --bbox "-10 40 10 60" \
  --start 2023-01-01 \
  --end   2023-06-01 \
  --output ./output
```

Credentials: `--username` / `--password` or `EARTHACCESS_USERNAME` / `EARTHACCESS_PASSWORD`.

#### ICESat-2

```sh
hydroeo fetch icesat2 \
  --bbox "-10 40 10 60" \
  --start 2023-01-01 \
  --end   2023-06-01 \
  --output ./output
```

No credentials required (uses SlideRule public API). Output is written as `atl13.parquet` inside the output directory.

#### Sentinel-3 / Sentinel-6

```sh
hydroeo fetch sentinel \
  --bbox "-10 40 10 60" \
  --start 2023-01-01 \
  --end   2023-06-01 \
  --product S3 \
  --output ./output \
  --creodias-username <user> \
  --creodias-password <pass>

# Or set credentials as environment variables:
export CREODIAS_USERNAME=...
export CREODIAS_PASSWORD=...
hydroeo fetch sentinel --bbox "-10 40 10 60" --start 2023-01-01 --end 2023-06-01 --product S6
```

Use `--product S3` for Sentinel-3 (default) or `--product S6` for Sentinel-6.


### Water body branches

A HydroEO project is built around a single **water body branch** — reservoirs/lakes, rivers, SWOT rasters, or SWOT Pixel Cloud for arbitrary areas. All four branches are mutually exclusive within one project config; a project targets one type only.

| Branch | Status | Products | Description |
| --- | --- | --- | --- |
| `reservoirs` | ✅ available | SWOT Lake, ICESat-2, Sentinel-3, Sentinel-6 | Closed water bodies: lakes, reservoirs. Defined by polygon shapefile + unique ID key. Supports multi-mission downloads and full timeseries processing. |
| `rivers` | 🧪 partial | SWOT Hydrocron | River projects support initialization plus SWOT Hydrocron download (public API, no credentials needed). With `rivers.aoi_path`, `initialize()` downloads SWORD v17b if needed, loads the continent-specific nodes/reaches gpkg, optionally buffers the AOI, subsets SWORD to the AOI, and maps each node/reach to the AOI `id_key` for output naming. Explicit `feature_numbers` (paired with `feature_type: nodes` or `feature_type: reaches`) are supported for SWOT downloads when `rivers.id` is provided as the output folder key. River preprocessing and plotting remain unimplemented. |
| `swot_raster` | ✅ available | SWOT Rasters (L2 HR/LR products) | Arbitrary areas of interest defined by bounding box or shapefile/geopackage geometries. Downloads SWOT raster products (SWOT_L2_HR_Raster_D, SWOT_L2_LR_SSH_2.0, SWOT_L2_HR_RIVERSP_2.0) via Earthdata. Extracts selected variables, applies quality filters, clips to AOI, and merges/reprojects tiles by date. Requires Earthdata credentials. |
| `swot_pixc` | ✅ available | SWOT Pixel Cloud (L2 PIXC) | Arbitrary areas of interest defined by bounding box or shapefile/geopackage geometries. Downloads SWOT Pixel Cloud data (SWOT_L2_HR_PIXC_D, SWOT_L2_HR_PIXC_2.0) via Earthdata. Extracts point data, filters by water class, computes derived heights (heightEGM), clips to AOI, and grids to regular rasters via binned statistics (median, mean, etc.) at specified resolution. Requires Earthdata credentials. |

### Project lifecycle

Each project exposes six top-level operations:

```
report() → initialize() → download() → update() → create_timeseries() → generate_summaries()
```

| Method | What it does |
| --- | --- |
| `project.report()` | Logs the project name, logs the number of loaded water bodies in the active branch, and returns a preview of the branch GeoDataFrame (`gdf.head()`). Useful as a quick inspection step after loading a project. |
| `project.initialize()` | Loads and validates config, resolves credentials, prepares output directories, and reads the water body input. For rivers with `aoi_path`, it ensures the SWORD v17b database exists under `project.main_dir`, then subsets the selected continent/layer to the AOI, using `rivers.buffer_meters` before the spatial filter when configured. For SWOT rasters and Pixel Cloud, validates AOI geometry and product selection. Reports all config issues in one message before doing any I/O. |
| `project.download()` | Downloads raw satellite data within the configured date range. **Reservoirs:** downloads from all enabled missions (SWOT, ICESat-2, Sentinel-3, Sentinel-6) keeping their existing raw-download layout. **Rivers:** SWOT Hydrocron API calls write filtered CSV outputs to `swot/rivers/<per_id>/<node_or_reach_id>/timeseries.csv`, where `<per_id>` comes from `rivers.id_key` value or fallback `rivers.id`. **SWOT Rasters:** downloads SWOT raster granules matching the AOI bounds and temporal range using Earthdata credentials, outputs to `swot_raster/<aoi_name>/raw/<product>/`. **SWOT Pixel Cloud:** downloads SWOT Pixel Cloud granules, filters by water class, computes derived heights, clips to AOI, grids via binned statistics, outputs to `swot_pixc/<aoi_name>/raster/`. |
| `project.update()` | Extends existing downloads from the latest observation up to today. Safe to run repeatedly — already-downloaded data is not re-fetched. |
| `project.create_timeseries()` | Extracts observations spatially matched to each water body, applies the configured cleaning filters, and writes per-body timeseries shapefiles. |
| `project.generate_summaries(show, save)` | Produces diagnostic plots for each water body and mission: raw crossings, filter effect, and merged multi-mission timeseries. |

## Products & Data Availability

This section consolidates all available data sources across HydroEO's four workflows.

### Reservoir & Lake Products

| Product | Satellite | Archive start | Credentials required |
| --- | --- | --- | --- | --- |
| `SWOT_L2_HR_LakeSP_D` | SWOT | 2023 | Earthdata + HydroWeb API key |
| `ATL13` | ICESat-2 | 2018 | None (SlideRule public API) |
| - | Sentinel-3A/B | 2016 | CREODIAS account |
| - | Sentinel-6 | 2020 | CREODIAS account |

### River Products

| Product | Satellite | Archive start | Credentials |
| --- | --- | --- | --- |
| SWOT Hydrocron (nodes/reaches) | SWOT | 2023 | None (public API) |

### SWOT Raster Products

| Product | Resolution | Archive start | Use case |
| --- | --- | --- | --- |
| `SWOT_L2_HR_Raster_D` | ~100m | 2024 | High-resolution grids: elevation, water surface height |
| `SWOT_L2_LR_SSH_2.0` | ~250m | 2023 | Low-resolution grids: sea surface height, open water |
| `SWOT_L2_HR_RIVERSP_2.0` | ~100m | 2024 | River surface: channel elevation, width, slope |

All raster products require Earthdata credentials.

### SWOT Pixel Cloud Products

| Product | Archive start | Use case |
| --- | --- | --- |
| `SWOT_L2_HR_PIXC_D` | 2024 | Point-based water heights, gridded to rasters via binning |
| `SWOT_L2_HR_PIXC_2.0` | 2024 | Point-based water heights (v2 format), gridded to rasters |

Both Pixel Cloud products require Earthdata credentials and support median/mean/max/min binning statistics.

### Output structure

All outputs are written under `project.main_dir`. Structure varies by branch:

#### Reservoirs
```
main_dir/
  <water_body_id>/
    raw_observations/      # extracted, spatially filtered shapefiles per mission
    cleaned_observations/  # filter-cleaned shapefiles per mission
    merged_progress/       # merged CSVs across all enabled missions per cleaning method
    ./                     # PNG diagnostics (crossings, cleaning, merging), single merged CSV across all enabled missions
  swot/                    # raw SWOT Lake SP downloads
  icesat2/                 # raw ICESat-2 downloads
  sentinel3/               # raw Sentinel-3 downloads
  sentinel6/               # raw Sentinel-6 downloads
```

#### Rivers
```
main_dir/
  swot/rivers/
    <aoi_name>/
      timeseries.csv       # SWOT Hydrocron data
```

#### SWOT Rasters
```
main_dir/
  swot_raster/<aoi_name>/
    raw/<product>/                # raw netCDF granule files
    processed/<product>/          # extracted/clipped GeoTIFFs per granule (temporary, deleted after merge)
    merged/                       # merged mosaics by date/variable
```

#### SWOT Pixel Cloud
```
main_dir/
  swot_pixc/<aoi_name>/
    raw/<product>/                # raw netCDF granule files
    trimmed/                      # preprocessed GeoJSON point data (filtered by water class)
    raster/                       # gridded rasters by date/variable
```

### Reservoir-specific: Cleaning filters

Applied during `create_timeseries()` for reservoir workflows only. Configured per mission under `processing_filters`:

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

### Credentials & Authentication Reference

This table consolidates all credential requirements across workflows and missions:

| Mission/Workflow | Product | Credentials | Environment variables | Where needed |
| --- | --- | --- | --- | --- |
| Reservoirs (SWOT) | `SWOT_L2_HR_LakeSP_D` | Earthdata account + HydroWeb API key | `EARTHDATA_USERNAME`, `EARTHDATA_PASSWORD`, `EODAG__HYDROWEB_NEXT__AUTH__CREDENTIALS__APIKEY` | PLD lake matching for reservoir outline correction |
| Reservoirs (ICESat-2) | `ATL13` | None required | — | SlideRule public API; optional local caching |
| Reservoirs (Sentinel-3) | `S3` | CREODIAS account | `CREODIAS_USERNAME`, `CREODIAS_PASSWORD` | CDSE/Copernicus download |
| Reservoirs (Sentinel-6) | `S6` | CREODIAS account | `CREODIAS_USERNAME`, `CREODIAS_PASSWORD` | CDSE/Copernicus download |
| Rivers (SWOT) | SWOT Hydrocron | None required | — | Public API; no auth needed |
| SWOT Raster | `SWOT_L2_HR_Raster_D`, `SWOT_L2_LR_SSH_2.0`, `SWOT_L2_HR_RIVERSP_2.0` | Earthdata account | `EARTHDATA_USERNAME`, `EARTHDATA_PASSWORD` | Earthdata/NASA download |
| SWOT Pixel Cloud | `SWOT_L2_HR_PIXC_D`, `SWOT_L2_HR_PIXC_2.0` | Earthdata account | `EARTHDATA_USERNAME`, `EARTHDATA_PASSWORD` | Earthdata/NASA download |

Note: Credentials can be provided in the config file under `earthaccess`, `hydroweb`, or `creodias` sections, or as environment variables. Environment variables take precedence and are recommended for automated/CI workflows.

## Logging

HydroEO logs to both console and file for debugging and monitoring:

- **Console**: INFO level (key operations and progress messages)
- **File**: DEBUG level (detailed execution trace, useful for troubleshooting)

Log files are automatically created in a `logs/` folder at the repository root with timestamped filenames (e.g., `HydroEO_2026-04-30_10-00-17.log`). The `logs/` directory is created automatically on the first run and excluded from git.

### Controlling logging output

By default, both console and file logging are enabled:

```python
from HydroEO.project import Project
from HydroEO.logging_config import setup_logging
import logging

setup_logging(logging.INFO)  # Console at INFO, file at DEBUG
```

To disable file logging (e.g., in tests or scripts where you only want console output):

```python
setup_logging(logging.INFO, enable_file_logging=False)
```

To enable DEBUG messages on the console:

```python
setup_logging(logging.DEBUG)  # Console and file both at DEBUG
```

### Log file location

All log files are stored in:
```
<project_root>/logs/
```

Log files follow the naming pattern: `HydroEO_YYYY-MM-DD_HH-MM-SS.log`

Multiple runs create separate log files, allowing you to review execution traces from previous runs without overwriting.

### Activating a use case
Each water body branch (`reservoirs`, `rivers`, `swot_raster`) supports an `enabled` flag (default: `true` when the section is present). Set `enabled: false` to keep a section in the file without activating it — useful in the unified template when switching between use cases without deleting sections.

### Shared date range
`project.startdate` and `project.enddate` act as a global fallback: any satellite section that omits its own `startdate`/`enddate` inherits these values. Per-satellite overrides are still supported by setting `startdate`/`enddate` inside the mission section.

### Incompatible satellite sources
ICESat-2, Sentinel-3, and Sentinel-6 require reservoir waterbody polygons for spatial filtering and have no effect in river, SWOT raster, or SWOT Pixel Cloud projects. Configuring them with `download: true` or `process: true` alongside a `rivers`, `swot_raster`, or `swot_pixc` branch (without a `reservoirs` section) emits a `UserWarning` at project load time.

### Reservoir-specific: Common per-mission processing options
Applied during `create_timeseries()` for reservoir workflows only. Available under each mission section (`swot`, `icesat2`, `sentinel3`, `sentinel6`):
- `processing_filters`: list of filters. Allowed values: `elevation`, `MAD`, `daily_mean`, `hampel`, `rolling_median`.
- `elevation_min_m`: lower bound used by elevation filter.
- `elevation_max_m`: upper bound used by elevation filter.
- `mad_threshold`: MAD outlier multiplier.

### SWOT
- `pld_match_max_distance_m`: max nearest-neighbour distance for PLD matching (Use Case A — reservoirs).
- `exclude_obs_id_values`: list of SWOT `obs_id` values to exclude during extraction (Use Case A).
- `hydrocron_fields.nodes` and `hydrocron_fields.reaches`: optional field lists for river Hydrocron requests (Use Case B — rivers).
- `quality_filters.nodes.max_q` and `quality_filters.reaches.max_q`: optional Hydrocron river quality thresholds (Use Case B). The default `2` keeps records where `node_q <= 2` or `reach_q <= 2`.
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
### SWOT Pixel Cloud

SWOT Pixel Cloud workflow configuration:

- `aoi.name`: identifier for this AOI (used in output directory paths and log files).
- `aoi.type`: one of `"bbox"`, `"shapefile"`, or `"geopackage"`. For non-bbox types, provide `aoi.path`.
- `aoi.bbox`: for type `"bbox"`, array of `[lon_min, lat_min, lon_max, lat_max]`.
- `aoi.path`: for type `"shapefile"` or `"geopackage"`, path to the geometry file.
- `product`: one of `"SWOT_L2_HR_PIXC_D"` (default) or `"SWOT_L2_HR_PIXC_2.0"`.
- `startdate` and `enddate`: temporal range as `[year, month, day]`.
- `classes`: optional water classification filtering (default: `open_water`, `water_near_land`). Available values: `land`, `land_near_water`, `water_near_land`, `open_water`, `dark_water`, `low_coh_water_near_land`, `open_low_coh_water`.
- `fields`: optional list of point fields to extract and grid (default: `heightEGM`). Available fields include: `heightEGM`, `height`, `geoid`, `water_frac`, `phase_noise_std`, `dheight_dphase`, `sig0`, etc.
- `grid_resolution`: output pixel size in meters (default: `100`). Must be specified; used for all gridding operations.
- `stat_method`: statistic method for binning (default: `"median"`). Choices: `median`, `mean`, `max`, `min`. Controls how multiple point observations in each grid cell are aggregated.
- `target_crs`: optional override — output EPSG code for the final raster (e.g., `"EPSG:32645"` for UTM 45N). By default uses `gis.global_crs`. Only set this when you need a different output CRS.

The download phase queries Earthdata for PIXC granules within the AOI bounding box and temporal range, tracking downloaded granules in `raw/<product>/downloaded.log`. Preprocessing extracts points matching the configured water classes, computes derived heights (e.g., `heightEGM` from geoid corrections), clips to AOI bounds, and saves to GeoJSON. The rasterization phase bins point data into regular grids using the specified statistic and resolution, outputting GeoTIFFs by date and variable. Raw netCDF files are retained for archival but can be deleted manually to save space after rasterization completes.
### Sentinel-3 and Sentinel-6
- `subset_file_id`: expected source netCDF filename inside each package.
- `sigma0_max`: maximum accepted `sigma0` during extraction.
- `download_threads`: thread count used for CREODIAS download batching.

### ICESat-2 ATL13 (SlideRule)

ICESat-2 ATL13 data is downloaded via [SlideRule](https://slideruleearth.io) using its
`atl13x` endpoint.  SlideRule streams results directly as a GeoDataFrame — no local HDF5
files are written.

SlideRule does **not** require Earthdata authentication — no credentials are needed for
the ICESat-2 download path.

If `download_dir` is provided, the returned GeoDataFrame is cached as `atl13.parquet`
inside that directory.  If `download_dir` is omitted, the directory defaults to
`{main_dir}/icesat2/` and no extra cache is written between download and process steps.

**Available ICESat-2 configuration options:**

- `atl13_fields`: optional list of ancillary beam-group fields from SlideRule (validated at config load time). Supported fields: `inland_water_body_id`, `inland_water_body_size`, `inland_water_body_type`, `segment_slope_trk_bdy`, `err_ht_water_surf`, `segment_quality`, `segment_geoid`, `segment_geoid_free2mean`, `segment_dem_ht`, `segment_near_sat_fract`.
- `atl13.pass_invalid`: allow invalid records to pass filtering (default: false).
- `atl13.beams`: optional list of beam names to select (e.g., `["gt1l", "gt2l", "gt3l"]` for left beams only).
- `atl13.spots`: optional list of spot numbers.

Common output columns: `height` (water surface elevation, EGM2008-corrected), `cycle_number`, `beam`, `rgt`, `ht_water_surf`, `stdev_water_surf`, `water_depth`, `spot`, `srcid`, and any requested ancillary fields.

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