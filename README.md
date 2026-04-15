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
Usage examples are available in [notebooks](./notebooks), and a complete sample configuration is in [notebooks/example_config.yaml](./notebooks/example_config.yaml).

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

Configuration is validated up-front. `Project.initialize()` now calls `Project.validate_config()` and reports common issues in one error message (missing keys, invalid dates, invalid optional parameter values, and missing input paths).

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

### Sentinel-3 and Sentinel-6
- `subset_file_id`: expected source netCDF filename inside each package.
- `sigma0_max`: maximum accepted `sigma0` during extraction.
- `download_threads`: thread count used for CREODIAS download batching.

### ICESat-2 ATL13 field selection
Use `icesat2.atl13_fields` to request ATL13 fields. Field names are validated at config load time.

Configuration uses canonical ATL13 field names. Extracted tabular and shapefile outputs use the column names shown below.

Supported ATL13 field names:
- `ht_ortho`
- `segment_lat`
- `segment_lon`
- `delta_time`
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
- `orbit_number`
- `rgt`
- `cycle_number`
- `beam`
- `file_name`

Field mapping from config key to output column:

| Config field name | Output column | Description |
| --- | --- | --- |
| `ht_ortho` | `height` | Orthometric water surface height for each segment. |
| `segment_lat` | `lat` | Segment latitude in decimal degrees. |
| `segment_lon` | `lon` | Segment longitude in decimal degrees. |
| `delta_time` | `date` | Seconds since ATLAS SDP epoch, converted to UTC datetime. |
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
| `orbit_number` | `orbit` | Orbit number for the granule. |
| `rgt` | `rgt` | Reference ground track identifier. |
| `cycle_number` | `cycle_number` | Repeat-cycle number for the granule. |
| `beam` | `beam` | ICESat-2 beam key used for extraction, for example `gt1l`. |
| `file_name` | `file_name` | Source ATL13 filename. |

Example:

```yaml
icesat2:
	atl13_fields:
		- ht_ortho
		- segment_lat
		- segment_lon
		- delta_time
		- inland_water_body_id
		- segment_quality
```

If you request `inland_water_body_type`, `inland_water_body_size`, and `inland_water_body_id` in config, the final output columns will appear as `wb_type`, `wb_size`, and `wb_id`.

## Available products
HydroEO currently supports 4 products for lake/reservoir workflows:

| Product key | Source mission | Typical coverage | Key retrievable parameters |
| --- | --- | --- | --- |
| `SWOT_LAKE` | SWOT | 2023-present | water surface elevation, water area/storage-related metrics, quality flags |
| `ATL13` | ICESat-2 | 2018-present | inland water surface elevation, along-track geometry, quality metadata |
| `S3` | Sentinel-3A/B | 2016-present | latitude, longitude, elevation, waveform, backscatter (`sig0`), time, altitude, tracker/range variables, geoid |
| `S6` | Sentinel-6 | 2020-present | latitude, longitude, altitude, `range_ocog`, wet/dry tropospheric corrections, backscatter, time, geoid |

Common across products: geolocation, observation time, water level/elevation-style measurements, and quality metadata.

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
- `HydroEO/system.py`: core system and reservoir handling.
- `HydroEO/satellites/swot.py`: SWOT query/download/subset handling.
- `HydroEO/satellites/icesat2.py`: ICESat-2 query/download/subset handling.
- `HydroEO/satellites/sentinel.py`: Sentinel-3 and Sentinel-6 handling.
- `HydroEO/downloaders/creodias.py`: CREODIAS/CDSE queries and download helpers.
- `HydroEO/downloaders/hydroweb.py`: HydroWeb integration and PLD support.

## More resources
- Notebooks: [notebooks](./notebooks)

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