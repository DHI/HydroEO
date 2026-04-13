![](images/HydroEO.png)
# Easy access Earth observation data for water resource applications

Repo to allow users with little EO (Earth Observation) knowledge to access and download altimetry over reservoirs and lakes for integration into larger water resource projects.

> [!CAUTION]
> HydroEO has had limited testing and further developments are likely to come. Please report any bugs or issues here: https://github.com/DHI/HydroEO/issues

## Installation
> Note: Ensure that `pip` and `git` are installed on your system before running the following command.
```sh
pip install git+https://github.com/DHI/HydroEO.git
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

For a full setup walkthrough, see [docs/getting_started.md](./docs/getting_started.md).

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
Install development dependencies first:

```sh
pip install -e ".[dev]"
```

Run test suites:

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
- API reference source page: [docs/reference.md](./docs/reference.md)