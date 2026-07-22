# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**HydroEO** is a Python package providing easy access to Earth observation (satellite altimetry) data for water resource applications. It automates download, processing, and analysis of altimetry observations from multiple satellite missions (SWOT, ICESat-2, Sentinel-3, Sentinel-6) over reservoirs, lakes, and rivers.

Four use cases, each driven by a scenario config in `configs/`:

1. **Reservoirs & lakes** (`configs/reservoirs.yaml`) — Multi-satellite download (SWOT Lake SP, ICESat-2 ATL13, Sentinel-3/6), Prior Lake Database (PLD) matching, full timeseries extraction/cleaning/multi-mission merge, diagnostic plots, optional dfs0 export. Status: stable.
2. **Rivers** (`configs/rivers.yaml`) — SWOT Hydrocron API download (public, no credentials), plus multi-satellite download, SWORD v17b database integration for node/reach matching. Timeseries filtering/cleaning is **not yet implemented** (preprocessing only). Status: test in progress.
3. **SWOT raster tiles** (`configs/swot_raster.yaml`) — Download SWOT L2 HR/LR SSH raster products, clip to arbitrary AOI, merge tiles by date/variable → georeferenced netCDF/GeoTIFF mosaics.
4. **SWOT pixel cloud** (`configs/swot_pixc.yaml`) — Download SWOT L2 PIXC point-cloud data, filter by water class, grid to regular rasters (median/mean/max/min).

Per-scenario docs live alongside each config: `configs/reservoirs.md`, `configs/rivers.md`, `configs/swot_raster.md`, `configs/swot_pixc.md`.

## Development Setup

Prerequisites: `uv` and `git`. Python 3.10–3.12 (pyproject currently pins `requires-python = ">=3.11"`).

```bash
git clone https://github.com/DHI/HydroEO.git
cd HydroEO
uv sync --all-extras   # installs dev + test extras
```

## Common Commands

```bash
make check      # lint + typecheck + test + doctest — run before considering work done
make lint       # ruff check HydroEO/
make format     # ruff format HydroEO/
make typecheck  # mypy HydroEO/ --config-file pyproject.toml
make test       # uv run pytest --disable-warnings
make doctest    # pytest --doctest-modules HydroEO
make coverage   # HTML coverage report (htmlcov/)
make docs       # mkdocs build
make clean      # remove dist/, htmlcov/, .coverage, site/
```

### Running tests directly

```bash
uv run pytest -m unit                                        # fast, fully mocked, no network
uv run pytest -m integration -v                               # live API calls, needs credentials
uv run pytest -m api_contract                                 # credential-free endpoint/schema checks
uv run pytest -m 'not sliderule'                               # skip SlideRule API tests if it's down
uv run pytest tests/unit/test_flows.py -v                      # one file
uv run pytest tests/unit/test_flows.py::test_download_pld_downloads_when_missing -v  # one test
uv run pytest tests/unit/test_flows.py -k pld -v               # keyword match
```

Markers are declared in `pyproject.toml` (`--strict-markers` enforced): `unit`, `integration`, `api_contract`, `sliderule`.

Test structure:
- `tests/unit/` — isolated, mocked, no network
- `tests/test_integration.py` — cross-module scenarios against live services
- `tests/test_e2e_run_pipeline.py` — full pipeline smoke tests
- `tests/test_api_contracts.py`, `tests/test_basic.py` — top-level contract/smoke checks
- `tests/conftest.py` and `tests/unit/test_flows.py` hold the shared fixtures (mock `Project` objects, GeoDataFrame builders) — check these first when a test's fixture setup is unclear

Test credentials (both naming schemes are used across tests — set all of them):

```bash
export EDL_USERNAME=... EDL_PASSWORD=...
export EARTHDATA_USERNAME="$EDL_USERNAME"
export EARTHDATA_PASSWORD="$EDL_PASSWORD"
export CREODIAS_USERNAME=... CREODIAS_PASSWORD=...
export HYDROWEB_API_KEY=...
```

CI (`.github/workflows/ci.yml`) runs ruff lint, then `pytest -m unit` on Ubuntu + Windows across Python 3.11/3.13. Integration tests are wired up but currently commented out in the workflow.

## Architecture

### Project lifecycle

```
Project(name, config)
  .report()             # log project + water bodies, return gdf.head()
  .initialize()         # validate config, prep dirs, download aux DBs (PLD/SWORD); collects ALL config errors before any I/O
  .download()           # fetch raw satellite data for the configured date range
  .update()             # extend existing downloads from latest observation to today; safe to re-run
  .create_timeseries()  # extract, clean, filter, merge — reservoirs only
  .generate_summaries() # diagnostic plots per water body/mission
```

`flows.py` is a **package** (`HydroEO/flows/`), split by concern, not a single module — despite what older docs may say. `HydroEO/flows/__init__.py` re-exports every public *and* private name (including `_`-prefixed helpers) from the submodules below, specifically so `from HydroEO import flows; flows._some_private_helper` and `patch.object(flows, "_name")` in tests keep working as if it were still one file:

- `_reservoir_init.py`, `_river_init.py` — PLD/SWORD setup and matching
- `_reservoir_download.py`, `_river_download.py` — per-mission download orchestration
- `_reservoir_pipeline.py`, `_river_pipeline.py` — end-to-end flow entry points
- `_river_common.py` — shared river/SWORD geometry helpers
- `_sentinel_shared.py` — Sentinel-3/6 download logic shared by both reservoirs and rivers
- `_clean_engine.py`, `_merge_engine.py` — timeseries cleaning/filtering and multi-mission merge
- `_summaries.py` — plot generation
- `_run_config.py` — CLI/run-level config plumbing
- `_constants.py` — flow-internal constants

When adding a new flow function, put it in the submodule matching its concern and add it to the `__init__.py` re-export list (both public and any private helpers a test might need to patch).

### Key top-level modules

| File | Purpose |
|------|---------|
| `HydroEO/project.py` | `Project` dataclass — config parsing and lifecycle orchestration |
| `HydroEO/flows/` | Standalone flow functions per lifecycle step (see above) |
| `HydroEO/waterbody.py` | `Reservoirs` / `Rivers` dataclasses — GeoDataFrame containers + config |
| `HydroEO/downloaders/` | Mission-agnostic downloaders: `creodias.py`, `dem.py`, `hydroweb.py` |
| `HydroEO/satellites/{swot,icesat2,sentinel}/` | Mission-specific download/preprocess/filter/extract logic |
| `HydroEO/utils/` | `general.py`, `geometry.py`, `timeseries.py`, `filters/` (elevation, MAD, hampel, rolling_median, Kalman merge) |
| `HydroEO/cli/` | `hydroeo` CLI entry point (Typer app) |
| `HydroEO/constants.py` | `MISSION_DEFAULTS`, field mappings, supported missions |
| `HydroEO/validation.py` | Config validation, aggregated error reporting |
| `HydroEO/plotting.py` | Diagnostic plots |
| `HydroEO/logging_config.py` | Logger setup (`setup_logging`) |

### Config system

YAML config → `Project.__init__()` → validation → `self.dirs`, `self.mission_options`, etc.

Key `Project` attributes:
- `prj.dirs` — output paths: `main`, `pld`, `sword`, `sword_subset`, `swot`, `icesat2`, `icesat2_processed`, `sentinel3`, `sentinel6`, `output`
- `prj.mission_options` — per-satellite settings (download/process flags, date overrides, filters, field lists)
- `prj.reservoirs` / `prj.rivers` — `Reservoirs`/`Rivers` instance (GeoDataFrame + config)
- `prj.keep_raw_pld`, `prj.keep_raw_sword` — booleans controlling raw-file cleanup after aux DB prep (default `False`)

Config validation happens in `project.initialize()` **before any I/O**: required sections (`project` + at least one of `reservoirs`/`rivers`/`swot_raster`/`swot_pixc`), credential fallback chain (config → env var → error), date-range tuple checks, mission-option type/value checks. All issues are collected and reported together rather than failing on the first one — preserve this pattern when adding new validation.

### PLD (Prior Lake Database) — reservoirs only

- Source: HydroWeb.next API (needs API key), SQLite layers merged and subset to AOI → `{main_dir}/aux/PLD/PLD_subset.gpkg`, with `present_in_pld.gpkg`/`missing_in_pld.gpkg` QA outputs alongside.
- `hydroweb.download_PLD(download_dir, bounds, raw_pld_path=None, keep_raw=True)` — downloads or reuses `raw_pld_path`; deletes temp files unless `keep_raw` or the path is outside `main_dir`.
- `flows._download_pld`, `flows._assign_pld_id` (spatial join → `prior_lake_id`), `flows._flag_missing_priors`.

### SWORD v17b — rivers only

- Source: Zenodo (public, no credentials) — per-continent GPKGs → `{main_dir}/aux/SWORD/gpkg/`; AOI subset → `{main_dir}/aux/SWORD/SWORD_subset.gpkg`.
- `flows._ensure_sword_database(prj)` resolves availability in priority order: full DB already present → user-provided zip (extract, honor `keep_raw_sword`) → user-provided directory (used as-is) → auto-download from Zenodo.
- `flows._prepare_rivers_from_sword(prj)`: if `sword_subset_path` is configured, reads it directly (skips DB + subsetting entirely — fastest path); otherwise calls `_ensure_sword_database` then spatially intersects with the AOI. Returns `target_features`, `target_id_col`, `target_ids` for Hydrocron queries.
- `rivers.continent_key` becomes optional once `sword_subset_path` is provided.

### Output layout

Everything under `{main_dir}` (from `project.main_dir` config):

- `aux/PLD/`, `aux/SWORD/` — auxiliary databases + QA/QC
- `raw/` — downloaded data per mission (`raw/swot/`, `raw/swot/{wb_id}/` for river Hydrocron CSVs, `raw/sentinel3/{id}/`, `raw/sentinel6/{id}/`, `raw/swot_raster/{aoi}/{product}/`, `raw/swot_pixc/{aoi}/{product}/`); `raw/icesat2/` is defined but stays empty since SlideRule writes straight to `processed/`
- `processed/` — `processed/icesat2/{id}/atl13.parquet`, `processed/swot_raster/{aoi}/{product}/` (per-variable GeoTIFFs), `processed/swot_pixc/{aoi}/` (trimmed GeoJSON)
- `results/` — flat namespace where reservoir IDs, river `wb_id`s, and AOI names are peers: `results/{reservoir_id}/` (raw/cleaned observations, merged timeseries CSVs, summary PNGs, optional dfs0), `results/{wb_id}/` (river diagnostic plots), `results/{aoi}/` (merged SWOT raster/pixc rasters)

Per-reservoir output detail (`results/{reservoir_id}/`):
```
raw_observations/           # raw extracted observations, SHP per mission
cleaned_observations/       # cleaned timeseries, CSV per mission
merged_progress/            # intermediate merged files
merged_timeseries.csv       # all missions merged (wide format)
all_cleaned_timeseries.csv  # long format with metadata
crossing_summary.png
cleaning_summary.png
merging_summary.png
[dfs0 files]                # if export_to_dfs0: true
```

### COP-DEM downloader (`HydroEO/downloaders/dem.py`)

Standalone CLI-only utility, **not** part of the project pipeline (`hydroeo fetch cop-dem`).

- Layers: `DEM30` (default), `DEM90`, `EDM`, `FLM`, `HEM`, `WBM` — comma-separated via `--dataset`.
- ZIPs download once; all requested layers extract in a single pass; each layer writes an independent `{output_basename}_{LAYER}.tif`.
- `DEM30` and `DEM90` cannot be combined (different CDSE datasets — raises `ValueError`). Aux-only requests default to GLO-30; if `DEM90` is requested, GLO-90 is used for all layers in that call.
- Key API: `VALID_LAYERS`, `LAYER_SUFFIXES`, `DEFAULT_LAYERS`, `_validate_layers`, `_resolve_cdse_dataset`, `download_cop_dem(minx, miny, maxx, maxy, output_dir, username, password, layers=["DEM30"], output_basename="cop_dem_merged") -> dict[str, Path]`.

## Credentials

| Workflow | Credentials | Env vars |
|---|---|---|
| Reservoirs — SWOT | Earthdata + HydroWeb API key | `EARTHDATA_USERNAME`, `EARTHDATA_PASSWORD`, `EODAG__HYDROWEB_NEXT__AUTH__CREDENTIALS__APIKEY` |
| Reservoirs — Sentinel-3/6 | CREODIAS | `CREODIAS_USERNAME`, `CREODIAS_PASSWORD` |
| Reservoirs — ICESat-2 | none | — |
| Rivers (Hydrocron) | none | — |
| SWOT raster/pixel cloud | Earthdata | `EARTHDATA_USERNAME`, `EARTHDATA_PASSWORD` |
| COP-DEM | CDSE (free, needs CCM data access enabled) | `CDSE_USERNAME`, `CDSE_PASSWORD` |

Env vars take precedence over config-file values. On Windows, credentials may need to be set via `os.environ` before calling `project.download()`.

## Style Conventions

- Type hints: Python 3.10+ syntax (`str | None`, not `Optional[str]`)
- Docstrings: NumPy style (Parameters/Returns/Raises)
- Logging: module-level `logger = logging.getLogger(__name__)`; `setup_logging(logging.DEBUG)` for full debug output. Logs also go to `logs/HydroEO_YYYY-MM-DD_HH-MM-SS.log` (gitignored, auto-created)
- Tests: mock `Project` via `SimpleNamespace` + geopandas fixtures (see `tests/conftest.py`, `tests/unit/test_flows.py`); mark every test `@pytest.mark.unit`/`integration`/`api_contract`/`sliderule`
- Keep `configs/*.md` and `configs/*.yaml` in sync with actual defaults when changing config keys or output structure

## Quick Debugging

- **Missing imports:** check `__all__` in the relevant package `__init__.py`
- **Config not loading:** validate YAML syntax; `project.initialize()` collects and reports all config errors together rather than failing fast
- **Credentials failing:** check env vars vs. config values (env wins); on Windows set via `os.environ[...]` before `project.download()`
- **Tests failing:** check fixture setup in `tests/conftest.py` / `tests/unit/test_flows.py`; run with `-v --tb=short`
- **Need full trace:** `setup_logging(logging.DEBUG)` for console + file debug output

## Known Limitations

- Rivers timeseries filtering/cleaning not yet implemented (preprocessing only)
- SWOT raster reprojection only supports UTM zones
- Sentinel API: only CREODIAS/CDSE supported (ESA Copernicus Hub discontinued)
- Windows: Spatialite may need manual PATH configuration for SWOT
