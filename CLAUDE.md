# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**HydroEO** is a Python package for accessing satellite altimetry EO data (SWOT, ICESat-2, Sentinel-3/6) for water resource applications — download, processing, and analysis over reservoirs, lakes, and rivers.

Four scenarios, each driven by a config in `configs/` (see `configs/*.md` for details):

1. **Reservoirs & lakes** (`reservoirs.yaml`) — multi-satellite download, PLD matching, full timeseries pipeline. Stable.
2. **Rivers** (`rivers.yaml`) — SWOT Hydrocron + multi-satellite download, SWORD v17b matching. Timeseries cleaning **not yet implemented** (preprocessing only). Test in progress.
3. **SWOT raster tiles** (`swot_raster.yaml`) — download SWOT L2 HR/LR SSH rasters, clip/merge to AOI mosaics.
4. **SWOT pixel cloud** (`swot_pixc.yaml`) — download SWOT L2 PIXC, filter by water class, grid to raster.

## Development Setup

Python 3.11+ (pyproject pins `>=3.11`), `uv` + `git`.

```bash
uv sync --all-extras
```

## Common Commands

```bash
make check      # lint + typecheck + test + doctest — run before considering work done
make lint / format / typecheck / test / doctest / coverage / docs / clean
```

```bash
uv run pytest -m unit                          # fast, mocked, no network
uv run pytest -m integration -v                # live APIs, needs credentials
uv run pytest -m api_contract                  # credential-free contract checks
uv run pytest tests/unit/test_flows.py -k pld -v
```

Markers (`unit`, `integration`, `api_contract`, `sliderule`) are declared in `pyproject.toml` (`--strict-markers`). CI runs ruff + `pytest -m unit` on Ubuntu/Windows, Python 3.11/3.13; integration tests are wired but commented out.

Test credentials (set both naming schemes — both are used across tests):

```bash
export EDL_USERNAME=... EDL_PASSWORD=...
export EARTHDATA_USERNAME="$EDL_USERNAME" EARTHDATA_PASSWORD="$EDL_PASSWORD"
export CREODIAS_USERNAME=... CREODIAS_PASSWORD=...
export HYDROWEB_API_KEY=...
```

`tests/conftest.py` and `tests/unit/test_flows.py` hold shared fixtures (mock `Project` via `SimpleNamespace`, GeoDataFrame builders) — check these first when fixture setup is unclear.

## Architecture

### Project lifecycle

```
Project(name, config)
  .report() / .initialize() / .download() / .update()
  .create_timeseries()   # reservoirs only
  .generate_summaries()
```

`initialize()` validates the whole config and collects **all** errors before any I/O — preserve this pattern when adding validation.

### `flows/` package

Split by concern, not one module. `flows/__init__.py` re-exports every public *and* private (`_`-prefixed) name from each submodule, so `flows._some_helper` and `patch.object(flows, "_name")` keep working. Submodules: `_reservoir_init`, `_river_init`, `_reservoir_download`, `_river_download`, `_reservoir_pipeline`, `_river_pipeline`, `_river_common`, `_sentinel_shared`, `_clean_engine`, `_merge_engine`, `_summaries`, `_run_config`, `_constants`. New flow functions go in the matching submodule and must be added to the `__init__.py` re-export list.

### Key modules

`project.py` (`Project` dataclass), `flows/` (see above), `waterbody.py` (`Reservoirs`/`Rivers`), `downloaders/` (`creodias.py`, `dem.py`, `hydroweb.py`), `satellites/{swot,icesat2,sentinel}/`, `utils/` (`general`, `geometry`, `timeseries`, `filters/`), `cli/` (Typer `hydroeo` entry point), `constants.py`, `validation.py`, `plotting.py`, `logging_config.py`.

### Config → outputs

YAML → `Project.__init__()` → validation → `prj.dirs`, `prj.mission_options`, `prj.reservoirs`/`prj.rivers`, `prj.keep_raw_pld`/`prj.keep_raw_sword`.

Everything lives under `{main_dir}`: `aux/PLD/`, `aux/SWORD/` (auxiliary DBs + QA), `raw/{mission}/`, `processed/{mission}/`, `results/{reservoir_id|wb_id|aoi}/` (results namespace is flat — reservoir IDs, river `wb_id`s, and AOI names are peers). Reservoir results include raw/cleaned observations, merged timeseries CSVs, summary PNGs, optional dfs0.

- **PLD** (reservoirs only): HydroWeb.next API → `aux/PLD/PLD_subset.gpkg`; see `hydroweb.download_PLD`, `flows._download_pld/_assign_pld_id/_flag_missing_priors`.
- **SWORD v17b** (rivers only): Zenodo (public) → `aux/SWORD/SWORD_subset.gpkg`; see `flows._ensure_sword_database`, `flows._prepare_rivers_from_sword` (skips DB entirely if `sword_subset_path` is configured).
- **COP-DEM downloader** (`downloaders/dem.py`) is a standalone CLI utility (`hydroeo fetch cop-dem`), not part of the project pipeline — see its docstrings for layer/dataset details.

## Credentials

| Workflow | Env vars |
|---|---|
| Reservoirs — SWOT | `EARTHDATA_USERNAME/PASSWORD`, `EODAG__HYDROWEB_NEXT__AUTH__CREDENTIALS__APIKEY` |
| Reservoirs — Sentinel-3/6 | `CREODIAS_USERNAME/PASSWORD` |
| Reservoirs — ICESat-2, Rivers | none |
| SWOT raster/pixel cloud | `EARTHDATA_USERNAME/PASSWORD` |
| COP-DEM | `CDSE_USERNAME/PASSWORD` |

Env vars take precedence over config-file values. On Windows, set credentials via `os.environ` before calling `project.download()`.

## Style Conventions

- Type hints: `str | None` (3.10+ syntax); docstrings: NumPy style
- Logging: module-level `logger = logging.getLogger(__name__)`; `setup_logging(logging.DEBUG)` for full trace. Logs also go to `logs/HydroEO_*.log` (gitignored)
- Tests: mock `Project` via `SimpleNamespace` + geopandas fixtures; mark every test `unit`/`integration`/`api_contract`/`sliderule`
- Keep `configs/*.md` and `configs/*.yaml` in sync when changing config keys or output structure

## Quick Debugging

- Missing imports → check `__all__` in the relevant package `__init__.py`
- Config not loading → validate YAML syntax; `initialize()` reports all errors together
- Credentials failing → check env vars vs. config (env wins); Windows needs `os.environ[...]` before `download()`
- Tests failing → check fixtures in `tests/conftest.py` / `tests/unit/test_flows.py`; run `-v --tb=short`

## Known Limitations

- Rivers timeseries filtering/cleaning not yet implemented
- SWOT raster reprojection only supports UTM zones
- Sentinel API: only CREODIAS/CDSE supported
- Windows: Spatialite may need manual PATH configuration for SWOT
