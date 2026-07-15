# HydroEO — AI Agent Guidelines

This document provides codebase knowledge for AI coding agents (Claude, Copilot, etc.) working on HydroEO.

## Project Overview

**HydroEO** is a Python package that provides easy access to Earth observation (satellite altimetry) data for water resource applications. It automates the download, processing, and analysis of altimetry observations from multiple satellite missions (SWOT, ICESat-2, Sentinel-3, Sentinel-6) over reservoirs, lakes, and rivers.

### Four Use Cases

1. **Reservoirs & Lakes** (`reservoirs` branch)
   - Multi-satellite data download (SWOT Lake SP, ICESat-2 ATL13, Sentinel-3, Sentinel-6)
   - Prior Lake Database (PLD) matching for reservoir refinement
   - Full timeseries extraction, cleaning, filtering, and multi-mission merge
   - Output: per-reservoir timeseries, diagnostic plots, optional dfs0 export

2. **Rivers** (`rivers` branch)
   - SWOT Hydrocron API download (public, no credentials needed)
   - Multi-satellite data download (SWOT Lake SP, ICESat-2 ATL13, Sentinel-3, Sentinel-6)
   - SWORD v17b database integration for river node/reach matching
   - Full timeseries extraction, cleaning, filtering, and multi-mission merge aligned to SWORD reaches or nodes
   - Output: per-river Hydrocron timeseries (CSV)

3. **SWOT Raster Tiles** (`swot_raster` branch)
   - Download SWOT L2 HR/LR SSH raster products
   - Clip to arbitrary AOI, merge tiles by date/variable
   - Output: georeferenced netCDF or GeoTIFF mosaics

4. **SWOT Pixel Cloud** (`swot_pixc` branch)
   - Download SWOT L2 PIXC (point cloud) data
   - Filter by water class, grid to regular rasters
   - Output: binned rasters (median/mean/max/min statistics)

---

## Development Setup

### Installation

```bash
# Prerequisites: uv and git installed
git clone https://github.com/DHI/HydroEO.git
cd HydroEO
uv sync  # Creates isolated venv and installs dependencies
```

### Python versions
HydroEO runs on Python 3.10–3.12.

### Common Commands (from Makefile)

```bash
make test          # Run pytest with disabled warnings
make lint          # Run ruff linter
make format        # Auto-format code with ruff
make typecheck     # mypy type checking
make check         # Run lint + typecheck + test (full validation)
make coverage      # Generate HTML coverage report
make docs          # Build MkDocs site (mkdocs.yml)
make clean         # Remove build/coverage artifacts
```

---

## Codebase Architecture

### Project Lifecycle

```
Project(config_file)
  └─ project.initialize()      # Load config, validate, setup directories, match PLD
  └─ project.download()         # Fetch raw satellite data
  └─ project.create_timeseries()  # Extract & clean observations per waterbody
  └─ project.generate_summaries() # Plots and diagnostics
```

### Key Files & Modules

| File | Purpose |
|------|---------|
| `HydroEO/project.py` | Main `Project` dataclass; config parsing and orchestration |
| `HydroEO/flows.py` | Standalone flow functions for each lifecycle step (initialize, download, timeseries, etc.) |
| `HydroEO/waterbody.py` | `Reservoirs` and `Rivers` dataclasses; GeoDataFrame containers and config |
| `HydroEO/downloaders/` | `creodias.py`, `dem.py`, `hydroweb.py` — mission-agnostic downloaders |
| `HydroEO/satellites/{swot, icesat2, sentinel}` | Mission-specific logic: download, preprocess, filter, extract |
| `HydroEO/utils/` | `general.py`, `geometry.py`, `timeseries.py`, `filters/` — utilities |
| `HydroEO/cli/` | CLI entry point (`hydroeo` command) |
| `HydroEO/constants.py` | MISSION_DEFAULTS, field mappings, supported missions |
| `HydroEO/validation.py` | Config validation, error reporting |

### Config System

**Flow:**  
YAML config file → `Project.__init__()` → validation → `self.dirs`, `self.mission_options`, `self.keep_raw_pld`, `self.keep_raw_sword`, etc.

**Key attributes:**
- `prj.dirs` — dictionary of output paths (main, pld, sword, sword_subset, swot, icesat2, sentinel3, sentinel6, output)
- `prj.mission_options` — per-satellite settings (download/process flags, date overrides, filters, field lists)
- `prj.reservoirs` / `prj.rivers` — waterbody `Reservoirs`/`Rivers` instances with GeoDataFrame + config
- `prj.keep_raw_pld` — boolean flag for PLD raw file cleanup (default False)
- `prj.keep_raw_sword` — boolean flag for SWORD raw file cleanup (default False)

### Flows Orchestration

All pipeline functions are standalone in `flows.py` and operate on the `Project` object:

```python
# Reservoirs workflow example
from HydroEO.project import Project
from HydroEO import flows

prj = Project(name="My Project", config="config.yaml")
flows.initialize_reservoirs(prj)  # PLD download + matching
flows.download_reservoirs(prj)     # Multi-satellite downloads
flows.create_reservoirs_timeseries(prj)  # Extraction + cleaning
flows.generate_reservoirs_summaries(prj)  # Plots
```

---

## PLD (Prior Lake Database)

### Download & Subset

- **Source:** HydroWeb.next API (requires API key)
- **Format:** SQLite files (from `SWOT_PRIOR_LAKE_DATABASE` collection)
- **Output:** `{main_dir}/aux/PLD/PLD_subset.gpkg` (geopackage)
- **QA/QC:** `present_in_pld.gpkg`, `missing_in_pld.gpkg` in same folder

### Configuration

```yaml
hydroweb:
  api_key: "your_api_key"          # Required; or set HYDROWEB_API_KEY env var
  raw_pld_path: "/path/to/PLD.zip" # Optional; if provided, skips download
  keep_raw_pld: false              # Delete raw zip after subset (default)
```

### Key Functions

- `hydroweb.download_PLD(download_dir, bounds, raw_pld_path=None, keep_raw=True)`
  - Downloads PLD from API (or reuses `raw_pld_path` if provided)
  - Merges SQLite layers, filters to bounds, exports GPKG
  - Cleanup: deletes temp folder + zip (unless `keep_raw=True` or `raw_pld_path` is outside `main_dir`)

- `flows._download_pld(prj)` — calls `download_PLD()` with project context
- `flows._assign_pld_id(prj)` — spatial join reservoirs to PLD lakes; sets `prior_lake_id` column
- `flows._flag_missing_priors(prj)` — exports matched/unmatched reservoirs to `aux/PLD/`

---

## SWORD v17b Database

### Download, Subsetting & Storage

- **Source:** Zenodo public repository (no credentials needed; [v17b dataset](https://zenodo.org/records/15299138))
- **Format:** GeoPackage (GPKG) files per continent/feature type (e.g., `eu_sword_nodes_v17b.gpkg`)
- **Storage:** `{main_dir}/aux/SWORD/gpkg/` (full database, if downloaded)
- **Subset Output:** `{main_dir}/aux/SWORD/SWORD_subset.gpkg` (spatial intersection with AOI)

### Configuration

```yaml
sword_db:
  raw_sword_path: "/path/to/SWORD_v17b_gpkg.zip"  # Optional; if provided, skip download
  sword_subset_path: "/path/to/SWORD_subset.gpkg" # Optional; if provided, skip download + subsetting
  keep_raw_sword: false                           # Delete downloaded zip after extraction (default)
```

### Key Functions

- `flows._ensure_sword_database(prj)` — orchestrates SWORD database availability
  - Checks if full SWORD (GPKGs) already in `prj.dirs["sword"]` → skips download
  - **User-provided zip:** extracts to `aux/SWORD/`; honours `keep_raw_sword`
  - **User-provided directory:** uses it directly (no extraction or cleanup)
  - **Auto-download:** Zenodo zip → `aux/SWORD/SWORD_v17b_gpkg.zip` → extract to `aux/SWORD/gpkg/`; honours `keep_raw_sword`

- `flows._prepare_rivers_from_sword(prj)` — prepares SWORD targets for rivers workflow
  - Gate 1: if `sword_subset_path` exists → reads it, skips database + subsetting
  - Gate 2: otherwise → calls `_ensure_sword_database()`, performs spatial intersection with AOI, saves to `SWORD_subset.gpkg`
  - Returns: `prj.rivers.target_features`, `target_id_col`, `target_ids` (for Hydrocron queries)

### Configuration Notes

- `sword_db` section is optional and only validated when `rivers` section is also present in `aoi_path` mode
- `rivers.continent_key` becomes optional when `sword_subset_path` is provided (subset already contains selected features)
- When `sword_subset_path` is provided, SWORD database download is skipped entirely (fastest path)

---

## Testing

### Test Structure

- **Unit tests:** `tests/unit/` — isolated, use mocks and fixtures
- **Integration tests:** `tests/test_integration.py` — cross-module scenarios
- **E2E tests:** `tests/test_e2e_run_pipeline.py` — full pipeline smoke tests
- **Fixtures:** `tests/conftest.py`, `tests/unit/test_flows.py` (project mocks, GeoDataFrames)

### Test Conventions

```python
import pytest
from types import SimpleNamespace
import geopandas as gpd

# Fixture pattern for mocking Project
@pytest.fixture
def mock_project_reservoirs(tmp_path):
    prj = SimpleNamespace()
    prj.dirs = {"main": str(tmp_path), "pld": str(tmp_path / "aux" / "PLD" / "PLD_subset.gpkg"), ...}
    prj.reservoirs = Reservoirs(gdf=gpd.GeoDataFrame(...), id_key="id", dirs=prj.dirs)
    return prj

# Test discovery: @pytest.mark.unit, @pytest.mark.integration, or @pytest.mark.e2e
@pytest.mark.unit
def test_my_function():
    ...
```

### Run Tests (with uv)

```bash
# All tests
uv run pytest

# Unit tests only
uv run pytest tests/unit/ -v

# By mark
uv run pytest -m unit
uv run pytest -m integration
uv run pytest -m e2e

# Specific file or test
uv run pytest tests/unit/test_flows.py::test_download_pld_downloads_when_missing -v

# PLD-specific tests
uv run pytest tests/unit/test_flows.py -k pld -v

# Coverage report
uv run pytest --cov-report html --cov=HydroEO tests/

# Or use Makefile (which wraps uv commands)
make test
make coverage
make check  # lint + typecheck + test (full validation)
```

---

## Directory Layout

```
HydroEO/
├── __init__.py
├── project.py              # Project dataclass
├── flows.py                # Standalone flow functions
├── waterbody.py            # Reservoirs, Rivers dataclasses
├── constants.py            # Defaults, mission config templates
├── validation.py           # Config validation
├── plotting.py             # Diagnostic plots
├── logging_config.py       # Logger setup
├── cli/                    # CLI entry point (hydroeo command)
├── downloaders/            # Non-mission-specific (creodias, dem, hydroweb)
├── satellites/             # Mission-specific logic
│   ├── swot/
│   ├── icesat2/
│   ├── sentinel/
│   └── {swot,icesat2,sentinel}/__init__.py  # __all__ exports
├── utils/
│   ├── general.py          # OS, file ops
│   ├── geometry.py         # CRS, spatial ops
│   ├── timeseries.py       # Temporal merge, statistics
│   └── filters/            # elevation, MAD, hampel, rolling_median
└── [tests/]                # Tests (sibling at repo root)
```

---

## Output Conventions

### Project Root (`{main_dir}`)

All outputs are placed under `{main_dir}` (from config `project.main_dir`) using a three-tier layout:

- `aux/PLD/` — PLD subset (GPKG) + QA/QC outputs (unchanged)
- `aux/SWORD/` — SWORD v17b subset (GPKG) + full database GPKGs (unchanged)
- `raw/` — downloaded satellite data
  - `raw/swot/` — SWOT Lake SP granule files; `raw/swot/{wb_id}/` — rivers Hydrocron CSVs
  - `raw/icesat2/` — defined but empty (SlideRule writes to `processed/`)
  - `raw/sentinel3/{id}/`, `raw/sentinel6/{id}/` — per-reservoir Sentinel subsets
  - `raw/swot_raster/{aoi}/{product}/` — SWOT raster netCDF granules
  - `raw/swot_pixc/{aoi}/{product}/` — SWOT pixel cloud netCDF granules
- `processed/` — intermediate/modified outputs
  - `processed/icesat2/{id}/atl13.parquet` — ICESat-2 SlideRule output
  - `processed/swot_raster/{aoi}/{product}/` — extracted per-variable GeoTIFFs
  - `processed/swot_pixc/{aoi}/` — trimmed GeoJSON point data
- `results/` — final outputs (flat namespace; reservoir IDs, river wb_ids, and AOI names are peers)
  - `results/{reservoir_id}/` — per-reservoir timeseries, observations, plots
  - `results/{wb_id}/` — per-river diagnostic plots
  - `results/{aoi}/` — merged SWOT raster TIFFs and SWOT pixc gridded rasters

### Per-Reservoir Outputs

```
{main_dir}/results/{reservoir_id}/
├── raw_observations/           # Raw extracted observations (SHP per mission)
├── cleaned_observations/       # Cleaned timeseries (CSV per mission)
├── merged_progress/            # Intermediate merged files
├── merged_timeseries.csv       # All missions merged (wide format)
├── all_cleaned_timeseries.csv  # Long format with metadata
├── crossing_summary.png
├── cleaning_summary.png
├── merging_summary.png
└── [dfs0 files]               # If export_to_dfs0: true
```

### Key `prj.dirs` keys

| Key | Path | Notes |
|-----|------|-------|
| `main` | `{main_dir}` | |
| `pld` | `{main}/aux/PLD/PLD_subset.gpkg` | |
| `sword` | `{main}/aux/SWORD/gpkg` | |
| `sword_subset` | `{main}/aux/SWORD/SWORD_subset.gpkg` | |
| `swot` | `{main}/raw/swot` | |
| `icesat2` | `{main}/raw/icesat2` | defined; empty |
| `icesat2_processed` | `{main}/processed/icesat2` | SlideRule parquet output |
| `sentinel3` | `{main}/raw/sentinel3` | |
| `sentinel6` | `{main}/raw/sentinel6` | |
| `output` | `{main}/results` | both Reservoirs and Rivers |

---

## Configuration Validation

Config validation occurs in `project.initialize()` (before any I/O):

- **Required sections:** `project`, at least one of `reservoirs/rivers/swot_raster/swot_pixc`
- **Credential fallbacks:** Config → environment variables → error
- **Date ranges:** Must be valid `[year, month, day]` tuples
- **Mission options:** Type/value checks on optional filters and field lists
- **Error reporting:** Collect all issues, report together before stopping

---

## Style & Conventions

### Code

- **Type hints:** Use Python 3.10+ syntax (e.g., `str | None` instead of `Optional[str]`)
- **Docstrings:** NumPy style (Parameters, Returns, Raises sections)
- **Logging:** Always use `logger = logging.getLogger(__name__)` at module level
- **Testing:** Use SimpleNamespace for mock projects, geopandas fixtures for test data
- **Imports:** Group std lib, external, local (isort-compatible)

### Documentation

- Update README.md when adding config keys or changing output structure
- Keep config examples in `configs/` in sync with defaults
- Document new CLI subcommands in README's "CLI" section

---

## Known Limitations & TODOs

- **Rivers timeseries filtering:** Not yet implemented (preprocessing only)
- **SWOT raster reprojection:** Only supports UTM zones; custom CRS support pending
- **Sentinel API:** Only CREODIAS/CDSE supported; ESA Copernicus Hub discontinued
- **Windows Spatialite:** May require manual PATH configuration for SWOT on Windows

---

## COP-DEM Downloader (`HydroEO/downloaders/dem.py`)

Standalone CLI-only utility — **not** integrated into the project pipeline.

### Available Layers

| Layer  | Description                | CDSE source     | Archive subfolder |
|--------|----------------------------|-----------------|-------------------|
| DEM30  | 30 m elevation (default)   | GLO-30-DGED     | `DEM/`            |
| DEM90  | 90 m elevation             | GLO-90-DGED     | `DEM/`            |
| EDM    | Editing Mask               | same as DEM res | `AUXFILES/`       |
| FLM    | Filling Mask               | same as DEM res | `AUXFILES/`       |
| HEM    | Height Error Mask          | same as DEM res | `AUXFILES/`       |
| WBM    | Water Body Mask            | same as DEM res | `AUXFILES/`       |

### Design

- `--dataset` accepts a comma-separated list of layer names (e.g. `DEM30,WBM,EDM`).
- ZIPs are downloaded **once**; all requested layers are extracted in a single pass.
- Each layer produces an independent output: `{output_basename}_{LAYER}.tif`.
- `DEM30` and `DEM90` cannot be combined (different CDSE datasets; raise `ValueError`).
- Aux-only requests (no DEM layer) default to the GLO-30 dataset.
- When `DEM90` is in the layer list, the GLO-90 dataset is used for all layers.
- Default layer when none specified: `["DEM30"]`.

### Key Constants / API

```python
VALID_LAYERS: set[str]           # {"DEM30","DEM90","EDM","FLM","HEM","WBM"}
LAYER_SUFFIXES: dict[str, str]   # e.g. {"WBM": "_WBM.tif"}
DEFAULT_LAYERS: list[str]        # ["DEM30"]
_CDSE_GLO30: str                 # "COP-DEM_GLO-30-DGED/2023_1"
_CDSE_GLO90: str                 # "COP-DEM_GLO-90-DGED/2023_1"

_validate_layers(layers)         # raises ValueError on bad input
_resolve_cdse_dataset(layers)    # returns CDSE dataset string

download_cop_dem(
    minx, miny, maxx, maxy,
    output_dir, username, password,
    layers=["DEM30"],            # list[str]
    output_basename="cop_dem_merged",
) -> dict[str, Path]             # {layer: output_path}
```

---

## References

- **Config Schema:** See `configs/` for per-scenario templates and `notebooks/example_config.yaml` for the unified all-in-one reference
- **Constants:** `HydroEO/constants.py` (MISSION_DEFAULTS, field mappings)
- **Validation Rules:** `HydroEO/validation.py` (config checks)
- **Build/Packaging:** `pyproject.toml` (dependencies, Python version)
- **Tests:** `tests/conftest.py`, `tests/unit/`, `tests/test_*.py`

---

## Quick Debugging

- **Missing imports:** Check `__all__` in package `__init__.py` files
- **Config not loading:** Validate YAML syntax and run `project.initialize()` which collects all errors
- **Credentials failing:** Check env vars or config values; use `os.environ["EARTHDATA_USERNAME"]` setters if needed
- **Tests failing:** Check fixture setup in `conftest.py` and `test_flows.py`; use `-v` and `--tb=short` flags
- **Logging:** Enable `setup_logging(logging.DEBUG)` for full debug output to console

