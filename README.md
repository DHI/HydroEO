![](images/HydroEO.png)
# Easy access Earth observation data for water resource applications

Python package for downloading, processing, and analysing satellite altimetry data over reservoirs, lakes, and rivers.

> [!CAUTION]
> HydroEO has had limited testing and further developments are likely to come. Please report any bugs or issues here: https://github.com/DHI/HydroEO/issues

## Installation

> Note: Ensure that `uv` and `git` are installed. For `uv`, follow [this link](https://docs.astral.sh/uv/getting-started/installation/) or on Windows PowerShell: `powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"`.

Then install HydroEO:

```sh
git clone https://github.com/DHI/HydroEO.git
cd HydroEO
uv sync
```

Python 3.10–3.12.

## Quick start

Pick a scenario config from [`configs/`](./configs/), fill in your paths and credentials, then:

```python
from HydroEO.project import Project

project = Project(name="my_project", config="configs/reservoirs.yaml")
project.initialize()
project.download()
project.create_timeseries()
project.generate_summaries()
```

Or via CLI: `hydroeo run configs/reservoirs.yaml`

## Use cases

| Branch | Config | Status | Satellites | Description |
| --- | --- | --- | --- | --- |
| `reservoirs` | [configs/reservoirs.yaml](configs/reservoirs.yaml) | ✅ | SWOT Lake SP, ICESat-2 ATL13, Sentinel-3/6 | Lakes and reservoirs from polygon input. Multi-mission timeseries with cleaning filters. → [Full docs](configs/reservoirs.md) |
| `rivers` | [configs/rivers.yaml](configs/rivers.yaml) | 🧪 partial | SWOT Hydrocron (public API) | River nodes/reaches from SWORD v17b. Download and diagnostic plots. → [Full docs](configs/rivers.md) |
| `swot_raster` | [configs/swot_raster.yaml](configs/swot_raster.yaml) | ✅ | SWOT L2 HR/LR Raster | Arbitrary AOI. Downloads, clips, and merges raster tiles by date. → [Full docs](configs/swot_raster.md) |
| `swot_pixc` | [configs/swot_pixc.yaml](configs/swot_pixc.yaml) | ✅ | SWOT L2 PIXC | Arbitrary AOI. Point cloud gridded to rasters via binned statistics. → [Full docs](configs/swot_pixc.md) |

## Project lifecycle

```
report() → initialize() → download() → update() → create_timeseries() → generate_summaries()
```

| Method | What it does |
| --- | --- |
| `project.report()` | Logs project name and loaded water bodies; returns `gdf.head()`. |
| `project.initialize()` | Validates config, prepares directories, downloads auxiliary databases (PLD, SWORD). Reports all config issues before doing any I/O. |
| `project.download()` | Downloads raw satellite data within the configured date range. |
| `project.update()` | Extends existing downloads from the latest observation to today. Safe to run repeatedly. |
| `project.create_timeseries()` | Extracts spatially matched observations, applies cleaning filters, writes per-body timeseries. *(Reservoirs only.)* |
| `project.generate_summaries()` | Produces diagnostic plots per water body and mission. |

## CLI

### Config-driven pipeline

```sh
hydroeo run        configs/reservoirs.yaml          # initialize + download
hydroeo initialize configs/reservoirs.yaml
hydroeo download   configs/reservoirs.yaml
hydroeo update     configs/reservoirs.yaml          # extend to today
hydroeo timeseries configs/reservoirs.yaml          # reservoirs only
hydroeo summaries  configs/reservoirs.yaml
hydroeo run        configs/reservoirs.yaml --name "Lake Geneva" --verbose
```

### Direct satellite downloads (no config file)

Download data using a bounding box and date range without a config file. See the per-scenario docs for full command examples:

- **SWOT raster** — see [configs/swot_raster.md](configs/swot_raster.md)
- **SWOT pixel cloud** — see [configs/swot_pixc.md](configs/swot_pixc.md)
- **ICESat-2** — `hydroeo fetch icesat2 --bbox "..." --start 2024-01-01 --end 2024-12-31` (no credentials)
- **Sentinel-3/6** — `hydroeo fetch sentinel --product S3 --bbox "..." --creodias-username <u> --creodias-password <p>`
- **COP-DEM** — `hydroeo fetch cop-dem --bbox "..." --cdse-username <u> --cdse-password <p>`

Use `--verbose` / `-v` for debug logging on any command.

## Credentials

| Workflow | Credentials | Environment variables |
| --- | --- | --- |
| Reservoirs — SWOT | Earthdata + HydroWeb API key | `EARTHDATA_USERNAME`, `EARTHDATA_PASSWORD`, `EODAG__HYDROWEB_NEXT__AUTH__CREDENTIALS__APIKEY` |
| Reservoirs — Sentinel-3/6 | CREODIAS | `CREODIAS_USERNAME`, `CREODIAS_PASSWORD` |
| Reservoirs — ICESat-2 | None | — |
| Rivers (Hydrocron) | None | — |
| SWOT Raster / Pixel Cloud | Earthdata | `EARTHDATA_USERNAME`, `EARTHDATA_PASSWORD` |
| COP-DEM | CDSE (free) | `CDSE_USERNAME`, `CDSE_PASSWORD` |

Credentials can be set in the config file or as environment variables; environment variables take precedence. On Windows, credentials may need to be set via `os.environ` before calling `project.download()`.

## Logging

HydroEO logs to console (INFO) and file (DEBUG). Log files are written to `logs/` (auto-created, gitignored) as `HydroEO_YYYY-MM-DD_HH-MM-SS.log`.

```python
from HydroEO.logging_config import setup_logging
import logging

setup_logging(logging.DEBUG)                            # both console and file at DEBUG
setup_logging(logging.INFO, enable_file_logging=False)  # console only
```

## Running tests

```sh
make test
```

See [tests/README.md](tests/README.md) for the full test suite, pytest markers, and CI/CD setup.
