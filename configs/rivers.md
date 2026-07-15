# Rivers — Configuration Reference

SWOT Hydrocron download for river nodes and reaches, matched to the SWORD v17b river network database, with optional ICESat-2 and Sentinel-3/6 as additional missions on the same targets. **SWOT itself needs no credentials** — Hydrocron is a public API — but ICESat-2/Sentinel-3/6 need the same credentials as the reservoirs workflow if enabled (see the root [README](../README.md)'s Credentials table).

> **Status:** Full pipeline supported — initialization, multi-mission download, extraction, cleaning, merge, and diagnostic plotting all work for rivers, sharing the same underlying engine reservoirs use. One gap: no `dfs0` export for rivers yet (reservoirs-only). See [Known limitations](#known-limitations) below.

## Quick start

At least one mission (`swot`, `icesat2`, `sentinel3`, or `sentinel6`) must be enabled (`download: true`) or nothing downloads — `download_rivers()` warns and returns early if none are, even with a valid `aoi_path`/SWORD setup (SWORD itself will still have been prepared by `initialize()`, so you may see SWORD files but no timeseries data if you hit this). SWOT via Hydrocron is the natural default/primary source for node/reach WSE, but it is not the *only* mission that can stand alone — ICESat-2/Sentinel-3/6 alone (without SWOT) will also download and process normally.

```python
from HydroEO.project import Project

project = Project(name="my_river_project", config="configs/rivers.yaml")
project.initialize()          # SWORD download + AOI subsetting
project.download()            # whichever of swot/icesat2/sentinel3/sentinel6 are enabled
project.create_timeseries()   # extraction + cleaning + merge
project.generate_summaries()  # plots
```

## Config reference

Start from [`configs/rivers.yaml`](rivers.yaml).

### `rivers`

| Key | Type | Required | Description |
| --- | --- | --- | --- |
| `enabled` | bool | ✅ | Must be `true` |
| **Option A — AOI file** | | | |
| `aoi_path` | string | ✅ (A) | Path to AOI `.shp` or `.gpkg` |
| `continent_key` | string | ✅ (A) | SWORD continent code (see table below) |
| `feature_type` | string | ✅ | `nodes` or `reaches` |
| `id_key` | string | ✅ (A) | AOI column used to name per-river output folders |
| `buffer_meters` | float | — | Optional AOI buffer before SWORD spatial subsetting. Also used as the extraction-corridor fallback if `extraction_buffer_meters` isn't set (see below). |
| **Option B — explicit IDs** | | | |
| `feature_numbers` | list | ✅ (B) | List of node or reach IDs |
| `feature_type` | string | ✅ | Must match IDs: `nodes` or `reaches` |
| `id` | string | ✅ (B) | Used as the output folder name |
| **ICESat-2/Sentinel-3/6 extraction (only relevant if those missions are enabled)** | | | |
| `extraction_buffer_meters` | float | — | Corridor half-width (m) around each SWORD target used to decide which raw ICESat-2/Sentinel points are plausibly real river water. Separate from `buffer_meters` above, which only decides which *targets* are in scope. If omitted: falls back to `buffer_meters`, then a flat `500.0` m default. |
| `width_buffer_factor` | float | — | Multiplier applied to SWORD's own `width` column when sizing the corridor per-target instead of using one flat distance (`distance = width/2 × width_buffer_factor`). Default `1.05`. Only used when SWORD `width` data is available; falls back to `extraction_buffer_meters` otherwise. |
| `max_node_assignment_meters` | float | — | Max distance (m) for assigning a raw ICESat-2/Sentinel point to its nearest node/reach. If omitted: falls back to the same resolved `extraction_buffer_meters` value above. |
| `overwrite_extraction` | bool | — | Default `false` (skip re-extracting a target whose output already exists, or — for ICESat-2 — whose downloaded data hasn't changed since the last extraction). `true` forces a full re-extraction from scratch. |
| `merging_options` | dict | — | Per-project override of any key in `flows.DEFAULT_RIVER_MERGING_OPTIONS` (Kalman/SVR/bias-correction/reach-slope-correction parameters). **Currently a direct copy of the reservoir-tuned defaults, not independently validated against real river data** — check kept/rejected observation counts on your own rivers before trusting them as-is. **If you use this pipeline, please also cite Schwatke et al. (2015)** — see the root [README](../README.md)'s Citation section. |

**SWORD continent codes:**

| Code | Region |
| --- | --- |
| `af` | Africa |
| `as` | Asia |
| `eu` | Europe |
| `na` | North America |
| `oc` | Oceania |
| `sa` | South America |

### `sword_db`

Optional shortcuts to avoid re-downloading SWORD on subsequent runs:

| Key | Default | Description |
| --- | --- | --- |
| `raw_sword_path` | — | Path to existing SWORD v17b zip or extracted folder; skips download |
| `sword_subset_path` | — | Path to a pre-made `SWORD_subset.gpkg`; skips download **and** spatial subsetting entirely |
| `keep_raw_sword` | `false` | Retain raw SWORD zip after extraction |

When `sword_subset_path` is provided:
- `rivers.continent_key` becomes optional (the subset already contains the right features).
- The full SWORD database download is skipped entirely — fastest option for repeated runs.

The SWORD subset is always written to `{main_dir}/aux/SWORD/SWORD_subset.gpkg` for future reuse.

### Additional missions: ICESat-2, Sentinel-3, Sentinel-6

Enabled the same way as for reservoirs (`download`/`process` flags per mission), clustered onto the same SWORD targets as SWOT via the corridor parameters above:

```yaml
icesat2:
  download: true
  process: true

sentinel3:
  download: true
  process: true
  sigma0_min: 0.0   # see note below

sentinel6:
  download: true
  process: true
  source: "earthdata"   # CREODIAS only distributes the Low Rate product;
                        # High Rate (20Hz Ku-band) needs "earthdata" and
                        # EARTHDATA_USERNAME/PASSWORD credentials
```

`sentinel3`/`sentinel6`'s `sigma0_min` is a water-only filter, off (`0.0`, a no-op) by default. Unlike a reservoir polygon, a buffered river corridor genuinely includes riverbank/vegetation, and unlike ICESat-2, Sentinel has no built-in water classification — tune this against real data for your rivers before relying on a nonzero value.

## Advanced: Hydrocron fields and quality filters

By default, a sensible set of fields is requested for each feature type. Add a `swot:` section to override:

```yaml
swot:
  hydrocron_fields:
    nodes:   ["node_id", "node_q", "reach_id", "time_str", "wse", "wse_u", "p_wse",
              "geoid_hght", "sword_version", "solid_tide", "load_tidef", "pole_tide",
              "width", "width_u", "p_width", "xovr_cal_q", "rdr_sig0", "xovr_cal_c", "dark_frac"]
    reaches: ["reach_id", "reach_q", "time_str", "wse", "wse_u", "slope", "slope_u",
              "slope2", "slope2_u", "width", "width_u", "geoid_hght", "solid_tide",
              "load_tidef", "pole_tide", "p_wse", "p_width"]
  quality_filters:
    nodes:   {max_q: 2}    # keep records where node_q <= 2
    reaches: {max_q: 2}    # keep records where reach_q <= 2
```

`reaches`' `slope`/`slope_u` fields feed the optional reach-level slope correction (`merging_options.use_reach_slope_correction`, off by default, reach-mode only — see `flows.DEFAULT_RIVER_MERGING_OPTIONS`); they're silently unused if you switch to node mode or don't enable that option.

## Output structure

```
{main_dir}/
  aux/SWORD/
    SWORD_subset.gpkg          # SWORD v17b subset (spatial intersection with AOI)
    gpkg/                      # full SWORD v17b database (if downloaded)
  raw/
    swot/
      <wb_id>/
        nodes_timeseries.csv    # or reaches_timeseries.csv
    icesat2/<wb_id>/            # if icesat2 enabled
    sentinel3/<wb_id>/          # if sentinel3 enabled
    sentinel6/<wb_id>/          # if sentinel6 enabled
  results/
    <wb_id>/                    # one map + one combined timeseries plot per waterbody
      <nodes|reaches>_map.png
      <nodes|reaches>_timeseries.png
    <target_id>/                # per node/reach (not nested under <wb_id>)
      raw_observations/
        swot.gpkg
        icesat2.gpkg             # if enabled
        sentinel3.gpkg           # if enabled
        sentinel6.gpkg           # if enabled
      cleaned_observations/
        swot.csv, icesat2.csv, ...
      merged_timeseries.csv
      merged_progress/          # per-processing-step diagnostic CSVs
      run_config.yaml           # per-target exclusions + merging-option overrides
      merging_summary.png       # (if save=True)
```

`<wb_id>` is the AOI feature's `id_key` value (Option A) or the configured `id` (Option B) — it groups targets for download batching and the two waterbody-level plots, but every target's actual data (`raw_observations/`, `cleaned_observations/`, `merged_timeseries.csv`, `run_config.yaml`) lives directly under its own `<target_id>` (node or reach ID), not nested inside `<wb_id>`.

## Per-target exclusion / merging-option overrides

Same mechanism as reservoirs — reachable via `Project` methods, not by importing `flows` directly:

```python
project.list_target_observations(23221000160051, target_type="rivers")
project.exclude_observations(23221000160051, platform="S3B", reason="...")
project.list_exclusions(23221000160051, target_type="rivers")
project.remove_exclusion(23221000160051, platform="S3B", orbit=1517)
project.set_merging_option(23221000160051, svr_radial_err=1.0)
```

`target_type` can be omitted — it's inferred automatically from whichever of reservoirs/rivers your project actually configures.

## Known limitations

- **River `merging_options` defaults are a direct copy of the reservoir-tuned defaults, not independently validated against real river data.** Check kept/rejected observation counts on your own rivers before trusting them as-is.
- **No `dfs0` export for rivers yet** (reservoirs-only).
- **Reach-level slope correction** (`use_reach_slope_correction`) has two assumptions not yet validated against real data: that SWOT's own reach-level WSE is approximately midpoint-referenced, and that the correction's sign convention actually reduces cross-mission scatter rather than increasing it. Confirm both against your own reach before relying on this in production.
- **Sentinel-6's `orbit_key` stability** hasn't been re-verified against real Sentinel-6 data the way Sentinel-3's was (see `flows.PRODUCT_TIMESERIES_KEYS`).