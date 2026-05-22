# Rivers — Configuration Reference

SWOT Hydrocron download for river nodes and reaches, matched to the SWORD v17b river network database. **No credentials required** — Hydrocron is a public API.

> **Status:** Partial. Initialization, SWOT Hydrocron download, and diagnostic plotting are supported. River timeseries filtering/preprocessing is not yet implemented.

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
| `buffer_meters` | float | — | Optional AOI buffer before SWORD spatial subsetting |
| **Option B — explicit IDs** | | | |
| `feature_numbers` | list | ✅ (B) | List of node or reach IDs |
| `feature_type` | string | ✅ | Must match IDs: `nodes` or `reaches` |
| `id` | string | ✅ (B) | Used as the output folder name |

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

## Output structure

```
{main_dir}/
  aux/SWORD/
    SWORD_subset.gpkg          # SWORD v17b subset (spatial intersection with AOI)
    gpkg/                      # full SWORD v17b database (if downloaded)
  raw/
    swot/
      <wb_id>/
        nodes_timeseries.csv   # or reaches_timeseries.csv
  results/
    <wb_id>/                   # diagnostic plots (map, timeseries)
```

Each node/reach gets its own subfolder within the AOI-feature folder so multiple nodes from the same AOI feature coexist safely.
