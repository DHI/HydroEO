# River Profile — Configuration Reference

Compute a cleaned, longitudinal water-surface-elevation (WSE) profile along a
river extending into the coastal zone of deltaic systems, from SWOT L2 HR
Raster tiles sampled at chainage points you supply. HydroEO downloads the
matching SWOT tiles for you (reusing the same download/preprocess pipeline as
[`swot_raster`](swot_raster.md), with `merge_tiles` forced off so each
acquisition date stays its own tile) and runs a configurable multi-stage
filter over the sampled values. Requires NASA Earthdata credentials.

## Limitations

The WSE interpolation/processing step assumes a continuous water-surface
profile along the river reach. It is therefore not valid across hydraulic
structures such as dams, weirs, and gates, where the water surface is
discontinuous.

**Workaround:** split the river network at hydraulic structures and process
each reach separately.

## Chainage input

`river_profile.chainage_path` must point to a point shapefile/geopackage
with:

- A defined CRS (any — reprojected internally as needed).
- A numeric column (`river_profile.chainage_field`, default `cngmeters`)
  giving each point's distance along the river, in metres. HydroEO does not
  compute this for you — it must already exist in the file (e.g. from a
  "generate points along line" + chainage-measurement GIS workflow).

`reverse_chainage: true` swaps which end of the river is distance 0 — HydroEO
computes `new_chainage = max(chainage) - chainage`, so whichever end
currently holds the maximum value becomes 0 and vice versa. Set it if your
file's 0 point is at the wrong end (e.g. at the downstream end when you want
0 to mark the upstream end).

## Config reference

Start from [`configs/river_profile.yaml`](river_profile.yaml).

| Key | Default | Description |
| --- | --- | --- |
| `name` | — | Identifier used as the output filename/directory prefix |
| `chainage_path` | — | Path to the chainage point shapefile/geopackage |
| `chainage_field` | `cngmeters` | Column with along-river distance (metres) |
| `reverse_chainage` | `false` | Flip chainage direction |
| `startdate` / `enddate` | project dates | Temporal range for the SWOT search |
| `aoi_buffer_meters` | `2000` | Buffer around the chainage extent used for the SWOT granule search/clip |
| `product` | `SWOT_L2_HR_Raster_D` | Only the HR raster product is supported |
| `variables` | `["wse"]` | Add `"geoid"` to also extract/save per-date geoid profiles |
| `granule_filter` | — | Optional glob pattern to filter granules (e.g. `"*100m*"`) |
| `zero_is_nodata` | `true` | Treat exact `0` as nodata in addition to the raster's own nodata value |
| `keep_intermediates` | `false` | `false` = only final profile shapefiles + final CSV; `true` = also write every intermediate stage |
| `plot_enable` | `true` | Write per-date + combined diagnostic PNGs |
| `plot_dpi` | `150` | Plot resolution |
| `plot_ylim` | — | Optional `[min, max]` y-axis override for all plots |
| `orbit_exclusions` | `[]` | Manual chainage-range exclusions per orbit (see below) |
| `filters` | see below | Per-stage filter settings |

### `orbit_exclusions`

For a river where a specific orbit/pass is known to produce bad data over a
specific chainage range (a manual correction, not something HydroEO detects
automatically), add entries like:

```yaml
orbit_exclusions:
  - orbit: "467"            # matched as a substring of the granule filename
    max_chainage_m: 50000   # exclude chainage <= 50 km for this orbit
  # min_chainage_m: ...     # optional lower bound; both bounds are optional
```

Points matching `min_chainage_m <= chainage <= max_chainage_m` for a granule
whose filename contains `orbit` are set to NaN before filtering. Leave this
empty unless you've identified a systematic artifact for your data.

## Filtering pipeline

Seven stages run in order, each individually toggleable via
`filters.<stage>.enabled`. All defaults below (other than `preclip`) match
field-tested values from the reference implementation.

| Stage | Purpose |
| --- | --- |
| `preclip` | Hard clip to a plausible elevation range (`min`/`max`, metres) before any statistical filtering. Defaults are intentionally broad (`-10` to `8000` m) — narrow this to your river's actual WSE range for better outlier rejection |
| `soft_clamp` | Detrend within along-river bins (`bin_width_m`) and softly clamp points away from the bin's dominant vertical mode — handles layover/multi-return without hard-masking |
| `hampel_1` | Distance-windowed (`win_m`) Hampel outlier filter around a robust local trend; `action: mask` drops flagged points, `action: replace` snaps them to the trend/window median |
| `rolling_quantile` | Local robust-trend regression + rolling quantile (`q`) of the residual, per window. **Its output becomes the working profile value carried into the remaining stages** (not merely a side reference) — this is what turns noisy per-point samples into a locally-smoothed profile; `density_cull` below decides which of these values to keep |
| `density_cull` | Drops points from stretches with too few nearby valid observations (`total_win_m` window; strictly below the `low_pct` percentile, or below the `abs_min` absolute count) |
| `hampel_2` | A second, typically milder Hampel pass after density culling |
| `spline_fill` | Fits an LSQ spline (degree `k`) to whatever values survived the stages above and pastes it **only** into remaining NaN gaps — never overwrites a value already produced by an earlier stage |

Disabling a stage passes its input straight through to the next stage
unchanged. Note that with `rolling_quantile` enabled (the default), disabling
`spline_fill` still yields a smoothed profile, not raw per-point
observations — to keep raw values, disable `rolling_quantile` too (and
inspect `profiles_raw`/`profiles_hampel1` via `keep_intermediates: true`).

## Quality report

After processing, a `quality_report.csv` is written to
`results/<name>/` with one row per acquisition date: how many points were
excluded by `orbit_exclusions`, how many were changed/flagged/dropped by
each filter stage, and the finite-point count before/after. Aggregate totals
are also logged.

## Nodes map

Every run writes `results/<name>/<name>_nodes_map.html`, a self-contained
interactive page (open it in a browser; needs internet for the map tiles and
the Leaflet/Plotly libraries):

- **Top** — every date's final profile along the river.
- **Bottom left** — the chainage points ("nodes") on a basemap, coloured by
  chainage, with node IDs drawn as you zoom in. Gray canvas and imagery
  basemaps are available.
- **Bottom right** — the final WSE timeseries of the selected node. Select a
  node by clicking it on the map, clicking a point on the profile, or typing
  its ID in the "Node" box. In the chart toolbar, the camera icon saves it as
  PNG and the disk icon downloads its data as CSV (`datetime_utc`, `label`,
  `node_id`, `distance_along_river_m`, `wse_m`).

Node IDs are `0..N-1` in chainage order (after `reverse_chainage`), and match
the `node_id` column of `<name>_profiles_final.csv` and of the
`profiles_final/*.shp` files. Any existing `node_id` column in the chainage
file is replaced. The timeseries shows final (spline-filled) values, so gaps
on a given date may be interpolated rather than observed.

## Output structure

```
{main_dir}/
  raw/
    swot_raster/<name>/<product>/       # raw SWOT granules (shared download step)
  processed/
    swot_raster/<name>/<product>/       # per-granule clipped wse (+geoid) tiles
  results/
    <name>/
      profiles_final/<name>_<date>_profile.shp   # always
      <name>_profiles_final.csv                  # always (node_id, distance, one column per date)
      <name>_nodes_map.html                      # always (interactive profile + nodes map + node timeseries)
      quality_report.csv                         # always
      plots/per_profile/*.png                    # if plot_enable
      plots/combined/*.png                       # if plot_enable
      profiles_geoid/                            # if "geoid" in variables (independent of keep_intermediates)
      profiles_raw/, profiles_prefilter/,        # only if keep_intermediates: true
      profiles_hampel1/
      <name>_profiles_raw.csv, _prefilter.csv, _hampel1.csv   # only if keep_intermediates
```

## Reference

For details on the processing workflow, see:

> Coppo Frias, M., Kittel, C. M. M., Nielsen, K., Shamsudduha, M., Hossain, S.,
> Musaeus, A. F., Toettrup, C., & Bauer-Gottwein, P. (2026). Resolving
> River-Coastal Water Surface Elevation Profiles with SWOT: Insights into
> Tide-Discharge Interactions in Large Deltas. *ESS Open Archive* (preprint).
> https://doi.org/10.22541/essoar.177170504.44018148/v2

If you use this workflow, please cite this publication:

```bibtex
@article{coppofrias2026swot,
  author    = {Coppo Frias, Monica and Kittel, Cécile M. M. and Nielsen, Karina and Shamsudduha, Mohammad and Hossain, Sazzad and Musaeus, Aske Folkmann and Toettrup, Christian and Bauer-Gottwein, Peter},
  title     = {Resolving River-Coastal Water Surface Elevation Profiles with SWOT: Insights into Tide-Discharge Interactions in Large Deltas},
  journal   = {ESS Open Archive},
  year      = {2026},
  publisher = {Wiley},
  note      = {Preprint},
  doi       = {10.22541/essoar.177170504.44018148/v2}
}
```

*Note: this reference points to a preprint and will be updated once the
peer-reviewed version is published.*
