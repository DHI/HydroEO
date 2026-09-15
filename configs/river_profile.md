# River Profile — Configuration Reference

Compute a cleaned, longitudinal water-surface-elevation (WSE) profile along a
river from SWOT L2 HR Raster tiles, sampled at chainage points you supply.
HydroEO downloads the matching SWOT tiles for you (reusing the same
download/preprocess pipeline as [`swot_raster`](swot_raster.md), with
`merge_tiles` forced off so each acquisition date stays its own tile) and
runs a configurable multi-stage filter over the sampled values. Requires
NASA Earthdata credentials.

## Chainage input

`river_profile.chainage_path` must point to a point shapefile/geopackage
with:

- A defined CRS (any — reprojected internally as needed).
- A numeric column (`river_profile.chainage_field`, default `cngmeters`)
  giving each point's distance along the river, in metres. HydroEO does not
  compute this for you — it must already exist in the file (e.g. from a
  "generate points along line" + chainage-measurement GIS workflow).

Set `reverse_chainage: true` if the column increases downstream but you want
distance 0 to be the upstream end (or vice versa) — HydroEO flips it via
`max(chainage) - chainage`.

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

Six stages run in order, each individually toggleable via
`filters.<stage>.enabled`. All defaults below (other than `preclip`) match
field-tested values from the reference implementation.

| Stage | Purpose |
| --- | --- |
| `preclip` | Hard clip to a plausible elevation range (`min`/`max`, metres) before any statistical filtering. Defaults are intentionally broad (`-10` to `8000` m) — narrow this to your river's actual WSE range for better outlier rejection |
| `soft_clamp` | Detrend within along-river bins (`bin_width_m`) and softly clamp points away from the bin's dominant vertical mode — handles layover/multi-return without hard-masking |
| `hampel_1` | Distance-windowed (`win_m`) Hampel outlier filter around a robust local trend; `action: mask` drops flagged points, `action: replace` snaps them to the trend/window median |
| `rolling_quantile` | Local robust-trend regression + rolling quantile (`q`) of the residual, per window — a smoothed reference used for density culling |
| `density_cull` | Drops points from stretches with too few nearby valid observations (`total_win_m` window; below the `low_pct` percentile or `abs_min` absolute count) |
| `hampel_2` | A second, typically milder Hampel pass after density culling |
| `spline_fill` | Fits an LSQ spline (degree `k`) to remaining valid points and pastes it **only** into remaining NaN gaps — never overwrites real data |

Disabling a stage passes its input straight through to the next stage
unchanged.

## Quality report

After processing, a `quality_report.csv` is written to
`results/<name>/` with one row per acquisition date: how many points were
excluded by `orbit_exclusions`, how many were changed/flagged/dropped by
each filter stage, and the finite-point count before/after. Aggregate totals
are also logged.

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
      <name>_profiles_final.csv                  # always
      quality_report.csv                         # always
      plots/per_profile/*.png                    # if plot_enable
      plots/combined/*.png                       # if plot_enable
      profiles_raw/, profiles_prefilter/,        # only if keep_intermediates: true
      profiles_hampel1/, profiles_geoid/
      <name>_profiles_raw.csv, _prefilter.csv, _hampel1.csv   # only if keep_intermediates
```
