"""Diagnostic plots for reservoirs and rivers, plus shared observation
loading helpers.

generate_reservoirs_summaries is tested with _load_product_timeseries
patched as a sibling; generate_rivers_summaries is tested with
_has_enough_observations_to_plot/_project_num_months/_river_target_corridor/
_load_merged_timeseries patched as siblings (tests/unit/test_flows.py) --
so all of these live in one module. _river_target_corridor is also called
directly (unpatched) by the download/extraction modules, which import it
from here.
"""

import logging
import os
import datetime

import geopandas as gpd
import pandas as pd

from HydroEO import plotting
from ._river_common import _group_river_targets_by_waterbody, _river_extraction_buffer_meters

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from HydroEO.project import Project

logger = logging.getLogger(__name__)


def _river_target_corridor(
    prj: "Project", target_ids, buffer_meters=None, width_buffer_factor=1.05,
):
    """
    Build one buffered, dissolved corridor polygon covering the given
    river targets (nodes or reaches), for use as the spatial AOI when
    downloading/extracting ICESat-2 and Sentinel-3/6 observations.

    This is deliberately a SEPARATE buffer distance from
    prj.rivers.buffer_meters (used earlier to decide which SWORD
    targets intersect the user's AOI at all) -- that question ("is this
    target in scope") and this one ("how far from the centerline could
    real river water still be, for a raw altimetry point to plausibly
    belong to this target") are different, and conflating them risks
    the same "one parameter doing two jobs badly" issue found elsewhere
    in this pipeline.

    Parameters
    ----------
    buffer_meters : float or None, optional
        Explicit, uniform buffer distance (meters), applied to every
        target regardless of its actual width. If None (default), uses
        each target's own SWORD "width" attribute instead: buffer
        distance = (width / 2) * width_buffer_factor. This is HALF the
        width, not the full width -- buffering a line expands it
        symmetrically by the given distance on EACH side, so a buffer
        of width/2 gives a corridor whose TOTAL span is approximately
        width * width_buffer_factor, matching the river's actual extent
        plus a margin, rather than doubling it. Falls back to
        _river_extraction_buffer_meters() (a flat scalar) if no usable
        "width" column is found -- e.g. if your SWORD data names it
        differently than assumed here, this degrades gracefully with a
        log message rather than failing.
    width_buffer_factor : float, optional
        Margin applied on top of each target's own width when using the
        width-based default. Default 1.05 -- 5% wider than the river's
        actual channel width. Only used when buffer_meters is None.

    NOTE: "width" is the expected SWORD column name per the standard
    SWORD data dictionary -- this has NOT been verified against a real
    downloaded SWORD file in this session (no sample data was
    available), unlike most other assumptions in this codebase. Check
    your actual target_features columns if width-based buffering
    doesn't seem to be kicking in.

    Returns a single-row GeoDataFrame in prj.global_crs (matching what
    icesat2.extract_observations/sentinel.extract_observations expect
    for their `features` argument, same as reservoirs), or None if no
    matching SWORD geometry is found for target_ids. Note the returned
    geometry may be a MultiPolygon if targets form disconnected pieces
    (e.g. separate reaches far enough apart that their buffers never
    touch) -- see _iter_geometry_pieces for how downloads handle this.
    """
    features = prj.rivers.target_features
    subset = features.loc[features[prj.rivers.target_id_col].isin(target_ids)]
    if subset.empty:
        return None

    local = subset.to_crs(prj.local_crs)

    if buffer_meters is not None:
        distances = buffer_meters
    elif "width" in local.columns and local["width"].notna().any():
        fallback_width = local["width"].median()
        distances = (local["width"].fillna(fallback_width) / 2) * width_buffer_factor
    else:
        logger.info(
            "No 'width' column found in SWORD target_features for this "
            "waterbody -- falling back to a flat extraction buffer "
            "instead of width-based sizing. Check your SWORD data's "
            "actual column names if this is unexpected."
        )
        distances = _river_extraction_buffer_meters(prj)

    buffered = local.buffer(distances)
    corridor = buffered.unary_union
    corridor_gdf = gpd.GeoDataFrame(
        geometry=[corridor], crs=prj.local_crs
    ).to_crs(prj.global_crs)
    return corridor_gdf


def _load_product_timeseries(data_dir, ext, products, reader_fn):
    """Load files of given extension from directory, optionally filtered to products."""
    if not os.path.exists(data_dir):
        return None
    df_list = []
    for file in os.listdir(data_dir):
        if file.endswith(ext):
            if not products or file.split(".")[0] in products:
                try:
                    df_list.append(reader_fn(os.path.join(data_dir, file)))
                except Exception as exc:
                    logger.warning(
                        "Failed to load %s from %s: %s",
                        file,
                        data_dir,
                        exc,
                    )
    return pd.concat(df_list) if df_list else None


def generate_reservoirs_summaries(
    prj: "Project", show: bool = False, save: bool = True
) -> None:
    """Generate per-reservoir plotting summaries.

    Parameters
    ----------
    prj : Project
        Project instance with reservoirs configuration
    show : bool
        Whether to display plots interactively
    save : bool
        Whether to save plots to disk
    """
    if not hasattr(prj, "reservoirs"):
        return

    for reservoir_id in prj.reservoirs.download_gdf[prj.reservoirs.id_key]:
        plotting.plot_crossings(
            gdf=prj.reservoirs.gdf,
            id_key=prj.reservoirs.id_key,
            reservoir_id=reservoir_id,
            output_dir=prj.dirs["output"],
            reservoir_type=prj.reservoirs.type,
            show=show,
            save=save,
        )

        plotting.plot_cleaning(
            reservoir_id=reservoir_id,
            output_dir=prj.dirs["output"],
            get_unfiltered_fn=lambda id, products: _load_product_timeseries(
                os.path.join(prj.dirs["output"], f"{id}", "raw_observations"),
                ".gpkg",
                products,
                lambda path: gpd.read_file(path).drop(columns=["geometry"]),
            ),
            get_cleaned_fn=lambda id, products: _load_and_parse_cleaned_timeseries(
                prj, id, products
            ),
            get_merged_fn=lambda id: _load_merged_timeseries(prj, id),
            reservoir_type=prj.reservoirs.type,
            show=show,
            save=save,
            products=getattr(prj, "to_process", None),
        )

        plotting.plot_merging(
            reservoir_id=reservoir_id,
            output_dir=prj.dirs["output"],
            reservoir_type=prj.reservoirs.type,
            show=show,
            save=save,
        )


def _load_and_parse_cleaned_timeseries(prj, id, products):
    """Load cleaned observations and parse dates."""
    df = _load_product_timeseries(
        os.path.join(prj.dirs["output"], f"{id}", "cleaned_observations"),
        ".csv",
        products,
        pd.read_csv,
    )
    if df is not None:
        df["date"] = pd.to_datetime(df.date, format="mixed", utc=True).dt.tz_convert(
            None
        )
        df = df.sort_values(by="date")
    return df


def _load_merged_timeseries(prj, id):
    """Load merged timeseries if it exists."""
    data_path = os.path.join(prj.dirs["output"], f"{id}", "merged_timeseries.csv")

    if os.path.exists(data_path):
        df = pd.read_csv(data_path)
        df["date"] = pd.to_datetime(df.date)
        df = df.sort_values(by="date")
        return df
    else:
        logger.warning(
            "%s does not exist, be sure to merge product timeseries first!",
            data_path,
        )
        return None


# ============================================================================
# RIVERS: Summaries & Visualization
# ============================================================================


def _project_num_months(prj: "Project") -> int:
    """
    Approximate number of months spanned by the project's configured
    date range -- used as a minimum-observation-count threshold for
    plotting (see _has_enough_observations_to_plot). Falls back to 1 if
    the project-level dates aren't resolvable for some reason.
    """
    project_cfg = prj.config.get("project", {})
    start = project_cfg.get("startdate")
    end = project_cfg.get("enddate")
    if not start or not end:
        return 1

    start_date = datetime.date(*start) if isinstance(start, list) else start
    end_date = datetime.date(*end) if isinstance(end, list) else end
    months = (
        (end_date.year - start_date.year) * 12
        + (end_date.month - start_date.month)
        + 1
    )
    return max(months, 1)


def _has_enough_observations_to_plot(prj: "Project", target_id, min_months: int) -> bool:
    """
    Whether a target has enough merged observations to be worth
    plotting -- more than min_months (the project's date range in
    months) or more than 2, whichever is larger. A reach/reservoir with
    only 1-2 points produces a plot that adds noise without telling you
    anything.
    """
    df = _load_merged_timeseries(prj, target_id)
    if df is None:
        return False
    threshold = max(min_months, 2)
    return len(df) > threshold


def generate_rivers_summaries(
    prj: "Project", show: bool = False, save: bool = True
) -> None:
    """Generate per-river plotting summaries.

    Parameters
    ----------
    prj : Project
        Project instance with rivers configuration
    show : bool
        Whether to display plots interactively
    save : bool
        Whether to save plots to disk
    """
    if not hasattr(prj, "rivers"):
        return

    waterbody_groups = _group_river_targets_by_waterbody(prj)
    min_months = _project_num_months(prj)

    for wb_id, target_ids in waterbody_groups.items():
        # Only plot targets with enough observations to be worth looking
        # at -- applies to all three plot types (map, time series, merge
        # progress) so a target excluded from one isn't confusingly still
        # shown in another.
        plottable_ids = [
            t for t in target_ids
            if _has_enough_observations_to_plot(prj, t, min_months)
        ]
        if not plottable_ids:
            logger.info(
                "Skipping plots for waterbody %s -- no targets with more "
                "than %d observations.", wb_id, max(min_months, 2),
            )
            continue

        # Compute the actual extraction corridor (same buffer resolution
        # used for real extraction, see _river_target_corridor) so the
        # shaded area shown is exactly what extraction uses, not an
        # approximation -- lets you visually judge whether width-based
        # buffering produced a reasonable corridor for this waterbody
        # (e.g. a lake-flagged reach whose SWORD width reflects a much
        # wider lake extent) without needing external knowledge of the
        # real river geometry.
        corridor_gdf = _river_target_corridor(
            prj, plottable_ids,
            buffer_meters=getattr(prj.rivers, "extraction_buffer_meters", None),
            width_buffer_factor=getattr(prj.rivers, "width_buffer_factor", 1.05),
        )
        corridor_geometry = corridor_gdf.geometry.iloc[0] if corridor_gdf is not None else None

        plotting.plot_river_crossings(
            prj, wb_id, plottable_ids, prj.dirs["output"], show=show, save=save,
            corridor_geometry=corridor_geometry,
        )

        plotting.plot_river_data(
            prj, wb_id, plottable_ids, prj.dirs["output"],
            get_merged_fn=lambda id: _load_merged_timeseries(prj, id),
            show=show, save=save,
        )

        for target_id in plottable_ids:
            plotting.plot_merging(
                reservoir_id=target_id,
                output_dir=prj.dirs["output"],
                reservoir_type="river",
                show=show,
                save=save,
            )


# ============================================================================
# MIKEIO
# ============================================================================


