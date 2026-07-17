"""HydroEO plotting utilities for summarizing altimetry observations and processing results."""

import logging
import os
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from cmcrameri import cm
import pandas as pd
import geopandas as gpd
import contextily as ctx
from typing import TYPE_CHECKING
from HydroEO.utils import general

if TYPE_CHECKING:
    from HydroEO.project import Project


logger = logging.getLogger(__name__)

# Shared platform color scheme, used by both plot_crossings (reservoir map)
# and plot_merging (progress plot) so a given platform's color is
# consistent across every plot type, not just within one function.
PLATFORM_COLORS_CMAP = cm.batlow.resampled(10)

# Mission-level colors: used by plot_crossings, which colors by
# raw_observations FILENAME (one file per mission, e.g. "sentinel3.gpkg"
# -- which can contain BOTH S3A and S3B rows together, since the file is
# named by mission, not by satellite instance).
PLATFORM_COLORS = {
    "icesat2": PLATFORM_COLORS_CMAP(0),
    "sentinel3": PLATFORM_COLORS_CMAP(3),
    "sentinel6": PLATFORM_COLORS_CMAP(6),
    "swot": PLATFORM_COLORS_CMAP(9),
}

# Satellite-instance colors: used by plot_merging/plot_river_data, which
# color by the actual "platform" COLUMN value in the data -- this is the
# specific satellite (e.g. "S3A"/"S3B", confirmed from
# sentinel.preprocess.extract_observations: platform = file_split[1] from
# the raw filename), not the generic mission name. S3A/S3B (and S6A/S6B,
# once Sentinel-6B exists) get distinct shades near their parent mission's
# PLATFORM_COLORS entry, so the two dicts stay visually related without
# being identical.
SATELLITE_COLORS = {
    "icesat2": PLATFORM_COLORS["icesat2"],
    "S3A": PLATFORM_COLORS_CMAP(2),
    "S3B": PLATFORM_COLORS_CMAP(4),
    "S6A": PLATFORM_COLORS_CMAP(5),
    "S6B": PLATFORM_COLORS_CMAP(7),
    "swot": PLATFORM_COLORS["swot"],
}

plt.rcParams["font.family"] = "serif"
plt.rcParams.update({"font.size": 10})


def plot_crossings(
    gdf,
    id_key,
    reservoir_id,
    output_dir,
    reservoir_type="reservoirs",
    show=True,
    save=False,
):
    """Plot raw observations crossing for a reservoir.

    Parameters
    ----------
    gdf : GeoDataFrame
        Reservoir geometries with id_key column.
    id_key : str
        Column name for reservoir identifier.
    reservoir_id : str/int
        The specific reservoir ID to plot.
    output_dir : str
        Base output directory containing raw_observations subdirectory.
    reservoir_type : str
        Type of reservoir (e.g., 'reservoirs').
    show : bool
        Whether to display the plot.
    save : bool
        Whether to save the plot to PNG.
    """
    colors = PLATFORM_COLORS

    # start figure
    fig, ax = plt.subplots()
    fig.suptitle(f"{reservoir_type}: {reservoir_id}")

    # extract item shape
    indx = gdf.loc[gdf[id_key] == reservoir_id].index[0]

    # set bounds
    xmin, ymin, xmax, ymax = gdf.loc[indx, "geometry"].bounds

    ax.set_xlim([xmin - 0.1, xmax + 0.1])
    ax.set_ylim([ymin - 0.1, ymax + 0.1])

    # plot reservoir
    gdf.loc[[indx]].plot(
        ax=ax,
        edgecolor="blue",
        facecolor="None",
        zorder=5,
        label="Reservoir outline",
        alpha=0.3,
    )

    # loop through each product in file and plot
    data_dir = os.path.join(output_dir, f"{reservoir_id}", "raw_observations")
    plotted_products = set()
    for file_name in os.listdir(data_dir):
        if file_name.endswith(".gpkg"):
            path_to_file = os.path.join(data_dir, file_name)
            product = file_name.split(".")[0]
            gdf_product = gpd.read_file(path_to_file)
            if product == "swot":
                zorder = 6
                alpha = 0.5
            elif product == "icesat2":
                zorder = 8
                alpha = 0.8
            else:
                zorder = 10
                alpha = 0.8
            gdf_product.plot(
                ax=ax,
                color=colors[product],
                edgecolor="none",
                alpha=alpha,
                zorder=zorder,
                label=product,
                markersize=8,
            )
            plotted_products.add(product)

    legend_handles = [
        Line2D([], [], color="blue", linewidth=2, label="Reservoir", alpha=0.3),
    ]

    for product in sorted(plotted_products):
        if product in colors:
            alpha = 0.5 if product == "swot" else 0.8
            legend_handles.append(
                Patch(facecolor=colors[product], alpha=alpha, label=product)
            )

    if legend_handles:
        ax.legend(
            handles=legend_handles,
            loc="lower center",
            bbox_to_anchor=(0.5, -0.2),
            ncol=len(legend_handles),
            frameon=False,
            # borderaxespad=0.0,
        )

    ax.set_xlabel("lon")
    ax.set_ylabel("lat")
    ax.tick_params(direction="in", width=1.5)
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)

    fig.tight_layout()
    if save:
        plt.savefig(
            os.path.join(output_dir, f"{reservoir_id}", "crossing_summary.png"), dpi=300
        )
    if show:
        plt.show()

    return None


def plot_cleaning(
    reservoir_id,
    output_dir,
    get_unfiltered_fn,
    get_cleaned_fn,
    get_merged_fn,
    reservoir_type="reservoirs",
    show=True,
    save=False,
    products=None,
):
    """Plot cleaning progression for a reservoir (unfiltered, cleaned, merged).

    Parameters
    ----------
    reservoir_id : str/int
        The specific reservoir ID to plot.
    output_dir : str
        Base output directory containing product timeseries.
    get_unfiltered_fn : callable
        Function that returns unfiltered timeseries DataFrame.
    get_cleaned_fn : callable
        Function that returns cleaned timeseries DataFrame.
    get_merged_fn : callable
        Function that returns merged timeseries DataFrame.
    reservoir_type : str
        Type of reservoir (e.g., 'reservoirs').
    show : bool
        Whether to display the plot.
    save : bool
        Whether to save the plot to PNG.
    products : list, optional
        List of products to filter on. Passed to get_*_fn functions.
    """
    cmap = cm.batlow.resampled(5)
    colors = {
        "icesat2": cmap(0),
        "S3A": cmap(1),
        "S3B": cmap(2),
        "S6A": cmap(3),
        "swot": cmap(4),
    }

    fig, main_ax = plt.subplots(3, 1, figsize=(10, 10))
    fig.suptitle(f"{reservoir_type}: {reservoir_id}")

    # plot unfiltered timeseries
    df = get_unfiltered_fn(reservoir_id, products)
    if df is not None:
        df = df[["date", "height", "platform", "product"]]
        df["date"] = pd.to_datetime(df["date"], format="mixed", utc=True)

        ax = main_ax[0]
        ax.set_title("Unfiltered Products")
        for platform in df.platform.unique():
            Q1 = df["height"].quantile(0.25)
            Q3 = df["height"].quantile(0.75)
            IQR = Q3 - Q1
            lower_bound = Q1 - 1.5 * IQR
            upper_bound = Q3 + 1.5 * IQR

            df_normal = df[
                (df["height"] >= lower_bound) & (df["height"] <= upper_bound)
            ]
            df_outliers = df[
                (df["height"] < lower_bound) | (df["height"] > upper_bound)
            ]
            df_normal.loc[df_normal.platform == platform].plot(
                ax=ax,
                x="date",
                y="height",
                c=colors[platform],
                kind="scatter",
                label=platform,
                s=5,
            )

        if not df_outliers.empty:
            for platform in df_outliers.platform.unique():
                platform_outliers = df_outliers[df_outliers.platform == platform]
                for __, row in platform_outliers.iterrows():
                    y_pos = upper_bound if row["height"] > upper_bound else lower_bound
                    marker = "^" if row["height"] > upper_bound else "v"
                    ax.scatter(
                        row["date"],
                        y_pos,
                        marker=marker,
                        color=colors[platform],
                        s=100,
                        zorder=5,
                    )

        handles, labels = ax.get_legend_handles_labels()
        if labels:
            ax.legend(
                loc="center left",
                bbox_to_anchor=(1.02, 0.5),
                borderaxespad=0.0,
            )
        ax.tick_params(axis="x", rotation=45)
        ax.xaxis.set_major_locator(mdates.AutoDateLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
        ax.tick_params(direction="in", width=1.5)
        for spine in ax.spines.values():
            spine.set_linewidth(1.5)

        # now plot cleaned timeseries
        df = get_cleaned_fn(reservoir_id, products)
        df = df[["date", "height", "platform", "product"]]
        df["date"] = pd.to_datetime(df["date"], format="mixed", utc=True)

        ax = main_ax[1]
        ax.set_title("Cleaned Products")
        for platform in df.platform.unique():
            df.loc[df.platform == platform].plot(
                ax=ax,
                x="date",
                y="height",
                c=colors[platform],
                kind="scatter",
                label=platform,
                s=5,
            )

        handles, labels = ax.get_legend_handles_labels()
        if labels:
            ax.legend(
                loc="center left",
                bbox_to_anchor=(1.02, 0.5),
                borderaxespad=0.0,
            )
        ax.tick_params(axis="x", rotation=45)
        ax.tick_params(direction="in")
        ax.xaxis.set_major_locator(mdates.AutoDateLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
        ax.tick_params(direction="in", width=1.5)
        for spine in ax.spines.values():
            spine.set_linewidth(1.5)

    # now plot merged timeseries
    df = get_merged_fn(reservoir_id)
    if df is not None:
        df = df[["date", "height"]]
        df["date"] = pd.to_datetime(df["date"], format="mixed", utc=True)

        ax = main_ax[2]
        ax.set_title("Merged Timeseries")
        df.plot(ax=ax, x="date", y="height", c="k", kind="scatter", label="merged", s=5)

        handles, labels = ax.get_legend_handles_labels()
        if labels:
            ax.legend(
                loc="center left",
                bbox_to_anchor=(1.02, 0.5),
                borderaxespad=0.0,
            )
        ax.tick_params(axis="x", rotation=45)
        ax.tick_params(direction="in")
        ax.xaxis.set_major_locator(mdates.AutoDateLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
        ax.tick_params(direction="in", width=1.5)
        for spine in ax.spines.values():
            spine.set_linewidth(1.5)

        fig.tight_layout()
        if save:
            plt.savefig(
                os.path.join(output_dir, f"{reservoir_id}", "cleaning_summary.png"),
                dpi=300,
            )
        if show:
            plt.show()


def plot_merging(
    reservoir_id,
    output_dir,
    reservoir_type="reservoirs",
    show=True,
    save=False,
):
    """Plot merging progression for a target (intermediate processing steps).

    Works identically for reservoirs and river targets (nodes/reaches) --
    reservoir_id is just used to build the output path
    ({output_dir}/{reservoir_id}/merged_progress), and reservoir_type is
    only a display label. For rivers, pass the node/reach id as
    reservoir_id and e.g. "river" as reservoir_type.

    Parameters
    ----------
    reservoir_id : str/int
        The specific reservoir or river target (node/reach) ID to plot.
    output_dir : str
        Base output directory containing merged_progress subdirectory.
    reservoir_type : str
        Display label (e.g. 'reservoirs' or 'river').
    show : bool
        Whether to display the plot.
    save : bool
        Whether to save the plot to PNG.
    """
    merged_dir = os.path.join(output_dir, str(reservoir_id), "merged_progress")

    if os.path.exists(merged_dir):
        file_names = os.listdir(merged_dir)

        # The actual, current save_progress filenames from
        # Timeseries.merge(), in pipeline order. Previously this
        # referenced "mad_filter.csv", which the pipeline has never
        # actually produced (bias_correct.csv is the real step there),
        # so bias_correct silently never appeared; spatial_correction.csv
        # and distance_penalty.csv (both optional, newer steps) were
        # missing entirely too.
        pipeline_steps = [
            "svr_linear.csv",
            "spatial_correction.csv",
            "bias_correct.csv",
            "daily_mad_error.csv",
            "distance_penalty.csv",
            "kalman.csv",
            "svr_radial.csv",
        ]
        present_steps = [f for f in pipeline_steps if f in file_names]

        if not present_steps:
            logger.warning(
                "No recognized merge-progress files found in %s -- nothing to plot.",
                merged_dir,
            )
            return

        # Base subplot count on recognized-and-present steps, not the raw
        # directory listing -- otherwise a mismatch between allocated
        # subplots and actually-plotted steps leaves blank axes.
        fig, main_ax = plt.subplots(len(present_steps), 1, figsize=(10, 10))
        fig.suptitle(f"{reservoir_type}: {reservoir_id}")

        # plt.subplots returns a bare Axes (not an array) when there's
        # only one subplot -- normalize so the loop below works either way.
        axes = np.atleast_1d(main_ax)

        for ax, file_name in zip(axes, present_steps):
            df = pd.read_csv(os.path.join(merged_dir, file_name))
            df["date"] = pd.to_datetime(df.date)

            if "platform" in df.columns:
                # Color by platform where the data still has it -- lost
                # after Kalman collapses to one value per date (kalman.csv,
                # svr_radial.csv), which fall through to the neutral case.
                for platform, group in df.groupby("platform"):
                    color = SATELLITE_COLORS.get(platform, "gray")
                    ax.scatter(
                        group["date"], group["height"],
                        color=color, label=platform, s=12,
                    )
                ax.legend(fontsize=8, loc="best", frameon=False)
            else:
                ax.scatter(df["date"], df["height"], color="k", s=12)

            ax.set_title(file_name)
            ax.tick_params(direction="in", width=1.5)
            for spine in ax.spines.values():
                spine.set_linewidth(1.5)

        fig.tight_layout()
        if save:
            plt.savefig(
                os.path.join(output_dir, f"{reservoir_id}", "merging_summary.png"),
                dpi=300,
            )
        if show:
            plt.show()

    return


def plot_river_crossings(
    prj: "Project",
    wb_id: str,
    target_ids: list,
    output_dir: str,
    save: bool = False,
    show: bool = False,
    zoom="auto",
    corridor_geometry=None,
) -> None:
    """Plot river target locations (nodes or reaches) on a basemap.

    Parameters
    ----------
    prj : Project
        Project instance with rivers configuration
    wb_id : str
        Waterbody identifier (river name or ID)
    target_ids : list
        List of target node or reach IDs to plot
    output_dir : str
        Base output directory where plots will be saved
    save : bool
        Whether to save the plot to PNG
    show : bool
        Whether to display the plot interactively
    zoom : int or "auto", optional
        Basemap tile zoom level. Default "auto" lets contextily pick a
        zoom appropriate to the actual plotted extent. A fixed zoom
        (e.g. 15, street-level detail) was previously hardcoded here --
        fine for a single reservoir's small extent (which is why the
        reservoir equivalent, plot_crossings, never needed this at all),
        but for a river AOI spanning a whole network, a fixed high zoom
        can mean fetching thousands of tiles over the network regardless
        of how small the local files are -- confirmed as the likely
        cause of a multi-minute stall on a real run. Pass an explicit
        int only if you specifically want more/less detail than the
        auto-selected level for a particular AOI size.
    corridor_geometry : shapely geometry, optional
        The actual buffered extraction corridor for this waterbody (see
        flows._river_target_corridor), drawn as a translucent shaded
        fill behind the reach/node lines. Whether width-based buffering
        (see _river_target_corridor) produces a reasonable corridor size
        -- e.g. for a lake-flagged reach whose SWORD "width" reflects a
        much wider lake extent rather than a narrow channel -- isn't
        something that can be judged from the numbers alone; seeing the
        actual shaded area against the real river geometry is what lets
        you decide whether an adjustment (an explicit
        extraction_buffer_meters, or a max_extraction_buffer_meters
        ceiling) is actually needed. None (default) skips the shading
        entirely -- purely additive, doesn't change any existing
        behavior if not provided.
    """
    subset_path = prj.dirs.get("sword_subset")
    if not subset_path or not os.path.exists(subset_path):
        raise FileNotFoundError(f"SWORD subset not found: {subset_path}")

    sword_gdf = gpd.read_file(subset_path)
    id_label = "nodes" if prj.rivers.target_id_col == "node_id" else "reaches"

    fig, ax = plt.subplots()
    fig.suptitle(f"{wb_id}: {id_label}")

    features = sword_gdf[sword_gdf[prj.rivers.target_id_col].isin(target_ids)]

    xmin, ymin, xmax, ymax = features["geometry"].total_bounds

    # If the corridor extends beyond the reaches' own bounds (it always
    # does, by definition -- it's a buffer around them), widen the axis
    # limits to show the full shaded area, not just the reach lines.
    if corridor_geometry is not None:
        cxmin, cymin, cxmax, cymax = corridor_geometry.bounds
        xmin, ymin = min(xmin, cxmin), min(ymin, cymin)
        xmax, ymax = max(xmax, cxmax), max(ymax, cymax)

    ax.set_xlim([xmin - 0.05, xmax + 0.05])
    ax.set_ylim([ymin - 0.05, ymax + 0.05])

    if corridor_geometry is not None:
        gpd.GeoDataFrame(geometry=[corridor_geometry], crs=features.crs).plot(
            ax=ax, facecolor="gray", edgecolor="gray", alpha=0.25,
            linewidth=0.5, zorder=1,
        )

    # Color each target individually rather than plotting everything in a
    # single black color -- otherwise a multi-reach map is nearly useless
    # for telling targets apart. Uses the same colormap/index scheme as
    # plot_river_data, so a given target's color is consistent between
    # the map and its time series line.
    cmap = cm.hawaii
    n = len(target_ids)
    for idx, target_id in enumerate(target_ids):
        target_features = features[features[prj.rivers.target_id_col] == target_id]
        if target_features.empty:
            continue
        color = cmap(idx / max(n - 1, 1))
        if id_label == "nodes":
            target_features.plot(
                ax=ax, color=color, markersize=5, edgecolor="none", zorder=3,
            )
        else:
            target_features.plot(ax=ax, color=color, zorder=3)
            # Reaches are lines -- unlike nodes, there's nothing marking
            # an individual reach's location distinctly. Add a marker at
            # each reach's midpoint (guaranteed to sit ON the line, unlike
            # .centroid which can fall off it for a curved geometry), in
            # the same color as its line.
            for geom in target_features.geometry:
                midpoint = geom.interpolate(0.5, normalized=True)
                ax.plot(
                    midpoint.x, midpoint.y,
                    marker="o", markersize=5, color=color, linestyle="none",
                    zorder=4,
                )

    ctx.add_basemap(
        ax,
        crs=features.crs,
        zoom=zoom,
        source=ctx.providers.CartoDB.Positron,  # ctx.providers.OpenStreetMap.Mapnik
        zorder=0,
    )

    ax.set_xlabel("lon")
    ax.set_ylabel("lat")
    ax.tick_params(direction="in", width=1.5)
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    fig.tight_layout()

    general.ifnotmakedirs(os.path.join(output_dir, f"{wb_id}"))

    if save:
        plt.savefig(
            os.path.join(output_dir, f"{wb_id}", f"{id_label}_map.png"), dpi=300
        )
    if show:
        plt.show()


def plot_river_data(
    prj: "Project",
    wb_id: str,
    target_ids: list,
    output_dir: str,
    get_merged_fn,
    show: bool = False,
    save: bool = False,
) -> None:
    """Plot water surface elevation (WSE) timeseries for river targets.

    Reads each target's actual merged/cleaned output (via get_merged_fn,
    e.g. flows._load_merged_timeseries) -- NOT the raw Hydrocron CSV.
    The raw per-waterbody CSV is single-mission (SWOT only) and
    unfiltered; the merged output reflects whatever missions were
    actually configured (ICESat-2/Sentinel-3/6 too, if enabled) and has
    already been through quality filtering, bias correction, and Kalman
    smoothing -- it's the thing actually worth looking at, and the two
    can otherwise look completely disconnected from each other.

    Parameters
    ----------
    prj : Project
        Project instance with rivers configuration
    wb_id : str
        Waterbody identifier (river name or ID)
    target_ids : list
        List of target node or reach IDs to plot
    output_dir : str
        Base output directory where plots will be saved
    get_merged_fn : callable
        Function taking a target_id and returning its merged timeseries
        DataFrame (with "date"/"height" columns) or None if unavailable.
    show : bool
        Whether to display the plot interactively
    save : bool
        Whether to save the plot to PNG
    """
    id_label = "nodes" if prj.rivers.target_id_col == "node_id" else "reaches"

    fig, ax = plt.subplots()
    fig.suptitle(f"{wb_id}: {id_label}")

    cmap = cm.hawaii
    n = len(target_ids)

    any_data = False
    for idx, target_id in enumerate(target_ids):
        df = get_merged_fn(target_id)
        if df is None or df.empty:
            continue
        any_data = True
        df = df.sort_values(by="date", ascending=True)
        df.plot(
            ax=ax,
            x="date",
            y="height",
            c=cmap(idx / max(n - 1, 1)),
            linewidth="0.5",
            linestyle="-.",
            marker="o",
            markersize=3,
            legend=False,
        )

    if not any_data:
        logger.warning(
            "No merged data available to plot for waterbody %s (%s) -- "
            "has create_rivers_timeseries() been run for these targets?",
            wb_id, id_label,
        )

    ax.set_xlabel("date")
    ax.set_ylabel("height [m]")
    ax.tick_params(direction="in", width=1.5)
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    fig.tight_layout()

    general.ifnotmakedirs(os.path.join(output_dir, f"{wb_id}"))

    if save:
        plt.savefig(
            os.path.join(output_dir, f"{wb_id}", f"{id_label}_timeseries.png"), dpi=300
        )
    if show:
        plt.show()

    return