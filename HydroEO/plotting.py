"""HydroEO plotting utilities for summarizing altimetry observations and processing results."""

import logging
import os
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from cmcrameri import cm
import pandas as pd
import geopandas as gpd
import contextily as ctx
from typing import TYPE_CHECKING
from HydroEO.utils import general

if TYPE_CHECKING:
    from HydroEO.project import Project


logger = logging.getLogger(__name__)

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
    cmap = cm.batlow.resampled(5)
    colors = {
        "icesat2": cmap(0),
        "sentinel3": cmap(2),
        "sentinel6": cmap(3),
        "swot": cmap(4),
    }

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
        if file_name.endswith(".shp"):
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
        ax.grid(True, linestyle="--", alpha=0.3)
        ax.tick_params(axis="x", rotation=45)
        ax.xaxis.set_major_locator(mdates.AutoDateLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))

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
        ax.grid(True, linestyle="--", alpha=0.3)
        ax.tick_params(axis="x", rotation=45)
        ax.tick_params(direction="in")
        ax.xaxis.set_major_locator(mdates.AutoDateLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))

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
        ax.grid(True, linestyle="--", alpha=0.3)
        ax.tick_params(axis="x", rotation=45)
        ax.tick_params(direction="in")
        ax.xaxis.set_major_locator(mdates.AutoDateLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))

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
    """Plot merging progression for a reservoir (intermediate processing steps).

    Parameters
    ----------
    reservoir_id : str/int
        The specific reservoir ID to plot.
    output_dir : str
        Base output directory containing merged_progress subdirectory.
    reservoir_type : str
        Type of reservoir (e.g., 'reservoirs').
    show : bool
        Whether to display the plot.
    save : bool
        Whether to save the plot to PNG.
    """
    merged_dir = os.path.join(output_dir, str(reservoir_id), "merged_progress")

    if os.path.exists(merged_dir):
        file_names = os.listdir(merged_dir)
        num_files = len(file_names)

        fig, main_ax = plt.subplots(num_files, 1, figsize=(10, 10))
        fig.suptitle(f"{reservoir_type}: {reservoir_id}")

        def _plot_file(i, fig, main_ax, file_name, file_names, dir):
            if file_name in file_names:
                i = i + 1
                ax = main_ax.flat[i]
                df = pd.read_csv(os.path.join(dir, file_name))
                df["date"] = pd.to_datetime(df.date)
                df.plot(ax=ax, x="date", y="height", c="k", kind="scatter")
                ax.set_title(file_name)
                ax.grid(True, linestyle="--", alpha=0.3)
            return i

        i = -1
        i = _plot_file(i, fig, main_ax, "svr_linear.csv", file_names, merged_dir)
        i = _plot_file(i, fig, main_ax, "mad_filter.csv", file_names, merged_dir)
        i = _plot_file(i, fig, main_ax, "daily_mad_error.csv", file_names, merged_dir)
        i = _plot_file(i, fig, main_ax, "kalman.csv", file_names, merged_dir)
        i = _plot_file(i, fig, main_ax, "svr_radial.csv", file_names, merged_dir)

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
    """
    sword_dir = prj.dirs["sword"]
    gpkg_name = f"{prj.rivers.continent_key}_sword_{prj.rivers.feature_type}_v17b.gpkg"
    gpkg_path = os.path.join(sword_dir, gpkg_name)

    if not os.path.exists(gpkg_path):
        raise FileNotFoundError(f"Expected SWORD file not found: {gpkg_path}")

    sword_gdf = gpd.read_file(gpkg_path)
    id_label = "nodes" if prj.rivers.target_id_col == "node_id" else "reaches"

    fig, ax = plt.subplots()
    fig.suptitle(f"{wb_id}: {id_label}")

    features = sword_gdf[sword_gdf[prj.rivers.target_id_col].isin(target_ids)]

    xmin, ymin, xmax, ymax = features["geometry"].total_bounds

    ax.set_xlim([xmin - 0.05, xmax + 0.05])
    ax.set_ylim([ymin - 0.05, ymax + 0.05])

    if id_label == "nodes":
        features.plot(ax=ax, color="black", markersize=5, edgecolor="none")
    else:
        features.plot(ax=ax, color="black")

    ctx.add_basemap(
        ax,
        crs=features.crs,
        zoom=15,
        source=ctx.providers.CartoDB.Positron,  # ctx.providers.OpenStreetMap.Mapnik
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
    show: bool = False,
    save: bool = False,
) -> None:
    """Plot water surface elevation (WSE) timeseries for river targets.

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
    show : bool
        Whether to display the plot interactively
    save : bool
        Whether to save the plot to PNG
    """
    id_label = "nodes" if prj.rivers.target_id_col == "node_id" else "reaches"
    csv_path = os.path.join(
        prj.dirs["swot"], "rivers", str(wb_id), f"{id_label}_timeseries.csv"
    )

    fig, ax = plt.subplots()
    fig.suptitle(f"{wb_id}: {id_label}")

    cmap = cm.hawaii
    n = len(target_ids)

    df = pd.read_csv(csv_path)
    df["date"] = pd.to_datetime(df.time_str)
    for idx, target_id in enumerate(target_ids):
        feature_df = df[df[prj.rivers.target_id_col] == target_id]
        feature_df = feature_df.sort_values(by="date", ascending=True)
        feature_df.plot(
            ax=ax,
            x="date",
            y="wse",
            c=cmap(idx / max(n - 1, 1)),
            linewidth="0.5",
            linestyle="-.",
            legend=False,
        )

    ax.set_xlabel("date")
    ax.set_ylabel("wse [m]")
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
