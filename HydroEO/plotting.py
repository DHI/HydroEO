"""HydroEO plotting utilities for summarizing altimetry observations and processing results."""

import logging
import os
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from cmcrameri import cm
import seaborn as sns
import pandas as pd
import geopandas as gpd


logger = logging.getLogger(__name__)


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
    sns.set()
    cmap = cm.batlow.resampled(5)
    colors = {
        "icesat2": cmap(0),
        "sentinel3": cmap(1),
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
        edgecolor="black",
        facecolor="None",
        zorder=5,
        label="Reservoir outline",
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
                zorder = 0
                alpha = 0.1
            else:
                zorder = 10
                alpha = 0.5
            gdf_product.plot(
                ax=ax,
                color=colors[product],
                edgecolor="none",
                alpha=alpha,
                zorder=zorder,
                label=product,
            )
            plotted_products.add(product)

    legend_handles = [
        Line2D(
            [],
            [],
            color="black",
            linewidth=1.0,
            label="Reservoir outline",
        )
    ]

    for product in sorted(plotted_products):
        if product in colors:
            legend_handles.append(
                Line2D(
                    [],
                    [],
                    marker="o",
                    linestyle="None",
                    color=colors[product],
                    label=product,
                    alpha=0.1 if product == "swot" else 0.5,
                )
            )

    if legend_handles:
        ax.legend(
            handles=legend_handles,
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            borderaxespad=0.0,
        )
    fig.tight_layout(rect=(0, 0, 0.82, 1))
    if save:
        plt.savefig(os.path.join(output_dir, f"{reservoir_id}", "crossing_summary.png"))
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
    sns.set()
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

        ax = main_ax[0]
        ax.set_title("Unfiltered Products")
        for platform in df.platform.unique():
            df.loc[df.platform == platform].plot(
                ax=ax,
                x="date",
                y="height",
                c=colors[platform],
                kind="scatter",
                label=platform,
            )

        handles, labels = ax.get_legend_handles_labels()
        if labels:
            ax.legend(
                loc="center left",
                bbox_to_anchor=(1.02, 0.5),
                borderaxespad=0.0,
            )
        ax.tick_params(axis="x", rotation=45)

        # now plot cleaned timeseries
        df = get_cleaned_fn(reservoir_id, products)
        df = df[["date", "height", "platform", "product"]]

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
            )

        handles, labels = ax.get_legend_handles_labels()
        if labels:
            ax.legend(
                loc="center left",
                bbox_to_anchor=(1.02, 0.5),
                borderaxespad=0.0,
            )
        ax.tick_params(axis="x", rotation=45)

    # now plot merged timeseries
    df = get_merged_fn(reservoir_id)
    if df is not None:
        df = df[["date", "height"]]

        ax = main_ax[2]
        ax.set_title("Merged Timeseries")
        df.plot(ax=ax, x="date", y="height", c="k", kind="scatter", label="merged")

        handles, labels = ax.get_legend_handles_labels()
        if labels:
            ax.legend(
                loc="center left",
                bbox_to_anchor=(1.02, 0.5),
                borderaxespad=0.0,
            )
        ax.tick_params(axis="x", rotation=45)

        fig.tight_layout(rect=(0, 0, 0.82, 1))
        if save:
            plt.savefig(
                os.path.join(output_dir, f"{reservoir_id}", "cleaning_summary.png")
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

        sns.set()
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
                os.path.join(output_dir, f"{reservoir_id}", "merging_summary.png")
            )
        if show:
            plt.show()

    return
