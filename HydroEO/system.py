from dataclasses import dataclass
import warnings

import os
import pandas as pd
import geopandas as gpd

from cmcrameri import cm
import seaborn as sns


import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import datetime

from HydroEO.satellites import swot, icesat2, sentinel
from HydroEO.utils import general, timeseries
from HydroEO.downloaders import hydroweb

from tqdm import tqdm


@dataclass
class System:
    gdf: gpd.GeoDataFrame
    id_key: str
    dirs: dict

    def report(self):
        print(f"Number of {self.type}: {len(self.gdf)}")
        return self.gdf.head()

    def download_altimetry(  # TODO: this function has gotten a bit complex and should maybe be refacotred and broken down
        self,
        product: str,
        startdate: tuple,
        enddate: tuple,
        credentials=None,
        grid: bool = False,
        update_existing: bool = False,
    ):
        supported_products = ["ATL13", "S3", "S6", "SWOT_LAKE"]

        # check if the requested product is supported
        product = product.upper()

        if product not in supported_products:
            raise ValueError(
                f'"{product}" is not accepted as a valid download product. Please provide a valid product.'
            )

        ### Set the download geometry for the search bounds
        # set grid as download bounds if needed TODO: add the setting of the download geometries outside of download function
        # if grid:
        #     if hasattr(self, "grid"):
        #         download_gdf = self.grid.copy()
        #     else:
        #         raise NameError(
        #             "Object has no atttribute: grid. Make sure to create grid before enabling this option"
        #         )
        # else:
        #     self.download_gdf = self.gdf

        # unpack and format the download dates
        startdate = datetime.date(*startdate)
        provided_start_date = startdate
        enddate = datetime.date(*enddate)

        ########################################################################################################################
        ##### We have to handle different downloads differently based on the way products are delivered in "granules"
        ### Swot data is provided in large granuals that likley cover most of the full AOI and downloading per reservoir or grid may be redundant
        ### Icesat-2 data can be easily subsetted within the product order so an individual download make sense
        ### Sentinel 3 and 6 data is provided in smaller granules so we can also request per small area and process one at a time
        ########################################################################################################################

        # First check if it is swot data that should be downloaded for the area
        if product == "SWOT_LAKE":
            # define and if needed create directory for each download geometry
            download_dir = os.path.join(self.dirs["swot"], rf"{self.type}")
            general.ifnotmakedirs(download_dir)

            # check if we need to update the start time for the download
            if update_existing:
                lastest_obs = swot.get_latest_obs_date(self.dirs["output"])
                lastest_obs = (
                    lastest_obs.year,
                    lastest_obs.month,
                    lastest_obs.day,
                )  # put back into list because download_altimetry expects this
                startdate = datetime.date(*lastest_obs)

            # grab coordinates of full area of interest
            coords = [
                (x, y)
                for x, y in self.download_gdf.unary_union.envelope.exterior.coords
            ]

            # query all available data — uses baseline D product (SWOT_L2_HR_LakeSP_D)
            print(
                f"Searching for {swot.SWOT_LAKE_SHORT_NAME} for aoi from {startdate} to {enddate}"
            )
            results = swot.query(
                aoi=coords,
                startdate=startdate,
                enddate=enddate,
            )

            # Filter to the prior-lake sub-collection granules.
            # Baseline D file names end with "_prior_*" (lower-case); baseline C
            # used "_Prior_*" (title-case). We match both for robustness.
            print(f"{len(results)} products returned from query")
            to_download = list()
            for result in results:
                link = result.data_links()[0]
                filename = link.split("/")[-1].lower()
                if "_prior_" in filename or "prior" in filename.split("_"):
                    to_download.append(result)
            print(f"{len(to_download)} prior-lake granules selected for download")

            # download the individual file
            _ = swot.download(to_download, download_directory=download_dir)

            # we want to subset the downloaded file to only include known waterbodies TODO: add functionality to skip subsetting if the subsetted file exists
            files_in_dir = [
                os.path.join(download_dir, f)
                for f in os.listdir(download_dir)
                if f.endswith(".zip")
            ]
            swot.subset_by_id(
                files_in_dir, self.download_gdf["prior_lake_id"].astype(int).values
            )

        elif product in ["ATL13"]:
            # loop through download geometry and download data
            for i in self.download_gdf.index:
                id = self.download_gdf.loc[i, self.id_key]
                print(f"\nDowloading data for id {id}:")

                ##### Check if we should edit search date parameters
                startdate = provided_start_date  # keep track of provided start date so we can refresh with each reservoir query
                # check if we need to update the start time for the download
                if update_existing:
                    raw_data_dir = os.path.join(
                        self.dirs["output"], str(id), "raw_observations"
                    )
                    latest_obs = icesat2.get_latest_obs_date(raw_data_dir)
                    if latest_obs is not None:
                        latest_obs = (
                            latest_obs.year,
                            latest_obs.month,
                            latest_obs.day,
                        )  # put back into list because download_altimetry expects this
                        startdate = datetime.date(*latest_obs)

                # grab coordinates of geometry
                coords = [
                    (x, y)
                    for x, y in self.download_gdf.loc[
                        i, "geometry"
                    ].envelope.exterior.coords
                ]

                # define and if needed create directory for each download geometry
                download_dir = os.path.join(
                    self.dirs["icesat2"], rf"{self.type}", rf"{id}"
                )
                general.ifnotmakedirs(download_dir)

                # query and download data
                print(
                    f"Searching for Icesat2 ATL13 for aoi from {startdate} to {enddate}"
                )
                _ = icesat2.query(
                    aoi=coords,
                    startdate=startdate,
                    enddate=enddate,
                    earthdata_credentials=credentials,
                    download_directory=download_dir,
                    product="ATL13",
                )

        elif product in ["S3", "S6"]:
            # set empty variables because no session has been started yet, this will occur at the first download and sessions will refresh automatically
            session_token = None
            session_start_time = None

            # loop through download geometry
            for i in self.download_gdf.index:
                # extract id for saving data
                id = self.download_gdf.loc[i, self.id_key]
                print(f"\nDowloading data for id {id}:")

                ##### Check if we should edit search date parameters
                startdate = provided_start_date  # keep track of provided start date so we can refresh with each reservoir query
                # check if we need to update the start time for the download
                if update_existing:
                    raw_data_dir = os.path.join(
                        self.dirs["output"], str(id), "raw_observations"
                    )
                    latest_obs = sentinel.get_latest_obs_date(raw_data_dir, product)
                    if latest_obs is not None:
                        latest_obs = (
                            latest_obs.year,
                            latest_obs.month,
                            latest_obs.day,
                        )  # put back into list because download_altimetry expects this
                        startdate = datetime.date(*latest_obs)

                # grab coordinates of geometry
                coords = [
                    (x, y)
                    for x, y in self.download_gdf.loc[
                        i, "geometry"
                    ].envelope.exterior.coords
                ]

                # make directiories for download, if needed
                if product == "S3":
                    download_dir = os.path.join(
                        self.dirs["sentinel3"], rf"{self.type}", rf"{id}"
                    )
                if product == "S6":
                    download_dir = os.path.join(
                        self.dirs["sentinel6"], rf"{self.type}", rf"{id}"
                    )
                general.ifnotmakedirs(download_dir)

                # query copernicus for download ids
                print(
                    f"Searching for Sentinel-{product} for aoi from {startdate} to {enddate}"
                )
                ids = sentinel.query(
                    aoi=coords,
                    startdate=startdate,
                    enddate=enddate,
                    product=product,
                    creodias_credentials=credentials,
                )

                # download granules that havent already been logged as downloaded
                session_token, session_start_time = sentinel.download(
                    ids,
                    download_directory=download_dir,
                    creodias_credentials=credentials,
                    token=session_token,
                    session_start_time=session_start_time,
                    threads=getattr(self, "mission_options", {})
                    .get("sentinel3" if product == "S3" else "sentinel6", {})
                    .get("download_threads", 1),
                )

                # once we have finished downloading all data for the aoi, we need to unzip the files keeping only .nc files
                general.unzip_dir_files_with_ext(
                    download_dir, download_dir, ".nc", show_progress=True
                )

                # subset the data so that we only keep what lies within the download geometry
                sentinel.subset(
                    aoi=coords,
                    download_dir=download_dir,
                    dest_dir=download_dir,
                    file_id=getattr(self, "mission_options", {})
                    .get("sentinel3" if product == "S3" else "sentinel6", {})
                    .get("subset_file_id", "enhanced_measurement.nc"),
                    product=product,
                    show_progress=True,
                )

                # clean up zip and unzipped folders keeping only the remaining subsetted data
                # general.remove_non_exts(download_dir, [".nc", ".log"]) # TODO: uncomment

    def get_unfiltered_product_timeseries(self, id, products: list = []):
        data_dir = os.path.join(self.dirs["output"], f"{id}", "raw_observations")

        # loop through raw directory and load in all raw observations
        df_list = list()
        for file in os.listdir(data_dir):
            if file.endswith(".shp"):
                # check if we process all products or only what is specified in the product list
                if (len(products) == 0) or (file.split(".")[0] in products):
                    gdf_path = os.path.join(data_dir, file)
                    gdf = gpd.read_file(gdf_path)
                    df = gdf.drop(columns=["geometry"])
                    df_list.append(df)

        # concatenate everything into one dataframe
        if len(df_list) > 0:
            df = pd.concat(df_list)
        else:
            df = None  # TODO: maybe raise an error instead?

        return df

    def clean_product_timeseries(
        self, products: list, filter_options_by_product: dict = None
    ):
        filter_options_by_product = filter_options_by_product or {}

        for id in tqdm(self.download_gdf[self.id_key]):
            for product in products:
                # get timeseries for id and each product to clean individually
                df = self.get_unfiltered_product_timeseries(id, [product])
                if df is not None:
                    product_options = filter_options_by_product.get(
                        product,
                        {
                            "processing_filters": ["elevation", "MAD"],
                            "elevation_min_m": 0.0,
                            "elevation_max_m": 8000.0,
                            "mad_threshold": 5.0,
                        },
                    )

                    # create a timeseries object
                    ts = timeseries.Timeseries(df, date_key="date", height_key="height")

                    # run all filters on timeseries
                    ts.clean(
                        product_options.get("processing_filters", ["elevation", "MAD"]),
                        filter_params={
                            "elevation_min_m": product_options.get(
                                "elevation_min_m", 0.0
                            ),
                            "elevation_max_m": product_options.get(
                                "elevation_max_m", 8000.0
                            ),
                            "mad_threshold": product_options.get("mad_threshold", 5.0),
                        },
                    )

                    # save filtered timeseries
                    export_dir = os.path.join(
                        self.dirs["output"], f"{id}", "cleaned_observations"
                    )
                    general.ifnotmakedirs(export_dir)
                    ts.export_csv(os.path.join(export_dir, f"{product}.csv"))

    def get_cleaned_product_timeseries(self, id, products: list = []):
        data_dir = os.path.join(self.dirs["output"], f"{id}", "cleaned_observations")

        if os.path.exists(data_dir):
            # loop through raw directory and load in all raw observations
            df_list = list()
            for file in os.listdir(data_dir):
                if file.endswith(".csv"):
                    # check if we process all products or only what is specified in the product list
                    if (len(products) == 0) or (file.split(".")[0] in products):
                        df_path = os.path.join(data_dir, file)
                        df = pd.read_csv(df_path)
                        df_list.append(df)

            # concatenate everything into one dataframe
            if len(df_list) > 0:
                df = pd.concat(df_list)
                df["date"] = pd.to_datetime(
                    df.date, format="mixed", utc=True
                ).dt.tz_convert(None)
                df = df.sort_values(by="date")
            else:
                df = None  # TODO: maybe raise an error instead?

        else:
            df = None

        return df

    def merge_product_timeseries(self, products: list):
        for id in tqdm(self.download_gdf[self.id_key]):
            if id not in []:
                # ["Lower Se San 2 + Lower Sre Pok 2"]:  # TODO: REMOVE THIS LINE!
                ts_list = list()
                for product in products:
                    # get timeseries for id and each product to clean individually
                    df = self.get_cleaned_product_timeseries(id, [product])
                    if df is not None:
                        ts_list.append(
                            timeseries.Timeseries(
                                df, date_key="date", height_key="height"
                            )
                        )

                # create a timeseries object of all available timeseries for id object
                if len(ts_list) > 0:
                    ts = timeseries.concat(ts_list)

                    data_dir = os.path.join(self.dirs["output"], f"{id}")
                    general.ifnotmakedirs(data_dir)
                    ts.export_csv(os.path.join(data_dir, "all_cleaned_timeseries.csv"))

                    # bias correct the timeseries
                    # ts.bias_correct()  # perhaps save intermediate step?

                    # run the merge function, runs linear SVR, MAD, Kalman filter and radial base SVR to get a final merged timeseries
                    # save the merged timeseries
                    data_dir = os.path.join(self.dirs["output"], f"{id}")
                    general.ifnotmakedirs(data_dir)
                    ts = ts.merge(
                        save_progress=True,
                        dir=os.path.join(data_dir, "merged_progress"),
                    )
                    ts.export_csv(os.path.join(data_dir, "merged_timeseries.csv"))

    def get_merged_timeseries(self, id):
        data_path = os.path.join(self.dirs["output"], f"{id}", "merged_timeseries.csv")

        if os.path.exists(data_path):
            df = pd.read_csv(data_path)
            df["date"] = pd.to_datetime(df.date)
            df = df.sort_values(by="date")
            return df
        else:
            warnings.warn(
                f"{data_path} does not exist, be sure to merge product timeseries first!"
            )
            return None

    def summarize_crossings_by_id(self, id, show=True, save=False):
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
        fig.suptitle(f"{self.type}: {id}")

        # extract item shape
        indx = self.gdf.loc[self.gdf[self.id_key] == id].index[0]

        # set boudns
        xmin, ymin, xmax, ymax = self.gdf.loc[indx, "geometry"].bounds

        ax.set_xlim([xmin - 0.1, xmax + 0.1])
        ax.set_ylim([ymin - 0.1, ymax + 0.1])

        # plot reservoir
        self.gdf.loc[[indx]].plot(
            ax=ax,
            edgecolor="black",
            facecolor="None",
            zorder=5,
            label="Reservoir outline",
        )

        # loop through each product in file and plot
        data_dir = os.path.join(self.dirs["output"], f"{id}", "raw_observations")
        plotted_products = set()
        for file_name in os.listdir(data_dir):
            if file_name.endswith(".shp"):
                path_to_file = os.path.join(data_dir, file_name)
                product = file_name.split(".")[0]
                gdf = gpd.read_file(path_to_file)
                if product == "swot":
                    zorder = 0
                    alpha = 0.1
                else:
                    zorder = 10
                    alpha = 0.5
                gdf.plot(
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
            plt.savefig(
                os.path.join(self.dirs["output"], f"{id}", "crossing_summary.png")
            )
        if show:
            plt.show()

        return None

    def summarize_cleaning_by_id(self, id, show=True, save=False):
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
        fig.suptitle(f"{self.type}: {id}")

        # plot unfiltered timeseries
        df = self.get_unfiltered_product_timeseries(id)
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
            df = self.get_cleaned_product_timeseries(id)
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
        df = self.get_merged_timeseries(id)
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
                    os.path.join(self.dirs["output"], f"{id}", "cleaning_summary.png")
                )
            if show:
                plt.show()

    def summarize_merging_by_id(self, id, show=True, save=False):
        merged_dir = os.path.join(self.dirs["output"], str(id), "merged_progress")

        if os.path.exists(merged_dir):
            file_names = os.listdir(merged_dir)
            num_files = len(file_names)

            sns.set()
            fig, main_ax = plt.subplots(num_files, 1, figsize=(10, 10))
            fig.suptitle(f"{self.type}: {id}")

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
            i = _plot_file(
                i, fig, main_ax, "daily_mad_error.csv", file_names, merged_dir
            )
            i = _plot_file(i, fig, main_ax, "kalman.csv", file_names, merged_dir)
            i = _plot_file(i, fig, main_ax, "svr_radial.csv", file_names, merged_dir)

            fig.tight_layout()
            if save:
                plt.savefig(
                    os.path.join(self.dirs["output"], f"{id}", "merging_summary.png")
                )
            if show:
                plt.show()

        return


@dataclass
class Reservoirs(System):
    def __post_init__(self):
        self.type = "reservoirs"

        self.dirs["output"] = os.path.join(self.dirs["main"], self.type)
        general.ifnotmakedirs(self.dirs["output"])

        self.geom_type = self.gdf.loc[0, "geometry"].geom_type

    def download_pld(self, overwrite=False):
        pld_path = self.dirs["pld"]

        # determine if we need to download or simply load the pld
        if (not os.path.exists(pld_path)) or (overwrite):
            print("Downloading PLD")
            download_dir = os.path.dirname(pld_path)
            file_name = os.path.basename(pld_path)
            bounds = list(self.gdf.unary_union.bounds)
            hydroweb.download_PLD(
                download_dir=download_dir, file_name=file_name, bounds=bounds
            )
        else:
            print("PLD located")

    def assign_pld_id(self, local_crs, max_distance):
        # load the pld
        pld = gpd.read_file(self.dirs["pld"])

        # perform the spatial join
        joined_gdf = gpd.sjoin_nearest(
            self.gdf.to_crs(local_crs),
            pld.to_crs(local_crs),
            how="left",
            max_distance=max_distance,
            distance_col="dist_to_pld",
        )
        joined_gdf = joined_gdf.to_crs(self.gdf.crs)

        # rename the joined columns to keep it clear
        joined_gdf = joined_gdf.rename(
            columns={"lake_id": "prior_lake_id", "res_id": "prior_res_id"}
        )

        # map null lake ids to a negative number -9999
        joined_gdf.loc[joined_gdf.prior_lake_id.isnull(), "prior_lake_id"] = -9999

        # reset the reservoir data frame to include the changes
        self.gdf = joined_gdf

    def flag_missing_priors(
        self,
    ):  # simple function to report what reservoirs do not have PLD shapes
        present = self.gdf.loc[self.gdf.prior_lake_id > 0].reset_index(drop=True)
        missing = self.gdf.loc[self.gdf.prior_lake_id < 0].reset_index(drop=True)

        # ESRI Shapefile limits field names to 10 chars. Rename explicitly so
        # writes are deterministic and avoid pyogrio truncation warnings.
        shp_field_map = {
            "index_right": "idx_right",
            "prior_lake_id": "prior_lake",
            "prior_res_id": "prior_res",
            "dist_to_pld": "dist_pld",
        }
        present_out = present.rename(columns=shp_field_map)
        missing_out = missing.rename(columns=shp_field_map)

        present_out.to_file(os.path.join(self.dirs["output"], "present_in_pld.shp"))
        missing_out.to_file(os.path.join(self.dirs["output"], "missing_in_pld.shp"))

        print(
            f"Out of the {len(self.gdf)} reservoirs, {len(present)} are present and {len(missing)} are missing from the PLD."
        )

    def extract_product_timeseries(self, products: list):
        if "icesat2" in products:
            for id in tqdm(
                self.download_gdf[self.id_key], desc="Extracting ICESat-2 ATL13 product"
            ):
                # filter gdf to only row with id
                sub_gdf = self.download_gdf.loc[self.download_gdf[self.id_key] == id]

                # prep export folder and paths
                download_dir = os.path.join(
                    self.dirs["icesat2"], rf"{self.type}", rf"{id}"
                )
                if os.path.exists(download_dir):
                    dst_dir = os.path.join(
                        self.dirs["output"], f"{id}", "raw_observations"
                    )
                    general.ifnotmakedirs(dst_dir)
                    dst_path = os.path.join(dst_dir, "icesat2.shp")

                    # extract observations within bounds and save as timeseries csv
                    icesat2.extract_observations(
                        src_dir=download_dir,
                        dst_path=dst_path,
                        features=sub_gdf,
                        atl13_fields=getattr(self, "mission_options", {})
                        .get("icesat2", {})
                        .get("atl13_fields"),
                        track_keys=getattr(self, "mission_options", {})
                        .get("icesat2", {})
                        .get("track_keys"),
                    )

        if "sentinel3" in products:
            for id in tqdm(
                self.download_gdf[self.id_key], desc="Extracting Sentinel-3 product"
            ):
                # filter gdf to only row with id
                sub_gdf = self.download_gdf.loc[self.download_gdf[self.id_key] == id]

                # prep export folder and paths
                download_dir = os.path.join(
                    self.dirs["sentinel3"], rf"{self.type}", rf"{id}"
                )
                if os.path.exists(download_dir):
                    dst_dir = os.path.join(
                        self.dirs["output"], f"{id}", "raw_observations"
                    )
                    general.ifnotmakedirs(dst_dir)
                    dst_path = os.path.join(dst_dir, "sentinel3.shp")

                    # extract observations within bounds and save as timeseries csv
                    sentinel.extract_observations(
                        src_dir=download_dir,
                        dst_path=dst_path,
                        features=sub_gdf,
                        sigma0_max=getattr(self, "mission_options", {})
                        .get("sentinel3", {})
                        .get("sigma0_max", 1e5),
                    )

        if "sentinel6" in products:
            for id in tqdm(
                self.download_gdf[self.id_key], desc="Extracting Sentinel-6 product"
            ):
                # filter gdf to only row with id
                sub_gdf = self.download_gdf.loc[self.download_gdf[self.id_key] == id]

                # prep export folder and paths
                download_dir = os.path.join(
                    self.dirs["sentinel6"], rf"{self.type}", rf"{id}"
                )
                if os.path.exists(download_dir):
                    dst_dir = os.path.join(
                        self.dirs["output"], f"{id}", "raw_observations"
                    )
                    general.ifnotmakedirs(dst_dir)
                    dst_path = os.path.join(dst_dir, "sentinel6.shp")

                    # extract observations within bounds and save as timeseries csv
                    sentinel.extract_observations(
                        src_dir=download_dir,
                        dst_path=dst_path,
                        features=sub_gdf,
                        sigma0_max=getattr(self, "mission_options", {})
                        .get("sentinel6", {})
                        .get("sigma0_max", 1e5),
                    )

        if "swot" in products:
            download_dir = os.path.join(self.dirs["swot"], rf"{self.type}")

            # extract observations within bounds and save as timeseries csv, slightly differnt format for downloading here due to the natrure of swot data
            swot.extract_observations(
                src_dir=download_dir,
                dst_dir=self.dirs["output"],
                dst_file_name="swot.shp",
                features=self.download_gdf,
                id_key=self.id_key,
                exclude_obs_id_values=getattr(self, "mission_options", {})
                .get("swot", {})
                .get("exclude_obs_id_values", ["no_data"]),
            )
