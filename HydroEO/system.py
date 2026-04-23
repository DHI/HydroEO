from dataclasses import dataclass
import logging
import warnings

import os
import pandas as pd
import geopandas as gpd
import datetime

from HydroEO.satellites import swot, icesat2, sentinel
from HydroEO.utils import general, timeseries
from HydroEO.downloaders import hydroweb
from HydroEO import plotting

from tqdm import tqdm

logger = logging.getLogger(__name__)


class HydroEODownloadError(RuntimeError):
    """Raised when a HydroEO download routine fails or returns no usable data."""


@dataclass
class System:
    gdf: gpd.GeoDataFrame
    id_key: str
    dirs: dict

    def report(self):
        logger.info("Number of %s: %s", self.type, len(self.gdf))
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
            logger.info(
                "Searching for %s for aoi from %s to %s",
                swot.SWOT_LAKE_SHORT_NAME,
                startdate,
                enddate,
            )
            results = swot.query(
                aoi=coords,
                startdate=startdate,
                enddate=enddate,
            )

            # Filter to the prior-lake sub-collection granules.
            # Baseline D file names end with "_prior_*" (lower-case); baseline C
            # used "_Prior_*" (title-case). We match both for robustness.
            logger.info("%s products returned from query", len(results))
            to_download = list()
            for result in results:
                link = result.data_links()[0]
                filename = link.split("/")[-1].lower()
                if "_prior_" in filename or "prior" in filename.split("_"):
                    to_download.append(result)
            logger.info(
                "%s prior-lake granules selected for download", len(to_download)
            )

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
                logger.info("Downloading data for id %s", id)

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

                # Grab coordinates of the actual reservoir geometry (not the envelope)
                # so the polygon centroid passed to the SlideRule AMS is accurate.
                geom = self.download_gdf.loc[i, "geometry"]
                if hasattr(geom, "geoms"):
                    geom = geom.geoms[0]  # MultiPolygon -> use first polygon
                coords = list(geom.exterior.coords)

                # define and if needed create directory for each download geometry
                download_dir = os.path.join(
                    self.dirs["icesat2"], rf"{self.type}", rf"{id}"
                )
                general.ifnotmakedirs(download_dir)

                # query and download data via SlideRule atl13x
                logger.info(
                    "Searching for Icesat2 ATL13 for aoi from %s to %s",
                    startdate,
                    enddate,
                )
                try:
                    _ = icesat2.query(
                        aoi=coords,
                        startdate=startdate,
                        enddate=enddate,
                        download_directory=download_dir,
                        atl13_options=getattr(self, "mission_options", {})
                        .get("icesat2", {})
                        .get("atl13", {}),
                        atl13_fields=getattr(self, "mission_options", {})
                        .get("icesat2", {})
                        .get("atl13_fields")
                        or None,
                    )
                except HydroEODownloadError as exc:
                    logger.warning("ICESat-2 download skipped for %s: %s", id, exc)

        elif product in ["S3", "S6"]:
            # set empty variables because no session has been started yet, this will occur at the first download and sessions will refresh automatically
            session_token = None
            session_start_time = None

            # loop through download geometry
            for i in self.download_gdf.index:
                # extract id for saving data
                id = self.download_gdf.loc[i, self.id_key]
                logger.info("Downloading data for id %s", id)

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
                logger.info(
                    "Searching for Sentinel-%s for aoi from %s to %s",
                    product,
                    startdate,
                    enddate,
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

    def get_unfiltered_product_timeseries(self, id, products: list = []):
        data_dir = os.path.join(self.dirs["output"], f"{id}", "raw_observations")

        if not os.path.exists(data_dir):
            return None

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

        ids_with_raw = [
            id
            for id in self.download_gdf[self.id_key]
            if os.path.exists(
                os.path.join(self.dirs["output"], f"{id}", "raw_observations")
            )
        ]
        if not ids_with_raw:
            logger.warning(
                "No raw observations found for any reservoir; skipping timeseries cleaning."
            )
            return

        for id in tqdm(ids_with_raw, desc="Cleaning product timeseries"):
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
        ids_with_cleaned = [
            id
            for id in self.download_gdf[self.id_key]
            if os.path.exists(
                os.path.join(self.dirs["output"], f"{id}", "cleaned_observations")
            )
        ]
        if not ids_with_cleaned:
            logger.warning(
                "No cleaned observations found for any reservoir; skipping timeseries merging."
            )
            return

        for id in tqdm(ids_with_cleaned, desc="Merging product timeseries"):
            if id not in []:
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
        """Plot raw observations crossing for a reservoir."""
        return plotting.plot_crossings(
            gdf=self.gdf,
            id_key=self.id_key,
            reservoir_id=id,
            output_dir=self.dirs["output"],
            reservoir_type=self.type,
            show=show,
            save=save,
        )

    def summarize_cleaning_by_id(self, id, show=True, save=False):
        """Plot cleaning progression (unfiltered, cleaned, merged timeseries)."""
        return plotting.plot_cleaning(
            reservoir_id=id,
            output_dir=self.dirs["output"],
            get_unfiltered_fn=self.get_unfiltered_product_timeseries,
            get_cleaned_fn=self.get_cleaned_product_timeseries,
            get_merged_fn=self.get_merged_timeseries,
            reservoir_type=self.type,
            show=show,
            save=save,
        )

    def summarize_merging_by_id(self, id, show=True, save=False):
        """Plot merging progression showing intermediate processing steps."""
        return plotting.plot_merging(
            reservoir_id=id,
            output_dir=self.dirs["output"],
            reservoir_type=self.type,
            show=show,
            save=save,
        )


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
            logger.info("Downloading PLD")
            download_dir = os.path.dirname(pld_path)
            file_name = os.path.basename(pld_path)
            bounds = list(self.gdf.unary_union.bounds)
            hydroweb.download_PLD(
                download_dir=download_dir, file_name=file_name, bounds=bounds
            )
        else:
            logger.info("PLD located")

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

        logger.info(
            "Out of the %s reservoirs, %s are present and %s are missing from the PLD.",
            len(self.gdf),
            len(present),
            len(missing),
        )

    def extract_product_timeseries(self, products: list):
        if "icesat2" in products:
            available_ids = [
                id
                for id in self.download_gdf[self.id_key]
                if os.path.exists(
                    os.path.join(
                        self.dirs["icesat2"], rf"{self.type}", rf"{id}", "atl13.parquet"
                    )
                )
            ]
            if not available_ids:
                logger.warning(
                    "No ICESat-2 downloads found; skipping timeseries extraction."
                )
            else:
                empty_ids = []
                for id in tqdm(available_ids, desc="Extracting ICESat-2 ATL13 product"):
                    sub_gdf = self.download_gdf.loc[
                        self.download_gdf[self.id_key] == id
                    ]
                    download_dir = os.path.join(
                        self.dirs["icesat2"], rf"{self.type}", rf"{id}"
                    )
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

                    if not os.path.exists(dst_path):
                        empty_ids.append(id)

                if empty_ids:
                    logger.warning(
                        "ICESat-2 timeseries empty for: %s (no observations passed the spatial filter or the download returned no data)",
                        ", ".join(str(i) for i in empty_ids),
                    )

        if "sentinel3" in products:
            available_ids = [
                id
                for id in self.download_gdf[self.id_key]
                if os.path.exists(
                    os.path.join(self.dirs["sentinel3"], rf"{self.type}", rf"{id}")
                )
            ]
            if not available_ids:
                logger.warning(
                    "No Sentinel-3 downloads found; skipping timeseries extraction."
                )
            else:
                empty_ids = []
                for id in tqdm(available_ids, desc="Extracting Sentinel-3 product"):
                    sub_gdf = self.download_gdf.loc[
                        self.download_gdf[self.id_key] == id
                    ]
                    download_dir = os.path.join(
                        self.dirs["sentinel3"], rf"{self.type}", rf"{id}"
                    )
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

                    if not os.path.exists(dst_path):
                        empty_ids.append(id)

                if empty_ids:
                    logger.warning(
                        "Sentinel-3 timeseries empty for: %s (no observations passed the spatial filter or the download returned no data)",
                        ", ".join(str(i) for i in empty_ids),
                    )

        if "sentinel6" in products:
            available_ids = [
                id
                for id in self.download_gdf[self.id_key]
                if os.path.exists(
                    os.path.join(self.dirs["sentinel6"], rf"{self.type}", rf"{id}")
                )
            ]
            if not available_ids:
                logger.warning(
                    "No Sentinel-6 downloads found; skipping timeseries extraction."
                )
            else:
                empty_ids = []
                for id in tqdm(available_ids, desc="Extracting Sentinel-6 product"):
                    sub_gdf = self.download_gdf.loc[
                        self.download_gdf[self.id_key] == id
                    ]
                    download_dir = os.path.join(
                        self.dirs["sentinel6"], rf"{self.type}", rf"{id}"
                    )
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

                    if not os.path.exists(dst_path):
                        empty_ids.append(id)

                if empty_ids:
                    logger.warning(
                        "Sentinel-6 timeseries empty for: %s (no observations passed the spatial filter or the download returned no data)",
                        ", ".join(str(i) for i in empty_ids),
                    )

        if "swot" in products:
            download_dir = os.path.join(self.dirs["swot"], rf"{self.type}")
            if not os.path.exists(download_dir):
                logger.warning(
                    "No SWOT downloads found; skipping timeseries extraction."
                )
            else:
                # extract observations within bounds and save as timeseries csv, slightly differnt format for downloading here due to the natrure of swot data
                empty_ids = swot.extract_observations(
                    src_dir=download_dir,
                    dst_dir=self.dirs["output"],
                    dst_file_name="swot.shp",
                    features=self.download_gdf,
                    id_key=self.id_key,
                    exclude_obs_id_values=getattr(self, "mission_options", {})
                    .get("swot", {})
                    .get("exclude_obs_id_values", ["no_data"]),
                )
                if empty_ids:
                    logger.warning(
                        "SWOT timeseries empty for: %s (no observations matched the prior lake ID or all were excluded)",
                        ", ".join(str(i) for i in empty_ids),
                    )
