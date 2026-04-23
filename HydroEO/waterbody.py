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
class WaterBody:
    gdf: gpd.GeoDataFrame
    id_key: str
    dirs: dict

    def report(self):
        logger.info("Number of %s: %s", self.type, len(self.gdf))
        return self.gdf.head()

    def _advance_startdate_if_existing(
        self, provided_start_date, update_existing, latest_obs
    ):
        """Return latest_obs date when update_existing is True and latest_obs is available."""
        if update_existing and latest_obs is not None:
            return datetime.date(latest_obs.year, latest_obs.month, latest_obs.day)
        return provided_start_date

    def _download_swot(self, startdate, enddate, update_existing):
        download_dir = os.path.join(self.dirs["swot"], rf"{self.type}")
        general.ifnotmakedirs(download_dir)

        latest_obs = (
            swot.get_latest_obs_date(self.dirs["output"]) if update_existing else None
        )
        startdate = self._advance_startdate_if_existing(
            startdate, update_existing, latest_obs
        )

        coords = [
            (x, y) for x, y in self.download_gdf.unary_union.envelope.exterior.coords
        ]

        logger.info(
            "Searching for %s for aoi from %s to %s",
            swot.SWOT_LAKE_SHORT_NAME,
            startdate,
            enddate,
        )
        results = swot.query(aoi=coords, startdate=startdate, enddate=enddate)

        # Filter to the prior-lake sub-collection granules.
        # Baseline D file names end with "_prior_*" (lower-case); baseline C
        # used "_Prior_*" (title-case). We match both for robustness.
        logger.info("%s products returned from query", len(results))
        to_download = []
        for result in results:
            link = result.data_links()[0]
            filename = link.split("/")[-1].lower()
            if "_prior_" in filename or "prior" in filename.split("_"):
                to_download.append(result)
        logger.info("%s prior-lake granules selected for download", len(to_download))

        _ = swot.download(to_download, download_directory=download_dir)

        files_in_dir = [
            os.path.join(download_dir, f)
            for f in os.listdir(download_dir)
            if f.endswith(".zip")
        ]
        swot.subset_by_id(
            files_in_dir, self.download_gdf["prior_lake_id"].astype(int).values
        )

    def _download_icesat2(self, startdate, enddate, update_existing):
        provided_start_date = startdate
        for i in self.download_gdf.index:
            id = self.download_gdf.loc[i, self.id_key]
            logger.info("Downloading data for id %s", id)

            raw_data_dir = os.path.join(
                self.dirs["output"], str(id), "raw_observations"
            )
            latest_obs = (
                icesat2.get_latest_obs_date(raw_data_dir) if update_existing else None
            )
            startdate = self._advance_startdate_if_existing(
                provided_start_date, update_existing, latest_obs
            )

            # Grab coordinates of the actual reservoir geometry (not the envelope)
            # so the polygon centroid passed to the SlideRule AMS is accurate.
            geom = self.download_gdf.loc[i, "geometry"]
            if hasattr(geom, "geoms"):
                geom = geom.geoms[0]  # MultiPolygon -> use first polygon
            coords = list(geom.exterior.coords)

            download_dir = os.path.join(self.dirs["icesat2"], rf"{self.type}", rf"{id}")
            general.ifnotmakedirs(download_dir)

            logger.info(
                "Searching for Icesat2 ATL13 for aoi from %s to %s", startdate, enddate
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

    def _download_sentinel(
        self, product, startdate, enddate, credentials, update_existing
    ):
        provided_start_date = startdate
        session_token = None
        session_start_time = None
        dir_key = "sentinel3" if product == "S3" else "sentinel6"

        for i in self.download_gdf.index:
            id = self.download_gdf.loc[i, self.id_key]
            logger.info("Downloading data for id %s", id)

            raw_data_dir = os.path.join(
                self.dirs["output"], str(id), "raw_observations"
            )
            latest_obs = (
                sentinel.get_latest_obs_date(raw_data_dir, product)
                if update_existing
                else None
            )
            startdate = self._advance_startdate_if_existing(
                provided_start_date, update_existing, latest_obs
            )

            coords = [
                (x, y)
                for x, y in self.download_gdf.loc[
                    i, "geometry"
                ].envelope.exterior.coords
            ]

            download_dir = os.path.join(self.dirs[dir_key], rf"{self.type}", rf"{id}")
            general.ifnotmakedirs(download_dir)

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

            session_token, session_start_time = sentinel.download(
                ids,
                download_directory=download_dir,
                creodias_credentials=credentials,
                token=session_token,
                session_start_time=session_start_time,
                threads=getattr(self, "mission_options", {})
                .get(dir_key, {})
                .get("download_threads", 1),
            )

            general.unzip_dir_files_with_ext(
                download_dir, download_dir, ".nc", show_progress=True
            )

            sentinel.subset(
                aoi=coords,
                download_dir=download_dir,
                dest_dir=download_dir,
                file_id=getattr(self, "mission_options", {})
                .get(dir_key, {})
                .get("subset_file_id", "enhanced_measurement.nc"),
                product=product,
                show_progress=True,
            )

    def download_altimetry(
        self,
        product: str,
        startdate: tuple,
        enddate: tuple,
        credentials=None,
        grid: bool = False,
        update_existing: bool = False,
    ):
        supported_products = ["ATL13", "S3", "S6", "SWOT_LAKE"]

        product = product.upper()
        if product not in supported_products:
            raise ValueError(
                f'"{product}" is not accepted as a valid download product. Please provide a valid product.'
            )

        startdate = datetime.date(*startdate)
        enddate = datetime.date(*enddate)

        if product == "SWOT_LAKE":
            self._download_swot(startdate, enddate, update_existing)
        elif product == "ATL13":
            self._download_icesat2(startdate, enddate, update_existing)
        elif product in ["S3", "S6"]:
            self._download_sentinel(
                product, startdate, enddate, credentials, update_existing
            )

    def _load_dataframes_from_dir(self, data_dir, ext, products, reader_fn):
        """Load files of a given extension from a directory, optionally filtered to named products."""
        if not os.path.exists(data_dir):
            return None
        df_list = []
        for file in os.listdir(data_dir):
            if file.endswith(ext):
                if not products or file.split(".")[0] in products:
                    df_list.append(reader_fn(os.path.join(data_dir, file)))
        return pd.concat(df_list) if df_list else None

    def get_unfiltered_product_timeseries(self, id, products: list = []):
        data_dir = os.path.join(self.dirs["output"], f"{id}", "raw_observations")
        return self._load_dataframes_from_dir(
            data_dir,
            ".shp",
            products,
            lambda path: gpd.read_file(path).drop(columns=["geometry"]),
        )

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
        df = self._load_dataframes_from_dir(data_dir, ".csv", products, pd.read_csv)
        if df is not None:
            df["date"] = pd.to_datetime(
                df.date, format="mixed", utc=True
            ).dt.tz_convert(None)
            df = df.sort_values(by="date")
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
            ts_list = []
            for product in products:
                df = self.get_cleaned_product_timeseries(id, [product])
                if df is not None:
                    ts_list.append(
                        timeseries.Timeseries(df, date_key="date", height_key="height")
                    )

            if len(ts_list) > 0:
                ts = timeseries.concat(ts_list)

                data_dir = os.path.join(self.dirs["output"], f"{id}")
                general.ifnotmakedirs(data_dir)
                ts.export_csv(os.path.join(data_dir, "all_cleaned_timeseries.csv"))

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
class Reservoirs(WaterBody):
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

    def _extract_per_reservoir_observations(
        self,
        dir_key,
        availability_check_fn,
        dst_filename,
        extract_fn,
        tqdm_desc,
        mission_label,
    ):
        """Find available IDs, loop with tqdm, extract observations, and warn about empty results."""
        available_ids = [
            id for id in self.download_gdf[self.id_key] if availability_check_fn(id)
        ]
        if not available_ids:
            logger.warning(
                "No %s downloads found; skipping timeseries extraction.", mission_label
            )
            return

        empty_ids = []
        for id in tqdm(available_ids, desc=tqdm_desc):
            sub_gdf = self.download_gdf.loc[self.download_gdf[self.id_key] == id]
            download_dir = os.path.join(self.dirs[dir_key], rf"{self.type}", rf"{id}")
            dst_dir = os.path.join(self.dirs["output"], f"{id}", "raw_observations")
            general.ifnotmakedirs(dst_dir)
            dst_path = os.path.join(dst_dir, dst_filename)
            extract_fn(download_dir=download_dir, dst_path=dst_path, features=sub_gdf)
            if not os.path.exists(dst_path):
                empty_ids.append(id)

        if empty_ids:
            logger.warning(
                "%s timeseries empty for: %s (no observations passed the spatial filter or the download returned no data)",
                mission_label,
                ", ".join(str(i) for i in empty_ids),
            )

    def extract_product_timeseries(self, products: list):
        if "icesat2" in products:
            self._extract_per_reservoir_observations(
                dir_key="icesat2",
                availability_check_fn=lambda id: os.path.exists(
                    os.path.join(
                        self.dirs["icesat2"], rf"{self.type}", rf"{id}", "atl13.parquet"
                    )
                ),
                dst_filename="icesat2.shp",
                extract_fn=lambda download_dir, dst_path, features: (
                    icesat2.extract_observations(
                        src_dir=download_dir,
                        dst_path=dst_path,
                        features=features,
                        atl13_fields=getattr(self, "mission_options", {})
                        .get("icesat2", {})
                        .get("atl13_fields"),
                        track_keys=getattr(self, "mission_options", {})
                        .get("icesat2", {})
                        .get("track_keys"),
                    )
                ),
                tqdm_desc="Extracting ICESat-2 ATL13 product",
                mission_label="ICESat-2",
            )

        if "sentinel3" in products:
            self._extract_per_reservoir_observations(
                dir_key="sentinel3",
                availability_check_fn=lambda id: os.path.exists(
                    os.path.join(self.dirs["sentinel3"], rf"{self.type}", rf"{id}")
                ),
                dst_filename="sentinel3.shp",
                extract_fn=lambda download_dir, dst_path, features: (
                    sentinel.extract_observations(
                        src_dir=download_dir,
                        dst_path=dst_path,
                        features=features,
                        sigma0_max=getattr(self, "mission_options", {})
                        .get("sentinel3", {})
                        .get("sigma0_max", 1e5),
                    )
                ),
                tqdm_desc="Extracting Sentinel-3 product",
                mission_label="Sentinel-3",
            )

        if "sentinel6" in products:
            self._extract_per_reservoir_observations(
                dir_key="sentinel6",
                availability_check_fn=lambda id: os.path.exists(
                    os.path.join(self.dirs["sentinel6"], rf"{self.type}", rf"{id}")
                ),
                dst_filename="sentinel6.shp",
                extract_fn=lambda download_dir, dst_path, features: (
                    sentinel.extract_observations(
                        src_dir=download_dir,
                        dst_path=dst_path,
                        features=features,
                        sigma0_max=getattr(self, "mission_options", {})
                        .get("sentinel6", {})
                        .get("sigma0_max", 1e5),
                    )
                ),
                tqdm_desc="Extracting Sentinel-6 product",
                mission_label="Sentinel-6",
            )

        if "swot" in products:
            download_dir = os.path.join(self.dirs["swot"], rf"{self.type}")
            if not os.path.exists(download_dir):
                logger.warning(
                    "No SWOT downloads found; skipping timeseries extraction."
                )
            else:
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
