"""The project module"""

from dataclasses import dataclass
from copy import deepcopy
import logging

import os
import yaml
import warnings

import geopandas as gpd
import datetime

from HydroEO.waterbody import Reservoirs, Rivers
from HydroEO import flows
from HydroEO.satellites.swot.raster import download_raster
from HydroEO.satellites.swot.pixc import download_pixc
from HydroEO.utils import general
from HydroEO.constants import MISSION_DEFAULTS
from HydroEO.validation import validate_config

logger = logging.getLogger(__name__)


@dataclass
class Project:
    """Create a new altimetry project from configuration file

    Parameters
    ----------
    name: str
        name of altimetry project
    config: str
        path to configuration file

    Examples
    --------
    >>> proj = Project(name="Lower Mekong", config='./data/lower_mekong.yaml')
    """

    name: str
    config: str

    def __post_init__(self):
        self.to_download = list()
        self.to_process = list()
        self.mission_options = dict()
        self.processing_options = dict()

        self.dirs = dict()
        self.startdates = dict()
        self.enddates = dict()

        ### Load in the config file and extract parameters
        with open(self.config, "rt") as f:
            self.config = yaml.safe_load(f.read())

        if self.config is None:
            self.config = {}

        self._apply_optional_defaults()
        validate_config(self.config)

        ### Define the project directory for saving outputs, etc
        if "project" in self.config.keys():
            self.dirs["main"] = self.config["project"]["main_dir"]
            general.ifnotmakedirs(self.dirs["main"])
        else:
            raise Warning("Project directory must be defined within configuration file")

        ##### Set credentials and access keys
        hydroweb_cfg = self.config.get("hydroweb", {})
        hydroweb_api_key = hydroweb_cfg.get("api_key") or os.environ.get(
            "HYDROWEB_API_KEY"
        )
        if hydroweb_api_key:
            os.environ["HYDROWEB_API_KEY"] = hydroweb_api_key

        # Set PLD paths and configuration
        self.dirs["pld"] = os.path.join(
            self.dirs["main"], "aux", "PLD", "PLD_subset.gpkg"
        )
        if "raw_pld_path" in hydroweb_cfg:
            self.dirs["pld_raw"] = hydroweb_cfg["raw_pld_path"]
        self.keep_raw_pld = hydroweb_cfg.get("keep_raw_pld", False)

        earthaccess_cfg = self.config.get("earthaccess", {})
        self.earthdata_user = (
            earthaccess_cfg.get("username")
            or os.environ.get("EARTHDATA_USERNAME")
            or os.environ.get("EDL_USERNAME")
        )
        self.earthdata_pass = (
            earthaccess_cfg.get("password")
            or os.environ.get("EARTHDATA_PASSWORD")
            or os.environ.get("EDL_PASSWORD")
        )

        if self.earthdata_user:
            os.environ["EARTHDATA_USERNAME"] = self.earthdata_user
            os.environ["EDL_USERNAME"] = self.earthdata_user
        if self.earthdata_pass:
            os.environ["EARTHDATA_PASSWORD"] = self.earthdata_pass
            os.environ["EDL_PASSWORD"] = self.earthdata_pass

        earthdata_token = (
            earthaccess_cfg.get("token")
            or os.environ.get("EARTHDATA_TOKEN")
            or os.environ.get("EDL_TOKEN")
        )
        if earthdata_token:
            os.environ["EARTHDATA_TOKEN"] = earthdata_token
            os.environ["EDL_TOKEN"] = earthdata_token

        creodias_cfg = self.config.get("creodias", {})
        self.creodias_user = creodias_cfg.get("username") or os.environ.get(
            "CREODIAS_USERNAME"
        )
        self.creodias_pass = creodias_cfg.get("password") or os.environ.get(
            "CREODIAS_PASSWORD"
        )

        if self.creodias_user:
            os.environ["CREODIAS_USERNAME"] = self.creodias_user
        if self.creodias_pass:
            os.environ["CREODIAS_PASSWORD"] = self.creodias_pass

        ##### Set spatialite folder in PATH if provided
        spatialite_folder = self.config.get("spatialite_folder")
        if spatialite_folder:
            os.environ["PATH"] = spatialite_folder + os.pathsep + os.environ["PATH"]

        ##### Set attributes for each satellite to be downloaded or processed
        self.__sat_init("swot")
        self.__sat_init("icesat2")
        self.__sat_init("sentinel3")
        self.__sat_init("sentinel6")

        ### Set the crs from the congiguration file if possible
        self.global_crs = "EPSG:4326"
        self.local_crs = None
        if "gis" in self.config.keys():
            if "global_crs" in self.config["gis"].keys():
                self.global_crs = self.config["gis"]["global_crs"]

            if "local_crs" in self.config["gis"].keys():
                self.local_crs = self.config["gis"]["local_crs"]

        ### load in elements for download and processing
        if "reservoirs" in self.config.keys() and self.config["reservoirs"].get(
            "enabled", True
        ):
            self.reservoirs = Reservoirs(
                gdf=gpd.read_file(self.config["reservoirs"]["path"]),
                id_key=self.config["reservoirs"]["id_key"],
                dirs=self.dirs,
            )

            self.reservoirs.gdf = self.reservoirs.gdf.to_crs(self.global_crs)
            self.reservoirs.mission_options = self.mission_options
            self.reservoirs.processing_options = self.processing_options
            self.reservoirs.export_to_dfs0 = self.config["reservoirs"].get(
                "export_to_dfs0", False
            )

        if "rivers" in self.config.keys() and self.config["rivers"].get(
            "enabled", True
        ):
            rivers_cfg = self.config["rivers"]

            rivers_aoi_gdf = None
            if rivers_cfg.get("aoi_path"):
                rivers_aoi_gdf = gpd.read_file(rivers_cfg["aoi_path"])
                if rivers_aoi_gdf.crs is None:
                    rivers_aoi_gdf = rivers_aoi_gdf.set_crs(self.global_crs)
                else:
                    rivers_aoi_gdf = rivers_aoi_gdf.to_crs(self.global_crs)

                rivers_id_key = rivers_cfg["id_key"]
                rivers_gdf = gpd.GeoDataFrame(
                    {"geometry": []}, geometry="geometry", crs=self.global_crs
                )
                input_mode = "aoi_path"
                target_ids = []
                target_id_col = (
                    "node_id"
                    if rivers_cfg.get("feature_type") == "nodes"
                    else "reach_id"
                )
            else:
                id_values = rivers_cfg.get("feature_numbers", [])
                rivers_id_key = rivers_cfg.get("id_key", "river_id")
                if rivers_cfg.get("feature_type") == "nodes":
                    input_mode = "feature_numbers"
                    target_id_col = "node_id"
                else:
                    input_mode = "feature_numbers"
                    target_id_col = "reach_id"

                target_ids = [int(value) for value in id_values]
                rivers_gdf = gpd.GeoDataFrame(
                    {"geometry": []}, geometry="geometry", crs=self.global_crs
                )

            self.rivers = Rivers(
                gdf=rivers_gdf,
                id_key=rivers_id_key,
                dirs=self.dirs,
            )

            if self.rivers.gdf.crs is None:
                self.rivers.gdf = self.rivers.gdf.set_crs(self.global_crs)
            else:
                self.rivers.gdf = self.rivers.gdf.to_crs(self.global_crs)

            self.rivers.mission_options = self.mission_options
            self.rivers.processing_options = self.processing_options
            self.rivers.input_mode = input_mode
            self.rivers.aoi_gdf = rivers_aoi_gdf
            self.rivers.continent_key = rivers_cfg.get("continent_key")
            self.rivers.feature_type = rivers_cfg.get("feature_type")
            self.rivers.buffer_meters = rivers_cfg.get("buffer_meters")
            self.rivers.configured_id = rivers_cfg.get("id")
            self.rivers.target_id_col = target_id_col
            self.rivers.target_ids = target_ids

        if "swot_raster" in self.config.keys() and self.config["swot_raster"].get(
            "enabled", True
        ):
            # Store the SWOT raster config for later use in download/preprocess
            # Will be instantiated in download() when needed
            self.swot_raster_config = self.config["swot_raster"]

        if "swot_pixc" in self.config.keys() and self.config["swot_pixc"].get(
            "enabled", True
        ):
            # Store the SWOT Pixel Cloud config for later use in download/preprocess
            self.swot_pixc_config = self.config["swot_pixc"]

        ### make sure we have a local crs (If we were not able to set it up from the config, grab it from one of the elements)
        if self.local_crs is None:
            if hasattr(self, "rivers"):
                rivers_crs_source = (
                    self.rivers.aoi_gdf
                    if self.rivers.aoi_gdf is not None
                    else self.rivers.gdf
                )
                if (
                    len(rivers_crs_source) > 0
                    and rivers_crs_source.geometry.notna().any()
                ):
                    self.local_crs = rivers_crs_source.estimate_utm_crs()
                elif (
                    len(self.rivers.gdf) > 0 and self.rivers.gdf.geometry.notna().any()
                ):
                    self.local_crs = self.rivers.gdf.estimate_utm_crs()
                else:
                    self.local_crs = self.global_crs

            elif hasattr(self, "reservoirs"):
                self.local_crs = self.reservoirs.gdf.estimate_utm_crs()
            elif hasattr(self, "swot_raster_config"):
                # For swot_raster, use global CRS as local if no local CRS specified
                self.local_crs = self.global_crs
            elif hasattr(self, "swot_pixc_config"):
                # For swot_pixc, use global CRS as local if no local CRS specified
                self.local_crs = self.global_crs
            else:
                raise UserWarning(
                    "Must provide a local crs or a river or reservoir shapefile to determine local crs"
                )

        ### Warn when lake/reservoir-only satellites are configured for rivers or raster modes
        if not hasattr(self, "reservoirs"):
            incompatible = [
                m
                for m in ["icesat2", "sentinel3", "sentinel6"]
                if m in self.to_download or m in self.to_process
            ]
            if incompatible:
                warnings.warn(
                    f"Satellite(s) {incompatible} are configured but have no effect without a "
                    "'reservoirs' section. ICESat-2, Sentinel-3, and Sentinel-6 require reservoir "
                    "waterbody polygons for spatial filtering. Remove these sections or add a "
                    "'reservoirs' section to silence this warning.",
                    UserWarning,
                    stacklevel=2,
                )

    # fucntion to set information for satellite downloads
    def __sat_init(self, name: str):
        if name in self.config.keys():
            if self.config[name]["download"]:
                self.to_download.append(name)

            if self.config[name]["process"]:
                self.to_process.append(name)

            if "download_dir" in self.config[name].keys():
                self.dirs[name] = self.config[name]["download_dir"]
            else:
                self.dirs[name] = os.path.join(self.dirs["main"], name)

            project_cfg = self.config.get("project", {})
            self.startdates[name] = self.config[name].get(
                "startdate"
            ) or project_cfg.get("startdate")
            self.enddates[name] = self.config[name].get("enddate") or project_cfg.get(
                "enddate"
            )

            self.mission_options[name] = {
                k: v
                for k, v in self.config[name].items()
                if k
                in [
                    "hydrocron_fields",
                    "quality_filters",
                    "atl13",
                    "atl13_fields",
                    "track_keys",
                    "subset_file_id",
                    "sigma0_max",
                    "download_threads",
                    "exclude_obs_id_values",
                    "pld_match_max_distance_m",
                ]
            }

            self.processing_options[name] = {
                "processing_filters": self.config[name].get(
                    "processing_filters", ["elevation", "MAD"]
                ),
                "elevation_min_m": self.config[name].get("elevation_min_m", 0.0),
                "elevation_max_m": self.config[name].get("elevation_max_m", 8000.0),
                "mad_threshold": self.config[name].get("mad_threshold", 5.0),
            }

    def _apply_optional_defaults(self):
        for mission, defaults in MISSION_DEFAULTS.items():
            if mission not in self.config:
                continue

            if not isinstance(self.config[mission], dict):
                continue

            for key, value in defaults.items():
                if key not in self.config[mission]:
                    self.config[mission][key] = deepcopy(value)

        # SlideRule returns core fields (height, lat/lon, date, rgt, cycle_number, beam)
        # by default — no forced field merging required.

    def report(self):
        logger.info("Project Name: %s", self.name)
        if hasattr(self, "rivers"):
            self.rivers.report()
        if hasattr(self, "reservoirs"):
            self.reservoirs.report()

    def _require_creodias_credentials(self):
        if not (self.creodias_user and self.creodias_pass):
            raise ValueError(
                "Missing CREODIAS credentials. Provide config['creodias']['username'/'password'] "
                "or set CREODIAS_USERNAME and CREODIAS_PASSWORD in the environment."
            )
        return (self.creodias_user, self.creodias_pass)

    def validate_config(self):
        """Validate loaded config and report all discovered issues at once.

        This is kept for backward compatibility with tests. Validation
        is performed automatically in __post_init__().
        """
        validate_config(self.config)
        return True

    def initialize(self):
        # Checks that we have all information needed for downloads
        if hasattr(self, "reservoirs"):
            flows.initialize_reservoirs(self)
        if hasattr(self, "rivers"):
            flows.initialize_rivers(self)

    def download(self):
        if hasattr(self, "reservoirs"):
            flows.download_reservoirs(self)
        if hasattr(self, "rivers"):
            flows.download_rivers(self)
        if hasattr(self, "swot_raster_config"):
            download_raster(
                config=self.swot_raster_config,
                project_dir=self.dirs["main"],
                global_crs=self.global_crs,
            )
        if hasattr(self, "swot_pixc_config"):
            download_pixc(
                config=self.swot_pixc_config,
                project_dir=self.dirs["main"],
            )

    def update(self):
        # get the current date of the system
        current_date = datetime.date.today()
        logger.info("Updating download archives up to %s", current_date)
        current_date = [current_date.year, current_date.month, current_date.day]

        if hasattr(self, "reservoirs"):
            if "swot" in self.to_download:
                logger.info("Updating SWOT Lake SP product")
            if "icesat2" in self.to_download:
                logger.info("Updating Icesat-2 ATL13 product")
            if "sentinel3" in self.to_download:
                logger.info("Updating Sentinel-3 Hydro product")
            if "sentinel6" in self.to_download:
                logger.info("Updating Sentinel-6 Hydro product")

            self.reservoirs.download(
                to_download=self.to_download,
                startdates=self.startdates,
                enddates=self.enddates,
                earthdata_credentials=(self.earthdata_user, self.earthdata_pass),
                creodias_credentials_provider=self._require_creodias_credentials,
                update_existing=True,
                enddate_overrides={
                    mission: current_date for mission in self.to_download
                },
            )

        if hasattr(self, "rivers"):
            self.rivers.download(
                to_download=self.to_download,
                startdates=self.startdates,
                enddates=self.enddates,
                update_existing=True,
                enddate_overrides={
                    mission: current_date for mission in self.to_download
                },
            )

    def create_timeseries(self):
        warnings.filterwarnings("ignore", module="pyogrio\\..*")
        if hasattr(self, "reservoirs"):
            flows.create_reservoirs_timeseries(self)
        if hasattr(self, "rivers"):
            logger.warning(
                "Rivers preprocessing is not implemented yet; skipping create_timeseries for rivers."
            )

    def generate_summaries(self, show=False, save=True):
        warnings.filterwarnings("ignore", module="pandas\\..*")
        logger.info("Plotting")
        if hasattr(self, "reservoirs"):
            flows.generate_reservoirs_summaries(self, show=show, save=save)
        if hasattr(self, "rivers"):
            flows.generate_rivers_summaries(self, show=show, save=save)
