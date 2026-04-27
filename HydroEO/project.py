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
from HydroEO.satellites.icesat2 import ATL13_DEFAULT_FIELDS
from HydroEO.flows import (
    ReservoirDownloadFlow,
    RiverDownloadFlow,
    PlottingFlow,
    PreprocessFlow,
)
from HydroEO.utils import general
from HydroEO.validation import (
    validate_config,
    SUPPORTED_TRACK_KEYS,
    DEFAULT_SWOT_HYDROCRON_FIELDS,
    DEFAULT_SWOT_QUALITY_FILTERS,
)

logger = logging.getLogger(__name__)

MISSION_DEFAULTS = {
    "swot": {
        "download": False,
        "process": False,
        "pld_match_max_distance_m": 100.0,
        "exclude_obs_id_values": ["no_data"],
        "hydrocron_fields": DEFAULT_SWOT_HYDROCRON_FIELDS,
        "quality_filters": DEFAULT_SWOT_QUALITY_FILTERS,
        "processing_filters": ["elevation", "MAD"],
        "elevation_min_m": 0.0,
        "elevation_max_m": 8000.0,
        "mad_threshold": 5.0,
    },
    "icesat2": {
        "download": False,
        "process": False,
        "atl13_fields": ATL13_DEFAULT_FIELDS,
        "atl13": {"pass_invalid": False, "beams": [], "spots": []},
        "track_keys": SUPPORTED_TRACK_KEYS,
        "processing_filters": ["elevation", "MAD"],
        "elevation_min_m": 0.0,
        "elevation_max_m": 8000.0,
        "mad_threshold": 5.0,
    },
    "sentinel3": {
        "download": False,
        "process": False,
        "subset_file_id": "enhanced_measurement.nc",
        "sigma0_max": 1e5,
        "download_threads": 1,
        "processing_filters": ["elevation", "MAD"],
        "elevation_min_m": 0.0,
        "elevation_max_m": 8000.0,
        "mad_threshold": 5.0,
    },
    "sentinel6": {
        "download": False,
        "process": False,
        "subset_file_id": "enhanced_measurement.nc",
        "sigma0_max": 1e5,
        "download_threads": 1,
        "processing_filters": ["elevation", "MAD"],
        "elevation_min_m": 0.0,
        "elevation_max_m": 8000.0,
        "mad_threshold": 5.0,
    },
}

# SlideRule's atl13x always returns the core fields (height, lat/lon in geometry,
# date index, beam, rgt, cycle_number) — no forced field merging is needed.
ICESAT2_REQUIRED_FIELDS: list[str] = []


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
        self.validate_config()

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

        if "PLD_path" in hydroweb_cfg:
            self.dirs["pld"] = hydroweb_cfg["PLD_path"]

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

        ##### Set attributes for each satellite to be downloaded or processed
        self.__sat_init("swot")
        self.__sat_init("icesat2")
        self.__sat_init("sentinel3")
        self.__sat_init("sentinel6")

        ### Set the crs from the congiguration file if possible
        if "gis" in self.config.keys():
            if "global_crs" in self.config["gis"].keys():
                self.global_crs = self.config["gis"]["global_crs"]
            else:
                self.global_crs = "EPSG:4326"

            if "local_crs" in self.config["gis"].keys():
                self.local_crs = self.config["gis"]["local_crs"]
            else:
                self.local_crs = None

        ### load in elements for download and processing
        if "reservoirs" in self.config.keys():
            self.reservoirs = Reservoirs(
                gdf=gpd.read_file(self.config["reservoirs"]["path"]),
                id_key=self.config["reservoirs"]["id_key"],
                dirs=self.dirs,
            )

            self.reservoirs.gdf = self.reservoirs.gdf.to_crs(self.global_crs)
            self.reservoirs.mission_options = self.mission_options
            self.reservoirs.processing_options = self.processing_options

        if "rivers" in self.config.keys():
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
                if "node_numbers" in rivers_cfg:
                    id_values = rivers_cfg["node_numbers"]
                    rivers_id_key = rivers_cfg.get("id_key", "river_id")
                    input_mode = "node_numbers"
                    target_id_col = "node_id"
                else:
                    id_values = rivers_cfg.get("reach_numbers", [])
                    rivers_id_key = rivers_cfg.get("id_key", "river_id")
                    input_mode = "reach_numbers"
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
            else:
                raise UserWarning(
                    "Must provide a local crs or a river or reservoir shapefile to determine local crs"
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

            self.startdates[name] = self.config[name]["startdate"]
            self.enddates[name] = self.config[name]["enddate"]

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

    def validate_config(self):
        """Validate loaded config and report all discovered issues at once."""
        validate_config(self.config)
        return True

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

    def initialize(self):
        self.validate_config()

        # Checks that we have all information needed for downloads
        if hasattr(self, "reservoirs"):
            # if we are processing swot data then we will need
            if "swot" in self.to_download or "swot" in self.to_process:
                # download PLD if needed
                self.reservoirs.download_pld(overwrite=False)

                # load the pld and associate the reservoirs with a "lake id"
                self.reservoirs.assign_pld_id(
                    local_crs=self.local_crs,
                    max_distance=self.mission_options.get("swot", {}).get(
                        "pld_match_max_distance_m", 100
                    ),
                )  # take the downloaded PLD and see where we have overlap with the input reservoirs
                self.reservoirs.flag_missing_priors()  # flag and export which reservoirs have entries in the PLD

            # assign download geometry (for reservoirs this is the same as the input boundaries)
            self.reservoirs.download_gdf = self.reservoirs.gdf

        if hasattr(self, "rivers"):
            self.rivers.prepare_download_targets(local_crs=self.local_crs)
            id_label = "node" if self.rivers.target_id_col == "node_id" else "reach"
            logger.info(
                "Initialized river %s ids: %s",
                id_label,
                ", ".join(str(target_id) for target_id in self.rivers.target_ids),
            )

    def download(self):
        if hasattr(self, "reservoirs"):
            ReservoirDownloadFlow(
                reservoirs=self.reservoirs,
                to_download=self.to_download,
                startdates=self.startdates,
                enddates=self.enddates,
                earthdata_credentials=(self.earthdata_user, self.earthdata_pass),
                creodias_credentials_provider=self._require_creodias_credentials,
            ).run(update_existing=False)
        if hasattr(self, "rivers"):
            RiverDownloadFlow(
                rivers=self.rivers,
                to_download=self.to_download,
                startdates=self.startdates,
                enddates=self.enddates,
                earthdata_credentials=(self.earthdata_user, self.earthdata_pass),
                creodias_credentials_provider=self._require_creodias_credentials,
            ).run(update_existing=False)

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

            ReservoirDownloadFlow(
                reservoirs=self.reservoirs,
                to_download=self.to_download,
                startdates=self.startdates,
                enddates=self.enddates,
                earthdata_credentials=(self.earthdata_user, self.earthdata_pass),
                creodias_credentials_provider=self._require_creodias_credentials,
            ).run(
                update_existing=True,
                enddate_overrides={
                    mission: current_date for mission in self.to_download
                },
            )

        if hasattr(self, "rivers"):
            RiverDownloadFlow(
                rivers=self.rivers,
                to_download=self.to_download,
                startdates=self.startdates,
                enddates=self.enddates,
                earthdata_credentials=(self.earthdata_user, self.earthdata_pass),
                creodias_credentials_provider=self._require_creodias_credentials,
            ).run(
                update_existing=True,
                enddate_overrides={
                    mission: current_date for mission in self.to_download
                },
            )

    def create_timeseries(self):
        warnings.filterwarnings("ignore", module="pyogrio\\..*")
        if hasattr(self, "reservoirs"):
            PreprocessFlow(
                reservoirs=self.reservoirs,
                to_process=self.to_process,
                processing_options=self.processing_options,
            ).run()
        if hasattr(self, "rivers"):
            logger.warning(
                "Rivers preprocessing is not implemented yet; skipping create_timeseries for rivers."
            )

    def generate_summaries(self, show=False, save=True):
        warnings.filterwarnings("ignore", module="pandas\\..*")
        if hasattr(self, "reservoirs"):
            PlottingFlow(self.reservoirs).run(show=show, save=save)
        if hasattr(self, "rivers"):
            logger.warning(
                "Rivers plotting is not implemented yet; skipping generate_summaries for rivers."
            )
