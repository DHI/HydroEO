"""The project module"""

from dataclasses import dataclass
from copy import deepcopy
import logging

import os
import yaml
import warnings

import geopandas as gpd
import datetime

from HydroEO.system import Reservoirs
from HydroEO.satellites.icesat2 import (
    ATL13_DEFAULT_FIELDS,
    SR_ATL13_VALID_ANCILLARY_FIELDS,
)
from HydroEO.flows import DownloadFlow, PlottingFlow, PreprocessFlow
from HydroEO.utils import general

logger = logging.getLogger(__name__)


SUPPORTED_TRACK_KEYS = ["gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"]
SUPPORTED_CLEAN_FILTERS = ["elevation", "MAD", "daily_mean", "hampel", "rolling_median"]

MISSION_DEFAULTS = {
    "swot": {
        "download": False,
        "process": False,
        "pld_match_max_distance_m": 100.0,
        "exclude_obs_id_values": ["no_data"],
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
        if "rivers" in self.config.keys():
            warnings.warn("Rivers system is not yet implemented, input will be ignored")

        if "reservoirs" in self.config.keys():
            self.reservoirs = Reservoirs(
                gdf=gpd.read_file(self.config["reservoirs"]["path"]),
                id_key=self.config["reservoirs"]["id_key"],
                dirs=self.dirs,
            )

            self.reservoirs.gdf = self.reservoirs.gdf.to_crs(self.global_crs)
            self.reservoirs.mission_options = self.mission_options
            self.reservoirs.processing_options = self.processing_options

        ### make sure we have a local crs (If we were not able to set it up from the config, grab it from one of the elements)
        if self.local_crs is None:
            if hasattr(self, "rivers"):
                self.local_crs = self.rivers.gdf.estimate_utm_crs()

            elif hasattr(self, "reservoirs"):
                self.local_crs = self.reservoirs.gdf.estimate_utm_crs()
            else:
                raise UserWarning(
                    "Must provide a local crs or a river or reservoir shapefile to determine local crs"
                )

        ### finally make sure we have initiated a buffered shape if we have rivers in the project, we must do this last so that we can work in the use defiend local crs if it exists
        if hasattr(self, "rivers"):
            self.rivers.make_buffer(local_crs=self.local_crs)

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

    @staticmethod
    def _is_valid_date_tuple(value):
        if not isinstance(value, list) or len(value) != 3:
            return False
        if not all(isinstance(v, int) for v in value):
            return False
        try:
            datetime.date(*value)
        except Exception:
            return False
        return True

    def validate_config(self):
        """Validate loaded config and report all discovered issues at once."""
        issues = []
        cfg = self.config or {}

        if "project" not in cfg or not isinstance(cfg["project"], dict):
            issues.append("Missing required section 'project'.")
        elif not cfg["project"].get("main_dir"):
            issues.append("Missing required key 'project.main_dir'.")

        if "reservoirs" not in cfg or not isinstance(cfg["reservoirs"], dict):
            issues.append("Missing required section 'reservoirs'.")
        else:
            if not cfg["reservoirs"].get("path"):
                issues.append("Missing required key 'reservoirs.path'.")
            elif not os.path.exists(cfg["reservoirs"]["path"]):
                issues.append(
                    f"Path in 'reservoirs.path' does not exist: {cfg['reservoirs']['path']}"
                )

            if not cfg["reservoirs"].get("id_key"):
                issues.append("Missing required key 'reservoirs.id_key'.")

        for mission in ["swot", "icesat2", "sentinel3", "sentinel6"]:
            if mission not in cfg:
                continue

            mission_cfg = cfg[mission]
            if not isinstance(mission_cfg, dict):
                issues.append(
                    f"Section '{mission}' must be a mapping of key/value pairs."
                )
                continue

            for key in ["download", "process"]:
                if key not in mission_cfg:
                    issues.append(f"Missing required key '{mission}.{key}'.")
                elif not isinstance(mission_cfg[key], bool):
                    issues.append(f"'{mission}.{key}' must be a boolean value.")

            for key in ["startdate", "enddate"]:
                if key not in mission_cfg:
                    issues.append(f"Missing required key '{mission}.{key}'.")
                elif not self._is_valid_date_tuple(mission_cfg[key]):
                    issues.append(
                        f"'{mission}.{key}' must be [year, month, day] with valid integer values."
                    )

            if self._is_valid_date_tuple(
                mission_cfg.get("startdate")
            ) and self._is_valid_date_tuple(mission_cfg.get("enddate")):
                if datetime.date(*mission_cfg["startdate"]) > datetime.date(
                    *mission_cfg["enddate"]
                ):
                    issues.append(
                        f"'{mission}.startdate' cannot be after '{mission}.enddate'."
                    )

            filters = mission_cfg.get("processing_filters", ["elevation", "MAD"])
            if not isinstance(filters, list) or any(
                not isinstance(v, str) for v in filters
            ):
                issues.append(
                    f"'{mission}.processing_filters' must be a list of strings."
                )
            else:
                invalid_filters = [
                    v for v in filters if v not in SUPPORTED_CLEAN_FILTERS
                ]
                if invalid_filters:
                    issues.append(
                        f"Invalid filters in '{mission}.processing_filters': {invalid_filters}. "
                        f"Valid filters are: {SUPPORTED_CLEAN_FILTERS}."
                    )

            min_h = mission_cfg.get("elevation_min_m", 0.0)
            max_h = mission_cfg.get("elevation_max_m", 8000.0)
            if not isinstance(min_h, (int, float)) or not isinstance(
                max_h, (int, float)
            ):
                issues.append(
                    f"'{mission}.elevation_min_m' and '{mission}.elevation_max_m' must be numeric."
                )
            elif min_h >= max_h:
                issues.append(
                    f"'{mission}.elevation_min_m' must be smaller than '{mission}.elevation_max_m'."
                )

            mad_threshold = mission_cfg.get("mad_threshold", 5.0)
            if not isinstance(mad_threshold, (int, float)) or mad_threshold <= 0:
                issues.append(f"'{mission}.mad_threshold' must be a positive number.")

            if mission == "swot":
                max_distance = mission_cfg.get("pld_match_max_distance_m", 100.0)
                if not isinstance(max_distance, (int, float)) or max_distance < 0:
                    issues.append(
                        "'swot.pld_match_max_distance_m' must be a non-negative number."
                    )

                excluded_obs = mission_cfg.get("exclude_obs_id_values", ["no_data"])
                if not isinstance(excluded_obs, list) or any(
                    not isinstance(v, str) for v in excluded_obs
                ):
                    issues.append(
                        "'swot.exclude_obs_id_values' must be a list of strings."
                    )

            if mission == "icesat2":
                track_keys = mission_cfg.get("track_keys", SUPPORTED_TRACK_KEYS)
                if not isinstance(track_keys, list) or any(
                    not isinstance(v, str) for v in track_keys
                ):
                    issues.append("'icesat2.track_keys' must be a list of strings.")
                else:
                    invalid_track_keys = [
                        k for k in track_keys if k not in SUPPORTED_TRACK_KEYS
                    ]
                    if invalid_track_keys:
                        issues.append(
                            f"Invalid entries in 'icesat2.track_keys': {invalid_track_keys}. "
                            f"Valid options are: {SUPPORTED_TRACK_KEYS}."
                        )

                fields = mission_cfg.get("atl13_fields", ATL13_DEFAULT_FIELDS)
                if not isinstance(fields, list) or any(
                    not isinstance(v, str) for v in fields
                ):
                    issues.append("'icesat2.atl13_fields' must be a list of strings.")
                else:
                    invalid_fields = [
                        v for v in fields if v not in SR_ATL13_VALID_ANCILLARY_FIELDS
                    ]
                    if invalid_fields:
                        valid_fields = sorted(SR_ATL13_VALID_ANCILLARY_FIELDS)
                        issues.append(
                            "Invalid ATL13 fields in 'icesat2.atl13_fields': "
                            f"{invalid_fields}. Valid options are: {valid_fields}."
                        )

                atl13_cfg = mission_cfg.get("atl13", {})
                if not isinstance(atl13_cfg, dict):
                    issues.append(
                        "'icesat2.atl13' must be a mapping of key/value pairs."
                    )
                else:
                    if "pass_invalid" in atl13_cfg and not isinstance(
                        atl13_cfg["pass_invalid"], bool
                    ):
                        issues.append("'icesat2.atl13.pass_invalid' must be a boolean.")
                    for _beam_key in ("beams", "spots"):
                        if _beam_key in atl13_cfg:
                            _v = atl13_cfg[_beam_key]
                            if not isinstance(_v, list) or any(
                                not isinstance(x, str) for x in _v
                            ):
                                issues.append(
                                    f"'icesat2.atl13.{_beam_key}' must be a list of strings."
                                )

            if mission in ["sentinel3", "sentinel6"]:
                sigma0_max = mission_cfg.get("sigma0_max", 1e5)
                if not isinstance(sigma0_max, (int, float)) or sigma0_max <= 0:
                    issues.append(f"'{mission}.sigma0_max' must be a positive number.")

                subset_file_id = mission_cfg.get(
                    "subset_file_id", "enhanced_measurement.nc"
                )
                if not isinstance(subset_file_id, str) or subset_file_id.strip() == "":
                    issues.append(
                        f"'{mission}.subset_file_id' must be a non-empty string."
                    )

                download_threads = mission_cfg.get("download_threads", 1)
                if not isinstance(download_threads, int) or download_threads <= 0:
                    issues.append(
                        f"'{mission}.download_threads' must be a positive integer."
                    )

        if issues:
            raise ValueError("Invalid configuration:\n - " + "\n - ".join(issues))

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

    def download(self):
        if hasattr(self, "reservoirs"):
            DownloadFlow(
                reservoirs=self.reservoirs,
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

            DownloadFlow(
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

    def create_timeseries(self):
        warnings.filterwarnings("ignore", module="pyogrio\\..*")
        if hasattr(self, "reservoirs"):
            PreprocessFlow(
                reservoirs=self.reservoirs,
                to_process=self.to_process,
                processing_options=self.processing_options,
            ).run()

    def generate_summaries(self, show=False, save=True):
        warnings.filterwarnings("ignore", module="pandas\\..*")
        if hasattr(self, "reservoirs"):
            PlottingFlow(self.reservoirs).run(show=show, save=save)
