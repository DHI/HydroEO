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
        self.merging_options = dict()

        self.dirs = dict()
        self.startdates = dict()
        self.enddates = dict()

        ### Load in the config file and extract parameters
        with open(self.config, "rt", encoding="utf-8") as f:
            self.config = yaml.safe_load(f.read())

        if self.config is None:
            self.config = {}

        self._apply_optional_defaults()
        validate_config(self.config)

        ### Define the project directory for saving outputs, etc
        if "project" in self.config.keys():
            self.dirs["main"] = general.normalize_path(self.config["project"]["main_dir"])
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
            self.dirs["pld_raw"] = general.normalize_path(hydroweb_cfg["raw_pld_path"])
        self.keep_raw_pld = hydroweb_cfg.get("keep_raw_pld", False)

        # Set SWORD database paths and configuration
        sword_db_cfg = self.config.get("sword_db", {})
        self.dirs["sword"] = os.path.join(self.dirs["main"], "aux", "SWORD", "gpkg")
        self.dirs["sword_subset"] = os.path.join(
            self.dirs["main"], "aux", "SWORD", "SWORD_subset.gpkg"
        )
        if "sword_subset_path" in sword_db_cfg:
            self.dirs["sword_subset"] = general.normalize_path(sword_db_cfg["sword_subset_path"])
        if "raw_sword_path" in sword_db_cfg:
            self.dirs["sword_raw"] = general.normalize_path(sword_db_cfg["raw_sword_path"])
            # Safety: if raw_sword_path is outside main_dir, force keep_raw_sword=True
            raw_path = os.path.abspath(sword_db_cfg["raw_sword_path"])
            main_path = os.path.abspath(self.dirs["main"])
            try:
                is_outside = not raw_path.startswith(main_path)
            except (TypeError, ValueError):
                is_outside = True
            if is_outside:
                self.keep_raw_sword = True
            else:
                self.keep_raw_sword = sword_db_cfg.get("keep_raw_sword", False)
        else:
            self.keep_raw_sword = sword_db_cfg.get("keep_raw_sword", False)

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
            # NOTE: this was previously read via getattr(prj.reservoirs,
            # "overwrite_extraction", False) in flows.py but never
            # actually wired to config -- meaning it silently always
            # defaulted to False regardless of what a user might have
            # tried to set. Fixed here.
            self.reservoirs.overwrite_extraction = self.config["reservoirs"].get(
                "overwrite_extraction", False
            )

            # User-configurable overrides for the merge()/Kalman/svr_radial
            # pipeline (see flows.DEFAULT_RESERVOIR_MERGING_OPTIONS for
            # every available key and its default). Only the keys the user
            # actually sets here are used to override the defaults -- any
            # keys not mentioned keep their default value, so a user only
            # needs to specify what they want to change.
            # NOTE: DEFAULT_RESERVOIR_MERGING_OPTIONS' svr_radial_err/gamma
            # trace back to DAHITI's own published values -- and DAHITI's
            # calibration is for LAKES specifically. A reservoir with
            # managed/operational water level changes (fill/drawdown
            # cycles) can have real dynamics on a much faster timescale
            # than a natural lake, making the lake-tuned defaults too
            # strict (rejecting genuine fast changes as if they were
            # noise). If your reservoirs behave more like this, override
            # svr_radial_err and svr_radial_gamma here rather than relying
            # on the lake-tuned defaults.
            self.merging_options = self.config["reservoirs"].get(
                "merging_options", {}
            )
            self.reservoirs.merging_options = self.merging_options

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

            # Corridor buffer for ICESat-2/Sentinel-3/6 extraction (see
            # flows._river_target_corridor). Deliberately separate from
            # buffer_meters above, which only decides which SWORD
            # targets count as "in the AOI" at all. If left unset
            # (None, the default), each target's own SWORD "width"
            # attribute is used instead of one flat value for every
            # target -- see width_buffer_factor below. Only falls back
            # to a flat value (prj.rivers.buffer_meters, then 500m) if
            # no usable "width" column is found.
            self.rivers.extraction_buffer_meters = rivers_cfg.get(
                "extraction_buffer_meters"
            )
            # Margin applied on top of each target's own SWORD width
            # when using the width-based default above (ignored if
            # extraction_buffer_meters is set explicitly). Default 1.05
            # = 5% wider than the target's actual channel width.
            self.rivers.width_buffer_factor = rivers_cfg.get(
                "width_buffer_factor", 1.05
            )
            # Max distance (m) for assigning a raw altimetry point to its
            # nearest SWORD target (see flows._assign_points_to_river_targets).
            # Falls back to the extraction buffer if not set.
            self.rivers.max_node_assignment_meters = rivers_cfg.get(
                "max_node_assignment_meters"
            )
            self.rivers.overwrite_extraction = rivers_cfg.get(
                "overwrite_extraction", False
            )

            # User-configurable overrides for the merge()/Kalman/svr_radial
            # pipeline (see flows.DEFAULT_RIVER_MERGING_OPTIONS for every
            # available key and its default). Set directly on
            # prj.rivers rather than routed through the shared
            # project-level self.merging_options reservoirs use, so a
            # project with both reservoirs and rivers configured doesn't
            # have them silently collide -- see flows._merge_timeseries's
            # per-target-type override lookup.
            # NOTE: DEFAULT_RIVER_MERGING_OPTIONS is currently a direct
            # copy of the reservoir defaults and has NOT been
            # independently validated against real river data -- unlike
            # the reservoir defaults, which were tuned this way against
            # real reservoirs. Treat it as a starting point to check
            # kept/rejected counts against, not a verified value.
            self.rivers.merging_options = rivers_cfg.get("merging_options", {})

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

        ### Warn when lake/reservoir-only satellites are configured for neither
        # reservoirs nor rivers mode. Previously this only checked for
        # 'reservoirs' since icesat2/sentinel3/sentinel6 had no river support
        # at all -- now that they work for rivers too (see
        # flows._download_rivers_icesat2/_download_rivers_sentinel), the
        # warning should only fire if NEITHER mode is configured.
        if not hasattr(self, "reservoirs") and not hasattr(self, "rivers"):
            incompatible = [
                m
                for m in ["icesat2", "sentinel3", "sentinel6"]
                if m in self.to_download or m in self.to_process
            ]
            if incompatible:
                warnings.warn(
                    f"Satellite(s) {incompatible} are configured but have no effect "
                    "without a 'reservoirs' or 'rivers' section. ICESat-2, Sentinel-3, "
                    "and Sentinel-6 require reservoir waterbody polygons or river "
                    "SWORD targets for spatial filtering. Remove these sections or add "
                    "a 'reservoirs'/'rivers' section to silence this warning.",
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
                self.dirs[name] = general.normalize_path(self.config[name]["download_dir"])
                #self.dirs[name] = Path(self.config[name]["download_dir"])
            else:
                self.dirs[name] = os.path.join(self.dirs["main"], "raw", name)

            if name == "icesat2":
                self.dirs["icesat2_processed"] = os.path.join(
                    self.dirs["main"], "processed", "icesat2"
                )

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
                    "sigma0_min",
                    "source",
                    "latency",
                    "short_name",
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

    def _require_earthdata_credentials(self):
        """
        Check upfront that EarthData credentials are available in some
        form earthaccess recognizes, before calling earthaccess.login().
        Without this check, earthaccess.login() -- called with no
        explicit strategy, same as this codebase's existing SWOT path --
        falls through to environment variables, then a .netrc file, then
        INTERACTIVE PROMPTING if neither is found (confirmed against
        earthaccess's own documentation). In a non-interactive run (a
        scheduled job, a CLI invocation) that prompt either hangs waiting
        for input that will never come, or raises a confusing low-level
        error deep inside earthaccess -- rather than a clear, immediate
        one here.
        """
        has_env = bool(
            os.environ.get("EARTHDATA_USERNAME") and os.environ.get("EARTHDATA_PASSWORD")
        )
        has_token = bool(os.environ.get("EARTHDATA_TOKEN"))
        netrc_path = os.environ.get(
            "NETRC",
            os.path.join(os.path.expanduser("~"), "_netrc" if os.name == "nt" else ".netrc"),
        )
        has_netrc = os.path.exists(netrc_path)

        if not (has_env or has_token or has_netrc):
            raise ValueError(
                "Missing EarthData credentials for the Sentinel-6 EarthData "
                "source (mission_options['sentinel6']['source'] = 'earthdata'). "
                "Set EARTHDATA_USERNAME and EARTHDATA_PASSWORD in the "
                "environment, or EARTHDATA_TOKEN, or create a .netrc file with "
                "your Earthdata Login credentials (register free at "
                "https://urs.earthdata.nasa.gov). Without one of these, "
                "earthaccess.login() falls through to an interactive prompt, "
                "which will hang in a non-interactive run rather than fail "
                "clearly."
            )

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
            flows.create_rivers_timeseries(self)

    def generate_summaries(self, show=False, save=True):
        warnings.filterwarnings("ignore", module="pandas\\..*")
        logger.info("Plotting")
        if hasattr(self, "reservoirs"):
            flows.generate_reservoirs_summaries(self, show=show, save=save)
        if hasattr(self, "rivers"):
            flows.generate_rivers_summaries(self, show=show, save=save)

    def _infer_target_type(self, target_type=None):
        """
        Resolve which target type (reservoirs/rivers) a per-target call
        applies to. If target_type is given explicitly, use it. Otherwise,
        infer it automatically when the project only has one of the two
        configured -- the common case -- and require an explicit choice
        only when both are configured, since there's no way to guess
        correctly between them.
        """
        if target_type is not None:
            return target_type
        has_reservoirs = hasattr(self, "reservoirs")
        has_rivers = hasattr(self, "rivers")
        if has_reservoirs and not has_rivers:
            return "reservoirs"
        if has_rivers and not has_reservoirs:
            return "rivers"
        if has_reservoirs and has_rivers:
            raise ValueError(
                "Both reservoirs and rivers are configured for this "
                "project -- specify target_type='reservoirs' or "
                "target_type='rivers' explicitly."
            )
        raise ValueError("Neither reservoirs nor rivers is configured for this project.")

    def list_target_observations(self, id, target_type=None):
        """
        Summarize what observations exist for a target (reservoir or
        river node/reach) at (platform, orbit) granularity -- the "what
        could I exclude" view. See generate_summaries()'s plots
        (platform-colored merge progress) for where a problem actually
        shows up in the data.
        """
        target_type = self._infer_target_type(target_type)
        return flows.list_target_observations(self, target_type, id)

    def exclude_observations(
        self, id, platform=None, orbit=None, date=None, reason=None, target_type=None,
    ):
        """
        Exclude observations from a target's merge, at whatever
        granularity is given (a whole platform, a specific orbit/pass
        value, a specific date, or any combination). Persisted to
        {output}/{id}/run_config.yaml -- survives across runs, and can be
        hand-edited or version-controlled directly. Re-run
        create_timeseries() (or just the merge step) afterward to apply it.

        Examples
        --------
        project.exclude_observations(my_id, platform="S3B", reason="bad calibration pass")
        project.exclude_observations(my_id, platform="S3B", orbit=1517)
        project.exclude_observations(my_id, date="2024-03-19")
        """
        target_type = self._infer_target_type(target_type)
        return flows.exclude_from_target(
            self, target_type, id, platform=platform, orbit=orbit, date=date, reason=reason,
        )

    def list_exclusions(self, id, target_type=None):
        """Current exclusion rules for a target, from its run_config.yaml."""
        target_type = self._infer_target_type(target_type)
        return flows.list_exclusions(self, target_type, id)

    def remove_exclusion(
        self, id, index=None, platform=None, orbit=None, date=None, target_type=None,
    ):
        """
        Remove one or more exclusion rules, either by position (index,
        from list_exclusions()) or by matching criteria -- the same
        platform/orbit/date fields used with exclude_observations().

        Examples
        --------
        project.remove_exclusion(my_id, platform="S3B", orbit=1517)
        project.remove_exclusion(my_id, index=0)
        """
        target_type = self._infer_target_type(target_type)
        return flows.remove_exclusion(
            self, target_type, id, index=index, platform=platform, orbit=orbit, date=date,
        )

    def set_merging_option(self, id, target_type=None, **kwargs):
        """
        Override one or more merging_options for just this one target.
        Example: project.set_merging_option(my_id, svr_radial_err=0.5)
        """
        target_type = self._infer_target_type(target_type)
        return flows.set_merging_option(self, target_type, id, **kwargs)