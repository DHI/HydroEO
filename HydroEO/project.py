"""The project module"""

from dataclasses import dataclass

import os
import yaml
import warnings

import geopandas as gpd
import datetime

from HydroEO.system import Reservoirs
from HydroEO.utils import general


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

        self.dirs = dict()
        self.startdates = dict()
        self.enddates = dict()

        ### Load in the config file and extract parameters
        with open(self.config, "rt") as f:
            self.config = yaml.safe_load(f.read())

        ### Define the project directory for saving outputs, etc
        if "project" in self.config.keys():
            self.dirs["main"] = self.config["project"]["main_dir"]
            general.ifnotmakedirs(self.dirs["main"])
        else:
            raise Warning("Project directory must be defined within configuration file")

        ##### Set credentials and access keys
        if "hydroweb" in self.config.keys():
            if "api_key" in self.config["hydroweb"]:
                # os.environ["EODAG__HYDROWEB_NEXT__AUTH__CREDENTIALS__APIKEY"] = (
                #     self.config["hydroweb"]["api_key"]
                # )
                os.environ["HYDROWEB_API_KEY"] = self.config["hydroweb"]["api_key"]

            if "PLD_path" in self.config["hydroweb"]:
                self.dirs["pld"] = self.config["hydroweb"]["PLD_path"]

        if "earthaccess" in self.config.keys():
            if "username" in self.config["earthaccess"].keys():
                self.earthdata_user = self.config["earthaccess"]["username"]
                self.earthdata_pass = self.config["earthaccess"]["password"]

                os.environ["EDL_USERNAME"] = self.earthdata_user
                os.environ["EDL_PASSWORD"] = self.earthdata_pass

            if "token" in self.config["earthaccess"].keys():
                os.environ["EDL_TOKEN"] = self.config["earthaccess"]["token"]

        if "creodias" in self.config.keys():
            self.creodias_user = self.config["creodias"]["username"]
            self.creodias_pass = self.config["creodias"]["password"]

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
            # self.rivers = Rivers(
            #     gdf=gpd.read_file(self.config["rivers"]["path"]),
            #     dirs=self.dirs,
            #     buffer_width=self.config["rivers"]["buffer_width"],
            #     grid_res=self.config["rivers"]["grid_res"],
            # )
            # self.rivers.gdf = self.rivers.gdf.to_crs(self.global_crs)

        if "reservoirs" in self.config.keys():
            self.reservoirs = Reservoirs(
                gdf=gpd.read_file(self.config["reservoirs"]["path"]),
                id_key=self.config["reservoirs"]["id_key"],
                dirs=self.dirs,
            )

            self.reservoirs.gdf = self.reservoirs.gdf.to_crs(self.global_crs)

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

    def report(self):
        print(f"Project Name: {self.name}")
        if hasattr(self, "rivers"):
            self.rivers.report()
        if hasattr(self, "reservoirs"):
            self.reservoirs.report()

    def initialize(self):
        # Checks that we have all information needed for downloads
        if hasattr(self, "reservoirs"):
            # if we are processing swot data then we will need
            if "swot" in self.to_download or "swot" in self.to_process:
                # download PLD if needed
                self.reservoirs.download_pld(overwrite=False)

                # load the pld and associate the reservoirs with a "lake id"
                self.reservoirs.assign_pld_id(
                    local_crs=self.local_crs, max_distance=100
                )  # take the downloaded PLD and see where we have overlap with the input reservoirs
                self.reservoirs.flag_missing_priors()  # flag and export which reservoirs have entries in the PLD

            # assign download geometry (for reservoirs this is the same as the input boundaries)
            self.reservoirs.download_gdf = self.reservoirs.gdf

    def download(self):
        if hasattr(self, "reservoirs"):
            if "swot" in self.to_download:
                self.reservoirs.download_altimetry(
                    product="SWOT_Lake",
                    startdate=self.startdates["swot"],
                    enddate=self.enddates["swot"],
                    update_existing=False,
                )

            if "icesat2" in self.to_download:
                self.reservoirs.download_altimetry(
                    product="ATL13",
                    startdate=self.startdates["icesat2"],
                    enddate=self.enddates["icesat2"],
                    update_existing=False,
                )

            if "sentinel3" in self.to_download:
                self.reservoirs.download_altimetry(
                    product="S3",
                    startdate=self.startdates["sentinel3"],
                    enddate=self.enddates["sentinel3"],
                    credentials=(self.creodias_user, self.creodias_pass),
                    update_existing=False,
                )
            if "sentinel6" in self.to_download:
                self.reservoirs.download_altimetry(
                    product="S6",
                    startdate=self.startdates["sentinel6"],
                    enddate=self.enddates["sentinel6"],
                    credentials=(self.creodias_user, self.creodias_pass),
                    update_existing=False,
                )

    def update(self):
        # get the current date of the system
        current_date = datetime.date.today()
        print(f"Updating download archives up to {current_date}")
        current_date = [current_date.year, current_date.month, current_date.day]

        if hasattr(self, "reservoirs"):
            if "swot" in self.to_download:
                print("Updating SWOT Lake SP product")
                self.reservoirs.download_altimetry(
                    product="SWOT_Lake",
                    startdate=self.startdates["swot"],
                    enddate=current_date,
                    update_existing=True,
                )
            if "icesat2" in self.to_download:
                print("Updating Icesat-2 ATL13 product")
                self.reservoirs.download_altimetry(
                    product="ATL13",
                    startdate=self.startdates["icesat2"],
                    enddate=current_date,
                    update_existing=True,
                )
            if "sentinel3":
                print("Updating Sentinel-3 Hydro product")
                self.reservoirs.download_altimetry(
                    product="S3",
                    startdate=self.startdates["sentinel3"],
                    enddate=current_date,
                    credentials=(self.creodias_user, self.creodias_pass),
                    update_existing=True,
                )
            if "sentinel6":
                print("Updating Sentinel-6 Hydro product")
                self.reservoirs.download_altimetry(
                    product="S6",
                    startdate=self.startdates["sentinel6"],
                    enddate=current_date,
                    credentials=(self.creodias_user, self.creodias_pass),
                    update_existing=True,
                )

    def create_timeseries(self):
        if hasattr(self, "reservoirs"):
            # extract the data
            self.reservoirs.extract_product_timeseries(self.to_process)

            # clean the individual timeseries
            self.reservoirs.clean_product_timeseries(
                products=self.to_process, filters=["elevation", "MAD"]
            )

            # merge the product timeseries
            self.reservoirs.merge_product_timeseries(products=self.to_process)

    def generate_summaries(self, show=False, save=True):
        if hasattr(self, "reservoirs"):
            for id in self.reservoirs.download_gdf[self.reservoirs.id_key]:
                self.reservoirs.summarize_crossings_by_id(id, show=show, save=save)
                self.reservoirs.summarize_cleaning_by_id(id, show=show, save=save)
                self.reservoirs.summarize_merging_by_id(id, show=show, save=save)
