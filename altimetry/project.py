"""The project module"""

from dataclasses import dataclass
import warnings
import typing

import os
import yaml

import pandas as pd
import geopandas as gpd
import shapely

from cmcrameri import cm
import matplotlib.pyplot as plt
from matplotlib.dates import date2num

from datetime import date, datetime

from altimetry.system import Rivers, Reservoirs
from altimetry.utils import utils


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
            utils.ifnotmakedirs(self.dirs["main"])
        else:
            raise Warning("Project directory must be defined within configuration file")

        ##### Set credentials and access keys
        if "hydroweb" in self.config.keys():
            if "api_key" in self.config["hydroweb"]:
                os.environ["EODAG__HYDROWEB_NEXT__AUTH__CREDENTIALS__APIKEY"] = (
                    self.config["hydroweb"]["api_key"]
                )

        if "earthaccess" in self.config.keys():
            self.earthdata_user = self.config["earthaccess"]["username"]
            self.earthdata_pass = self.config["earthaccess"]["password"]

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
            self.rivers = Rivers(
                gdf=gpd.read_file(self.config["rivers"]["path"]),
                dirs=self.dirs,
                buffer_width=self.config["rivers"]["buffer_width"],
                grid_res=self.config["rivers"]["grid_res"],
            )
            self.rivers.gdf = self.rivers.gdf.to_crs(self.global_crs)

        if "reservoirs" in self.config.keys():
            self.dirs["pld"] = self.config["reservoirs"]["prior_path"]

            self.reservoirs = Reservoirs(
                gdf=gpd.read_file(self.config["reservoirs"]["path"]), dirs=self.dirs
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
