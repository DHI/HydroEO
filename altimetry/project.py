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


from altimetry.sat_utils import icesat2, sentinel
from altimetry.system import Rivers, Reservoirs
from altimetry import geometry, utils

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
    name : str
    config : str

    def __post_init__(self):

        ### Load in the config file and extract parameters
        with open(self.config, 'rt') as f:
            self.config = yaml.safe_load(f.read())

        ### Define the project directory for saving outputs, etc
        if 'project' in self.config.keys():
            self.main_dir = self.config['project']["main_dir"]
            utils.ifnotmakedirs(self.main_dir)
        else:
            raise Warning("Project directory must be defined within configuration file")


        ### Parameters associated with icesat2 downloads
        if 'earthdata' in self.config.keys():
            self.earthdata_user = self.config['earthdata']["username"]
            self.earthdata_pass = self.config['earthdata']["password"]
        
        if 'icesat2' in self.config.keys():

            if 'download_dir' in self.config['icesat2'].keys():
                self.icesat2_dir = self.config['icesat2']["download_dir"]
            else:
                self.icesat2_dir = os.path.join(self.main_dir, 'icesat2')

            self.icesat2_startdate = self.config['icesat2']["startdate"]
            self.icesat2_enddate = self.config['icesat2']["enddate"]
        
        ### Parameters associated with Sentinel downloads
        if 'creodias' in self.config.keys():
            self.creodias_user = self.config['creodias']["username"]
            self.creodias_pass = self.config['creodias']["password"]
        
        if 'sentinel3' in self.config.keys():

            if 'download_dir' in self.config['sentinel3'].keys():
                self.sentinel3_dir = self.config['sentinel3']["download_dir"]
            else:
                self.sentinel3_dir = os.path.join(self.main_dir, 'sentinel3')

            self.sentinel3_startdate = self.config['sentinel3']["startdate"]
            self.sentinel3_enddate = self.config['sentinel3']["enddate"]

        if 'sentinel6' in self.config.keys():

            if 'download_dir' in self.config['sentinel6'].keys():
                self.sentinel6_dir = self.config['sentinel6']["download_dir"]
            else:
                self.sentinel6_dir = os.path.join(self.main_dir, 'sentinel6')

            self.sentinel6_startdate = self.config['sentinel6']["startdate"]
            self.sentinel6_enddate = self.config['sentinel6']["enddate"]


        ### Set the crs from the congiguration file is possible
        if 'gis' in self.config.keys():

            if 'global_crs' in self.config['gis'].keys():
                self.global_crs = self.config['gis']['global_crs']
            else:
                self.global_crs = 'EPSG:4326'

            if 'local_crs' in self.config['gis'].keys():
                self.local_crs = self.config['gis']['local_crs']
            else:
                self.local_crs = None


        ### load in elements for download and processing
        if 'rivers' in self.config.keys():

            self.rivers = Rivers(gdf = gpd.read_file( self.config['rivers']['path']),
                                 icesat2_dir = self.icesat2_dir,
                                 sentinel3_dir = self.sentinel3_dir,
                                 sentinel6_dir = self.sentinel6_dir,
                                 output_dir = os.path.join(self.main_dir, "rivers"),
                                 buffer_width = self.config['rivers']['buffer_width'],
                                 grid_res     = self.config['rivers']['grid_res']
                                 )
            self.rivers.gdf = self.rivers.gdf.to_crs(self.global_crs)


        if 'reservoirs' in self.config.keys():

            self.reservoirs = Reservoirs(gdf = gpd.read_file(self.config['reservoirs']['path']),
                                         icesat2_dir = self.icesat2_dir,
                                         sentinel3_dir = self.sentinel3_dir,
                                         sentinel6_dir = self.sentinel6_dir,
                                         output_dir = os.path.join(self.main_dir, "reservoirs")
                                         )
            
            self.reservoirs.gdf = self.reservoirs.gdf.to_crs(self.global_crs)


        ### make sure we have a local crs (If we were not able to set it up from the config, grab it from one of the elements)
        if self.local_crs is None:
            if hasattr(self, 'rivers'):
                self.local_crs = self.rivers.gdf.estimate_utm_crs()

            elif hasattr(self, 'reservoirs'):
                self.local_crs = self.reservoirs.gdf.estimate_utm_crs()
            else:
                raise UserWarning("Must provide a local crs or a river or reservoir shapefile to determine local crs")
            

        ### finally make sure we have initiated a buffered shape if we have rivers in the project, we must do this last so that we can work in the use defiend local crs if it exists
        if hasattr(self, 'rivers'):
            self.rivers.make_buffer(local_crs=self.local_crs)




    def report(self):
        print(f"Project Name: {self.name}")
        if hasattr(self, 'rivers'):
            self.rivers.report()
        if hasattr(self, 'reservoirs'):
            self.reservoirs.report()