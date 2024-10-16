"""The project module"""
from dataclasses import dataclass
import warnings
import typing

import os
import yaml

import pandas as pd
import geopandas as gpd

import matplotlib.pyplot as plt
from matplotlib.dates import date2num

from datetime import date, datetime

from altimetry import icesat2, geometry, utils

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

        # Load in the config file and extract parameters
        with open(self.config, 'rt') as f:
            self.config = yaml.safe_load(f.read())

        if 'earthdata' in self.config.keys():
            self.earthdata_user = self.config['earthdata']["username"]
            self.earthdata_pass = self.config['earthdata']["password"]
        
        if 'icesat2' in self.config.keys():
            self.icesat2_dir = self.config['icesat2']["download_dir"]
            self.icesat2_startdate = self.config['icesat2']["startdate"]
            self.icesat2_enddate = self.config['icesat2']["enddate"]

        if 'gis' in self.config.keys():

            if 'global_crs' in self.config['gis'].keys():
                self.global_crs = self.config['gis']['global_crs']
            else:
                self.global_crs = 'EPSG:4326'

            if 'local_crs' in self.config['gis'].keys():
                self.local_crs = self.config['gis']['local_crs']
            else:
                self.local_crs = None

        # load in elements
        if 'rivers' in self.config.keys():
            self.rivers = gpd.read_file(self.config['rivers']['path'])
            self.rivers = self.rivers.to_crs(self.global_crs)
            self.buffer_width = self.config['rivers']['buffer_width']
            self.grid_res = self.config['rivers']['grid_res']
        else:
            self.rivers = None

        if 'reservoirs' in self.config.keys():
            #self.reservoirs = gpd.read_file(self.config['reservoirs']['path'])
            #self.reservoirs = self.reservoirs.to_crs(self.global_crs)
            self.reservoirs = None
            warnings.warn("Support for reservoirs is not yet implemented. The reservoir shapefile will not be read into memory")
        else:
            self.reservoirs = None

        # make sure we have a local crs
        if self.local_crs is None:
            if self.rivers is not None:
                self.local_crs = self.rivers.estimate_utm_crs()
            elif self.reservoirs is not None:
                self.local_crs = self.reservoirs.estimate_utm_crs()
            else:
                raise UserWarning("Must provide a local crs or a river or reservoir shapefile to determine local crs")

        # make a buffered version of the rivers for use later on
        if self.rivers is not None:

            # buffer rivers within local crs
            self.buffered_rivers = self.rivers.copy()
            self.buffered_rivers.geometry = self.buffered_rivers.to_crs(self.local_crs).buffer(self.buffer_width)
            self.buffered_rivers = self.buffered_rivers.to_crs(self.global_crs)

    def report(self):
        print(f"Project Name: {self.name}")
        self.river_report()
        self.reservoir_report()

    def river_report(self):
        if self.rivers is not None:
            print(f"Number of rivers: {len(self.rivers)}")
            return self.rivers.head()
        else:
            print(f"Number of rivers: 0")
            return None
        
    def reservoir_report(self):
        if self.reservoirs is not None:
            print(f"Number of reservoirs: {len(self.reservoirs)}")
            return self.reservoirs.head()
        else:
            print(f"Number of reservoirs: 0")
            return None
        
    def grid_rivers(self, visualize=False):
            
        if self.rivers is not None:

            # break down large river system into grid for download and search
            area_bounds = self.buffered_rivers.unary_union.bounds
            grid_gdf = gpd.GeoDataFrame(geometry=geometry.grid_bounds(area_bounds, self.grid_res), crs=self.buffered_rivers.crs)
            grid_gdf = grid_gdf.loc[grid_gdf.intersects(self.buffered_rivers.unary_union)].reset_index(drop=True)
            grid_gdf['id'] = grid_gdf.index

            self.river_grid = grid_gdf

            if visualize:
                fig, ax = plt.subplots()
                self.rivers.plot(ax=ax)
                self.river_grid.plot(ax=ax, edgecolor='black', facecolor='None')
                ax.set_title('Full River System')
                fig.tight_layout()
                plt.show()

        else:
            warnings.warn("Cannot create grid because not rivers have been loaded.")


    def download_by_grid(self): # TODO: consider option to overwrite data already in folder or to start where download left off

        for i in self.river_grid.index:

            # grab coordinates of river
            coords = [(x, y) for x, y in self.river_grid.loc[i, 'geometry'].exterior.coords]
            id = self.river_grid.loc[i, 'id']

            # define bounds of data search
            startdate = date(*self.icesat2_startdate)
            enddate   = date(*self.icesat2_enddate)

            # define and if needed create directory for each river stretch
            download_directory = os.path.join(self.icesat2_dir, f"{id}")
            utils.ifnotmakedirs(download_directory)

            # make a simple aoi based on coordinates
            order_ids = icesat2.query(aoi=coords, startdate=startdate, enddate=enddate, earthdata_credentials=None, download_directory=download_directory, product='ATL13')


    def plot_data_by_grid(self, id):

        download_directory = os.path.join(self.icesat2_dir, f"{id}")
        files = list(os.listdir(download_directory))

        # start figure
        fig, ax = plt.subplots()

        # set boudns
        xmin, ymin, xmax, ymax = self.river_grid.loc[id, 'geometry'].bounds
        ax.set_xlim([xmin-0.1, xmax+0.1])
        ax.set_ylim([ymin-0.1, ymax+0.1])

        # plot rid and river
        self.river_grid.loc[[id]].plot(ax=ax, edgecolor='black', facecolor='None')
        self.rivers.plot(ax=ax)

        # load and plot all files for all tracks and crossings
        for file in files:
            infile= os.path.join(download_directory, file)
            for key in ['gt1l', 'gt1r', 'gt2l', 'gt2r', 'gt3l', 'gt3r']:

                data = icesat2.ATL13(infile, key)

                if data.check_height_data():
                    data_df = data.read()
                    data_gdf = gpd.GeoDataFrame(data_df, geometry=gpd.points_from_xy(data_df.lon, data_df.lat))
                    data_gdf = data_gdf.loc[data_gdf.within(self.buffered_rivers.unary_union)].reset_index(drop=True)

                    if len(data_gdf) > 0:
                        data_gdf['color'] = [int(date2num(i)) for i in data_gdf["date"].values]
                        data_gdf.plot(ax=ax, column='color', vmin=int(date2num(datetime(2019, 1, 1))), vmax=int(date2num(datetime(2024, 12, 31))), cmap='jet', alpha=0.5)

        fig.tight_layout()
        plt.show()


    def find_river_crossings(self):

        center_points = list()

        for id in self.river_grid.index:


            download_directory = os.path.join(self.icesat2_dir, f"{id}")

            if os.path.exists(download_directory):

                # load and plot all files for all tracks and crossings
                files = list(os.listdir(download_directory))        
                for file in files:
                    infile= os.path.join(download_directory, file)

                    for key in ['gt1l', 'gt1r', 'gt2l', 'gt2r', 'gt3l', 'gt3r']:

                        data = icesat2.ATL13(infile, key)
                        
                        if data.check_height_data():
                            data_df = data.read()
                            data_gdf = gpd.GeoDataFrame(data_df, geometry=gpd.points_from_xy(data_df.lon, data_df.lat))
                            data_gdf = data_gdf.loc[data_gdf.within(self.buffered_rivers.unary_union)].reset_index(drop=True)

                            if len(data_gdf) > 0:
                                center_points.append(data_gdf.loc[[0]])
                                
        center_points = pd.concat(center_points)
        
        self.crossings = center_points


    def plot_river_crossings(self):

        fig, ax = plt.subplots()

        # plot river
        self.rivers.plot(ax=ax)
        self.crossings['color'] = [int(date2num(i)) for i in self.crossings["date"].values]
        self.crossings.plot(ax=ax, column="color", vmin=int(date2num(datetime(2019, 1, 1))), vmax=int(date2num(datetime(2024, 12, 31))), cmap='Spectral')

        fig.tight_layout()
        plt.show()
