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

        if 'project' in self.config.keys():
            self.main_dir = self.config['project']["main_dir"]
            utils.ifnotmakedirs(self.main_dir)

            self.river_dir = os.path.join(self.main_dir, "rivers")
            utils.ifnotmakedirs(self.river_dir)

            self.reservoir_dir = os.path.join(self.main_dir, "reservoirs")
            utils.ifnotmakedirs(self.reservoir_dir)

        else:
            raise Warning("Project directory must be defined within configuration file")
            

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
            self.reservoirs = gpd.read_file(self.config['reservoirs']['path'])
            self.reservoirs = self.reservoirs.to_crs(self.global_crs)
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
            self.river_grid.to_file(os.path.join(self.river_dir, 'river_grid.shp'))

            if visualize:
                fig, ax = plt.subplots()
                self.rivers.plot(ax=ax)
                self.river_grid.plot(ax=ax, edgecolor='black', facecolor='None')
                ax.set_title('Full River System')
                fig.tight_layout()
                plt.show()

        else:
            warnings.warn("Cannot create grid because no rivers have been loaded.")

    def load_river_grid(self):

        self.river_grid = gpd.read_file(os.path.join(self.river_dir, 'river_grid.shp'))


    def download_by_grid(self): # TODO: consider option to overwrite data already in folder or to start where download left off

        for i in self.river_grid.index:

            # grab coordinates of river
            coords = [(x, y) for x, y in self.river_grid.loc[i, 'geometry'].exterior.coords]
            id = self.river_grid.loc[i, 'id']

            # define bounds of data search
            startdate = date(*self.icesat2_startdate)
            enddate   = date(*self.icesat2_enddate)

            # define and if needed create directory for each river stretch
            download_directory = os.path.join(self.icesat2_dir, rf"rivers\{id}")
            utils.ifnotmakedirs(download_directory)

            # make a simple aoi based on coordinates
            order_ids = icesat2.query(aoi=coords, startdate=startdate, enddate=enddate, earthdata_credentials=None, download_directory=download_directory, product='ATL13')


    def plot_data_by_grid(self, id):

        download_directory = os.path.join(self.icesat2_dir, rf"rivers\{id}")
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
                        data_gdf.plot(ax=ax, column='color', vmin=int(date2num(datetime(2019, 1, 1))), vmax=int(date2num(datetime(2024, 12, 31))), cmap=cm.batlow, alpha=0.5)

        fig.tight_layout()
        plt.show()

    
    def plot_timeseries_by_grid(self, id):

        # lists to hold plotting data
        height_list = list()
        date_list = list()
        lat_list = list()
        lon_list = list()

        # load and find average height for all values in river
        download_directory = os.path.join(self.icesat2_dir, rf"rivers\{id}")
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
                        
                        # find average height and date
                        height_list.append(data_gdf.height.mean())
                        date_list.append(data_gdf.date.mean())
                        lat_list.append(data_gdf.lat.mean())
                        lon_list.append(data_gdf.lon.mean())

        plot_df = pd.DataFrame({'date':date_list , 'height':height_list, 'lat':lat_list, 'lon':lon_list})

        # start figure
        fig, ax = plt.subplots()

        plot_df.plot(ax=ax, x='date', y='height', kind='scatter', c=plot_df['lon'], cmap=cm.batlow)

        fig.tight_layout()
        plt.show()



    def find_river_crossings(self):
        """Method to take all available icesat2 data and find the singular value to use for the river crossing.
        Currently takes the closest ATL13 value to the river centerline. River crossings are saved as a geopandas geodataframe.

        Raises
        ------
        Warning
            _description_
        """

        # list to hold the filtered center points and heights
        center_points = list()

        # process crossings by grid id
        for id in self.river_grid.index:

            # make sure that we have downloaded data for this grid cell, if not skip and process only what we have
            download_directory = os.path.join(self.icesat2_dir, rf"rivers\{id}")
            if os.path.exists(download_directory):

                # load date for all tracks and crossings
                files = list(os.listdir(download_directory))        
                for file in files:
                    for key in ['gt1l', 'gt1r', 'gt2l', 'gt2r', 'gt3l', 'gt3r']:
                        infile= os.path.join(download_directory, file)
                        data = icesat2.ATL13(infile, key)
                        
                        # make sure we have height data within the icesat file
                        if data.check_height_data():

                            # read data, turn into dataframe and filter by whats in the buffered river region
                            data_df = data.read()
                            data_gdf = gpd.GeoDataFrame(data_df, geometry=gpd.points_from_xy(data_df.lon, data_df.lat))
                            data_gdf = data_gdf.loc[data_gdf.within(self.buffered_rivers.unary_union)].reset_index(drop=True)

                            # ensure that we have at least two points in the river zone to process and centerline value
                            if len(data_gdf) > 1:

                                # take a line of the crossing and calculate its intersection with the river centerline
                                representative_line = shapely.LineString([(data_gdf['geometry'].iloc[0].x, data_gdf['geometry'].iloc[0].y), (data_gdf['geometry'].iloc[-1].x, data_gdf['geometry'].iloc[-1].y)])
                                crossing_point = shapely.intersection(representative_line, self.rivers.unary_union)

                                # there are a number of things that can happen when runnign the intersection so try to account for these
                                if crossing_point.is_empty:
                                    # then there was not a clean intersection, and we should do nothing with the data
                                    pass

                                elif crossing_point.geom_type == 'Point':
                                    # proceed as expected, one clear crossing means we can grab the data value at the point closest to the centerline or take the mean of the crossing
                                    index, _, _ = geometry.find_closest_geom(crossing_point, data_gdf.geometry.values)
                                    center_points.append(data_gdf.loc[[index]])

                                elif crossing_point.geom_type == 'MultiPoint':
                                    # It is possible that a passing crosses multiple points in a bending rivers so we grab each point clossest to the centerline at each respective crossing
                                    for point in crossing_point.geoms:
                                        index, _, _ = geometry.find_closest_geom(point, data_gdf.geometry.values)
                                        center_points.append(data_gdf.loc[[index]])

                                else:
                                    # something is not as expected to raise a warning to be investigated
                                    raise Warning("crossing_point is not a point or multilinestring warrenting investigation")

                            # if only one valid observation then take this (TODO: Consider dropping this, we could have one value on the edge of the river thats not really representative of the centerline)
                            elif len(data_gdf) > 0:
                                center_points.append(data_gdf.iloc[[0]])

        # push all of the crossing centerpoints into one geodataframe and reset the index
        self.crossings = pd.concat(center_points).reset_index(drop=True)
        self.crossings = self.crossings.set_crs(self.global_crs)

        # save copy of river crossings so we dont have to repeat the process
        self.crossings.to_file(os.path.join(self.river_dir, rf"icesat2_crossings.shp"))

    def load_river_crossings(self):
        self.crossings = gpd.read_file(os.path.join(self.river_dir, rf"icesat2_crossings.shp"))

    def plot_river_crossings(self):

        fig, ax = plt.subplots()

        # plot river
        self.rivers.plot(ax=ax)
        self.crossings['color'] = [int(date2num(i)) for i in self.crossings["date"].values]
        self.crossings.plot(ax=ax, column="color", vmin=int(date2num(datetime(2019, 1, 1))), vmax=int(date2num(datetime(2024, 12, 31))), cmap=cm.batlow)

        fig.tight_layout()
        plt.show()



    ##### Functions associated with the reservoirs (the download process is handeled a bit differently, perhaps there is a better way to make an elements class in which rivers and reservoirs inherit from?)
    def download_by_reservoir(self, start_id=0):

        for id in self.reservoirs.index[start_id:]:

            # grab coordinates of river
            coords = [(x, y) for x, y in self.reservoirs.loc[id, 'geometry'].envelope.exterior.coords]
            #id = self.reservoirs.loc[i, 'id']

            # define bounds of data search
            startdate = date(*self.icesat2_startdate)
            enddate   = date(*self.icesat2_enddate)

            # define and if needed create directory for each river stretch
            download_directory = os.path.join(self.icesat2_dir, rf"reservoirs\{id}")
            utils.ifnotmakedirs(download_directory)

            # make a simple aoi based on coordinates
            order_ids = icesat2.query(aoi=coords, startdate=startdate, enddate=enddate, earthdata_credentials=None, download_directory=download_directory, product='ATL13')

    def plot_data_by_reservoir(self, id):

        download_directory = os.path.join(self.icesat2_dir, rf"reservoirs\{id}")
        files = list(os.listdir(download_directory))

        # start figure
        fig, ax = plt.subplots()

        # set boudns
        xmin, ymin, xmax, ymax = self.reservoirs.loc[id, 'geometry'].bounds
        ax.set_xlim([xmin-0.1, xmax+0.1])
        ax.set_ylim([ymin-0.1, ymax+0.1])

        # plot reservoir
        self.reservoirs.loc[[id]].plot(ax=ax, edgecolor='black', facecolor='None')

        # load and plot all files for all tracks and crossings
        for file in files:
            infile= os.path.join(download_directory, file)
            for key in ['gt1l', 'gt1r', 'gt2l', 'gt2r', 'gt3l', 'gt3r']:

                data = icesat2.ATL13(infile, key)

                if data.check_height_data():
                    data_df = data.read()
                    data_gdf = gpd.GeoDataFrame(data_df, geometry=gpd.points_from_xy(data_df.lon, data_df.lat))
                    data_gdf = data_gdf.loc[data_gdf.within(self.reservoirs.unary_union)].reset_index(drop=True)

                    if len(data_gdf) > 0:
                        data_gdf['color'] = [int(date2num(i)) for i in data_gdf["date"].values]
                        data_gdf.plot(ax=ax, column='color', vmin=int(date2num(datetime(2019, 1, 1))), vmax=int(date2num(datetime(2024, 12, 31))), cmap=cm.batlow, alpha=0.5)

        fig.tight_layout()
        plt.show()


    def plot_timeseries_by_reservoir(self, id):

        # lists to hold plotting data
        height_list = list()
        date_list = list()
        lat_list = list()
        lon_list = list()

        # load and find average height for all values in river
        download_directory = os.path.join(self.icesat2_dir, rf"reservoirs\{id}")
        files = list(os.listdir(download_directory))
        for file in files:
            infile= os.path.join(download_directory, file)
            for key in ['gt1l', 'gt1r', 'gt2l', 'gt2r', 'gt3l', 'gt3r']:

                data = icesat2.ATL13(infile, key)

                if data.check_height_data():
                    data_df = data.read()
                    data_gdf = gpd.GeoDataFrame(data_df, geometry=gpd.points_from_xy(data_df.lon, data_df.lat))
                    data_gdf = data_gdf.loc[data_gdf.within(self.reservoirs.unary_union)].reset_index(drop=True)

                    if len(data_gdf) > 0:
                        
                        # find average height and date
                        height_list.append(data_gdf.height.mean())
                        date_list.append(data_gdf.date.mean())
                        lat_list.append(data_gdf.lat.mean())
                        lon_list.append(data_gdf.lon.mean())

        plot_df = pd.DataFrame({'date':date_list , 'height':height_list, 'lat':lat_list, 'lon':lon_list})

        # start figure
        fig, ax = plt.subplots()

        plot_df.plot(ax=ax, x='date', y='height', kind='scatter', c=plot_df['lon'], cmap=cm.batlow)

        fig.tight_layout()
        plt.show()