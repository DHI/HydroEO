from dataclasses import dataclass
import warnings

import os
import pandas as pd
import geopandas as gpd
import shapely

from cmcrameri import cm
import matplotlib.pyplot as plt
from matplotlib.dates import date2num

from datetime import date, datetime

from altimetry import icesat2, geometry, utils

from tqdm import tqdm



@dataclass
class System:

    gdf: gpd.GeoDataFrame
    icesat2_dir: str
    output_dir: str

    def __post_init__(self):

        utils.ifnotmakedirs(self.output_dir)

    def report(self):
        print(f"Number of {self.type}: {len(self.gdf)}")
        return self.gdf.head()
    
    def download_icesat2(self, startdate: tuple, enddate:tuple, grid:bool=False, start_index=0):

        # set grid as download bounds if needed
        if grid:
            if hasattr(self, 'grid'):
                download_gdf = self.grid
            else:
                raise NameError('Object has no atttribute: grid. Make sure to create grid before enabling this option')
        else:
            download_gdf = self.gdf

        # unpack and format the download dates
        startdate = date(*startdate)
        enddate   = date(*enddate)

        # loop through download geometry and download data
        for i in download_gdf.index[start_index:]:

            # grab coordinates of geometry
            coords = [(x, y) for x, y in download_gdf.loc[i, 'geometry'].envelope.exterior.coords]
            id = download_gdf.loc[i, 'dl_id']

            # define and if needed create directory for each download geometry
            download_dir = os.path.join(self.icesat2_dir, rf"{self.type}\{id}")
            utils.ifnotmakedirs(download_dir)

            # make a simple aoi based on coordinates
            order_ids = icesat2.query(aoi=coords, startdate=startdate, enddate=enddate, earthdata_credentials=None, download_directory=download_dir, product='ATL13')


    def load_crossings(self):
        self.crossings = gpd.read_file(os.path.join(self.output_dir, rf"icesat2_crossings.shp"))


    def map_all_crossings(self):

        fig, ax = plt.subplots()

        # plot river
        self.gdf.plot(ax=ax)
        self.crossings.plot(ax=ax, column="height", cmap=cm.batlow, legend=True, legend_kwds={'label': 'Height (m)'})

        fig.tight_layout()
        plt.show()


    
    def plot_timeseries_by_id(self, id):

        # start figure
        fig, ax = plt.subplots()

        # extract this grids crossings from the main river crossing dataframe
        data_gdf = self.crossings.loc[self.crossings.dl_id == id]

        # make a scatter plot with this data
        data_gdf.plot(ax=ax, x='date', y='height', kind='scatter', c=data_gdf['lon'], cmap=cm.batlow)

        fig.tight_layout()
        plt.show()



@dataclass
class Rivers(System):

    buffer_width : float
    grid_res : float

    def __post_init__(self):

        self.type = 'rivers'

    def make_buffer(self, local_crs):

        # make a buffered version of the rivers for use later on
        self.buffered_gdf = self.gdf.copy()
        self.buffered_gdf.geometry = self.buffered_gdf.to_crs(local_crs).buffer(self.buffer_width)
        self.buffered_gdf = self.buffered_gdf.to_crs(self.gdf.crs) # return the crs to that of the rivers file


    def make_grid(self):

        # break down large river system into grid for download and search
        area_bounds = self.buffered_gdf.unary_union.bounds
        
        # make a dataframe from a list of grid polygons
        grids = geometry.grid_bounds(area_bounds, self.grid_res)
        grid_gdf = gpd.GeoDataFrame(geometry=grids, crs=self.buffered_gdf.crs)

        # keep only the grids that intersect with the river line
        grid_gdf = grid_gdf.loc[grid_gdf.intersects(self.buffered_gdf.unary_union)].reset_index(drop=True)
        grid_gdf['dl_id'] = grid_gdf.index # set a grid id based on the remaing grids

        # save grid and set attribute
        grid_gdf.to_file(os.path.join(self.output_dir, f'grid.shp'))
        self.grid = grid_gdf


    def load_grid(self):

        self.grid = gpd.read_file(os.path.join(self.output_dir, f'grid.shp'))
            

    def visualize_grid(self):

        fig, ax = plt.subplots()
        self.gdf.plot(ax=ax)
        self.grid.plot(ax=ax, edgecolor='black', facecolor='None')
        ax.set_title('Full River System')
        fig.tight_layout()
        plt.show()


    def extract_crossings(self):
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
        for id in tqdm(self.grid.index):

            # make sure that we have downloaded data for this grid cell, if not skip and process only what we have
            download_directory = os.path.join(self.icesat2_dir, rf"{self.type}\{id}")
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
                            data_gdf = data_gdf.loc[data_gdf.within(self.buffered_gdf.unary_union)].reset_index(drop=True)
                            data_gdf['dl_id'] = id

                            # ensure that we have at least two points in the river zone to process and centerline value
                            if len(data_gdf) > 1:

                                # take a line of the crossing and calculate its intersection with the river centerline
                                representative_line = shapely.LineString([(data_gdf['geometry'].iloc[0].x, data_gdf['geometry'].iloc[0].y), (data_gdf['geometry'].iloc[-1].x, data_gdf['geometry'].iloc[-1].y)])
                                crossing_point = shapely.intersection(representative_line, self.gdf.unary_union)

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
        self.crossings = self.crossings.set_crs(self.gdf.crs)

        # save copy of river crossings so we dont have to repeat the process
        self.crossings.to_file(os.path.join(self.output_dir, rf"icesat2_crossings.shp"))


    def map_crossings_by_grid(self, id):

        # start figure
        fig, ax = plt.subplots()

        # set bounds
        xmin, ymin, xmax, ymax = self.grid.loc[id, 'geometry'].bounds
        ax.set_xlim([xmin-0.1, xmax+0.1])
        ax.set_ylim([ymin-0.1, ymax+0.1])

        # plot grid and river
        self.grid.loc[[id]].plot(ax=ax, edgecolor='black', facecolor='None')
        self.gdf.plot(ax=ax)

        # extract this grids crossings from the main river crossing dataframe
        data_gdf = self.crossings.loc[self.crossings.dl_id == id]

        if len(data_gdf) > 0:
            data_gdf['color'] = [int(date2num(i)) for i in data_gdf["date"].values]
            data_gdf.plot(ax=ax, column='color', vmin=int(date2num(datetime(2019, 1, 1))), vmax=int(date2num(datetime(2024, 12, 31))), cmap=cm.batlow, alpha=0.5, legend=True)

        fig.tight_layout()
        plt.show()

    
    def crossings_to_rivers(self, riv_key):

        # add the river name (or id) to the closest crossing point
        self.crossings = gpd.sjoin_nearest(self.crossings, self.gdf[[riv_key, 'geometry']], how='left')


    def plot_river_profile(self, riv_key, name, delta):

        # extract the river associated with the name
        river = self.gdf.loc[self.gdf[riv_key] == name]

        # extract the crossings matching the river key name
        crossings = self.crossings.loc[self.crossings[riv_key] == name] # this should be a gdf
        
        #### Start the figure
        # now plot the river profile by distance downriver
        fig, main_ax = plt.subplots(1, 2, figsize=(10, 5))

        ax = main_ax[0]
        river.plot(ax=ax)
        crossings.plot(ax=ax, column='height', cmap=cm.batlow, legend=True, legend_kwds={'label': 'Height (m)'})

        # Now we convert to local crs and do the calculateions for where on the river the crossing is
        river = river.to_crs(river.estimate_utm_crs())
        river_line = river.geometry.values[0]

        crossings = crossings.to_crs(crossings.estimate_utm_crs())

        # interpolate the river linestring into points
        if river_line.geom_type == 'MultiLinesString':

            # TODO: convert to linestring somehow?
            pass

        # now assume we are dealing with a linestring
        # use shapely interpolate to get the points at a known interval downstream
        river_points, point_dist = geometry.line_to_points(river_line, delta)

        # find the nearest linepoint to the crossings and record the distance along the river
        for i in crossings.index:
            point = crossings.loc[i, 'geometry']
            index, _, _ = geometry.find_closest_geom(point, river_points)
            crossings.loc[i, 'dist'] = point_dist[index]

        # make the 1d profile
        ax = main_ax[1]
        plot = ax.scatter(crossings.dist/1000, crossings.height, c=[d.month for d in crossings.date], cmap=cm.brocO)
        #plot = crossings.plot(ax=ax, x='dist', y='height', kind='scatter', c=[d.month for d in crossings.date], cmap=cm.batlow, legend=True)
        cbar = fig.colorbar(plot, label='Month')
        ax.set_xlabel("Distance downriver (km)")
        ax.set_ylabel("Height (m)")

        fig.tight_layout()
        plt.show()


@dataclass
class Reservoirs(System):

    def __post_init__(self):

        self.type = 'reservoirs'
        self.gdf['dl_id'] = self.gdf.index


    def extract_crossings(self):

        gdf_list = list()

        for id in tqdm(self.gdf.index):

            # load and plot all files for all tracks and crossings
            download_directory = os.path.join(self.icesat2_dir, rf"{self.type}\{id}")
            if os.path.exists(download_directory):

                files = list(os.listdir(download_directory))
                for file in files:
                    for key in ['gt1l', 'gt1r', 'gt2l', 'gt2r', 'gt3l', 'gt3r']:

                        infile= os.path.join(download_directory, file)
                        data = icesat2.ATL13(infile, key)

                        if data.check_height_data():
                            data_df = data.read()
                            data_gdf = gpd.GeoDataFrame(data_df, geometry=gpd.points_from_xy(data_df.lon, data_df.lat))
                            data_gdf = data_gdf.loc[data_gdf.within(self.gdf.unary_union)].reset_index(drop=True)

                            if len(data_gdf) > 0:

                                # if we have data for the reservoir add it to the full reservoir dataframe
                                # add a column keeping track of the reservoir id
                                data_gdf['dl_id'] = id
                                gdf_list.append(data_gdf)

        # once all reservoirs are cycled, concatenate them all into a main data frame of valid information
        self.crossings = pd.concat(gdf_list).reset_index(drop=True)
        self.crossings = self.crossings.set_crs(self.gdf.crs)
        self.crossings.to_file(os.path.join(self.output_dir, f'icesat2_crossings.shp'))


    def map_crossings_by_id(self, id):

        # start figure
        fig, ax = plt.subplots()

        # set boudns
        xmin, ymin, xmax, ymax = self.gdf.loc[id, 'geometry'].bounds
        ax.set_xlim([xmin-0.1, xmax+0.1])
        ax.set_ylim([ymin-0.1, ymax+0.1])

        # plot reservoir
        self.gdf.loc[[id]].plot(ax=ax, edgecolor='black', facecolor='None')

        data_gdf = self.crossings.loc[self.crossings.dl_id == id]
        
        if len(data_gdf) > 0:
            data_gdf['color'] = [int(date2num(i)) for i in data_gdf["date"].values]
            data_gdf.plot(ax=ax, column='color', vmin=int(date2num(datetime(2019, 1, 1))), vmax=int(date2num(datetime(2024, 12, 31))), cmap=cm.batlow, alpha=0.5, legend=True)

        fig.tight_layout()
        plt.show()