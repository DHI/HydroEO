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

from altimetry.sat_utils import swot, icesat2, sentinel
from altimetry import geometry, utils

from tqdm import tqdm



@dataclass
class System:

    gdf: gpd.GeoDataFrame
    dirs: dict

    def report(self):
        print(f"Number of {self.type}: {len(self.gdf)}")
        return self.gdf.head()
    
    def download_altimetry(self, product:str, startdate: tuple, enddate:tuple, credentials=None, grid:bool=False, start_index=0):

        supported_products = ['ATL13', 'S3', 'S6', 'SWOT_LAKE']

        # check if the requested product is supported
        product = product.upper()
        if product not in supported_products:
            raise ValueError(f'"{product}" is not accepted as a valid download product. Please provide a valid product.')


        ### Set the download geometry for the search bounds
        # set grid as download bounds if needed
        if grid:
            if hasattr(self, 'grid'):
                download_gdf = self.grid.copy()
            else:
                raise NameError('Object has no atttribute: grid. Make sure to create grid before enabling this option')
        else:
            download_gdf = self.download_gdf.copy()

        # unpack and format the download dates
        startdate = date(*startdate)
        enddate   = date(*enddate)

        ########################################################################################################################
        ##### We have to handle different downloads differently based on the way products are delivered in "granules"
        ### Swot data is provided in large granuals that likley cover most of the full AOI and downloading per reservoir or grid may be redundant
        ### Icesat-2 data can be easily subsetted within the product order so an individual download make sense
        ### Sentinel 3 and 6 data is provided in smaller granules so we can also request per small area and process one at a time

        # First check if it is swot data that should be downloaded for the area
        if product == 'SWOT_LAKE':

            # define and if needed create directory for each download geometry
            download_dir = os.path.join(self.dirs['swot'], rf"{self.type}")
            utils.ifnotmakedirs(download_dir)

            # grab coordinates of full area of interest
            coords = [(x, y) for x, y in download_gdf.unary_union.envelope.exterior.coords]

            # query all available data
            results = swot.query(aoi=coords, startdate=startdate, enddate=enddate, earthdata_credentials=credentials, product='SWOT_L2_HR_LakeSP_2.0')

            # loop through results and check if it is the product we want (right now it must include "Prior"), later we can try and get the observed water bodies and match them
            for result in results:

                # just based on that we know where to find the long product name and download link, may need to change if naming conventions change!
                product_pieces = result.data_links()[0].split('/')[-1].split('_')
                if 'Prior' in product_pieces:

                    # download the individual file
                    files = swot.download(result, download_directory=download_dir)

                    # we want to subset the downloaded file to only include known waterbodies
                    swot.subset_by_id(files, download_gdf.dl_id.astype(int).values)


        elif product in ['S3', 'S6', 'ATL13']:

            # loop through download geometry and download data
            for i in download_gdf.index[start_index:]:

                # grab coordinates of geometry
                coords = [(x, y) for x, y in download_gdf.loc[i, 'geometry'].envelope.exterior.coords]
                id = download_gdf.loc[i, 'dl_id']

                print(id)

                if product == 'ATL13':

                    # define and if needed create directory for each download geometry
                    download_dir = os.path.join(self.dirs['icesat2'], rf"{self.type}\{id}")
                    utils.ifnotmakedirs(download_dir)

                    # query and download data
                    _ = icesat2.query(aoi=coords, startdate=startdate, enddate=enddate, earthdata_credentials=credentials, download_directory=download_dir, product='ATL13')

                elif product == 'S3':
                
                    # define and if needed create directory for each download geometry
                    download_dir = os.path.join(self.dirs['sentinel3'], rf"{self.type}\{id}")
                    utils.ifnotmakedirs(download_dir)

                    # query and download data
                    sentinel.query(aoi=coords, startdate=startdate, enddate=enddate, creodias_credentials=credentials, download_directory=download_dir, product='S3')

                    # once we have finished downloading all data for the aoi, we need to unzip everything and then subset it
                    utils.unzip_dir_files(download_dir, download_dir)
                    sentinel.subset(aoi = coords, download_dir = download_dir, dest_dir = download_dir, product = product)

                    # clean up zip and unzipped folders keeping only the remaining subsetted data
                    utils.remove_non_exts(download_dir, ['.nc', '.log'])

                elif product == 'S6':
                
                    # define and if needed create directory for each download geometry
                    download_dir = os.path.join(self.dirs['sentinel6'], rf"{self.type}\{id}")
                    utils.ifnotmakedirs(download_dir)

                    # query and download data
                    sentinel.query(aoi=coords, startdate=startdate, enddate=enddate, creodias_credentials=credentials, download_directory=download_dir, product='S6')

                    # once we have finished downloading all data for the aoi, we need to unzip everything and then subset it
                    utils.unzip_dir_files_with_ext(download_dir, download_dir, '.nc')
                    sentinel.subset(aoi = coords, download_dir = download_dir, dest_dir = download_dir, product = product)

                    # clean up zip and unzipped folders keeping only the remaining subsetted data
                    utils.remove_non_exts(download_dir, ['.nc', '.log'])

    def get_unfiltered_timeseries_by_id(self, id):

        data_dir = os.path.join(self.dirs['output'], f"{id}", "raw_observations")

        # loop through raw directory and load in all raw observations
        df_list = list()
        for file in os.listdir(data_dir):
            if file.endswith('.shp'):

                gdf_path = os.path.join(data_dir, file)
                gdf = gpd.read_file(gdf_path)
                df = gdf.drop(columns = ['geometry'])
                df_list.append(df)

        # concatenate everything into one dataframe for plotting
        df = pd.concat(df_list)

        return df


    # def load_crossings(self):
    #     self.crossings = gpd.read_file(os.path.join(self.dirs['output'], rf"icesat2_crossings.shp"))


    # def map_all_crossings(self):

    #     fig, ax = plt.subplots()

    #     # plot river
    #     self.gdf.plot(ax=ax)
    #     self.crossings.plot(ax=ax, column="height", cmap=cm.batlow, legend=True, legend_kwds={'label': 'Height (m)'})

    #     fig.tight_layout()
    #     plt.show()


    
    # def plot_timeseries_by_id(self, id, summarize=True):

    #     # start figure
    #     fig, ax = plt.subplots()

    #     # extract this grids crossings from the main river crossing dataframe
    #     data_gdf = self.crossings.loc[self.crossings.dl_id == id]

    #     if summarize:
    #         # apply simple mean groupby to make daily estimate (There is no adjustment for bias here just an example of how to get to a single daily value quickly)
    #         data_gdf['date'] = data_gdf.date.dt.floor('d')
    #         data_gdf = data_gdf[['date', 'height', 'lon']].groupby(by='date').mean().reset_index()

    #     # make a scatter plot with this data
    #     data_gdf.plot(ax=ax, x='date', y='height', kind='scatter', c=data_gdf['lon'], cmap=cm.batlow)

    #     fig.tight_layout()
    #     plt.show()



@dataclass
class Reservoirs(System):

    def __post_init__(self):

        self.type = 'reservoirs'

        self.dirs['output'] = os.path.join(self.dirs['main'], self.type)
        utils.ifnotmakedirs(self.dirs['output'])

        self.geom_type = self.gdf.loc[0, 'geometry'].geom_type
        # self.gdf['dl_id'] = self.gdf.index


    def assign_pld_id(self, local_crs, max_distance):

        # load the pld
        pld = gpd.read_file(self.dirs['pld'])

        # perform the spatial join
        joined_gdf = gpd.sjoin_nearest(self.gdf.to_crs(local_crs), pld.to_crs(local_crs), how='left', max_distance=max_distance, distance_col='dist_to_pld')
        joined_gdf = joined_gdf.to_crs(self.gdf.crs)

        # rename the joined columns to keep it clear
        joined_gdf = joined_gdf.rename(columns={'lake_id':'prior_lake_id', 'res_id':'prior_res_id'})

        # reset the reservoir data frame to include the changes
        self.gdf = joined_gdf


    def flag_missing_priors(self): # simple function to report what reservoirs do not have PLD shapes

        present = self.gdf.loc[self.gdf.prior_lake_id > 0].reset_index(drop=True)
        missing = self.gdf.loc[self.gdf.prior_lake_id.isnull()].reset_index(drop=True)

        present.to_file(os.path.join(self.dirs['output'], 'present_in_pld.shp'))
        missing.to_file(os.path.join(self.dirs['output'], 'missing_in_pld.shp'))

        print(f"Out of the {len(self.gdf)} reservoirs, {len(present)} area present and {len(missing)} are missing from the PLD.")


    def assign_reservoir_polygons(self): # function to set a reservoir polygon to the reservoirs in the system based on a pld id TODO: need to account for getting shapes for non hits

        # load the pld
        pld = gpd.read_file(self.dirs['pld'])

        download_gdf = self.gdf.loc[self.gdf.prior_lake_id > 0].reset_index(drop=True)

        geometries = [pld.loc[pld.lake_id == lake_id].geometry.values[0] for lake_id in download_gdf.prior_lake_id]

        download_gdf['geometry'] = geometries

        download_gdf['dl_id'] = download_gdf.prior_lake_id.astype(int)

        self.download_gdf = download_gdf



    def extract_product_timeseries(self, products: list):

        if 'icesat2' in products:

            for id in tqdm(self.download_gdf.dl_id, desc='Extracting ICESat-2 ATL13 product'):

                # filter gdf to only row with id
                sub_gdf = self.download_gdf.loc[self.download_gdf.dl_id == id]
                    
                # prep export folder and paths
                download_dir = os.path.join(self.dirs['icesat2'], rf"{self.type}\{id}")
                if os.path.exists(download_dir):

                    dst_dir = os.path.join(self.dirs['output'], f"{id}", "raw_observations")
                    utils.ifnotmakedirs(dst_dir)
                    dst_path = os.path.join(dst_dir, 'icesat2.shp')

                    # extract observations within bounds and save as timeseries csv
                    icesat2.extract_observations(src_dir=download_dir, dst_path=dst_path, features=sub_gdf)

        if 'sentinel3' in products:

            for id in tqdm(self.download_gdf.dl_id, desc='Extracting Sentinel-3 product'):

                # filter gdf to only row with id
                sub_gdf = self.download_gdf.loc[self.download_gdf.dl_id == id]

                # prep export folder and paths
                download_dir = os.path.join(self.dirs['sentinel3'], rf"{self.type}\{id}")
                if os.path.exists(download_dir):

                    dst_dir = os.path.join(self.dirs['output'], f"{id}", "raw_observations")
                    utils.ifnotmakedirs(dst_dir)
                    dst_path = os.path.join(dst_dir, 'sentinel3.shp')

                    # extract observations within bounds and save as timeseries csv
                    sentinel.extract_observations(src_dir=download_dir, dst_path=dst_path, features=sub_gdf)

        if 'sentinel6' in products:

            for id in tqdm(self.download_gdf.dl_id, desc='Extracting Sentinel-6 product'):

                # filter gdf to only row with id
                sub_gdf = self.download_gdf.loc[self.download_gdf.dl_id == id]

                # prep export folder and paths
                download_dir = os.path.join(self.dirs['sentinel6'], rf"{self.type}\{id}")
                if os.path.exists(download_dir):

                    dst_dir = os.path.join(self.dirs['output'], f"{id}", "raw_observations")
                    utils.ifnotmakedirs(dst_dir)
                    dst_path = os.path.join(dst_dir, 'sentinel6.shp')

                    # extract observations within bounds and save as timeseries csv
                    sentinel.extract_observations(src_dir=download_dir, dst_path=dst_path, features=sub_gdf)

        if 'swot' in products:

            download_dir = os.path.join(self.dirs['swot'], rf"{self.type}")

            # extract observations within bounds and save as timeseries csv, slightly differnt format for downloading here due to the natrure of swot data
            swot.extract_observations(src_dir=download_dir, dst_dir=self.dirs['output'], dst_file_name='swot.shp', features=self.download_gdf)


    # def map_crossings_by_id(self, id):

    #     # start figure
    #     fig, ax = plt.subplots()

    #     # set boudns
    #     xmin, ymin, xmax, ymax = self.gdf.loc[id, 'geometry'].bounds
    #     ax.set_xlim([xmin-0.1, xmax+0.1])
    #     ax.set_ylim([ymin-0.1, ymax+0.1])

    #     # plot reservoir
    #     self.gdf.loc[[id]].plot(ax=ax, edgecolor='black', facecolor='None')

    #     data_gdf = self.crossings.loc[self.crossings.dl_id == id]
        
    #     if len(data_gdf) > 0:
    #         data_gdf['color'] = [int(date2num(i)) for i in data_gdf["date"].values]
    #         data_gdf.plot(ax=ax, column='color', vmin=int(date2num(datetime(2019, 1, 1))), vmax=int(date2num(datetime(2024, 12, 31))), cmap=cm.batlow, alpha=0.5, legend=True)

    #     fig.tight_layout()
    #     plt.show()




@dataclass
class Rivers(System):

    buffer_width : float
    grid_res : float

    def __post_init__(self):

        self.type = 'rivers'

        self.dirs['output'] = os.path.join(self.dirs['main'], self.type)
        utils.ifnotmakedirs(self.dirs['output'])

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
        grid_gdf.to_file(os.path.join(self.dirs['output'], f'grid.shp'))
        self.grid = grid_gdf


    def load_grid(self):

        self.grid = gpd.read_file(os.path.join(self.dirs['output'], f'grid.shp'))
            

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
            download_directory = os.path.join(self.dirs['icesat2'], rf"{self.type}\{id}")
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
        self.crossings.to_file(os.path.join(self.dirs['output'], rf"icesat2_crossings.shp"))


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
