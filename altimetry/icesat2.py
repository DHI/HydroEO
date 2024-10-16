
from dataclasses import dataclass
import os
import warnings
import typing
import datetime
import icepyx as ipx

import h5py
import numpy as np
import pandas as pd

"""
# import necessary packages


import netCDF4



import rasterio
import tables as tb

from datetime import datetime, timedelta
from typing import List, Tuple
"""





def query(aoi: list, startdate: datetime.date, enddate: datetime.date, earthdata_credentials: tuple, download_directory: str, product: str  ='ATL13') -> object:
    """Query and download ICESat-2 products through NASA NSIDC DAAC
    Reference: Scheick, J. et al., (2019). icepyx: Python tools for obtaining and working with ICESat-2 data. https://github.com/icesat2py/icepyx.
    TODO:   Lookout for icepyx v1.x is being deprecations late 2024; the back-end systems on which it relies will be shut down as of late 2024.
            At that time, upgrade to icepyx v2.x, which uses the new NASA Harmony back-end, will be required.

    Parameters
    ----------
    aoi : list
        List of Polygon coordinates for the area to search within. Should be provided as coordinate pairs in decimal degrees as 
        [(longitude1, latitude1), (longitude2, latitude2), … (longitude_n,latitude_n), (longitude1,latitude1)] or 
        [longitude1, latitude1, longitude2, latitude2, … longitude_n,latitude_n, longitude1,latitude1].
        Your list must contain at least four points, where the first and last are identical.
    startdate : datetime.date
        The earliest date for which to search for data.
    enddate : datetime.date
        The latest date for which to search for data.
    earthdata_credentials : tuple
        Login credentials for earth data user account. First entry should be username while last entry should be password.
    download_directory : str
        Path to the folder where icesat data should be downloaded
    product : str, optional
        ICESat-2 data product ID, also known as “short name” (e.g. ATL03). Available data products can be found at: https://nsidc.org/data/icesat-2/data-sets, by default 'ATL13'

    Returns
    -------
    object
        icepyx Query object containing the order ids associated with the request.
    """
    
    # define query object
    region_a = ipx.Query(product = product, 
                         spatial_extent = aoi,
                         date_range = [startdate.strftime('%Y-%m-%d'), enddate.strftime('%Y-%m-%d')],
                         start_time='00:00:00',
                         end_time='23:59:59'
                        )
    
    # provide credentials
    if earthdata_credentials is not None:
        un, pw = earthdata_credentials
        region_a.earthdata_login(un, pw)

    # order granules
    region_a.order_granules()
    region_a.download_granules(download_directory)

    return region_a.granules.orderIDs


@dataclass
class ATL13:
    """
    ICESat-2 utility class for reading ATL13 data from .h5 files.
    
    Arguments:
    ----------
        infile: ATL03 file path (.h5)
        track_key: ICESat-2 ground track key
    """
    infile : str
    track_key : str

    def __post_init__(self):
    
        self.file_name = os.path.basename(self.infile)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            with h5py.File(self.infile, 'r') as src:

                # Segment track heights
                self.setattr_from_file(src, 'height_seg', 'ht_ortho')
                self.setattr_from_file(src, 'lat_seg', 'segment_lat')
                self.setattr_from_file(src, 'lon_seg', 'segment_lon')
                
                # Water body information
                self.setattr_from_file(src, 'wb_id', 'inland_water_body_id')
                self.setattr_from_file(src, 'wb_size', 'inland_water_body_size')
                self.setattr_from_file(src, 'wb_type', 'inland_water_body_type')

                # Segment information - error and quality
                self.setattr_from_file(src, 'wb_slope', 'segment_slope_trk_bdy')
                self.setattr_from_file(src, 'height_err', 'err_ht_water_surf')
                self.setattr_from_file(src, 'quality_seg', 'segment_quality')
                                
                # Tide/DEM corrections
                self.setattr_from_file(src, 'geoid_track', 'segment_geoid')
                self.setattr_from_file(src, 'geoid_corr_track', 'segment_geoid_free2mean')
                self.setattr_from_file(src, 'dem_track', 'segment_dem_ht')

                # Saturation fraction - The fraction of pulses within the short segment determined to be nearly saturated based on ATL03 geosegment rate input.
                self.setattr_from_file(src, 'sat_frac_track', 'segment_near_sat_fract')

                # time and date attributes
                if 'delta_time' in src[self.track_key].keys():
                    delta_time_seg = np.asarray(src[self.track_key]['delta_time'], float)
                    atlas_offset   = np.asarray(src['ancillary_data']['atlas_sdp_gps_epoch'])[0]
                    self.date = [(datetime.datetime(1980, 1, 6) + datetime.timedelta(seconds=c2_time+atlas_offset)) for c2_time in delta_time_seg]
                else:
                    self.date=None

    def setattr_from_file(self, src: h5py.File, name: str, attribute:str):

        if attribute in src[self.track_key].keys():
            self.__setattr__(name, np.asarray(src[self.track_key][attribute]))
        else:
            self.__setattr__(name, None)

    def check_height_data(self) -> bool:
            
        return self.height_seg is not None

    def read(self) -> pd.DataFrame:
        """
        Reads data (height, coordinates and geophysical corrections) from ATL03 file,
        and removes points outside buffered DEM elevation bands.
                
        Returns:
        ----------
            track_df: Dataframe with requested variables
        """
                
        track_df = pd.DataFrame(data={
            'height' : self.height_seg,
            'lat' : self.lat_seg,
            'lon' : self.lon_seg,
            'date' : self.date,
            'wb_type' : self.wb_type,
            'wb_size' : self.wb_size,
            'wb_id' : self.wb_id,
            'dem' : self.dem_track,
            'sat_frac_track': self.sat_frac_track,
            'beam' : self.track_key,
            'file_name' : self.file_name
            
        }, index = np.arange(len(self.height_seg)))
                
        return track_df


# def read_ice_ATL13_nc(ice_file, water_mask_file, dem_file):
#     """

#     Parameters
#     ----------
#     s3_folder : STRING
#         Foldername and path for Level 2 netcdf files
#     source : STRING
#              "GPOD" or "SciHub" - source of Sentinel-3 netcdf file

#     Returns
#     -------
#     height : Array
#         Heights in track referenced to geoid (EGM2008)

#     """

#     beams = ['gt1l', 'gt1r', 'gt2l', 'gt2r', 'gt3l', 'gt3r']
#     nc = netCDF4.Dataset(ice_file)
#     lat = 'segment_lat'
#     lon = 'segment_lon'
#     # height = 'ht_water_surf'
#     tai = 'delta_time'
#     geoid = 'ht_ortho'
#     atlas_offset = nc['ancillary_data']['atlas_sdp_gps_epoch'][:].filled()

#     vs = pd.DataFrame()
#     for beam in beams:
#         if geoid in nc[beam].variables:
#             _vs = pd.DataFrame({'height': nc[beam][geoid][:].filled(),
#                                 'wb_id': nc[beam]['inland_water_body_id'][:].filled(),
#                                 'wb_size': nc[beam]['inland_water_body_size'][:].filled(),
#                                 'wb_type': nc[beam]['inland_water_body_type'][:].filled(),
#                                 'lat': nc[beam][lat][:].filled(),
#                                 'lon': nc[beam][lon][:].filled(),
#                                 'date': [(datetime(1980, 1, 6) +
#                                           timedelta(seconds=c2_time+atlas_offset)).date() for
#                                          c2_time in nc[beam][tai][:].filled()],
#                                 'beam': np.repeat(beam, len(nc[beam][lat][:].filled()))})
    
#             vs = vs.append(_vs)
#     vs['file_name'] = np.repeat(os.path.split(ice_file)[1], len(vs['height']))

#     # Add water occurence value if possible
#     if os.path.isfile(water_mask_file):
#         with rasterio.open(water_mask_file) as dataset:
#             band1 = dataset.read(1)
#             # extract position of the Sentinel-3 acquisition
#             row, col = dataset.index(vs['lon'], vs['lat'])
#             water_mask_list = []
#             for i in range(len(col)):
#                 # if the water mask includes the position of the acqusitions
#                 if (row[i] < band1.shape[0]) and (col[i] < band1.shape[1]) and (
#                         row[i] > 0) and (col[i] > 0):
#                     water_mask_list.append(band1[row[i], col[i]])
#                 else:
#                     water_mask_list.append(np.nan)
#         vs['water_mask'] = water_mask_list

#     # Add elevation value if possible
#     if os.path.isfile(dem_file):
#         with rasterio.open(dem_file) as dataset:
#             band1 = dataset.read(1)
#             # extract position of the Sentinel-3 acquisition
#             row, col = dataset.index(vs['lon'], vs['lat'])
#             # Do the same thing as for water mask but for DEM
#             dem_list = []
#             for i in range(len(col)):
#                 # if the water mask includes the position of the acqusitions
#                 if (row[i] < band1.shape[0]) and (col[i] < band1.shape[1]) and (
#                         row[i] > 0) and (col[i] > 0):
#                     dem_list.append(band1[row[i], col[i]])
#                 else:
#                     dem_list.append(np.nan)
#         vs['dem'] = dem_list
    
#     return vs


# def read_ice_ATL03_nc(ice_file, water_mask_file, dem_file, dem_thresh=30,
#                       wm_thresh=0):

#     """
#     Parameters
#     ----------
#     s3_folder : STRING
#         Foldername and path for Level 2 netcdf files
#     source : STRING
#         "GPOD" or "SciHub" - source of Sentinel-3 netcdf file        

#     Returns
#     -------
#     height : Array
#         Heights in track referenced to geoid (EGM2008)
        
#     """
    
#     beams = ['gt1l', 'gt1r', 'gt2l', 'gt2r', 'gt3l', 'gt3r']
#     nc = netCDF4.Dataset(ice_file)
#     lat = 'lat_ph'
#     lon = 'lon_ph'
#     tai = 'delta_time'
#     atlas_offset = nc['ancillary_data']['atlas_sdp_gps_epoch'][:].filled()

#     vs = pd.DataFrame()
#     for beam in beams:
#         try:
#             conf_water = nc[beam]['heights']['signal_conf_ph'][:].filled()[:,-1]
#             conf_land = nc[beam]['heights']['signal_conf_ph'][:].filled()[:,0]
            
            
            
#             heights = nc[beam]['heights']['h_ph'][:].filled()[(conf_water > 1) | (conf_land > 1)]
#             geoid_corr = nc[beam]['geophys_corr']['geoid'][:].filled()
#             geoid_corr[geoid_corr == nc[beam]['geophys_corr']['geoid'][:].fill_value] = np.nan
#             geoid = np.interp(nc[beam]['heights']['delta_time'][:].filled(), 
#                               nc[beam]['geophys_corr']['delta_time'][:].filled(),
#                               geoid_corr)[(conf_water > 1) | (conf_land > 1)]

#             _vs = pd.DataFrame({'height'   : heights - geoid,
#                                 'lat'      : nc[beam]['heights'][lat][:].filled()[(conf_water > 1) | (conf_land > 1)],
#                                 'lon'      : nc[beam]['heights'][lon][:].filled()[(conf_water > 1) | (conf_land > 1)],
#                                 'date'     : [(datetime(1980,1,6) + timedelta(seconds = c2_time+atlas_offset)).date() 
#                                          for c2_time in nc[beam]['heights'][tai][:].filled()[(conf_water > 1) | (conf_land > 1)]],
#                                 'beam'     : np.repeat(beam, len(heights)),
#                                 'conf_land': conf_land[(conf_water > 1) | (conf_land > 1)],
#                                 'conf_water': conf_water[(conf_water > 1) | (conf_land > 1)]})
    
#             vs = pd.concat([vs, _vs])
#         except IndexError:
#             pass

#     if len(vs) > 0:
#         # Add water occurence value if possible
#         if os.path.isfile(water_mask_file):
#             with rasterio.open(water_mask_file) as dataset:
#                 band1 = dataset.read(1)
#                 # extract position of the Sentinel-3 acquisition
#                 row, col = dataset.index(vs['lon'], vs['lat'])
#                 water_mask_list = []
#                 for i in range(len(col)):
#                     # if the water mask includes the position of the acqusitions
#                     if (row[i] < band1.shape[0]) and (col[i] < band1.shape[1]) and (
#                             row[i] > 0) and (col[i] > 0):
#                         water_mask_list.append(band1[row[i], col[i]])
#                     else:
#                         water_mask_list.append(np.nan)
#             vs['water_mask'] = water_mask_list
#             vs = vs.loc[vs.water_mask > wm_thresh]
    
#         # Add elevation value if possible
#         if os.path.isfile(dem_file):
#             with rasterio.open(dem_file) as dataset:
#                 band1 = dataset.read(1)
#                 # extract position of the Sentinel-3 acquisition
#                 row, col = dataset.index(vs['lon'], vs['lat'])
#                 # Do the same thing as for water mask but for DEM
#                 dem_list = []
#                 for i in range(len(col)):
#                     # if the water mask includes the position of the acqusitions
#                     if (row[i] < band1.shape[0]) and (col[i] < band1.shape[1]) and (
#                             row[i] > 0) and (col[i] > 0):
#                         dem_list.append(band1[row[i], col[i]])
#                     else:
#                         dem_list.append(np.nan)
#             vs['dem'] = dem_list
    
#     return vs


# def set_weak_beam(infile : str) -> list[str]:
    
#     """
#     Orders track key list according to strong/weak beams, with strong beams
#     first for each pair
    
#     Arguments:
#     ----------
#         infile: path to ATL03 file
    
#     Returns:
#     ----------
#         track_keys: list of ordered track keys with strong beams first in each pair
#     """
    
#     with h5py.File(infile, 'r') as src:
#         try: 
#             sc_orient = np.asarray(list(src['orbit_info']['sc_orient']))
#         except:
#             print('sc_orient not found - estimate direction by photon count')
            
#             if check_height_data(infile, 'gt1l') and check_height_data(infile, 'gt1r'):
#                 gt1l = np.asarray(list(src['gt1l']['ht_ortho']))
#                 gt1r = np.asarray(list(src['gt1r']['ht_ortho']))
                
#                 if np.divide(gt1l.shape[0], gt1r.shape[0]) > 2:
#                     sc_orient = [0]
#                 elif np.divide(gt1l.shape[0], gt1r.shape[0]) < 0.5:
#                     sc_orient = [1]
#                 else:
#                     print('Inconclussive assesment of weak track')
#                     sc_orient = [0]   
#             else:
#                 print('Inconclussive assesment - no height data found')
#                 sc_orient = [0]
            
#     # set the order of the track keys, strong beam is the first in each pair
#     if sc_orient[0] == 1:     
#         track_keys = ['gt1r', 'gt1l', 'gt2r', 'gt2l', 'gt3r', 'gt3l']
#     elif sc_orient[0] == 0:
#         track_keys = ['gt1l', 'gt1r', 'gt2l', 'gt2r', 'gt3l', 'gt3r']
    
#     return track_keys


# def check_height_data(infile : str, track_key : str) -> bool:
    
#     """
#     Check if height data is present in ATL03 file for given track key
    
#     Arguments:
#     -----------
#         infile: path to ATL03 file
#         track_key: ICESat-2 beam key (gt1l, gt1r, gt2l, gt2r, gt3l, gt3r)
    
#     Returns:
#     ----------
#         height_data_present: boolean indicating presence of height data in file
#     """
    
#     with h5py.File(infile, 'r') as src:
#         if track_key in list(src.keys()):
#             if 'ht_ortho' in list(src[track_key].keys()):
#                 height_data_present = True 
#                 print('height data found')
#             else:
#                 height_data_present = False
#                 print('height data not found')        
#         else:
#             height_data_present = False
#             print('track data not found')
        
#     return height_data_present 