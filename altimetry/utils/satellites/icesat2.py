from dataclasses import dataclass
import os
import warnings
import typing
import datetime
import icepyx as ipx

import h5py
import numpy as np
import pandas as pd
import geopandas as gpd

from tqdm import tqdm


def query(
    aoi: list,
    startdate: datetime.date,
    enddate: datetime.date,
    earthdata_credentials: tuple,
    download_directory: str,
    product: str = "ATL13",
) -> object:
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
    query = ipx.Query(
        product=product,
        spatial_extent=aoi,
        date_range=[startdate.strftime("%Y-%m-%d"), enddate.strftime("%Y-%m-%d")],
        start_time="00:00:00",
        end_time="23:59:59",
    )

    # provide credentials
    if earthdata_credentials is not None:
        un, pw = earthdata_credentials
        query.earthdata_login(un, pw)

    # TODO: add in some kind of log to not redownload data?
    # This is difficult to do because as of now there is no way with icepyx to select sprcific granuals returned from a query,
    # so even if we were to log the downloaded granules, there is no way to ensure that we dont download them again later if
    # the same granuale appears in a new set of search parameters. Thus we may be downloading a granule with no actual crossing
    # over a reservoir many times over. Consider pull request to be able to add this functionality later

    # order granules
    query.order_granules()
    query.download_granules(download_directory)


@dataclass
class ATL13:
    """
    ICESat-2 utility class for reading ATL13 data from .h5 files.

    Arguments:
    ----------
        infile: ATL03 file path (.h5)
        track_key: ICESat-2 ground track key
    """

    infile: str
    track_key: str

    def __post_init__(self):
        self.file_name = os.path.basename(self.infile)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            with h5py.File(self.infile, "r") as src:
                # Segment track heights
                self.setattr_from_file(src, "height_seg", "ht_ortho")
                self.setattr_from_file(src, "lat_seg", "segment_lat")
                self.setattr_from_file(src, "lon_seg", "segment_lon")

                # Water body information
                self.setattr_from_file(src, "wb_id", "inland_water_body_id")
                self.setattr_from_file(src, "wb_size", "inland_water_body_size")
                self.setattr_from_file(src, "wb_type", "inland_water_body_type")

                # Segment information - error and quality
                self.setattr_from_file(src, "wb_slope", "segment_slope_trk_bdy")
                self.setattr_from_file(src, "height_err", "err_ht_water_surf")
                self.setattr_from_file(src, "quality_seg", "segment_quality")

                # Tide/DEM corrections
                self.setattr_from_file(src, "geoid_track", "segment_geoid")
                self.setattr_from_file(
                    src, "geoid_corr_track", "segment_geoid_free2mean"
                )
                self.setattr_from_file(src, "dem_track", "segment_dem_ht")

                # Saturation fraction - The fraction of pulses within the short segment determined to be nearly saturated based on ATL03 geosegment rate input.
                self.setattr_from_file(src, "sat_frac_track", "segment_near_sat_fract")

                # time and date attributes
                if "delta_time" in src[self.track_key].keys():
                    delta_time_seg = np.asarray(
                        src[self.track_key]["delta_time"], float
                    )
                    atlas_offset = np.asarray(
                        src["ancillary_data"]["atlas_sdp_gps_epoch"]
                    )[0]
                    self.date = [
                        (
                            datetime.datetime(1980, 1, 6)
                            + datetime.timedelta(seconds=c2_time + atlas_offset)
                        )
                        for c2_time in delta_time_seg
                    ]
                else:
                    self.date = None

    def setattr_from_file(self, src: h5py.File, name: str, attribute: str):
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

        track_df = pd.DataFrame(
            data={
                "height": self.height_seg,
                "lat": self.lat_seg,
                "lon": self.lon_seg,
                "date": self.date,
                "wb_type": self.wb_type,
                "wb_size": self.wb_size,
                "wb_id": self.wb_id,
                "dem": self.dem_track,
                "sat_frac_track": self.sat_frac_track,
                "beam": self.track_key,
                "file_name": self.file_name,
            },
            index=np.arange(len(self.height_seg)),
        )

        return track_df


def extract_observations(src_dir, dst_path, features):
    # read icesat data for each availble option in directory
    gdf_list = list()
    files = list(os.listdir(src_dir))
    for file in files:
        for key in ["gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"]:
            infile = os.path.join(src_dir, file)
            data = ATL13(infile, key)

            # check if the height data is actually present
            if data.check_height_data():
                data_df = data.read()
                data_gdf = gpd.GeoDataFrame(
                    data_df, geometry=gpd.points_from_xy(data_df.lon, data_df.lat)
                )

                # filter observations to ensure they fall within geometry
                data_gdf = data_gdf.loc[
                    data_gdf.within(features.unary_union)
                ].reset_index(drop=True)

                # if we have data for the reservoir add it to the reservoir specific dataframe
                if len(data_gdf) > 0:
                    gdf_list.append(data_gdf)

    # once all tracks are processed combine them and save in the destination dir
    if len(gdf_list) > 0:
        observations = pd.concat(gdf_list).reset_index(drop=True)
        observations["platform"] = "icesat2"
        observations["product"] = "ATL13"
        observations = observations.set_crs(features.crs)
        observations.to_file(dst_path)


def get_latest_obs_date(data_dir):
    shp_path = os.path.join(data_dir, "icesat2.shp")
    if os.path.exists(shp_path):
        gdf = gpd.read_file(shp_path)
        last_obs_date = max(gdf.date.values).astype(datetime.date)
        return last_obs_date
