from dataclasses import dataclass
import warnings
import os
from tqdm import tqdm

import netCDF4

import datetime

import numpy as np
import pandas as pd
import geopandas as gpd
import shapely

from HydroEO.downloaders import creodias
from HydroEO.utils import geometry
from HydroEO.utils.general import center_longitude


def query(
    aoi: list,
    startdate: datetime.date,
    enddate: datetime.date,
    product: str = "S3",
) -> object:
    """
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
            Sentinel Satellite product we are interested in. As of now, only Sentinel3 (SR_2_LAN___) and Sentinel6 (P4_2__LR_____) are supported
    """
    # ensure product label is in all caps to standardize
    product = product.upper()

    # We want to expect a similar formatting for aoi queries as with icesat2 so we will need to convert the list of coords into a polygon wihtin the WKT format
    # first check if the coordinates are given as a list of tuples or as a list of adjacent points
    aoi = geometry.format_coord_list(aoi)

    # build the query parameters based on the requeste product
    if product == "S3":
        params = {
            "collection": "Sentinel3",
            "start_date": startdate,
            "end_date": enddate,
            "geometry": shapely.Polygon(aoi).wkt,
            "productType": "SR_2_LAN_HY",
            "instrumentShortName": "SRAL",
        }

    elif product == "S6":
        params = {
            "collection": "Sentinel6",
            "start_date": startdate,
            "end_date": enddate,
            "geometry": shapely.Polygon(aoi).wkt,
            "productType": "P4_2__LR_____",
            "instrumentShortName": "P4",
        }

    else:
        raise ValueError(
            f'"{product}" is unrecognized as a valid sentinel product for query'
        )

    results = creodias.query(**params)
    ids = [result["id"] for result in results.values()]

    return ids


def download(
    ids: list,
    download_directory: str,
    creodias_credentials: tuple,
    token: str = None,
    session_start_time: str = None,
):
    # Check if we have a progress log file in this directory, if not make it
    log_path = os.path.join(download_directory, "downloaded.log")
    if not os.path.exists(log_path):
        with open(log_path, "w") as log:
            pass

    # open the log file with reading and writing access
    with open(log_path, "r") as log:
        # first read all of the downloaded ids
        downloaded_ids = [line.rstrip() for line in log]
        already_downloaded = [id for id in ids if id in downloaded_ids]
        ids_to_download = [id for id in ids if id not in downloaded_ids]

    # download single product by product ID, if not already downloaded
    print(f"{len(ids)} products returned by query.")
    print(f"{len(already_downloaded)} products already downloaded.")
    print(f"Downloading {len(ids_to_download)} files.")

    # Download ids that are not listed as downloaded already in the logfile
    return creodias.download_list(
        ids_to_download,
        username=creodias_credentials[0],
        password=creodias_credentials[1],
        token=token,
        session_start_time=session_start_time,
        outdir=download_directory,
        log_file=log_path,
        threads=1,
    )

    """ # Code using the single download function which generates a token for each download and hits the session limit
    # prefer the above code

    for i in tqdm(range(0, len(ids))):
        if ids[i] not in downloaded_ids:
            outfile = os.path.join(
                download_directory, file_prefix + str(ids[i]) + ".zip"
            )
            creodias.download(
                ids[i],
                outfile=outfile,
                username=creodias_credentials[0],
                password=creodias_credentials[1],
            )
            with open(log_path, "a") as log:
                log.write(ids[i] + "\n")  # add the id to the downloaded log
    """


def subset(
    aoi: list,
    download_dir: str,
    dest_dir: str,
    file_id: str = "enhanced_measurement.nc",
    product: str = "S3",
    show_progress=False,
):
    """
    Subset sentinel-3 netcdf files

    Parameters
    ----------
    download_dir : String
        Directory with all .SENX folders.
    dest_dir : String
        Directory to save cropped .nc files.
    extent : List
        Corner coordinates of area of interest.
    file_id : string, optional
        Original netcdf to crop. The default is r'enhanced_measurement.nc'.
    extension : string, optional
        Folder extension. The default is '.SEN3' for Sentinel-3
    """

    # ensure correct formatting of coordinate list and extract corner coordinates
    # Upper left corner = lon min, lat max / Lower right corner = lon max, lat min
    ulx, lry, lrx, uly = shapely.Polygon(geometry.format_coord_list(aoi)).bounds
    # print(ulx, lry, lrx, uly)
    if ulx > lrx or lry > uly:
        raise ValueError("Study area extent conflicts - please check coordinates.")

    # define product specific variables
    product = product.upper()

    EXTENSION = {"S3": ".SEN3", "S6": ".SEN6"}

    VARIABLES = {
        "S3": [
            "lat_20_ku",
            "lon_20_ku",
            "elevation_ocog_20_ku",
            "waveform_20_ku",
            "sig0_ocog_20_ku",
            "time_20_ku",
            "alt_20_ku",
            "tracker_range_20_ku",
            "range_ocog_20_ku",
            "lat_01",
            "geoid_01",
        ],
        "S6": [
            "latitude",
            "longitude",
            "altitude",
            "range_ocog",
            "model_wet_tropo_cor_measurement_altitude",
            "model_dry_tropo_cor_measurement_altitude",
            "sig0_ocog",
            "time_tai",
            "geoid",
        ],
    }

    ATTRIBUTES = {
        "S3": ["first_meas_lat", "last_meas_lat", "pass_number"],
        "S6": [
            "first_measurement_latitude",
            "last_measurement_latitude",
            "pass_number",
        ],
    }

    # loop through the downloads folder and subset any products with the right extention
    pbar = tqdm(
        total=int(len(os.listdir(download_dir)) / 2),
        desc=f"Subsetting data in {os.path.basename(download_dir)}",
        unit="file",
        disable=not show_progress,
    )
    for folder in os.listdir(download_dir):
        # print(EXTENSION[product], folder)
        if folder.endswith(EXTENSION[product]):
            pbar.update(1)
            # within the sentinel file, select the correct type of mearurements with the file_id key
            if product == "S3":
                file = os.path.join(download_dir, folder, file_id)

            elif product == "S6":
                for potential in os.listdir(os.path.join(download_dir, folder)):
                    if "STD" in potential.split("_"):
                        if potential.endswith(".nc"):
                            file = os.path.join(download_dir, folder, potential)

            if os.path.isfile(file):
                # we will save a new file with the same name as the sentinel folder but only keeping the nc file we choose
                nc_cropped = os.path.join(
                    dest_dir,
                    "sub_"
                    + os.path.split(folder)[-1].split(EXTENSION[product])[0]
                    + ".nc",
                )

                # now we read the file
                # print(file)
                nc = netCDF4.Dataset(file)

                # SciHub contains 1Hz and 20Hz variables
                # Get index of values within AOI at both frequencies
                # The way that variabels are accesed/called within sentinel 3 and 6 are different, so we extract the indices as needed
                if product == "S3":
                    freq_20 = "_20_ku"
                    freq_01 = "_01"

                    lat20 = nc["lat" + freq_20][:]
                    lon20 = nc["lon" + freq_20][:]

                    lat01 = nc["lat" + freq_01][:]
                    lon01 = nc["lon" + freq_01][:]

                elif product == "S6":
                    freq_20 = "data_20"
                    freq_01 = "data_01"

                    lat = "latitude"
                    lon = "longitude"

                    lat20 = nc[freq_20]["ku"][lat][:]
                    lon20 = nc[freq_20]["ku"][lon][:]

                    lat01 = nc[freq_01][lat][:]
                    lon01 = nc[freq_01][lon][:]

                # Adjust the longitude
                lon20 = center_longitude(lon20)
                # print(min(lon20), max(lon20))

                min_index20 = len(lat20) + 1
                max_index20 = 0
                min_index01 = len(lat01) + 1
                max_index01 = 0

                # Get the indices within the bounds for 20Hz
                selected = np.where(
                    (lon20 <= lrx) & (lon20 >= ulx) & (lat20 >= lry) & (lat20 <= uly)
                )[0]
                # print(len(selected))

                if len(selected) > 0:
                    if selected[0] < min_index20:
                        min_index20 = selected[0]

                    if selected[-1] > max_index20:
                        max_index20 = selected[-1]

                # Get the indices within the bounds for 1Hz
                selected01 = np.where(
                    (lon01 <= lrx) & (lon01 >= ulx) & (lat01 >= lry) & (lat01 <= uly)
                )[0]

                if len(selected01) > 0:
                    if selected01[0] < min_index01:
                        min_index01 = selected01[0]

                    if selected01[-1] > max_index01:
                        max_index01 = selected01[-1]

                nc.close()

                # if we have data within the bounds, start copying all extra data to new file
                # print(max_index20)
                # print(lon20)
                # print(lat20)
                if max_index20 > 0:
                    # print("Copying data to new file")
                    with (
                        netCDF4.Dataset(file) as src,
                        netCDF4.Dataset(nc_cropped, "w") as dst,
                    ):
                        if product == "S3":
                            __crop_s3(
                                src,
                                dst,
                                (min_index20, max_index20),
                                (min_index01, max_index01),
                            )

                        elif product == "S6":
                            __crop_s6(
                                src,
                                dst,
                                (min_index20, max_index20),
                                (min_index01, max_index01),
                            )

                    # Check that everything crucial is there:
                    with (
                        netCDF4.Dataset(file) as src,
                        netCDF4.Dataset(nc_cropped, "r") as dst,
                    ):
                        if product == "S3":
                            for var in VARIABLES[product]:
                                if var not in src.variables.keys():
                                    print(f"{var} missing from copied file.")

                            for att in ATTRIBUTES[product]:
                                if att not in src.ncattrs():
                                    print(f"{att} missing from copied file.")

                        if product == "S6":
                            for var in VARIABLES[product]:
                                if var not in src["data_20"]["ku"].variables.keys():
                                    print(f"{var} missing from copied file.")

                            for att in ATTRIBUTES[product]:
                                if att not in src.ncattrs():
                                    print(f"{att} missing from copied file.")

                else:
                    # print("No data within bounds, no subsetted file created")
                    pass


def __crop_s3(src, dst, index_bounds_20, index_bounds_01):
    min_index20, max_index20 = index_bounds_20
    min_index01, max_index01 = index_bounds_01

    # copy global attributes all at once via dictionary
    dst.setncatts(src.__dict__)

    # copy dimensions
    for name, dimension in src.dimensions.items():
        if "20_" in name:
            dst.createDimension(
                name,
                (
                    len(np.arange(min_index20, max_index20 + 1, 1))
                    if not dimension.isunlimited()
                    else None
                ),
            )
        elif "01" in name:
            dst.createDimension(
                name,
                (
                    len(np.arange(min_index01, max_index01 + 1, 1))
                    if not dimension.isunlimited()
                    else None
                ),
            )
        else:
            dst.createDimension(
                name, (len(dimension) if not dimension.isunlimited() else None)
            )

    # copy all file data except for the excluded
    for name, variable in src.variables.items():
        if "20_" in name:
            dst.createVariable(name, variable.datatype, variable.dimensions)
            dst[name].setncatts({k: variable.getncattr(k) for k in variable.ncattrs()})
            # print(min_index20, max_index20)
            # print(src[name][min_index20 : max_index20 + 1])
            dst[name][:] = src[name][min_index20 : max_index20 + 1]

        elif "01" in name:
            dst.createVariable(name, variable.datatype, variable.dimensions)
            dst[name].setncatts({k: variable.getncattr(k) for k in variable.ncattrs()})
            dst[name][:] = src[name][min_index01 : max_index01 + 1]


def __crop_s6(src, dst, index_bounds_20, index_bounds_01):
    min_index20, max_index20 = index_bounds_20
    min_index01, max_index01 = index_bounds_01

    # copy global attributes all at once via dictionary
    dst.setncatts(src.__dict__)

    # go through all groups within file
    for gr_name, group in src.groups.items():
        # create a duplicate group in the destination
        dst.createGroup(gr_name)
        dst[gr_name].setncatts(src[gr_name].__dict__)

        if "_01" in gr_name:
            for name, dimension in src[gr_name].dimensions.items():
                dst[gr_name].createDimension(
                    name,
                    (
                        len(np.arange(min_index01, max_index01 + 1, 1))
                        if not dimension.isunlimited()
                        else None
                    ),
                )
            for name, variable in src[gr_name].variables.items():
                dst[gr_name].createVariable(
                    name, variable.datatype, variable.dimensions
                )
                dst[gr_name][name].setncatts(
                    {k: variable.getncattr(k) for k in variable.ncattrs()}
                )
                dst[gr_name][name][:] = src[gr_name][name][
                    min_index01 : max_index01 + 1
                ]

            for gr_name_bw, group in src[gr_name].groups.items():
                dst[gr_name].createGroup(gr_name_bw)
                dst[gr_name][gr_name_bw].setncatts(src[gr_name][gr_name_bw].__dict__)

                for name, dimension in src[gr_name][gr_name_bw].dimensions.items():
                    dst[gr_name][gr_name_bw].createDimension(
                        name,
                        (
                            len(np.arange(min_index01, max_index01 + 1, 1))
                            if not dimension.isunlimited()
                            else None
                        ),
                    )

                # Create all variables for each subgroup
                for name, variable in src[gr_name][gr_name_bw].variables.items():
                    dst[gr_name][gr_name_bw].createVariable(
                        name, variable.datatype, variable.dimensions
                    )
                    dst[gr_name][gr_name_bw][name].setncatts(
                        {k: variable.getncattr(k) for k in variable.ncattrs()}
                    )
                    dst[gr_name][gr_name_bw][name][:] = src[gr_name][gr_name_bw][name][
                        min_index01 : max_index01 + 1
                    ]

        elif "_20" in gr_name:
            for gr_name_bw, group in src[gr_name].groups.items():
                dst[gr_name].createGroup(gr_name_bw)
                dst[gr_name][gr_name_bw].setncatts(src[gr_name][gr_name_bw].__dict__)
                for name, dimension in src[gr_name][gr_name_bw].dimensions.items():
                    dst[gr_name][gr_name_bw].createDimension(
                        name,
                        (
                            len(np.arange(min_index20, max_index20 + 1, 1))
                            if not dimension.isunlimited()
                            else None
                        ),
                    )

                # Create all variables for each subgroup
                for name, variable in src[gr_name][gr_name_bw].variables.items():
                    dst[gr_name][gr_name_bw].createVariable(
                        name, variable.datatype, variable.dimensions
                    )
                    dst[gr_name][gr_name_bw][name].setncatts(
                        {k: variable.getncattr(k) for k in variable.ncattrs()}
                    )
                    dst[gr_name][gr_name_bw][name][:] = src[gr_name][gr_name_bw][name][
                        min_index20 : max_index20 + 1
                    ]

        else:
            for name, dimension in src[gr_name].dimensions.items():
                dst[gr_name].createDimension(
                    name, (len(dimension) if not dimension.isunlimited() else None)
                )


@dataclass
class Sentinel:
    infile: str

    def check_height_data(self) -> bool:
        return hasattr(self, "height")


@dataclass
class Sentinel3(Sentinel):
    infile: str

    def __post_init__(self):
        self.file_name = os.path.basename(self.infile)

        # TODO: some kind of check to make sure we can open the data

        # proceed with opening the data
        with netCDF4.Dataset(self.infile, "r") as src:
            # TODO: check if 20Hz elevation is available
            self.__setattr_from_file(src, "elev", "elevation_ocog_20_ku")
            self.__setattr_from_file(src, "lat", "lat_20_ku")
            self.__setattr_from_file(src, "lon", "lon_20_ku")
            self.__setattr_from_file(src, "wf", "waveform_20_ku")
            self.__setattr_from_file(src, "sig0", "sig0_ocog_20_ku")
            self.__setattr_from_file(src, "alt", "alt_20_ku")
            self.__setattr_from_file(src, "tracker_range", "tracker_range_20_ku")
            self.__setattr_from_file(src, "range_OCOG", "range_ocog_20_ku")
            self.epoch_OCOG = (
                -self.tracker_range + self.range_OCOG - 0.55
            ) * 128 / 60 + 43  # 128 bins for 60m

            # set date attribute
            self.__setattr_from_file(src, "tai", "time_20_ku")
            self.date = [
                (datetime.datetime(2000, 1, 1) + datetime.timedelta(seconds=c2_time))
                for c2_time in self.tai
            ]

            # Get geoid elevation at 20Hz resolution from 1Hz dataset
            self.__setattr_from_file(src, "lat_01", "lat_01")
            self.__setattr_from_file(src, "geoid_01", "geoid_01")

            # interpolate the geoid
            self.__interpolate_geoid()

            # Get retracked WSE if any valid observations over land
            self.height = self.elev - self.geoid

            # Set some auxiliary data
            self.first_meas_lat = src.getncattr("first_meas_lat")
            self.last_meas_lat = src.getncattr("last_meas_lat")
            self.pass_number = src.getncattr("pass_number")

    def __setattr_from_file(self, src: netCDF4.Dataset, name: str, attribute: str):
        if attribute in src.variables.keys():
            self.__setattr__(name, src[attribute][:].filled())
        else:
            self.__setattr__(name, None)

    def __interpolate_geoid(self):
        if len(self.geoid_01) > 1:
            # perform the interpolation
            geoid_interp = np.interp(self.lat, self.lat_01, self.geoid_01)

            # np.interp sorts to ascending values - if descending, this has to be corrected.
            if self.lat_01[0] > self.lat_01[1]:
                geoid_interp = np.interp(
                    self.lat,
                    np.flip(self.lat_01, axis=0),
                    np.flip(self.geoid_01, axis=0),
                )

            elif len(self.lat_01) == len(self.elev):
                geoid_interp = self.geoid_01

        else:
            geoid_interp = self.elev * 0

        # now reset the geoid to be the new interpolated geoid
        self.geoid = geoid_interp

    def read(self) -> pd.DataFrame:
        # Create a Data Frame with coordinates, OCOG height, date
        vs = pd.DataFrame(
            data={
                "height": self.height,
                "sigma0": self.sig0,
                "lat": self.lat,
                "lon": self.lon,
                "date": self.date,
                "alt": self.alt,
                "geoid": self.geoid,
                "tracker_range": self.tracker_range,
                "range_OCOG": self.range_OCOG,
                "epoch_OCOG": self.epoch_OCOG,
            },
            index=np.arange(len(self.height)),
        )

        # Apply the same geophysical corrections as for OCOG
        vs["geo_cor"] = (vs["alt"] - vs["range_OCOG"] - vs["geoid"]) - vs["height"]

        # Store waveform
        wf_vector = pd.DataFrame(self.wf)
        vs["wf"] = wf_vector.apply(lambda r: tuple(r), axis=1).apply(np.array)

        # Determine if the path is ascending or ascending as well as its path, orbit and file name
        vs["sat_path"] = np.repeat(
            "descending"
            if self.first_meas_lat - self.last_meas_lat > 0
            else "ascending",
            len(self.height),
        )
        vs["pass"] = np.repeat(self.pass_number, len(self.height))

        # Sentinel-3 naming convention:
        # OBS we add "sub" in front when subsetting
        # MMM_OL_L_TTTTTT_yyyymmddThhmmss_YYYYMMDDTHHMMSS_YYYYMMDDTHHMMSS_[instance ID]_GGG_/
        # [class ID].SEN3
        # Instance ID: DDDD_CCC_LLL____: Duration, cycle number, relative orbit number,
        # 4 underscores
        # GGG: Center which generated the file
        # Class ID: P_XX_NNN - P = platform, XX - timeliness, NNN - version
        vs["file_name"] = np.repeat(self.file_name, len(self.height))

        file_name_split = self.file_name.split("_")
        vs["orbit"] = np.repeat(
            list(filter(None, file_name_split))[10], len(self.height)
        )
        vs["relative_orbit"] = np.repeat(
            list(filter(None, file_name_split))[9], len(self.height)
        )
        vs["platform"] = np.repeat(
            list(filter(None, file_name_split))[1], len(self.height)
        )

        # There are sometimes issues with the time in variable time_20_ku so we also store
        # the acquisition date
        str_date = vs["file_name"][0][20:35]
        date = datetime.datetime(
            int(str_date[0:4]),
            int(str_date[4:6]),
            int(str_date[6:8]),
            int(str_date[9:11]),
            int(str_date[11:13]),
        )
        vs["date_acquisition"] = np.repeat(date, len(self.height))

        return vs


@dataclass
class Sentinel6(Sentinel):
    infile: str

    def __post_init__(self):
        self.file_name = os.path.basename(self.infile)

        # TODO: some kind of check to make sure we can open the data

        # proceed with opening the data
        with netCDF4.Dataset(self.infile, "r") as src:
            # subset into the file to extract 20 and 1hz data
            self.data_20 = src["data_20"]["ku"]
            self.data_01 = src["data_01"]

            # pull data from the 20Hz TODO: could add a set from file attribute as with sentinel 3 class, this would check if the attribute exists before loading
            self.lat = self.data_20["latitude"][:].filled()
            self.lon = self.data_20["longitude"][:].filled()
            self.sig0 = self.data_20.variables["sig0_ocog"][:].filled()
            self.peakiness = self.data_20["peakiness"][:].filled()
            self.tai = self.data_20.variables["time_tai"][:].filled()
            self.geoid = self.data_20["geoid"][:].filled()
            self.altitude = self.data_20.variables["altitude"][:].filled()
            self.range_ocog = self.data_20.variables["range_ocog"][:].filled()
            self.wet_tropo = self.data_20.variables[
                "model_wet_tropo_cor_measurement_altitude"
            ][:].filled()
            self.dry_tropo = self.data_20.variables[
                "model_dry_tropo_cor_measurement_altitude"
            ][:].filled()

            # set the length of data points
            self.length = len(self.altitude)

            # calculate the height
            try:
                self.__interpolate_tides()
                self.height = (
                    self.altitude
                    - self.range_ocog
                    + self.wet_tropo
                    + self.solid_earth_tide
                    + self.pole_tide
                    + self.dry_tropo
                    - self.geoid
                )
            except IndexError:
                warnings.warn(f"No height data found within file: {self.file_name}")

            # Set some auxiliary data
            self.first_meas_lat = src.getncattr("first_measurement_latitude")
            self.last_meas_lat = src.getncattr("last_measurement_latitude")
            self.pass_number = src.getncattr("pass_number")
            self.absolute_rev_number = src.getncattr("absolute_rev_number")

    def __interpolate_tides(self):
        # temperary variables from 1Hz frequency to use for interpolation
        lat_01 = self.data_01["latitude"][:].filled()
        solid = self.data_01["solid_earth_tide"][
            :
        ].filled()  # check if we can open solid variable?
        pole = self.data_01["pole_tide"][:].filled()

        ### interpolate tides
        if len(solid) > 1:
            self.solid_earth_tide = np.interp(self.lat, lat_01, solid)
            self.pole_tide = np.interp(self.lat, lat_01, pole)

        # np.interp sorts to ascending values - if descending, this has to be corrected.
        if lat_01[0] > lat_01[1]:
            self.solid_earth_tide = np.interp(self.lat, np.flip(lat_01), np.flip(solid))
            self.pole_tide = np.interp(self.lat, np.flip(lat_01), np.flip(pole))

    def read(self):
        # create initial dataframe
        vs = pd.DataFrame(
            {
                "height": self.height,
                "sigma0": self.sig0,
                "lat": self.lat,
                "lon": self.lon,
                "date": [
                    (
                        datetime.datetime(2000, 1, 1)
                        + datetime.timedelta(seconds=c2_time)
                    )  # .date()
                    for c2_time in self.tai
                ],
                "wf_peakiness": self.peakiness,
                "geoid": self.geoid,
            }
        )

        # calculate additional data columns
        vs["sat_path"] = np.repeat(
            "descending"
            if self.first_meas_lat - self.last_meas_lat > 0
            else "ascending",
            self.length,
        )
        vs["pass"] = np.repeat(self.pass_number, self.length)
        vs["orbit"] = np.repeat(self.absolute_rev_number, self.length)
        vs["file_name"] = np.repeat(self.file_name, self.length)

        return vs


def extract_observations(src_dir, dst_path, features):
    # read data for each availble option in directory
    gdf_list = list()
    files = list(os.listdir(src_dir))

    for file in files:
        try:
            file_split = file.split("_")
            if file_split[0] == "sub":  # assume that we have subsetted the data already
                # use differnt read funciton base don which sentinel we use
                platform = file_split[1]
                if "S3" in platform:
                    infile = os.path.join(src_dir, file)
                    data = Sentinel3(infile)
                    # hard code product name, will need to change if incorperating other products in future
                    product = "SR_2_LAN_HY"

                elif "S6" in platform:
                    infile = os.path.join(src_dir, file)
                    data = Sentinel6(infile)
                    # hard code product name, will need to change if incorperating other products in future
                    product = "P4_2__LR_____"

                else:
                    raise Exception(
                        "unsure which satellite product file is associated with"
                    )

                # read the data, if we have height data
                if data.check_height_data():
                    data_df = data.read()
                    data_gdf = gpd.GeoDataFrame(
                        data_df, geometry=gpd.points_from_xy(data_df.lon, data_df.lat)
                    )

                    # filter observations to ensure they fall within geometry
                    data_gdf = data_gdf.loc[
                        data_gdf.within(features.unary_union)
                    ].reset_index(drop=True)
                    if len(data_gdf) > 0:
                        # filter by sigma0
                        data_gdf = data_gdf.loc[data_gdf.sigma0 < 1e5].reset_index(
                            drop=True
                        )

                        # add platform by satellite type
                        data_gdf["platform"] = platform
                        data_gdf["product"] = product

                        # if we have data for the reservoir add it to the reservoir specific dataframe
                        gdf_list.append(data_gdf)
        except Exception:
            print("Unable to open sentinel file")

    # once all tracks are processed combine them and save in the destination dir
    if len(gdf_list) > 0:
        observations = pd.concat(gdf_list).reset_index(drop=True)
        observations = observations.set_crs(features.crs)
        observations.to_file(dst_path)


def get_latest_obs_date(data_dir, product):
    if product.upper() == "S3":
        shp_path = os.path.join(data_dir, "sentinel3.shp")
    elif product.upper() == "S6":
        shp_path = os.path.join(data_dir, "sentinel6.shp")

    if os.path.exists(shp_path):
        gdf = gpd.read_file(shp_path)
        last_obs_date = max(gdf.date.values).astype(datetime.date)
        return last_obs_date
