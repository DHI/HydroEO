
import os
import datetime


import shapely

from creodias_finder.query import query as query_creodias
from creodias_finder.download import download as download_creodias



def query(aoi: list, startdate: datetime.date, enddate: datetime.date, creodias_credentials: tuple, download_directory: str, product: str  ='S3') -> object:
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

    # We want to expect a similar formatting for aoi queries as with icesat2 so we will need to convert the list of coords into a polygon wihtin the WKT format
    # first check if the coordinates are given as a list of tuples or as a list of adjacent points
    if not hasattr(aoi[0], '__iter__'):

        print("Formatting coordinates")

        # make sure we have an even number of elements in the list
        if len(aoi)%2 > 0:
            raise ValueError("The inputed list of elements is not even. Please retry the query with an even list of coordinates or a list of coordinate pairs")

        # unpack the 1d list into coordinate pairs
        new_list = list()
        while (len(aoi) > 0):
            new_list.append((aoi.pop[0], aoi.pop(0)))
        aoi = new_list
    
    # now we assume that we have a list of tuples that we can make a shapely polygon out of
    aoi = shapely.Polygon(aoi)

    # build the query parameters based on the requeste product
    if product.upper() == 'S3':

        file_prefix = 'S3_'

        params = {
            'collection' : 'Sentinel3',
            'start_date' : startdate,
            'end_date' : enddate,
            'geometry' : aoi.wkt,
            'productType' : 'SR_2_LAN___',
            'instrumentShortName' : 'SRAL'
        }
        
    elif product.upper() == 'S6':

        file_prefix = 'S6_'

        params = {
            'collection' : 'Sentinel6',
            'start_date' : startdate,
            'end_date' : enddate,
            'geometry' : aoi.wkt,
            'productType' : 'P4_2__LR_____',
            'instrumentShortName' : 'P4'
        }

    else:
        raise ValueError(f'"{product}" is unrecognized as a valid sentinel product for query')
    
    results = query_creodias(**params)
    ids = [result['id'] for result in results.values()]

    # download single product by product ID
    for i in range(0, len(ids)):
        outfile = os.path.join(download_directory, file_prefix + str(i) + '.zip')
        download_creodias(ids[i], outfile=outfile, username=creodias_credentials[0], password=creodias_credentials[1])
    
    return