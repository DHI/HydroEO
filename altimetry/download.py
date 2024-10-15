


# def icepyx_query(aoi_bbox, startdate, enddate, earthdata_credentials, download_directory, product='ATL13'):
#     """
#     Download through NASA NSIDC DAAC
#     Reference: Scheick, J. et al., (2019). icepyx: Python tools for obtaining and working with ICESat-2 data. https://github.com/icesat2py/icepyx.

#     Parameters
#     ----------
#     aoi_bbox : DataFrame
#         AOI bounding box (from geopandas.geometry.bounds)
#         DataFrame with columns minx, miny, maxx, maxy values containing the bounds for the AOI geometry.
#     startdate : String
#         YYYY-MM-DD - startdate for download.
#     enddate : String
#         YYYY-MM-DD - enddate for download.
#     earthdata_credentials : Filename path
#         Textfile with username and password for EarthData on two separate lines
#     download_directory : Path
#         Directory to save files.
#     product : string, optional
#         Default is ATL13
#         ICESat-2 product to be downloaded.

#     Returns
#     -------
#     None.

#     """
#     un, pw = read_credentials(earthdata_credentials)
    
#     region_a = ipx.Query(product, aoi_bbox ,[startdate.strftime('%Y-%m-%d'), enddate.strftime('%Y-%m-%d')], \
#                                start_time='00:00:00', end_time='23:59:59')
#     region_a.earthdata_login(un, pw)
#     region_a.order_granules()

#     return region_a.granules.orderIDs