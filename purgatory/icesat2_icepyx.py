import datetime

import icepyx as ipx


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

    # provide credentials TODO: assess if we can remove this authentification step now that the environment variables are set
    if earthdata_credentials is not None:
        un, pw = earthdata_credentials
        query.earthdata_login(un, pw)

    # TODO: add in some kind of log to not redownload data?
    # This is difficult to do because as of now there is no way with icepyx to select sprcific granuals returned from a query,
    # so even if we were to log the downloaded granules, there is no way to ensure that we dont download them again later if
    # the same granuale appears in a new set of search parameters. Thus we may be downloading a granule with no actual crossing
    # over a reservoir many times over. Consider pull request to be able to add this functionality later

    if hasattr(query.granules, "avail"):
        # order granules
        query.order_granules()
        query.download_granules(download_directory)

    else:
        print("No available granules found")
