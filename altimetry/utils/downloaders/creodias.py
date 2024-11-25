import shutil
from pathlib import Path
import requests
from tqdm import tqdm
import datetime
from six.moves.urllib.parse import urlencode
from six import string_types
import dateutil.parser
from shapely.geometry import shape
import time
import os

##### Global variables
API_URL = (
    "https://catalogue.dataspace.copernicus.eu/resto/api/collections/{collection}"
    "/search.json?maxRecords=1000"
)
DOWNLOAD_URL = "https://zipper.creodias.eu/download"
TOKEN_URL = "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token"


##### Functions associated with queries
def query(
    collection,
    start_date=None,
    end_date=None,
    geometry=None,
    status="ONLINE",
    **kwargs,
):
    """Query the EOData Finder API

    Parameters
    ----------
    collection: str, optional
        the data collection, corresponding to various satellites
    start_date: str or datetime
        the start date of the observations, either in iso formatted string or datetime object
    end_date: str or datetime
        the end date of the observations, either in iso formatted string or datetime object
        if no time is specified, time 23:59:59 is added.
    geometry: WKT polygon or object impementing __geo_interface__
        area of interest as well-known text string
    status : str
        allowed online/offline/all status (ONLINE || OFFLINE || ALL)
    **kwargs
        Additional arguments can be used to specify other query parameters,
        e.g. productType=L1GT
        See https://documentation.dataspace.copernicus.eu/APIs/OpenSearch.html for details

    Returns
    -------
    dict[string, dict]
        Products returned by the query as a dictionary with the product ID as the key and
        the product's attributes (a dictionary) as the value.
    """
    query_url = _build_query(
        API_URL.format(collection=collection),
        start_date,
        end_date,
        geometry,
        status,
        **kwargs,
    )

    query_response = {}
    while query_url:
        response = requests.get(query_url)
        response.raise_for_status()
        data = response.json()
        for feature in data["features"]:
            query_response[feature["id"]] = feature
        query_url = _get_next_page(data["properties"]["links"])
    return query_response


def _build_query(
    base_url, start_date=None, end_date=None, geometry=None, status=None, **kwargs
):
    query_params = {}

    if start_date is not None:
        start_date = _parse_date(start_date)
        query_params["startDate"] = start_date.isoformat()
    if end_date is not None:
        end_date = _parse_date(end_date)
        end_date = _add_time(end_date)
        query_params["completionDate"] = end_date.isoformat()

    if geometry is not None:
        query_params["geometry"] = _parse_geometry(geometry)

    if status is not None:
        query_params["status"] = status

    for attr, value in sorted(kwargs.items()):
        value = _parse_argvalue(value)
        query_params[attr] = value

    url = base_url
    if query_params:
        url += f"&{urlencode(query_params)}"

    return url


def _get_next_page(links):
    for link in links:
        if link["rel"] == "next":
            return link["href"]
    return False


def _parse_date(date):
    if isinstance(date, datetime.datetime):
        return date
    elif isinstance(date, datetime.date):
        return datetime.datetime.combine(date, datetime.time())
    try:
        return dateutil.parser.parse(date)
    except ValueError:
        raise ValueError(
            "Date {date} is not in a valid format. Use Datetime object or iso string"
        )


def _add_time(date):
    if date.hour == 0 and date.minute == 0 and date.second == 0:
        date = date + datetime.timedelta(hours=23, minutes=59, seconds=59)
        return date
    return date


def _tastes_like_wkt_polygon(geometry):
    try:
        return geometry.replace(", ", ",").replace(" ", "", 1).replace(" ", "+")
    except Exception:
        raise ValueError("Geometry must be in well-known text format")


def _parse_geometry(geom):
    try:
        # If geom has a __geo_interface__
        return shape(geom).wkt
    except AttributeError:
        if _tastes_like_wkt_polygon(geom):
            return geom
        raise ValueError(
            "geometry must be a WKT polygon str or have a __geo_interface__"
        )


def _parse_argvalue(value):
    if isinstance(value, string_types):
        value = value.strip()
        if not any(
            value.startswith(s[0]) and value.endswith(s[1])
            for s in ["[]", "{}", "//", "()"]
        ):
            value.replace(" ", "+")
        return value
    elif isinstance(value, (list, tuple)):
        # Handle value ranges
        if len(value) == 2:
            value = "[{},{}]".format(*value)
            return value
        else:
            raise ValueError(
                "Invalid number of elements in list. Expected 2, received " "{}".format(
                    len(value)
                )
            )
    else:
        raise ValueError(
            "Additional arguments can be either string or tuple/list of 2 values"
        )


##### Functions associated with downloads
def _get_token(username, password):
    token_data = {
        "client_id": "cdse-public",
        "username": username,
        "password": password,
        "grant_type": "password",
    }
    response = requests.post(TOKEN_URL, data=token_data).json()
    try:
        return response["access_token"]
    except KeyError:
        raise RuntimeError(f"Unable to get token. Response was {response}")


def _download_raw_data(url, outfile, show_progress):
    """Downloads data from url to outfile.incomplete and then moves to outfile"""
    outfile_temp = str(outfile) + ".incomplete"
    try:
        downloaded_bytes = 0
        with requests.get(url, stream=True, timeout=100) as req:
            # analyze status code
            if req.status_code == 200:
                with tqdm(
                    unit="B", unit_scale=True, disable=not show_progress
                ) as progress:
                    chunk_size = 2**20  # download in 1 MB chunks
                    with open(outfile_temp, "wb") as fout:
                        for chunk in req.iter_content(chunk_size=chunk_size):
                            if chunk:  # filter out keep-alive new chunks
                                fout.write(chunk)
                                progress.update(len(chunk))
                                downloaded_bytes += len(chunk)

                shutil.move(outfile_temp, outfile)

            else:
                print(f"Download failed: response was {req.status_code}")

    finally:
        try:
            Path(outfile_temp).unlink()
        except OSError:
            pass


def download(uid, token, outfile, show_progress=True):
    """Download a file from CreoDIAS to the given location

    Parameters
    ----------
    uid:
        CreoDIAS UID to download
    username:
        Username
    password:
        Password
    outfile:
        Path where incomplete downloads are stored
    """
    url = f"{DOWNLOAD_URL}/{uid}?token={token}"
    _download_raw_data(url, outfile, show_progress)


def download_list(
    uids,
    username,
    password,
    token,
    session_start_time,
    outdir,
    threads=1,
    show_progress=True,
    log_file=False,
):
    """Downloads a list of UIDS

    Parameters
    ----------
    uids:
        A list of UIDs
    username:
        Username
    password:
        Password
    outdir:
        Output direcotry


    Returns
    -------
    dict
        mapping uids to paths to downloaded files
    """

    def _token_age(session_start_time):
        return (time.time() - session_start_time) / 60

    if token is None:
        print("Generating session token")
        token = _get_token(username, password)
        session_start_time = time.time()
    print(f"Session token age: {_token_age(session_start_time):.2f} minutes")

    if show_progress:
        if len(uids) > 0:
            pbar = tqdm(
                total=len(uids),
                desc=f"Downloading files to {os.path.basename(outdir)}",
                unit="file",
            )

    for uid in uids:
        # assess age of token, (Expires every ten minutes so we refresh every 9 minutes)
        if _token_age(session_start_time) > 9:
            print(
                f"\nSession token age: {_token_age(session_start_time):.2f} minutes. Refreshing session token now."
            )
            session_start_time = time.time()
            token = _get_token(username, password)

        # download file
        outfile = Path(outdir) / f"{uid}.zip"
        download(uid, token=token, outfile=outfile, show_progress=False)

        if log_file:
            with open(log_file, "a") as log:
                log.write(uid + "\n")  # add the id to the downloaded log

        if show_progress:
            pbar.update(1)

    return token, session_start_time
