help_message = """
Download products from your Hydroweb.next project (https://hydroweb.next.theia-land.fr) using EODAG (https://github.com/CS-SI/eodag)
This script is an example tuned for your last Hydroweb.next project but feel free to adapt it for future requests.
Follow these steps:
1. If not already done, install EODAG and packaging latest version using `pip install -U eodag packaging` or `conda update eodag packaging`
2a. Generate an API-Key from Hydroweb.next portal in your user settings
2b. Carefully store your API-Key
- either in your eodag configuration file (usually ~/.config/eodag/eodag.yml, automatically generated the first time you use eodag) in auth/credentials/apikey="PLEASE_CHANGE_ME"
- or in an environment variable `export EODAG__HYDROWEB_NEXT__AUTH__CREDENTIALS__APIKEY="PLEASE_CHANGE_ME"`
3. You can change download directory by modifying the variable path_out. By default, current path is used.
4. You are all set, run this script `python download_SWOT_Prior_Lake_Database.py`

For more information, please refer to EODAG Documentation https://eodag.readthedocs.io
"""
import os
from eodag import EODataAccessGateway, SearchResult
from eodag import __version__ as eodag_version
from eodag import setup_logging

from altimetry import utils

def download(download_dir : str, bounds : list):

    # create download directory if needed
    utils.ifnotmakedirs(download_dir)

    # Set timeout to 30s
    os.environ["EODAG__HYDROWEB_NEXT__SEARCH__TIMEOUT"] = "30"

    dag = EODataAccessGateway()

    # ---------------------------------------------------------------------------------------------------------------------
    # This command will perform a search using provided query arguments.
    # - specify a collection using the `productType` key
    # - add time restrictions using the `start` and `end` keys (e.g. "start": "2020-05-01" , "end": "2020-05-10T00:00:00Z",
    #   UTC ISO8601 format)
    # - add spatial restrictions using the "geom" key (e.g. "geom": "POLYGON ((1 43, 2 43, 2 44, 1 44, 1 43))" WKT string,
    #   a bounding-box list [lonmin, latmin, lonmax, latmax] can also be passed )
    # - more query arguments can be used, see
    #   https://eodag.readthedocs.io/en/stable/notebooks/api_user_guide/4_search.html?#Search-parameters
    query_args = {
        "productType": "SWOT_PRIOR_LAKE_DATABASE",
        "geoms" : bounds,
        "items_per_page": 2000, # Default search criteria when iterating over collection pages
    }

    # Iterate over all pages to find all products
    search_results = ([])
    for i, page_results in enumerate(dag.search_iter_page(**query_args)):
        search_results.extend(page_results)

    # This command actually downloads the matching products
    downloaded_paths = dag.download_all(search_results, outputs_prefix=download_dir)

