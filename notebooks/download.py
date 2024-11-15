import os

from altimetry.project import Project
from altimetry.utils.downloads import hydroweb


##### Initialize project and reservoir system from config file
mekong = Project(name="Mekong", config=r".\data\config.yaml")
reservoirs = mekong.reservoirs
print(mekong.report())


##### Download or load PLD and assign as a download data frame
# determine if we need to download or simply load the pld
print("Assessing SWOT Prior Lake Database (PLD) for download geometries and ids")
pld_path = mekong.reservoirs.dirs["pld"]
if not os.path.exists(pld_path):
    download_dir = os.path.dirname(pld_path)
    file_name = os.path.basename(pld_path)
    bounds = list(mekong.reservoirs.gdf.unary_union.bounds)

    hydroweb.download_PLD(download_dir=download_dir, file_name=file_name, bounds=bounds)

# load the pld and associate the reservoirs with a "lake id"
reservoirs.assign_pld_id(
    local_crs=mekong.local_crs, max_distance=100
)  # take the downloaded PLD and see where we have overlap with the input reservoirs
reservoirs.flag_missing_priors()  # flag and export which reservoirs have entries in the PLD
reservoirs.assign_reservoir_polygons()  # set the download polygons from the swot database


##### Download data
print("Downloading data")

download = ["S3", "S6"]

# reservoirs.download_altimetry(product='SWOT_Lake', startdate=mekong.swot_startdate, enddate=mekong.swot_enddate)
# reservoirs.download_altimetry(product='ATL13', startdate=mekong.startdates['icesat2'], enddate=mekong.enddates['icesat2'], start_index=0)

if "S3" in download:
    reservoirs.download_altimetry(
        product="S3",
        startdate=mekong.startdates["sentinel3"],
        enddate=mekong.enddates["sentinel3"],
        credentials=(mekong.creodias_user, mekong.creodias_pass),
    )

if "S6" in download:
    reservoirs.download_altimetry(
        product="S6",
        startdate=mekong.startdates["sentinel6"],
        enddate=mekong.enddates["sentinel6"],
        credentials=(mekong.creodias_user, mekong.creodias_pass),
    )
