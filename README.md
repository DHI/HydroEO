# rk-altimetry: Easy access altimetry for water resource applications

Repo to allow users with little EO (Earth Observation) knowledge to access and download altimetry over reservoirs and lakes for integration into larger water resource projects.

> [!CAUTION]
> rk-altimetry is experimental and under development.
> * The package is expected to be ready for limited use by Janurary 2025

## What is currently included in the library?
- (Needs testing) - Provide shape file of reservoir or reservoirs of interest to initiate search and dowload of available data
- (Needs testing) - Supported Satellite products include
    - SWOT Lake SP product
    - ICESat-2 ATL13 inland water product
    - Sentinel 3A and 3B inland hydrology product
    - Sentinel 6 inland hydrology product
- (Needs testing) - Get report of all available crossings over reservoir(s) of interest
- (Under development) - Bias correct between satellite orbits and simple cleaning per product timeseries
- (Under development) - Estimation of reservoir state based on kalman filter for multi product timeseries over individual reservoirs
- (Needs testing) - Allow for near real time updates of water levels of alreadey initialized reservoirs

## What could come?
- (not started) - Report of data gaps or questionable timeseries
- (not started) - Implementation of downloads for virtual stations along a provided river shapefile

## How to get started
1) Familiarize yourself with the config file
The altimetry project is initialized entirely from information within the provided config file. The configuration file includes information on the project directory, reservoir shape file location, gis info as well as credentials for downloading data from various locations. Importantly, this is also where you will specify which satellite products you wish to download and process and for which dates to download.

Project information within the config file:
```
project :
  main_dir : "C:\\Users\\username\\altimetry_project" # main directory in which to store processed outputs

gis :
  global_crs : 'EPSG:4326'

reservoirs :
  path : "C:\\Users\\username\\altimetry_project\\reservoirs.shp" # path to the shapefile holding one or more resevoirs per feature
  id_key : 'project' # the key within the shapefile to the column that holds the unique reservoir ids
```


Example of how to provide credentials for ICESat-2 download
```
earthaccess: # create an Earth Data account at https://urs.earthdata.nasa.gov/ 
  username : ""
  password : ""
```

Example of specifying download criteria for ICESat-2
```
icesat2:
  download : True
  process : True
  download_dir : "C:\\Users\\username\\altimetry_project\\data\\icesat2" # directory in which to store raw ICESat-2 files (if not provided, a directory will be made within the project directory)
  startdate : [2024, 1, 1] # [year, month, day] format
  enddate   : [2024, 11, 01]
```

Please see the example configuration file in the notebooks folder for a complete example to download all support products

2) Create free accounts for data downloads
- Data provided by NASA (SWOT and ICESat-2) is accessed through the earthaccess python package which allows for easy downloads using your Earth Data Login credentials. You can register for a free Earth Data Login account at https://urs.earthdata.nasa.gov/. By default, earthaccess will look for your Earth Data Login credentials in a .netrc file, or in environment variables EARTHDATA_USERNAME and EARTHDATA_PASSWORD. If you do not set one of these before running, your credentials must be provided within the config file. 

- If wishing to download SWOT data, reference to the Prior Lake Database must be made. If you do not already have a downloaded version, a subset of the PLD will be downloaded from Hydroweb. If you do not have an account already, you should create an account here: https://hydroweb.next.theia-land.fr/. Once your account is made, navigate to the User settings and create an API key. Copy this key and add it to the corresponding spot within the config file.

- Data provided by ESA (Sentinel 3 and Sentinel 6) is accessed through the copernicus data space. Create a free account here: https://dataspace.copernicus.eu/. Credentials must be provided in the config file directly.

3) Gather waterbody polygons within a single shapefile
- reservoir information must be provided within a single shapefile in which the path is specified within the config file. A column must be specified in which a waterbodies unique id is provided. The project will automatically process data for all provided polygons.

4) Make your first downloads
- Write your config file
- Import the package
- Initialize your project
- Run your downloads
- See the example notebook for a complete example

```
from altimetry.project import Project
altimetry_project = Project(name="my_altimetry_project", config="config.yaml")
altimetry_project.initialize()
altimetry_project.download()
```

## More details into available satellite products

### Surface Water Ocean and Topography (SWOT) mission

### ATLAS/ICESat-2 L3A Inland Water Surface Height (ATL13)
![](images/icesat2-hqprint.jpg)
https://science.nasa.gov/wp-content/uploads/2023/06/icesat2-hqprint-print.jpg?w=4096&format=jpeg

ICESat-2 is the second of NASA's "Ice, Cloud, and Land Elevation Satellites" carrying a photon-counting laser altimeter called "ATLAS". While the satellite is orignially intended to monitor changes in the cryosphere, it has proven extremely usefull for inland water applications. The ATLAS instrument has six beams organized into 3 pairs. The distance between two beams within a pair is 90 meters and 3.3km between pairs. The 6 beam design increases the chances of capturing waterbodies and allows a user to infer the slope of a river reach from the 6 nearby crossing which proves usefull in 1D hydrodynamic modelling (https://www.nature.com/articles/s41597-023-02215-x). Given its intended use, the orbit design for ICESat-2 prioritizes latitudes closer to the poles, meaning that if your area of interest is closer to the equator you can expect a more temporally sparse timeseries for a given crossing. The satellite operates with a repeat cycle of 91 days meaning that you can expect a maximum of 91 days between water surface height observations, though if two orbits cross nearby the target, you can have more frequent observations. 

The ATLAS/ICESat-2 L3A Inland Water Surface Height (ATL13 - https://nsidc.org/data/atl13/versions/6) product provides along-track water surface heights and descriptive statistics for inland water bodies. Inland water bodies include lakes, reservoirs, rivers, bays, estuaries and a 7km near-shore buffer. This product is a good start for the development of an interal altimetry utility because the product is already processed and provides reliable levels without much expert inspection. The downside being the sparse temporal series making its use for creating virtual stations along rivers limited. For more information on this satellite see: https://nsidc.org/data/icesat-2 and https://icesat-2.gsfc.nasa.gov/.

Summary of mission details:
- Dates: October 2018 to present
- Spatial bounds: [N:90 S:-90 E:180 W:-180 degrees]
- Repeat orbit: 91-days
- 6 beams per pass

Example of ground track and beam setup for ICESat-2:
![](images/ICESat2BeamPattern.png)
https://icesat-2.gsfc.nasa.gov/science/specs


### Sentinel-3A, 3B
The sentinel series can provide more temporally dense inland water observations. Currently there is no support for downloading and including this data within virtual stations.


### Sentinel-6


## Python development resources
If you're interested in learning more about best practices for developing Python packages, check out the following resources:

- [Python Package Development at DHI](https://dhi.github.io/python-package-development/)
- [Scientific Python Library Development Guide](https://learn.scientific-python.org/development/)
