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
Information on what accounts must be created to dowload what products.
Include link to example notebook
Inlcue short example of simplist way to get started


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
