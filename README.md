# rk-altimetry: Easy access altimetry for water resource applications

Repo to allow users with little EO (Earth Observation) knowledge to access and download altimetry over rivers and reservoirs for integration into larger water resource projects.

## What altimetry missions are available?

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


### (Development not begun) Sentinel-3A, 3B, 6
The sentinel series can provide more temporally dense inland water observations. Currently there is no support for downloading and including this data within virtual stations.

## What is currently included in the library?
- (Needs testing) - Provide shape file of river or rivers of interest to initiate search and dowload of available
- (Needs testing) - Get report of all available crossings over river of interest
- (Needs testing) - provide similar support for reservoirs
- (Needs testing) - provide means of grouping reservoirs and rivers within a project class that can easily process a basin or mike project area to provide usefull data to mike cloud
- (Under development) - Implement S3A, S3B and S6 data downloads

## What could come?
- (not started) - Allow for near real time updates of water levels along river
- (not started) - Suggest grouping of crossings into virtual stations along river by reach?
- (not started) - Group crossings over river of interest into selected virtual stations
- (not started) - Organize historical observations for virtual stations and provide reports for each virtual station

## Python development resources
If you're interested in learning more about best practices for developing Python packages, check out the following resources:

- [Python Package Development at DHI](https://dhi.github.io/python-package-development/)
- [Scientific Python Library Development Guide](https://learn.scientific-python.org/development/)
