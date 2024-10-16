# rk-altimetry: Easy access altimetry for water resource applications

Repo to allow users with little EO (Earth Observation) knowledge to access and download altimetry over rivers and reservoirs for integration into larger water resource projects.

## What altimetry missions are available?

### ATLAS/ICESat-2 L3A Inland Water Surface Height (ATL13)
ATLAS/ICESat-2 L3A Inland Water Surface Height (ATL13 - https://nsidc.org/data/atl13/versions/6) can be accessed from 13 October 2018 to 11 November 2020 within [N:90S:-90E:180W:-180]. The satellite is orignially intended to monitor changes in the cryosphere but has proven extremely usefull for inland water applications. The ATL13 product contains along-track water surface heights and descriptive statistics for inland water bodies. Inland water bodies include lakes, reservoirs, rivers, bays, estuaries and a 7km near-shore buffer. Given its intended use, the orbit design for ICESat-2 prioritizes latitudes closer to the poles, meaning that if your area of interest is closer to the equator you can expect a more temporally sparse timeseries for a given virtual station. The satellite operates with a repeat cycle of 91 days meaning that you can expect a maximum of 91 days between water surface height observations. This product is a good start for the development of an interal altimetry utility because the product is already processed and provides reliable levels without much expert inspection. For more information on this satellite see: https://nsidc.org/data/icesat-2 and https://icesat-2.gsfc.nasa.gov/.

### (Development not begun) Sentinel-3A, 3B, 6
The sentinel series can provide more temporally dense inland water observations. Currently there is no support for downloading and including this data within virtual stations.

## What is currently included in the library?

## What could come?
- (Under development) - Provide shape file of river or rivers of interest to initiate search and dowload of available
- (Under development) - Get report of all available crossings over river of interest
- (Under development) - Suggest grouping of crossings into virtual stations along river
- (Under development) - Group observations of ICESat-2 crossings over river of interest into selected virtual stations
- (Under development) - Organize historical observations for virtual stations and provide reports for each virtual station
- (Under developmemt) - Allow for near real time updates of water levels along river

- (Under development) - provide similar support for reservoirs
. (Under development) - provide means of grouping reservoirs and rivers within a project class that can easily process a basin or mike project area to provide usefull data to mike cloud

## Python development resources

If you're interested in learning more about best practices for developing Python packages, check out the following resources:

- [Python Package Development at DHI](https://dhi.github.io/python-package-development/)
- [Scientific Python Library Development Guide](https://learn.scientific-python.org/development/)
