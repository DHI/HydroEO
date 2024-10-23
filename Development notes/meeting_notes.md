# Meeting Notes and project design 

## Meeting on 23-10-2024

### Main goal
- provide water level timeseries from multimission altimetry datasets
- Focus on reservoirs first
- include rivers as bonus [keep in mind when developing]

### Datasets to include
- SWORD database
- SWOT lake and river reach averaged products
- Icesat-ATL13 [consider ATL22]
- Sentinel-3
- Sentinel-6

### Clear documentation of what the project is intended for
- explain through easy to understand readmes and a summary document or slide deck for presentations
- what can we use it for?
- what can you not use it for?
- What are the next steps

### Start to answer questions on filling the gaps
- how to identify what is an expected timeseries - raise flag
- point out when reservoir is not in sword database
- point out where there is no data available as part of output

### Main steps
1) Start from sword database
    - make sure that user can download this database and it will serve as a reference for downloads, processing and how data is provided

2) allow user to add own shapefile of reservoirs (and rivers)
    - match these shape files with sword database

3) download data products for each reservoir (and reach)

4) subset data, group and save raw data in table format for each reservoir (and reach)

5) combine individual sensor timeseries into a single reach and reservoir level timeseries
    - Provide an average timeseries but also individual satellite data point to user to assess what they want
    - will likely need to bias correct timeseries
    - perhaps an option to omit specific data based on what they see

6) provide user with flags for the downloads and timeseries aggregation to identify if a dataproduct is suspect for this element or missing satellite data

7) implement methods for updating timeseries files with new data