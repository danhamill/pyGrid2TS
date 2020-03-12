# pyGrid2TS
Python tools for time series analyses of gridded data

Welcome to pyGrid2TS, the python package that can process gridded datasets in to DSS time series. The goal of this package is to leverage python geoprocessing tools and parallel processing to develop time series for a variety of climatic gridded datasets.

## Dependencies
- rasterio
- rasterstats
- geopandas
- joblib
- pydsstools


## Useage

### Input

1. Directory of daily SWE surfaces
2. Shapefile of watershed boundary
3. Shapefile of subbasin boundaries

### Output

1. Daily clipped raster for watershed and subbasin boundaries
2. DSS file


```python
(base) C:\workspace\git_clones\pyGrid2TS>python pygrid2ts.py -h
usage: python pygrid2ts.py --basin-shp total_watershed_dissolved.shp --sbasin-shp total_watershed.shp --ds SONDAS --raster-dir snodas_dir --dss output_timeseries.dss --oroot output_dir --bname RIRIE --dtfmt %Y%m%d --mconv 1000.0

Zonal Statistic Software to Create DSS Time Series

optional arguments:
  -h, --help            show this help message and exit
  --basin-shp BASIN_SHP
                        Dissolved Basin Shapefile
  --sbasin-shp SBASIN_SHP
                        Sub Basin Shapefile
  --ds DS               Raster Datasource
  --raster-dir FPATH    Directory for rasters to be processed
  --dss DSS_FILE        Path to Output DSS file
  --oroot OROOT         Output Directory for Basin Rasters
  --bname BASIN         Basin Name
  --dtfmt DTFMT         strftime Date Format
  --mconv M_CONV        Converstion Factor for Input Rasters to Meters
 ```

