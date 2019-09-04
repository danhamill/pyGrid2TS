# -*- coding: utf-8 -*-
#author          :Daniel Hamill
#email           :Daniel.D.Hamill@usace.army.mil


import os
os.environ['GDAL_DATA'] = r"C:\Anaconda3\envs\gal_env\Library\share\gdal"
from features.grid import get_grids
from features.timeseries import ParseTS
from glob import glob
from rasterstats import zonal_stats
from joblib import Parallel, delayed
import geopandas as gpd
from features.utils import check_crs
import rasterio
import numpy as np
import time


class Zonal_Stat(object):
    def __init__(self, basin_gdf, sbasin_gdf, grid, oRoot):
        self.crs = None
        self.zstat = None
        self.file = None
        self.basin_file = None
        self.basin_gdf = basin_gdf
        self.sbasin_gdf = sbasin_gdf
        self.grid = grid
        self.clip_rast_path = oRoot
        self.basin_avg = None
        self.sbasin_avg = []
        self.basin_vol = None
        self.sbasin_vol = []
        self.sbasin_names = []

    @staticmethod
    def test(gdf, grid):
        if grid.crs == gdf.crs:
            return True
        else:
            gdf = gdf.to_crs(grid.crs)
            print('Converted shapefile crs to grid crs: {0}'.format(grid.crs))

    def get_basin_raster(self):
        '''
        Function to extract basin raster and save to
        '''
        assert len(self.basin_gdf) == 1, "Basin shapefile has more than one feature!!"

        zs = zonal_stats(self.basin_gdf, self.grid.fpath, raster_out=True)
        stuff = zs[0]
        try:
            with rasterio.open(self.clip_rast_path + os.sep + 'Total-basin-' + self.grid.date.strftime('%Y-%m-%d') + '.tif', 'w',driver='GTiff',
                                   height = stuff['mini_raster_array'].shape[0], width = stuff['mini_raster_array'].shape[1],
                                   dtype=rasterio.dtypes.int16,
                                   crs=self.grid.crs.to_proj4(), count = 1, transform=stuff['mini_raster_affine']) as dst:
                        dst.write(stuff['mini_raster_array'].astype('int16').filled(-9999),1)
                        dst.nodata = -9999
        except:
            time.sleep(1)
            with rasterio.open(self.clip_rast_path + os.sep + 'Total-basin-' + self.grid.date.strftime('%Y-%m-%d') + '.tif', 'w',driver='GTiff',
                                   height = stuff['mini_raster_array'].shape[0], width = stuff['mini_raster_array'].shape[1],
                                   dtype=rasterio.dtypes.int16,
                                   crs=self.grid.crs.to_proj4(), count = 1, transform=stuff['mini_raster_affine']) as dst:
                        dst.write(stuff['mini_raster_array'].astype('int16').filled(-9999),1)
                        dst.nodata = -9999
        self.grid.data = stuff['mini_raster_array'].astype('int16')
        self.grid.affine = stuff['mini_raster_affine']
        self.basin_avg = self.grid.data.mean()/1000
        self.basin_vol = self.grid.data.filled(0).sum()*self.grid.affine[0]*self.grid.affine[0]/1000

    def get_sbasin_stats(self):
        assert self.sbasin_gdf.columns.contains('name'), "No Name field in subbasin shapefile"
        zs = zonal_stats(self.sbasin_gdf, self.grid.data, affine = self.grid.affine, raster_out=True,nodata=-9999, all_touched=False)
        for item, name  in zip(zs,self.sbasin_gdf.name.tolist()):
            self.sbasin_names.append(name)
            self.sbasin_avg.append(item['mean']/1000)
            vol = item['mini_raster_array'].astype('int16').filled(0).sum()*self.grid.affine[0]*self.grid.affine[0]/1000
            self.sbasin_vol.append(vol)
            with rasterio.open(self.clip_rast_path + os.sep + name +'-' + self.grid.date.strftime('%Y-%m-%d') + '.tif', 'w',driver='GTiff',
                               height = item['mini_raster_array'].shape[0], width = item['mini_raster_array'].shape[1],
                               dtype=rasterio.dtypes.int16,
                               crs=self.grid.crs.to_proj4(), count = 1, transform=item['mini_raster_affine']) as dst:
                    dst.write(item['mini_raster_array'].astype('int16').filled(-9999),1)
                    dst.nodata = -9999

def get_zs(basin_gdf, sbasin_gdf, grid, oRoot):
    zs = Zonal_Stat(basin_gdf, sbasin_gdf, grid, oRoot)
    zs.get_basin_raster()
    zs.get_sbasin_stats()
    return zs



def main():
#    files1 = glob(r"\\rsgis-base.crrel.usace.army.mil\study\snow\nohrsc_gdal\conus_tiffs\us_ssmv11034tS__T0001TTNATS*.tif")
#    files2 = glob(r"\\rsgis-base.crrel.usace.army.mil\study\snow\nohrsc_gdal\conus_tiffs\zz_ssmv11034tS__T0001TTNATS*.tif")
##
#    files = files1 + files2

    files = glob(r"D:\snodas\*.tif")

    ts = ParseTS(files,'%Y%m%d', month_start=9, month_end=6)

    glist = Parallel(n_jobs=-1, verbose=10)(delayed(get_grids)(fname, date) for fname, date in zip(ts.ts.flist, ts.ts.dates))


    basin_gdf = gpd.read_file(r"D:\souris\shp\Darling_dissolved.shp")
    basin_gdf = check_crs(basin_gdf)
    basin_gdf.columns = basin_gdf.columns.str.lower()

    sbasin_gdf = gpd.read_file(r"D:\souris\shp\Darling_DAGrouping.shp")
    sbasin_gdf = check_crs(sbasin_gdf)
    sbasin_gdf.columns = sbasin_gdf.columns.str.lower()


    oRoot = r"D:\souris\rasters\Total_Watershed"
    zs_list = Parallel(n_jobs=-1, verbose = 10)(delayed(get_zs)(basin_gdf, sbasin_gdf, grid, oRoot) for grid in glist)

#    zs_list = []
#    for grid in glist[30:31]:
#        print(grid.date)
#        zs = Zonal_Stat(basin_gdf, sbasin_gdf, grid, oRoot)
#        zs.get_basin_raster()
#        zs.get_sbasin_stats()
#        zs_list.append(zs)



if __name__ == '__main__':
    main()
