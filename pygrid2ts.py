# -*- coding: utf-8 -*-
import os
os.environ['GDAL_DATA'] = r"C:\Anaconda3\envs\gal_env\Library\share\gdal"
from features.grid import get_grids
from features.timeseries import ParseTS, dateutil
from glob import glob
from rasterstats import zonal_stats
from joblib import Parallel, delayed
import geopandas as gpd
from features.utils import check_crs
import rasterio


class Zonal_Stat(object):
    def __init__(self, gdf, grid, oRoot):
        self.shp = gdf
        self.crs = None
        self.zstat = None
        self.file = None
        self.gdf = gdf
        self.grid = grid
        self.clip_rast_path = oRoot
        
    @staticmethod
    def test(gdf, grid):
        if grid.crs == gdf.crs:
            return True
        else:
            gdf = gdf.to_crs(grid.crs)
            print('Converted shapefile crs to grid crs: {0}'.format(grid.crs))
            
    def get_basin_raster(self):
        tmp_gdf = None
        if len(self.gdf)>1:
            self.gdf.loc[:,'Dissove'] = 'dissolve_me'
            tmp_gdf = self.gdf.dissolve(by='disolve_me')
            print('!!! Input shapefile had multiple features, dissolved for total basin')
        if not tmp_gdf:
            zs = zonal_stats(self.gdf, self.grid.fpath, raster_out=True)
        else:
            zs = zonal_stats(tmp_gdf.gdf, self.grid.fpath, raster_out=True)
        stuff = zs[0]
        with rasterio.open(self.clip_rast_path + os.sep + 'Total_basin_' + self.grid.date.strftime('%Y-%m-%d') + '.tif', 'w',driver='GTiff', 
                               height = stuff['mini_raster_array'].shape[0], width = stuff['mini_raster_array'].shape[1], 
                               dtype=rasterio.dtypes.int16, 
                               crs=self.grid.crs.to_proj4(), count = 1, transform=stuff['mini_raster_affine']) as dst:
                    dst.write(stuff['mini_raster_array'].astype('int16').filled(-9999),1)
                    dst.nodata = -9999
        self.grid.data = stuff['mini_raster_array'].astype('int16')
        self.grid.affine = stuff['mini_raster_affine']
        
    def get_stats(self):
        zs = zonal_stats(self.gdf, self.grid.data, affine = self.grid.affine, raster_out=True,nodata=-9999)
        return zs

def get_zs(gdf, grid):
    zs = Zonal_Stat(gdf, grid)
    return zs.get_stats()
    
    
    

def main():
    files1 = glob(r"\\rsgis-base.crrel.usace.army.mil\study\snow\nohrsc_gdal\conus_tiffs\us_ssmv11034tS__T0001TTNATS*.tif")
    files2 = glob(r"\\rsgis-base.crrel.usace.army.mil\study\snow\nohrsc_gdal\conus_tiffs\zz_ssmv11034tS__T0001TTNATS*.tif")
    
    files = files1 + files2

    #files = glob(r"C:\Users\u4rreddh\Documents\temp\*.tif")

    ts = ParseTS(files,'%Y%m%d', month_start=9, month_end=6)
        
    glist = Parallel(n_jobs=-1, verbose=10)(delayed(get_grids)(fname, date) for fname, date in zip(ts.ts.flist, ts.ts.dates))
    
    
    gdf = gpd.read_file(r"D:\souris\shp\Darling_DAGrouping.shp")
    gdf = check_crs(gdf)    
    
    zs_list = Parallel(n_jobs=-1, verbose = 10)(delayed(get_zs)(gdf, grid) for grid in glist)
    

        

if __name__ == '__main__':
    main()






        