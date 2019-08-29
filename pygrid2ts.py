# -*- coding: utf-8 -*-

from features.grid import Grid, get_grids
from features.timeseries import TimeSeries
from glob import glob
from rasterstats import zonal_stats
from joblib import Parallel, delayed
import geopandas as gpd


class Zonal_Stat(object):
    def __init__(self, gdf, grid):
        self.shp = gdf
        self.crs = None
        self.zstat = None
        self.file = None
        self.gdf = gdf
        self.grid = grid
        
    @staticmethod
    def test(gdf, grid):
        if grid.crs == gdf.crs:
            return True
        else:
            gdf = gdf.to_crs(grid.crs)
            print('Converted shapefile crs to grid crs: {0}'.format(grid.crs))
        
    def get_stats(self):
        zs = zonal_stats(self.gdf, self.grid.fpath, raster_out=True)
        return zs

def get_zs(gdf, grid):
    zs = Zonal_Stat(gdf, grid)
    zs.test(gdf, grid)
    return zs.get_stats()
    
    
    
class ParseTS(object):
    def __init__(self, files):
        ts = TimeSeries(files)
        ts.import_ts()
        self.ts = ts

def main():
    files = glob(r"D:\ririe\rasters\dissolved_watershed\*alb.tif")
    
    ts = ParseTS(files)
        
    glist = Parallel(n_jobs=-1, verbose=10)(delayed(get_grids)(i) for i in ts.ts.flist)
    
    gdf = gpd.read_file(r"D:\ririe\shp\total_watershed_dissolved.shp")
    zs_list = []
    
    for grid in glist:
        zs_list.append(get_zs(gdf,grid))
        



    

    


if __name__ == '__main__':
    main()






        