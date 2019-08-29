# -*- coding: utf-8 -*-
import os
import rasterio
import dateutil.parser as dparse

class Grid(object):
    def __init__(self, filepath):
        self.name = None
        self.crs = None
        self.date = None
        self.date_fmt = None
        self.driver = None
        self.fpath = filepath

        
    @staticmethod
    def test(filepath):
        try:
            ds = rasterio.open(filepath)
            driver = ds.driver
        except:
            driver = "Not Found"
        return driver
    
    def import_grid_info(self, dtfmt = '%Y%m%d'):
        self.name = os.path.basename(self.fpath)
        bname = os.path.basename(self.fpath)
        self.date = dparse.parse(bname.split('.')[0],fuzzy=True)
        with rasterio.open(self.fpath) as ds:
            self.crs = ds.crs
            self.driver = ds.driver
            
            
def get_grids(fpath):
    g = Grid(fpath)
    g.import_grid_info()
    return g
