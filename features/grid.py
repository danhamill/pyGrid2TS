# -*- coding: utf-8 -*-
import os
import rasterio


class Grid(object):
    def __init__(self, filepath, date):
        self.name = None
        self.crs = None
        self.date = date
        self.date_fmt = None
        self.driver = None
        self.fpath = filepath
        self.data = None
        self.affine = None

        
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
        with rasterio.open(self.fpath) as ds:
            self.crs = ds.crs
            self.driver = ds.driver           
            
            
def get_grids(fname, date):
    g = Grid(fname, date)
    g.import_grid_info()
    return g
