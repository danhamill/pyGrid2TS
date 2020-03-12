# -*- coding: utf-8 -*-
#author          :Daniel Hamill
#email           :Daniel.D.Hamill@usace.army.mil
import os
import rasterio


class Grid(object):
    def __init__(self, filepath, date, dtfmt = '%Y%m%d'):
        self.name = None
        self.crs = None
        self.date = date
        self.date_fmt = None
        self.driver = None
        self.fpath = filepath
        self.data = None
        self.affine = None
        self.dtfmt = dtfmt


    @staticmethod
    def test(filepath):
        try:
            ds = rasterio.open(filepath)
            driver = ds.driver
        except:
            driver = "Not Found"
        return driver

    def import_grid_info(self):
        self.name = os.path.basename(self.fpath)
        with rasterio.open(self.fpath) as ds:
            self.crs = ds.crs
            self.driver = ds.driver
            self.data = ds.read(1)


def get_grids(fname, date, dtfmt = '%Y%m%d'):
    g = Grid(fname, date, dtfmt)
    g.import_grid_info()
    return g
