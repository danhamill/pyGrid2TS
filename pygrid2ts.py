# -*- coding: utf-8 -*-
#author          :Daniel Hamill
#email           :Daniel.D.Hamill@usace.army.mil


import os
os.environ['GDAL_DATA'] = r"C:\Anaconda3\envs\gal_env\Library\share\gdal"
from features.grid import get_grids
from features.timeseries import ParseTS
from features.utils import date_util
from glob import glob
from rasterstats import zonal_stats
from joblib import Parallel, delayed
import geopandas as gpd
from features.utils import check_crs
import rasterio
import numpy as np
import time
from features.utils import zstat2dss
import pandas as pd
from argparse import ArgumentParser
import datetime
import argparse

class Zonal_Stat(object):
    def __init__(self, basin_gdf, sbasin_gdf, grid, oRoot, ds, basin, m_conv):
        """Short summary.

        Parameters
        ----------
        basin_gdf : type
            Description of parameter `basin_gdf`.
        sbasin_gdf : type
            Description of parameter `sbasin_gdf`.
        grid : type
            Description of parameter `grid`.
        oRoot : type
            Description of parameter `oRoot`.
        ds : type
            Description of parameter `ds`.
        basin : type
            Description of parameter `basin`.
        m_conv : type
            Description of parameter `m_conv`.

        Returns
        -------
        type
            Description of returned object.

        """
        self.basin = basin
        self.ds = ds
        self.crs = None
        self.zstat = None
        self.file = None
        self.basin_file = None
        self.basin_gdf = basin_gdf
        self.sbasin_gdf = sbasin_gdf
        self.grid = grid
        self.clip_rast_path = oRoot
        self.basin_avg = None
        self.basin_dates = []
        self.sbasin_avg = []
        self.basin_vol = None
        self.sbasin_vol = []
        self.sbasin_names = []
        self.sbasin_dates = []
        self.m_conv = float(m_conv)

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

        raster_dtype = stuff['mini_raster_array'].dtype

        if np.float32  == stuff['mini_raster_array'].dtype:
            raster_dtype = rasterio.dtypes.float32
            array_dtype = 'float32'
        elif np.int16 == stuff['mini_raster_array'].dtype:
            raster_dtype = rasterio.dtypes.int16
            array_dtype = 'int16'
        elif np.int64 == stuff['mini_raster_array'].dtype:
            raster_dtype = rasterio.dtypes.int16
            array_dtype = 'int16'
        elif np.float64 == stuff['mini_raster_array'].dtype:
            raster_dtype = rasterio.dtypes.float32
            array_dtype = 'float32'

        try:
            os.makedirs(self.clip_rast_path + os.sep + 'Total_Watershed' + os.sep + self.ds)
        except:
            pass

        try:
            with rasterio.open(self.clip_rast_path + os.sep + 'Total_Watershed' + os.sep + self.ds + os.sep + 'Total-basin-' + self.grid.date.strftime('%Y-%m-%d') + '.tif', 'w',driver='GTiff',
                                   height = stuff['mini_raster_array'].shape[0], width = stuff['mini_raster_array'].shape[1],
                                   dtype=raster_dtype,
                                   crs=self.grid.crs.to_proj4(), count = 1, transform=stuff['mini_raster_affine']) as dst:
                        dst.write(stuff['mini_raster_array'].astype(array_dtype).filled(-9999),1)
                        dst.nodata = -9999
        except:
            time.sleep(0.5)
            with rasterio.open(self.clip_rast_path + os.sep + 'Total_Watershed' + os.sep + self.ds + os.sep + 'Total-basin-' + self.grid.date.strftime('%Y-%m-%d') + '.tif', 'w',driver='GTiff',
                                   height = stuff['mini_raster_array'].shape[0], width = stuff['mini_raster_array'].shape[1],
                                   dtype=raster_dtype,
                                   crs=self.grid.crs.to_proj4(), count = 1, transform=stuff['mini_raster_affine']) as dst:
                        dst.write(stuff['mini_raster_array'].astype(array_dtype).filled(-9999),1)
                        dst.nodata = -9999
        self.grid.data = stuff['mini_raster_array'].astype(array_dtype)
        self.grid.affine = stuff['mini_raster_affine']
        self.basin_avg = stuff['mini_raster_array'].astype(array_dtype).mean()/self.m_conv
        self.basin_vol = self.grid.data.filled(0).sum()*self.grid.affine[0]*self.grid.affine[0]/self.m_conv
        self.basin_dates.append(self.grid.date)

    def get_sbasin_stats(self):
        assert self.sbasin_gdf.columns.contains('name'), "No Name field in subbasin shapefile"
        zs = zonal_stats(self.sbasin_gdf, self.grid.data, affine = self.grid.affine, raster_out=True,nodata=-9999, all_touched=False)

        stuff = zs[0]
        if np.float32  == stuff['mini_raster_array'].dtype:
            raster_dtype = rasterio.dtypes.float32
            array_dtype = 'float32'
        elif np.int16 == stuff['mini_raster_array'].dtype:
            raster_dtype = rasterio.dtypes.int16
            array_dtype = 'int16'
        elif np.int64 == stuff['mini_raster_array'].dtype:
            raster_dtype = rasterio.dtypes.int16
            array_dtype = 'int16'
        elif np.float64 == stuff['mini_raster_array'].dtype:
            raster_dtype = rasterio.dtypes.float32
            array_dtype = 'float32'
        del stuff

        for item, name  in zip(zs,self.sbasin_gdf.name.str.replace(' ' , '-').tolist()):
            self.sbasin_names.append(name)
            if item['mean'] is None :
                self.sbasin_avg.append(0)
            else:
                self.sbasin_avg.append(item['mean']/self.m_conv)
            vol = item['mini_raster_array'].astype(array_dtype).filled(0).sum()*self.grid.affine[0]*self.grid.affine[0]/self.m_conv
            self.sbasin_vol.append(vol)
            self.sbasin_dates.append(self.grid.date)

            try:
                os.makedirs(self.clip_rast_path + os.sep + name.replace(' ' , '-') + os.sep + self.ds)
            except:
                pass
            try:

                with rasterio.open(self.clip_rast_path + os.sep + name + os.sep + self.ds + os.sep + name +'-' + self.grid.date.strftime('%Y-%m-%d') + '.tif', 'w',driver='GTiff',
                                   height = item['mini_raster_array'].shape[0], width = item['mini_raster_array'].shape[1],
                                   dtype=raster_dtype,
                                   crs=self.grid.crs.to_proj4(), count = 1, transform=item['mini_raster_affine']) as dst:
                        dst.write(item['mini_raster_array'].astype(array_dtype).filled(-9999),1)
                        dst.nodata = -9999
            except:
                time.sleep(0.5)
                with rasterio.open(self.clip_rast_path + os.sep + name + os.sep + self.ds + os.sep + name +'-' + self.grid.date.strftime('%Y-%m-%d') + '.tif', 'w',driver='GTiff',
                                   height = item['mini_raster_array'].shape[0], width = item['mini_raster_array'].shape[1],
                                   dtype=raster_dtype,
                                   crs=self.grid.crs.to_proj4(), count = 1, transform=item['mini_raster_affine']) as dst:
                        dst.write(item['mini_raster_array'].astype(array_dtype).filled(-9999),1)
                        dst.nodata = -9999

def get_zs(basin_gdf, sbasin_gdf, grid, oRoot,ds, basin, m_conv=1):
    zs = Zonal_Stat(basin_gdf, sbasin_gdf, grid, oRoot,ds, basin, m_conv)
    zs.get_basin_raster()
    zs.get_sbasin_stats()
    return zs



def main(a):
    basin_shp,sbasin_shp,ds,dss_file,oRoot,basin, dtfmt, m_conv, fpath = a.basin_shp,a.sbasin_shp,a.ds,a.dss_file,a.oRoot,a.basin, a.dtfmt, a.m_conv, a.fpath
    basin_gdf = gpd.read_file(basin_shp)
    basin_gdf = check_crs(basin_gdf)
    basin_gdf.columns = basin_gdf.columns.str.lower()

    sbasin_gdf = gpd.read_file(sbasin_shp)
    sbasin_gdf = check_crs(sbasin_gdf)
    sbasin_gdf.columns = sbasin_gdf.columns.str.lower()

    files = glob(fpath + os.sep + '*.tif')

    assert len(files)> 1, "Input Raster Directory {0} is empty".format(fpath)
    assert isinstance(date_util(files[0], dtfmt), datetime.datetime), "Date format {0} could not be found in raster {1}".format(dtfmt, files[0])

    ts = ParseTS(files,dtfmt, month_start=9, month_end=6)
    glist = Parallel(n_jobs=-1, verbose=10)(delayed(get_grids)(fname, date) for fname, date in zip(ts.ts.flist, ts.ts.dates))
    zs_list = Parallel(n_jobs=-1, verbose = 10)(delayed(get_zs)(basin_gdf, sbasin_gdf, grid, oRoot,ds, basin, 1e3) for grid in glist)
    sbasin, tbasin =zstat2dss(zs_list, basin, ds,dss_file)
    sbasin.to_csv(oRoot + os.sep + ds + "_sbasin_zonal_stats.csv")
    tbasin.to_csv(oRoot + os.sep + ds + "_tbasin_zonal_stats.csv")
#    #checking area
#    zs = zs_list[2324]
#    shp1 = gpd.read_file(r"E:\ririe\shp\total_watershed.shp")
#    shp1 = shp1.to_crs(sbasin_gdf.crs)
#    all_shp = gpd.GeoDataFrame(pd.concat([sbasin_gdf, shp1], ignore_index=True))
#    tmp = zonal_stats(all_shp, zs.grid.data, affine =zs.grid.affine, raster_out=True, nodata=-9999, stats=['mean','count','sum'])
#
#
    # ds = 'SNODAS'
    # files = glob(r"E:\ririe\rasters\testing\Total_Watershed\SNODAS\*.tif")
    #
    # ts = ParseTS(files,'%Y-%m-%d', month_start=9, month_end=6)
    # glist = Parallel(n_jobs=-1, verbose=10)(delayed(get_grids)(fname, date) for fname, date in zip(ts.ts.flist, ts.ts.dates))
    # zs_list = Parallel(n_jobs=-1, verbose = 10)(delayed(get_zs)(basin_gdf, sbasin_gdf, grid, oRoot,ds, basin, 1e3) for grid in glist)
    # snodas_sbasin, snodas_tbasin = zstat2dss(zs_list, basin, ds, dss_file)
    #
    # ds = 'SNOTEL_IDW'
    # files = glob(r'E:\ririe\rasters\testing\Total_Watershed\SNOTEL_IDW\*.tif')
    # ts = ParseTS(files,'%Y-%m-%d', month_start=9, month_end=6)
    # glist = Parallel(n_jobs=-1, verbose=10)(delayed(get_grids)(fname, date) for fname, date in zip(ts.ts.flist, ts.ts.dates))
    # zs_list = Parallel(n_jobs=-1, verbose = 10)(delayed(get_zs)(basin_gdf, sbasin_gdf, grid, oRoot,ds, basin, 1) for grid in glist)
    # snotel_sbasin, snotel_tbasin = zstat2dss(zs_list, basin, ds+'_grid_res', dss_file)
    #
    # ds = 'SNOTEL_NEW_ELEV_IDW'
    # files = glob(r"E:\ririe\rasters\SNOTEL_GRIDS\IDW_NEW_ELEV\*.tif")
    # ts = ParseTS(files,'%Y-%m-%d', month_start=9, month_end=6)
    # glist = Parallel(n_jobs=-1, verbose=10)(delayed(get_grids)(fname, date) for fname, date in zip(ts.ts.flist, ts.ts.dates))
    # zs_list = Parallel(n_jobs=-1, verbose = 10)(delayed(get_zs)(basin_gdf, sbasin_gdf, grid, oRoot,ds, basin, 1) for grid in glist)
    # snotel_new_sbasin, snotel_new_tbasin = zstat2dss(zs_list, basin, ds+'_grid_res', dss_file)
#
#    ds = 'SNOTEL_IDW_5000'
#    files = glob(r'E:\ririe\rasters\SNOTEL_GRIDS\IDW_5000\*.tif')
#    ts = ParseTS(files,'%Y-%m-%d', month_start=9, month_end=6)
#    glist = Parallel(n_jobs=-1, verbose=10)(delayed(get_grids)(fname, date) for fname, date in zip(ts.ts.flist, ts.ts.dates))
#    zs_list = Parallel(n_jobs=-1, verbose = 10)(delayed(get_zs)(basin_gdf, sbasin_gdf, grid, oRoot,ds, basin, 1) for grid in glist)
#    snodas_sbasin, snodas_tbasin = zstat2dss(zs_list, basin, ds+'_grid_res', dss_file)


#    for grid in glist:
#        break
#        a = get_zs(basin_gdf, sbasin_gdf, grid, oRoot,ds, basin, 1e-3)
#


#    ua_sbasin_ann_max = uofa_sbasin.reset_index(level=1).groupby([pd.Grouper(level=0, freq='Y'),'name']).max()
#    ua_tbasin_ann_max = uofa_tbasin.reset_index(level=1).groupby([pd.Grouper(level=0, freq='Y'),'name']).max()
#
#    sdas_sbasin_ann_max = snodas_sbasin.reset_index(level=1).groupby([pd.Grouper(level=0, freq='Y'),'name']).max()
#    sdas_tbasin_ann_max = snodas_tbasin.reset_index(level=1).groupby([pd.Grouper(level=0, freq='Y'),'name']).max()
#
#
#    sdas_sbasin_cum_max = sdas_sbasin_ann_max.groupby(level=0).sum()
#    ua_sbasin_cum_max = ua_sbasin_ann_max.groupby(level=0).sum()
#
#    ua_sbasin_cum_max.sort_values('vol')
#
#
#    def nse(simulation_s, evaluation):
#        nse_ = 1 - (np.sum((evaluation - simulation_s) ** 2, axis=0, dtype=np.float64) /
#                    np.sum((evaluation - np.mean(evaluation)) ** 2, dtype=np.float64))
#
#        return nse_
#
#    import matplotlib.pyplot as plt
#
#
#    fig,ax = plt.subplots()
#    ax.plot_date(x=ua_tbasin_ann_max.index.get_level_values(0), y=ua_tbasin_ann_max.vol.values, label = 'Total Basin', ls = '-')
#    ax.plot_date(x=ua_sbasin_cum_max.index, y=ua_sbasin_cum_max.vol.values, label = 'Subbasin Summation', ls = '--')
#    ax.set_title('UofA Annual Maximum SWE')
#    ax.legend(bbox_to_anchor=(0.8,-0.1),ncol=2)
#    ax.set_ylabel("SWE Volume [m$^3$]")
#    fig.tight_layout()
#
#    fig,ax = plt.subplots()
#    ax.plot_date(x=sdas_tbasin_ann_max.index.get_level_values(0), y=sdas_tbasin_ann_max.vol.values, label = 'Total Basin', ls = '-')
#    ax.plot_date(x=sdas_sbasin_cum_max.index, y=sdas_sbasin_cum_max.vol.values, label = 'Subbasin Summation', ls = '--')
#    ax.set_title('SNODAS Annual Maximum SWE')
#    ax.set_ylabel("SWE Volume [m$^3$]")
#    ax.legend(bbox_to_anchor=(0.8,-0.1),ncol=2)
#    fig.tight_layout()


#    #checking area
#    zs = zs_list[40]
#    shp1 = gpd.read_file(r"E:\ririe\shp\total_watershed.shp")
#    shp1 = shp1.to_crs(sbasin_gdf.crs)
#    all_shp = gpd.GeoDataFrame(pd.concat([sbasin_gdf, shp1], ignore_index=True))
#    tmp = zonal_stats(all_shp, zs.grid.data, affine =zs.grid.affine, raster_out=True, nodata=-9999, stats=['mean','count','sum'])



    #    zs_list = []
#    for grid in glist[4864:4865]:
#        print(grid.date)
#        zs = Zonal_Stat(basin_gdf, sbasin_gdf, grid, oRoot,ds)
#        zs.get_basin_raster()
#        zs.get_sbasin_stats()
#        zs_list.append(zs)



if __name__ == '__main__':

    #python pygrid2ts.py --basin-shp "E:\ririe\shp\total_watershed_dissolved.shp" --sbasin-shp "E:\ririe\shp\total_watershed.shp" --ds SONDAS --raster-dir "E:\SNODAS" --dss "E:\ritie\output_timeseries.dss" --oroot "E:\ririe\beta" --bname RIRIE --dtfmt %Y%m%d --mconv 1000.0

    parser = ArgumentParser(description='Zonal Statistic Software to Create DSS Time Series',
                            usage = "python pygrid2ts.py --basin-shp {0} --sbasin-shp {1} --ds {2} --raster-dir {3} --dss {4} --oroot {5} --bname RIRIE --dtfmt {6} --mconv {7}".format("total_watershed_dissolved.shp","total_watershed.shp", "SONDAS","snodas_dir","output_timeseries.dss", "output_dir", "%%Y%%m%%d", 1000.0))
    parser.add_argument('--basin-shp', dest='basin_shp', help='Dissolved Basin Shapefile', type=str )
    parser.add_argument('--sbasin-shp', dest='sbasin_shp', help='Sub Basin Shapefile', type=str )
    parser.add_argument('--ds', dest='ds', help='Raster Datasource', type=str )
    parser.add_argument('--raster-dir', dest='fpath', help='Directory for rasters to be processed', type=str )
    parser.add_argument('--dss', dest='dss_file', help='Path to Output DSS file', type=str )
    parser.add_argument('--oroot', dest='oRoot', help='Output Directory for Basin Rasters', type=str )
    parser.add_argument('--bname', dest='basin', help='Basin Name', type=str )
    parser.add_argument('--dtfmt', dest='dtfmt', help='strftime Date Format', type=str )
    parser.add_argument('--mconv', dest='m_conv', help='Converstion Factor for Input Rasters to Meters ', type=float )

    a = parser.parse_args()

    #a = argparse.Namespace(basin='RIRIE', basin_shp='E:\\ririe\\shp\\total_watershed_dissolved.shp', ds='SONDAS', dss_file='E:\\ritie\\output_timeseries.dss', dtfmt='%Y%m%d', fpath='E:\\SNODAS', m_conv=1000.0, oRoot='E:\\ririe\\beta', sbasin_shp='E:\\ririe\\shp\\total_watershed.shp')
    print(a)
    print(a.basin_shp)
    assert os.path.exists(a.basin_shp), "Basin Shapefile {1} not found".format(a.basin_shp)
    print(a.sbasin_shp)
    assert os.path.exists(a.sbasin_shp), "Subbasin Shapefile {1} not found".format(a.sbasin_shp)
    print(a.ds)
    print(a.dss_file)
    if os.path.exists(a.dss_file):
        print('Output DSS file {0} exists'.format(a.dss_file))
    else:
        print("Will Create new DSS file: {0}".format(a.dss_file))
    print(a.oRoot)
    print(a.basin)
    print(a.dtfmt)
    print(a.m_conv)

    main(a)
