# -*- coding: utf-8 -*-
#author          :Daniel Hamill
#email           :Daniel.D.Hamill@usace.army.mil

from rasterio import crs
from datetime import datetime
import re
import os
import numpy as np
import pandas as pd
from pydsstools.heclib.dss import HecDss
from pydsstools.core import TimeSeriesContainer


def test_func():
    return 'test_return'

def get_albers():
    return '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'

def check_crs(gdf):

    '''
    Simple function to compare projections
    '''

    s2 = crs.CRS.from_dict(gdf.crs)
    s1 = crs.CRS.from_proj4(get_albers())

    if s1 == s2:
        pass
    else:
        gdf = gdf.to_crs(s1.to_dict())
    return gdf


def date_util(name,dtfmt):
    '''
    Utility to find dates from complex strings:
    Requires self.dtfmt to be correcly specified
    '''
    parts = dtfmt.split('%')[1:]

    
    fname = ''.join(os.path.split(name)[1].split('.')[0].split('_'))
    str2regex = {'Y':'[1-2][0-9][0-9][0-9]',
        'm':'[0-1][0-9]',
        'd':'[0-3][0-9]',
        'Y-':'[1-2][0-9][0-9][0-9]-',
        'm-':'[0-1][0-9]-',
        'd-':'[0-3][0-9]-'}

    regex_string = ''.join([str2regex[i] for i in parts])
    date_str = re.search(regex_string,fname)
    return datetime.strptime(date_str.group(), dtfmt)


def zstat2dss(zs_list, basin, ds, dss_file):
    dates = [i.grid.date for i in zs_list]
    sbasin_avg = [np.round(i.sbasin_avg,4) for i in zs_list]
    sbasin_vol = [i.sbasin_vol for i in zs_list]
    basin_avg = [np.round(i.basin_avg,4) for i in zs_list]
    basin_vol = [i.basin_vol for i in zs_list]
    names = [i.sbasin_names for i in zs_list]


    date_rav = np.ravel(np.repeat(dates,len(zs_list[0].sbasin_names)))
    sbasin_avg_rav = np.ravel(sbasin_avg)
    sbasin_vol_rav = np.ravel(sbasin_vol)
    names_rav = np.ravel(names)


    idx = pd.MultiIndex.from_tuples(zip(date_rav, names_rav), names=['date','name'])
    sbasin = pd.DataFrame(index=idx,data={'mean_swe':sbasin_avg_rav, 'vol':sbasin_vol_rav} )
    sbasin = sbasin.sort_index(level=0)

    idx = pd.MultiIndex.from_tuples(zip(dates, np.ravel(np.repeat('Total_Basin',len(dates)))), names=['date','name'])
    tbasin = pd.DataFrame(index=idx,data={'mean_swe':basin_avg, 'vol':basin_vol} )
    tbasin = tbasin.sort_index(level=0)

    idx = pd.date_range(sbasin.index.get_level_values(0).min(), sbasin.index.get_level_values(0).max())
    
    
    fid = HecDss.Open(dss_file,version=6)
    fid.close()

    for name, group in sbasin.groupby(level=1):
        #group.loc[:, 'wy'] = np.where(group.index.get_level_values(0).month>9,group.index.get_level_values(0).year+1,group.index.get_level_values(0).year)
        group.index = group.index.droplevel(1)
        group.index = group.index.sort_values()

        group=group.reindex(idx, fill_value=0)

        
        start_date =group.index.min().strftime('%d%b%Y %H:%M:%S')
        pname = '/{0}/{1}/AVG_SWE//1DAY/{2}/'.format(basin,name.upper().replace(' ', '_'), ds)

        print(pname)
        tsc = TimeSeriesContainer()
        tsc.granularity = 60 #seconds i.e. minute granularity
        tsc.numberValues = group.mean_swe.size
        tsc.startDateTime=start_date
        tsc.pathname = pname
        tsc.units = "M"
        tsc.type = "INST-VAL"
        tsc.interval = 1
        #must a +ve integer for regular time-series
        #actual interval implied from E part of pathname
        tsc.values =group.mean_swe.values
        #values may be list,array, numpy array

        fid = HecDss.Open(dss_file)
        fid.deletePathname(tsc.pathname)
        status = fid.put(tsc)
        fid.close()

        pname = '/{0}/{1}/VOL//1DAY/{2}/'.format(basin,name.upper().replace(' ', '_'), ds)
        print(pname)

        tsc = TimeSeriesContainer()
        tsc.granularity = 60 #seconds i.e. minute granularity
        tsc.numberValues = group.index.size
        tsc.startDateTime=start_date
        tsc.pathname = pname
        tsc.units = "CUBIC_METERS"
        tsc.type = "INST-VAL"
        tsc.interval = 1
        #must a +ve integer for regular time-series
        #actual interval implied from E part of pathname
        tsc.values =group.vol.values
        #values may be list,array, numpy array

        fid = HecDss.Open(dss_file)
        fid.deletePathname(tsc.pathname)
        status = fid.put(tsc)
        fid.close()
    for name, group in tbasin.groupby(level=1):
        
        #group.loc[:, 'wy'] = np.where(group.index.get_level_values(0).month>9,group.index.get_level_values(0).year+1,group.index.get_level_values(0).year)

        group.index = group.index.droplevel(1)
        group.index = group.index.sort_values()

        group=group.reindex(idx, fill_value=0)

        start_date =group.index.min().strftime('%d%b%Y %H:%M:%S')
        pname = '/{0}/{1}/AVG_SWE//1DAY/{2}/'.format(basin,name.upper().replace(' ', '_'), ds)


        print(pname)
        tsc = TimeSeriesContainer()
        tsc.granularity = 60 #seconds i.e. minute granularity
        tsc.numberValues = group.mean_swe.size
        tsc.startDateTime=start_date
        tsc.pathname = pname
        tsc.units = "M"
        tsc.type = "INST-VAL"
        tsc.interval = 1
        #must a +ve integer for regular time-series
        #actual interval implied from E part of pathname
        tsc.values =group.mean_swe.values
        #values may be list,array, numpy array

        fid = HecDss.Open(dss_file)
        fid.deletePathname(tsc.pathname)
        status = fid.put(tsc)
        fid.close()

        pname = '/{0}/{1}/VOL//1DAY/{2}/'.format(basin,name.upper().replace(' ', '_'), ds)
        print(pname)

        tsc = TimeSeriesContainer()
        tsc.granularity = 60 #seconds i.e. minute granularity
        tsc.numberValues = group.index.size
        tsc.startDateTime=start_date
        tsc.pathname = pname
        tsc.units = "CUBIC_METERS"
        tsc.type = "INST-VAL"
        tsc.interval = 1
        #must a +ve integer for regular time-series
        #actual interval implied from E part of pathname
        tsc.values =group.vol.values
        #values may be list,array, numpy array

        fid = HecDss.Open(dss_file)
        fid.deletePathname(tsc.pathname)
        status = fid.put(tsc)
        fid.close()
    return sbasin, tbasin