# -*- coding: utf-8 -*-

from rasterio import crs
from datetime import datetime
import re
import os

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
    fname = '-'.join(os.path.split(name)[1].split('.')[0].split('_'))
    str2regex = {'Y':'[1-2][0-9][0-9][0-9]',
        'm':'[0-1][0-9]',
        'd':'[0-3][0-9]'}
    
    regex_string = ''.join([str2regex[i] for i in parts])
    date_str = re.search(regex_string,fname)
    return datetime.strptime(date_str.group(), dtfmt)