# -*- coding: utf-8 -*-

import datetime as dt
import pandas as pd
import os
import numpy as np
from .utils import date_util



class TimeSeries(object):
    def __init__(self, files,dtfmt,month_start=None, month_end=None):
        self.start_date = None
        self.end_date = None
        self.flist = files
        self.f_df = None
        self.dates = None
        self.month_start = month_start
        self.month_end = month_end
        self.dtfmt = dtfmt
        
    @staticmethod
    def test(files):
        a = [0 for a  in files if os.path.isfile(a)]
        if len(a) == len(files): 
            return True
        else:
            return [i for i in files if not os.path.isfile(i)]
     

    
    def import_ts(self):
        names = ['-'.join(os.path.split(i)[1].split('.')[0].split('_')) for i in self.flist]
        dates =[date_util(i,self.dtfmt) for i in names]
        
        #no subseting of data
        if not self.month_start and not self.month_end:
            print('No grids excluded')
            pass
        
        #subset the data using month_start and month_end
        else:
            b =[i.strftime("%m") for i in dates]
            b = np.array(b).astype(int)
            if self.month_end < self.month_start:
                idx = np.argwhere((b<self.month_end) | (b>self.month_start))
                d_a = np.array(dates)
                dates = d_a[idx].ravel().tolist()
                self.flist = np.array(self.flist)[idx].ravel().tolist()    
            else:
                idx = np.argwhere((b<self.month_end) | (b<self.month_start))
                d_a = np.array(dates)
                dates = d_a[idx].ravel().tolist()
                self.flist = np.array(self.flist)[idx].ravel().tolist()       
            
        self.dates = dates
        self.start_date = min(dates)
        self.end_date = max(dates)
        
    def get_ts_df(self):
        return pd.DataFrame({'fname':self.flist, 'dates':self.dates})
        
class ParseTS(object):
    def __init__(self, files, dtfmt,month_start=None,month_end=None):
        ts = TimeSeries(files,dtfmt, month_start, month_end)
        ts.import_ts()
        self.ts = ts
    
        
        
        