# -*- coding: utf-8 -*-

import datetime as dt
import pandas as pd
import os
import dateutil.parser as dparse



class TimeSeries(object):
    def __init__(self, files):
        self.start_date = None
        self.end_date = None
        self.flist = files
        self.f_df = None
        self.dates = None
        
    @staticmethod
    def test(files):
        a = [0 for a  in files if os.path.isfile(a)]
        if len(a) == len(files): 
            return True
        else:
            return [i for i in files if not os.path.isfile(i)]
    
    def import_ts(self,dtfmt = '%Y%m%d'):
        a = [dparse.parse(os.path.split(i)[1],fuzzy=True) for i in self.flist]
        self.dates = a
        self.start_date = min(a)
        self.end_date = max(a)
        
    def get_ts_df(self):
        return pd.DataFrame({'fname':self.flist, 'dates':self.dates})
        
        
        
        
        