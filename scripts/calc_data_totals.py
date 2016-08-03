# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 10:07:43 2016

@author: rstreet
"""

from os import path
from sys import argv
import glob
import data_counter

def calc_data_totals():
    """Function to calculate the total number of frames of different types
    acquired within a given time frame"""
    
    params = parse_args_data()

    data = data_counter.DataCounter()
    data.night_dir = params['date_search_string']
    nframes = 0
    frame_types = [ 'nbiases','ndarks','nflats','nscience_raw',\
                    'nscience_quicklook', 'nscience_reduced' ]

    dir_sort = params['dir_list']
    dir_sort.sort()
    
    for night_dir in dir_sort:
        night_counter = data_counter.DataCounter()
        night_counter.night_dir = path.basename(night_dir)
        night_counter.calc_nightly_data_totals(night_dir,params['out_dir'])
        data.sum_nightly_data_totals(night_counter)
        print night+' '+night_data.summary()
        
    print('\nData holdings: ')
    print data.summary()
    print('\nTotal number of frames: '+str(nframes))

def parse_args_data():
    """Function to harvest the parameters required for the data summary code
    from the commandline or prompts"""
    
    params = {}
    if len(argv) != 3:
        params['top_data_dir'] = raw_input('Please enter the path to the data directory: ')
                params['out_dir'] = raw_input('Please enter the output directory path: ')
        params['date_search_string'] = raw_input('Please enter the date search string [e.g. yyyymmd?]: ')
    else:
        params['top_data_dir'] = argv[1]
        params['out_dir'] = argv[2]
        params['date_search_string'] = argv[3]
    
    params['dir_list'] = glob.glob(path.join(params['top_data_dir'], \
                                params['date_search_string']))
        
    return params

if __name__ == '__main__':
    calc_data_totals()