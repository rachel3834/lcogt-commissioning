# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 10:07:43 2016

@author: rstreet
"""

from os import path
from sys import argv
import glob

def calc_nightly_data_totals(night_dir):
    """Function to calculate the number of frames of different types available
    within a single nights data directory
    """
    data = {'nbiases': 0, 'ndarks': 0, 'nflats': 0, 'nscience': 0}
    frame_list = glob.glob( path.join( night_dir, '*.fits.fz' ) )
    for frame in frame_list:
        if '-b00' in frame:
            data['nbiases'] = data['nbiases'] + 1
        elif '-d00' in frame:
            data['ndarks'] = data['ndarks'] + 1
        elif '-f00' in frame:
            data['nflats'] = data['nflats'] + 1
        elif '-e00' in frame:
            data['nscience'] = data['nscience'] + 1
    return data
    
def calc_data_totals(params):
    """Function to calculate the total number of frames of different types
    acquired within a given time frame"""

    data = {'nbiases': 0, 'ndarks': 0, 'nflats': 0, 'nscience': 0}
    nframes = 0
    for night_dir in params['dir_list']:
        night_data = calc_nightly_data_totals(night_dir)
        for ftype, fcount in night_data.items():
            data[ftype] = data[ftype] + fcount
            nframes = nframes + fcount

    print('Data holdings: ')
    print('N bias frames: ' + str(data['nbiases']))
    print('N dark frames: ' + str(data['ndarks']))
    print('N flat fields: ' + str(data['nflats']))
    print('N science frames: ' + str(data['nscience']))
    print('\nTotal number of frames: '+str(nframes))

def parse_args_data():
    """Function to harvest the parameters required for the data summary code
    from the commandline or prompts"""
    
    params = {}
    if len(argv) != 3:
        params['top_data_dir'] = raw_input('Please enter the path to the data directory: ')
        params['date_search_string'] = raw_input('Please enter the date search string [e.g. yyyymmd?]: ')
    else:
        params['top_data_dir'] = argv[1]
        params['date_search_string'] = argv[2]
    
    params['dir_list'] = glob.glob(path.join(params['top_data_dir'], params['date_search_string']))
    
    return params

if __name__ == '__main__':
    parse_args_data()