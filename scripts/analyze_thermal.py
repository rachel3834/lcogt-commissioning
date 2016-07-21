# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 13:30:32 2016

@author: rstreet
"""

import analyze_calibs
import glob
from os import argv

def analyze_thermal_stability():
    """Function to analyze the thermal stability of a camera"""
    
    params = parse_args_thermal()
    
    for night_dir in params['dir_list']:
        params['data_dir'] = night_dir
        frames = analyze_calibs.FrameSet(params)
        frames.make_frame_listings()
    
        (ts, currents) = frames.measure_dark_current()
        print ts
        print currents
        
def parse_args_thermal():
    """Function to harvest the parameters required for thermal properties 
    analysis from the commandline or prompts"""
    
    params = {}
    if len(argv) != 4:
        params['top_data_dir'] = raw_input('Please enter the path to the data directory: ')
        params['out_dir'] = raw_input('Please enter the path to the output directory: ')
        params['date_search_string'] = raw_input('Please enter the date search string [e.g. yyyymmd?]: ')
    else:
        params['top_data_dir'] = argv[2]
        params['out_dir'] = argv[3]
        params['date_search_string'] = argv[4]
    
    params['dir_list'] = glob.glob(path.join(params['data_dir'], params['date_search_string']))
    
    return params
    
if __name__ == '__main__':
    analyze_thermal_stability()