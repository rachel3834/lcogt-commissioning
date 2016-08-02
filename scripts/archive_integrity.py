# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 20:13:49 2016

@author: rstreet
"""
from os import path
import glob

def chk_archive_integrity(params):
    """Function to verify that all available raw science frames have 
    corresponding reduced output
    """
    
    params = parse_args()

    for night_dir in params['dir_list']:
        night = night_dir.split('/')[-2]
        raw_dir = path.join( night_dir, 'raw' )
        red_dir= path.join( night_dir, 'banzai' )
        raw_frames = glob.glob( raw_dir, '*.fits.fz' )
        red_frames = glob.glob( red_dir, '*.fits.fz' )
        
        if len(raw_frames) != len(red_frames):
            print('Warning: '+str(night)+' has '+str(len(raw_frames))+\
                    ' raw frames but '+str(len(red_frames))+' reduced frames')
            
    
def parse_args():
    """Function to harvest the parameters required for the data summary code
    from the commandline or prompts"""
    
    params = {}
    if len(argv) != 3:
        params['top_data_dir'] = raw_input('Please enter the path to the data directory: ')        
        params['date_search_string'] = raw_input('Please enter the date search string [e.g. yyyymmd?]: ')
    else:
        params['top_data_dir'] = argv[1]
        params['date_search_string'] = argv[2]
    
    params['dir_list'] = glob.glob(path.join(params['top_data_dir'], \
                                params['date_search_string']))
                                
    return params

if __name__ == '__main__':
    chk_archive_integrity()