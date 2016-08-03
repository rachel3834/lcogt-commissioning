# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 20:35:12 2016

@author: rstreet
"""

from os import path
from sys import argv
import glob
import archive_access
from astropy.io import fits
import data_counter

def chk_archive_for_target():
    """Function to search the archive for frames on a specific object"""
    
    params = parse_args()
    key_list = [ 'OBJECT' ]
    nframes = 0
    log = open( path.join(params['out_dir'], 'data_listing.log'), 'w' )
    for night_dir in params['dir_list']:
        night_counter = data_counter.DataCounter()
        night_counter.night_dir = path.basename(night_dir)
        night_counter.calc_nightly_object_total(night_dir,params['out_dir'],\
                                                params['object_name'])
        output = night_counter.object_data_summary(params['object_name'])
        log.write(output+'\n')
        
        nframes = nframes + night_counter.nscience_raw
        
    log.write('\nTotal frames for '+params['object_name']+': '+str(nframes)+'\n')
    log.close()
    

def parse_args():
    """Function to harvest the parameters required for the data summary code
    from the commandline or prompts"""
    
    params = {}
    if len(argv) != 5:
        params['top_data_dir'] = raw_input('Please enter the path to the data directory: ')
        params['out_dir'] = raw_input('Please enter the output directory path: ')        
        params['date_search_string'] = raw_input('Please enter the date search string [e.g. yyyymmd?]: ')
        params['object_name'] = raw_input('Please enter the targets OBJECT keyword: ')
    else:
        params['top_data_dir'] = argv[1]
        params['out_dir'] = argv[2]
        params['date_search_string'] = argv[3]
        params['object_name'] = argv[4]
    
    params['dir_list'] = glob.glob(path.join(params['top_data_dir'], \
                                params['date_search_string']))
                                
    return params

if __name__ == '__main__':
    chk_archive_for_target()