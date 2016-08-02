# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 20:35:12 2016

@author: rstreet
"""

from os import path
import glob
import archive_access
from astropy.io import fits
from shutil import remove

def chk_archive_for_target():
    """Function to search the archive for frames on a specific object"""
    
    params = parse_args()
    key_list = [ 'OBJECT' ]
    
    log = open( path.join(params['out_dir'], 'data_listing.log'), 'r' )
    total_frames = 0
    for night_dir in params['dir_list']:
        night = night_dir.split('/')[-2]
        frames = glob.glob( path.join(night_dir,'raw') )
        nframes = 0
        for frame in frames:
            frame_pars = get_frame_params(params, frame, keywords)
            if frame_pars['OBJECT'] == params['target']:
                nframes = nframes + 1
                red_frame = path.join(night_dir,'banzai')
                if path.isfile(red_frame) == False:
                    log.write('Warning: No reduced frame for '+\
                            path.basename(frame)+'\n')
        log.write(str(night)+' '+str(nframes)+' for target\n')
        total_frames = total_frames + nframes
    log.write('\nTotal frames for '+params['target']+': '+str(total_frames)+'\n')
    log.close()
    

def get_frame_params(params, frame, key_list):
    """Function to object the requested header keyword parameters from the
    FITS header of a frame"""
    
    keywords = {}
    uframe = archive_access.fetch_frame(frame,params['out_dir'])
    hdr = fits.getheader(uframe)
    for key in key_list:
        keywords[key] = hdr[key]
    remove(uframes)
    
    return keywords

def parse_args():
    """Function to harvest the parameters required for the data summary code
    from the commandline or prompts"""
    
    params = {}
    if len(argv) != 3:
        params['top_data_dir'] = raw_input('Please enter the path to the data directory: ')
        params['out_dir'] = raw_input('Please enter the output directory path: ')        
        params['date_search_string'] = raw_input('Please enter the date search string [e.g. yyyymmd?]: ')
        params['target'] = raw_input('Please enter the targets OBJECT keyword: ')
    else:
        params['top_data_dir'] = argv[1]
        params['out_dir'] = argv[2]
        params['date_search_string'] = argv[3]
        params['target'] = argv[4]
    
    params['dir_list'] = glob.glob(path.join(params['top_data_dir'], \
                                params['date_search_string']))
                                
    return params

if __name__ == '__main__':
    chk_archive_for_target()