# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 15:39:37 2016

@author: rstreet
"""

from os import path, remove
from sys import exit, argv
from shutil import copy
import compression_handler

def fetch_frame(frame_path,working_dir):
    """Function to fetch a copy of an input frame from the data archive
    to the working analysis directory, and uncompress the file for analysis
    Returns the uframe image path
    """
    frame_path = frame_path.replace('\n','')
    if path.isdir(frame_path) == False:
        wframe = path.join(working_dir, path.basename(frame_path))
        if path.isfile(frame_path) == True and path.isfile(wframe) == False:            
            copy(frame_path, wframe)
        elif path.isfile(frame_path) == False:
            print('ERROR: Missing archive data '+frame_path)
            exit()
        uframe = wframe.replace('.fz','')
        if path.isfile(uframe) == False:
            compression_handler.uncompress( wframe )
        remove(wframe)
    else:
        print 'Error: missing file name (trailing blank line in input?)'
        exit()
        
    return uframe

if __name__ == '__main__':
    
    if len(argv) == 3:
        frame_path = argv[1]
        working_dir = argv[2]
    else:
        frame_path = raw_input('Please enter the path to the archive frame: ')
        working_dir = raw_input('Please enter the output directory path: ')
    uframe = fetch_frame(frame_path,working_dir)
    