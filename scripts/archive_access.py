# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 15:39:37 2016

@author: rstreet
"""

from os import path, remove
from shutil import copy
import compression_handler

def fetch_frame(frame_path,working_dir):
    """Function to fetch a copy of an input frame from the data archive
    to the working analysis directory, and uncompress the file for analysis
    Returns the uframe image path
    """
    frame_path = frame_path.replace('\n','')
    
    wframe = path.join(working_dir, path.basename(frame_path))
    if path.isfile(wframe) == False:            
        copy(frame_path, wframe)
    uframe = wframe.replace('.fz','')
    if path.isfile(uframe) == False:
        compression_handler.uncompress( wframe )
    remove(wframe)
    
    return uframe
