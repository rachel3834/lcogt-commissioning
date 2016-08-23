# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 15:10:03 2016

@author: rstreet
"""

import numpy as np
from astropy.io import fits
from sys import argv

# Bad Pixel Mask Data
# Remember, coordinates must be given in Python array indices, i.e. numbered
# from zero rather than in pixel coordinates.
bpm_data = {
            'fl11': [ (1201, 1202, 0, 2047) ]
            }

class BadPixelMask:
    def __init__(self):
        self.camera = None
        self.mask = np.zeros(1)
        
    def make_bpm(self, camera_id, frame_size):
        self.camera = camera_id
        self.mask = np.zeros(frame_size)
        
        if self.camera in bpm_data:
            
            for c in bpm_data[self.camera]:
                self.mask[c[2]:c[3],c[0]:c[1]] = 1.0
        
        else:
            print('Warning: No BPM available for '+self.camera)
        
    def output_bpm(self):
        hdu = fits.PrimaryHDU(self.mask,header=fits.Header())
        hdulist = fits.HDUList([hdu])
        file_path = 'bad_pixel_mask_'+self.camera+'.fits'
        hdulist.writeto(file_path, clobber=True)

if __name__ == '__main__':
    if len(argv) != 4:
        camera_id = raw_input('Please enter a camera ID [e.g. fl01]: ')
        frame_size_x = raw_input('Please enter the size of the final (munged) frame in X dimension: ')
        frame_size_y = raw_input('Please enter the size of the final (munged) frame in Y dimension: ')
    else:
        camera_id = argv[1]
        frame_size_x = argv[2]
        frame_size_y = argv[3]
        
    camera_id = str(camera_id).lower()
    frame_size = ( int(frame_size_x), int(frame_size_y) )
    
    bpm = BadPixelMask()
    bpm.make_bpm(camera_id,frame_size)
    bpm.output_bpm()
    