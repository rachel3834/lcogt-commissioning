# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 14:08:21 2016

@author: rstreet
"""

class SinistroFrame:
    
    def __init__(self, image_file):
        self.image_file = image_file
        self.raw_quad_dims = ( 2080, 2048 )
        self.bias_sec = ( 1, 2048, 2049, 2080 )
        # Pixel range: xmin, xmax, ymin, ymax
        self.image_sec = ( 1, 2048, 1, 2048 )
        self.NAXIS1 = (self.image_sec[3] - self.image_sec[2] + 1)*2
        self.NAXIS2 = (self.image_sec[1] - self.image_sec[0] + 1)*2

        # Transformations to switch the quadrants so that the image sections line up
        self.transformations =  {
                                1: [],
                                2: ['fliplr'],
                                3: ['fliplr','flipud'],
                                4: ['flipud']
                                }
 
        self.quad_dimensions = {
                                1: [ImageSec[0]-1,ImageSec[1],ImageSec[2]-1,ImageSec[3]],
                                2: [ImageSec[1],NAXIS1,ImageSec[2]-1,ImageSec[3]],
                                3: [ImageSec[1],NAXIS1,ImageSec[3],NAXIS2],
                                4: [ImageSec[0]-1,ImageSec[1],ImageSec[3],NAXIS2]
                                }
    
    def read_3d_data(self):
        
        if path.isfile(self.image_file) == False:
            print 'Error: Cannot find frame '+self.image_file
            exit()

        imageobj = fits.open(self.image_file)
        self.header = imageobj[0].header
        data_cube = imageobj[0].data
        self.data_cube = data_cube.astype('float')
        imageobj.close()
