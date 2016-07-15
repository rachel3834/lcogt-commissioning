# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 11:12:18 2016

@author: rstreet
"""

import numpy as np
from os import path
from sys import argv, exit
from astropy.io import fits
from matplotlib import pyplot

def image_fft(image):
    """Function to compute an FFT of an image data array."""
    
    fft_image = abs(np.fft.rfft(image)).mean(axis=0)
    
    expected_noise = np.random.normal(0, image.std(), size=image.shape)
    
    fft_noise = abs(np.fft.rfft(expected_noise)).mean(axis=0)
    
    return fft_image, fft_noise


def plot_image_fft():
    
    params = parse_args()
    
    image = fits.getdata(params['image_path'])
    
    (fft_image, fft_noise) = image_fft(image)
    
    log_fft_image = np.log10(fft_image)
    log_fft_noise = np.log10(fft_noise)
    
    pyplot.plot(log_fft_image,'m-')
    pyplot.plot(log_fft_noise,'k-')
    pyplot.xlabel('Cycles per line')
    pyplot.ylabel('log10(FFT)')
    (xmin,xmax,ymin,ymax) = pyplot.axis()
    ymax = ymax * 1.05
    if xmax > len(fft_image):
        xmax = len(fft_image)
    pyplot.axis([xmin,xmax,ymin,ymax])
    plotname = path.join( params['output_path'], \
                        path.splitext(path.basename(params['image_path']))[0]+'_fft_log.png' )
    pyplot.savefig(plotname)
    pyplot.close(1)
    
    
def parse_args():
    """Function to harvest and parse the commandline arguments for noise 
    analysis"""

    params = {}
    if len(argv) == 1:
        params['image_path'] = raw_input('Please enter the path to an image FITS file: ')
        params['output_path'] = raw_input('Please enter the output directory path: ')
    else:
        params['image_path'] = argv[1]
        params['output_path'] = argv[2]

    if path.isfile(params['image_path']) == False:
        print('ERROR: Cannot find input image '+params['image_path'])
        exit()
    if path.isdir(params['output_path']) == False:
        print('ERROR: No such output directory '+params['output_path'])
        
    return params
    
if __name__ == '__main__':
    plot_image_fft()