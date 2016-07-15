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


def plot_image_fft(fig,image):
    """Function to create the plot of an FFT, given a data array of an 
    image region"""
    
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

    return fig
    
def plot_whole_frame_fft(params):
    """Function to plot the FFT of a whole frame"""
    
    image = fits.getdata(params['image_path'])
    fig = pyplot.figure(1)
    fig = plot_image_fft(fig,image)
    plotname = path.join( params['out_dir'], \
         path.splitext(path.basename(params['image_path']))[0]+'_fft_log.png' )
    pyplot.savefig(plotname)
    pyplot.close(1)

def plot_quadrant_ffts(params):
    """Function to plot separate FFTs of each quadrant of a Sinistro frame"""
    
    image = fits.getdata(params['image_path'])

    header = fits.getheader(params['image_path'])    
    NAXIS1 = float(header['NAXIS1'])
    NAXIS2 = float(header['NAXIS2'])
    
    params['regions'] = {
                        0: [1.0,(NAXIS1/2.0),1.0,(NAXIS2/2.0)],
                        1: [(NAXIS1/2.0),NAXIS1,1.0,(NAXIS2/2.0)],
                        2: [(NAXIS1/2.0),NAXIS1,(NAXIS2/2.0),NAXIS2],
                        3: [1.0,(NAXIS1/2.0),(NAXIS2/2.0),NAXIS2]
                        }

    plotfmt = [ 'r', 'b', 'm', 'g' ]
    plotord = [ 2, 3, 1, 4 ]
    fig = pyplot.figure(2)
    pyplot.rcParams['font.size'] = 10.0
    
    for q in range(0,4,1):
        region = params['regions'][q]
        
        quad_image = image[region[0]:region[1],region[2]:region[3]]
        
        ax = pyplot.subplot(2,2,q+1)
        pyplot.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9,\
            wspace=0.4,hspace=0.4)
            
        fig = plot_image_fft(fig,quad_image)
        
    plotname = path.join( params['out_dir'], \
         path.splitext(path.basename(params['image_path']))[0]+'_fft_log.png' )
    pyplot.savefig(plotname)
    pyplot.close(2)


def parse_args():
    """Function to harvest and parse the commandline arguments for noise 
    analysis"""

    params = {}
    if len(argv) == 1:
        params['image_path'] = raw_input('Please enter the path to an image FITS file: ')
        params['out_dir'] = raw_input('Please enter the output directory path: ')
    else:
        params['image_path'] = argv[1]
        params['out_dir'] = argv[2]

    if path.isfile(params['image_path']) == False:
        print('ERROR: Cannot find input image '+params['image_path'])
        exit()
    if path.isdir(params['out_dir']) == False:
        print('ERROR: No such output directory '+params['out_dir'])
        
    return params
    
if __name__ == '__main__':
    params = parse_args()
    plot_quadrant_ffts(params)