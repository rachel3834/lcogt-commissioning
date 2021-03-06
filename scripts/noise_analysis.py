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

def image_fft(image,sumaxis=0):
    """Function to compute an FFT of an image data array."""
    
    axis = 1
    if sumaxis==1: axis=0
    
    fft_image = abs(np.fft.rfft(image,axis=axis)).mean(axis=sumaxis)
    
    gaussian_noise = np.random.normal(0, image.std(), size=image.shape)
    
    fft_noise = abs(np.fft.rfft(gaussian_noise,axis=axis)).mean(axis=axis)
    
    return fft_image, fft_noise


def plot_image_fft(fig,image,axis=0):
    """Function to create the plot of an FFT, given a data array of an 
    image region"""
    
    (fft_image, fft_noise) = image_fft(image,sumaxis=axis)
    
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
    
    if 'image_data' in params.keys():
        image = params['image_data']
    else:
        image = fits.getdata(params['image_path'])
    fig = pyplot.figure(1)
    fig = plot_image_fft(fig,image)
    plotname = path.join( params['out_dir'], \
         path.splitext(path.basename(params['image_path']))[0]+'_fft_log.png' )
    pyplot.savefig(plotname)
    pyplot.close(1)

def plot_quadrant_ffts(params):
    """Function to plot separate FFTs of each quadrant of a Sinistro frame"""
    
    (image, params) = get_image_data(params)
    axes = { 0: 'x', 1: 'y' }
    for axis in axes.keys():
        plotfmt = [ 'r', 'b', 'm', 'g' ]
        plotord = [ 2, 3, 1, 4 ]
        fig = pyplot.figure(2)
        pyplot.rcParams['font.size'] = 10.0
        
        for q,qid in enumerate(plotord):
            region = params['regions'][qid-1]
            
            quad_image = image[region[0]:region[1],region[2]:region[3]]
            
            ax = pyplot.subplot(2,2,q+1)
            pyplot.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9,\
                wspace=0.4,hspace=0.4)
            
            fig = plot_image_fft(fig,quad_image,axis=axis)
            pyplot.title('Quadrant '+str(q+1))
        plotname = path.join( params['out_dir'], \
             path.splitext(path.basename(params['image_path']))[0]+'_fft_log_'+axes[axis]+'.png' )
        pyplot.savefig(plotname)
        pyplot.close(2)

def get_image_data(params):
    
    if 'image_data' in params.keys():
        image = params['image_data']
    else:
        image = fits.getdata(params['image_path'])

    header = fits.getheader(params['image_path'])  
    NAXIS1 = float(header['NAXIS1'])
    NAXIS2 = float(header['NAXIS2'])
    
    params['regions'] = {
                        0: [1,int(NAXIS1/2.0),1,int(NAXIS2/2.0)],
                        1: [int(NAXIS1/2.0),int(NAXIS1),1,int(NAXIS2/2.0)],
                        2: [int(NAXIS1/2.0),int(NAXIS1),int(NAXIS2/2.0),int(NAXIS2)],
                        3: [1,int(NAXIS1/2.0),int(NAXIS2/2.0),int(NAXIS2)]
                        }
                        
    return image, params
    
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
        exit()
        
    return params


def plot_quadrant_hist(params,logy=True,xrange=None):
    """Function to plot histograms of the pixel values of an image, separated
    into separate quadrants"""
    
    def get_hist_axis_range(data):
        idx = np.where(data < 15.0)
        (hist_data,bins) = np.histogram(data[idx],bins=200)
        idx = hist_data > 10
        return (bins[idx].min(), bins[idx].max())
        
    (image, params) = get_image_data(params)
    plotfmt = [ 'r', 'b', 'm', 'g' ]
    plotord = [ 2, 3, 1, 4 ]
    fig = pyplot.figure(3)
    pyplot.rcParams['font.size'] = 10.0
    
    for q,qid in enumerate(plotord):
        region = params['regions'][qid-1]
        quad_image = image[region[0]:region[1],region[2]:region[3]]
        if 'nbins' not in params.keys():
            nbins = int(((quad_image.max()-quad_image.min())/7.5)*2.0)
        else:
            nbins = int(params['nbins'])
            
        ax = pyplot.subplot(2,2,q+1)
        pyplot.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9,\
            wspace=0.4,hspace=0.4)
        
        data = quad_image.flatten()
        
        xaxis_range = get_hist_axis_range(data)
        
        idx = np.where(data >= xaxis_range[0])
        data = data[idx]
        idx = np.where(data <= xaxis_range[1])
        data = data[idx]
        pyplot.hist(data,bins=nbins,color='w',\
                range=(data.min(),data.max()),log=logy)
        stats = basic_stats(quad_image)
        
        pyplot.title('Quadrant '+str(qid)+\
                    '\nMedian = '+str(round(stats['median'],3))+'e-, '+\
                    'Mean = '+str(round(stats['mean'],3))+'e-\n'+\
                    'St.D = '+str(round(stats['std'],3))+'e-',\
                    fontsize=10)
        if q in [2,3]: pyplot.xlabel('Pixel value [e-]')
        pyplot.ylabel('Frequency')
        
        (xr1,xr2,yr1,yr2) = pyplot.axis()
        if xrange != None:
            xr1 = xrange[0]
            xr2 = xrange[1]
        dx = (xr2 - xr1)/10.0
        dy = (yr2 - yr1)/10.0
        yr2 = yr2 + 3.0*dy
        pyplot.axis([xr1,xr2,yr1,yr2])
        yval = yr2 - 1.5*dy
        pyplot.plot(np.array([0.0,stats['std']]),np.array([yval]*2),'k-')
    
    if xrange != None:
        plotname = path.join( params['out_dir'], \
         path.splitext(path.basename(params['image_path']))[0]+'_hist_zoom.png' )
    else:
        plotname = path.join( params['out_dir'], \
         path.splitext(path.basename(params['image_path']))[0]+'_hist.png' )
    pyplot.savefig(plotname)
    pyplot.close(3)

def basic_stats(image):
    
    stats = {}
    stats['mean'] = image.mean()
    stats['median'] = np.median(image)
    stats['std'] = image.std()
    return stats

def quadrant_stats(params):
    """Function to calculate basic statistics for each quadrant of a munged, 
    2D image"""
    
    (image, params) = get_image_data(params)
    image_stats = {}
    for q in range(0,4,1):
        region = params['regions'][q]
        quad_image = image[region[0]:region[1],region[2]:region[3]]
        image_stats[q] = basic_stats(quad_image)
    return image_stats
    
if __name__ == '__main__':
    params = parse_args()
    plot_quadrant_ffts(params)
    plot_quadrant_hist(params)