# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 11:12:18 2016

@author: rstreet
"""
import argparse
import logging
import numpy as np
from os import path
import matplotlib.pyplot as plt
from scipy.signal import medfilt

from common.Image import Image
log = logging.getLogger(__name__)


def image_fft(image):
    """Function to compute an FFT of an image data array."""
    

    fft_image_0 = abs(np.fft.rfft(image,axis=0)).mean(axis=1)
    fft_image_1 = abs(np.fft.rfft(image,axis=1)).mean(axis=0)
    gaussian_noise = np.random.normal(0, image[5:-5,5:-5].std(), size=image.shape)
    fft_noise = abs(np.fft.rfft(gaussian_noise,axis=0)).mean(axis=1)
    return fft_image_0, fft_image_1, fft_noise

def findpeaks (data):
    basedata = data - medfilt (data, 11)
    basedata[0] = 0
    peakidx= np.where (basedata > np.std (basedata[2:]) * 4)[0]

    goodpeaks = []
    log.debug("{} promising peaks {}".format (len(peakidx), str(peakidx)))

    for ii in range (len(peakidx)):
        if peakidx[ii]>5:
            goodpeaks.append (peakidx[ii])
    return goodpeaks


def plot_image_fft(fig,image):
    """Function to create the plot of an FFT, given a data array of an 
    image region"""
    
    (fft_image_0, fft_image_1, fft_noise) = image_fft(image)
    log_fft_image_0 = np.log10(fft_image_0)
    log_fft_image_1 = np.log10(fft_image_1)
    log_fft_noise = np.log10(fft_noise)
    

    plt.plot(log_fft_noise,color='grey', label="gauss equivalent")
    plt.plot(log_fft_image_1, label="FFT X")

    indexes = findpeaks(log_fft_image_1)
    values = log_fft_image_1[indexes]
    print (indexes, values)
    if len (values) > 0:
        plt.plot (indexes, values, '.')


    plt.plot(log_fft_image_0, label="FFT Y")

    plt.xlabel('Cycles per line')
    plt.ylabel('log10(FFT)')
    plt.legend()
    (xmin,xmax,ymin,ymax) = plt.axis()
    ymax = ymax * 1.05
    if xmax > len(fft_image_0):
        xmax = len(fft_image_0)
    plt.axis([xmin,xmax,ymin,ymax])

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

    image = Image(args.fitsfile[0], gaincorrect=False, overscancorrect=False, skycorrect=False)
    log.info ("Working on image {}".format (args.fitsfile))

    plotord = [ 1, 3, 4, 2 ]
    fig = plt.figure(2)
    for q,qid in enumerate(plotord):

        quad_image = (image.data[q] - np.mean (image.data[q]))
        ax = plt.subplot(2,2,q+1)
        plt.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9,\
            wspace=0.4,hspace=0.4)
            
        fig = plot_image_fft(fig,quad_image)
        plt.title('Quadrant '+str(q+1))
    plotname = path.join( "fft.png")
    plt.savefig(plotname, dpi=200)
    plt.close(2)

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

    parser = argparse.ArgumentParser(description='Noise analysis of image.')

    parser.add_argument('fitsfile', type=str, nargs=1,
                    help='Fits files for cross talk measurement')

    parser.add_argument('--loglevel', dest='loglevel', default='INFO', choices=['DEBUG', 'INFO'],
                        help='Set the debug level')
    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.loglevel.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')
    return args


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
    args = parse_args()
    plot_quadrant_ffts(args)
