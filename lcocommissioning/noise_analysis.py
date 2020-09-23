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


def image_fft(image, samplefreq = 0):
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


def plot_image_fft(fig,image, samplerate = 0):
    """Function to create the plot of an FFT, given a data array of an 
    image region"""


    (fft_image_0, fft_image_1, fft_noise) = image_fft(image)
    log_fft_image_0 = np.log10(fft_image_0)
    log_fft_image_1 = np.log10(fft_image_1)
    log_fft_noise = np.log10(fft_noise)

    if samplerate >0:
        freqfactor = samplerate / (2* len(fft_image_0))
    else:
        freqfactor = 1

    freq_0 = np.arange (len(fft_image_0)) * freqfactor

    plt.plot(freq_0, log_fft_noise,color='grey', label="gauss equivalent")
    plt.plot(freq_0, log_fft_image_1, label="FFT X")

    indexes = findpeaks(log_fft_image_1)
    values = log_fft_image_1[indexes]
    peaklist = zip (freq_0[indexes], values)
    for (f,v) in peaklist:
        print (f,v)

    if len (values) > 0:
        plt.plot (freq_0[indexes], values, '.')

    plt.plot(freq_0, log_fft_image_0, label="FFT Y")

    if samplerate >0:
        plt.xlabel('Frequency [Hz]')
    else:
        plt.xlabel('Cycles per line')

    plt.ylabel('log10(FFT)')
    #plt.legend()
    (xmin,xmax,ymin,ymax) = plt.axis()
    ymax = ymax * 1.05
    if xmax > np.max (freq_0):
        xmax = np.max (freq_0)
    plt.axis([xmin,xmax,ymin,ymax])

    return fig
    


def plot_quadrant_ffts(params):
    """Function to plot separate FFTs of each quadrant of a Sinistro frame"""

    image = Image(args.fitsfile[0], gaincorrect=False, overscancorrect=False, skycorrect=False)
    log.info ("Working on image {}".format (args.fitsfile))

    plotord = [ 1, 3, 4, 2 ]
    plotord = [1]
    fig = plt.figure(2)
    for q,qid in enumerate(plotord):

        quad_image = (image.data[q] - np.mean (image.data[q]))
        ax = plt.subplot(2,2,q+1)
        plt.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9,\
            wspace=0.4,hspace=0.4)

        fig = plot_image_fft(fig,quad_image, samplerate=args.pr)
        plt.title('Quadrant '+str(q+1))
    plotname = path.join( "fft.png")
    plt.savefig(plotname, dpi=200)
    plt.close(2)


    
def parse_args():
    """Function to harvest and parse the commandline arguments for noise 
    analysis"""

    parser = argparse.ArgumentParser(description='Noise analysis of image.')
    parser.add_argument ('--pr', type=float, default = 0, help="Pixel rate in kHz")

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
    for q in range(0,1,4):
        region = params['regions'][q]
        quad_image = image[region[0]:region[1],region[2]:region[3]]
        image_stats[q] = basic_stats(quad_image)
    return image_stats





if __name__ == '__main__':
    args = parse_args()
    plot_quadrant_ffts(args)
