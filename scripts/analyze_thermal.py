# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 13:30:32 2016

@author: rstreet
"""

import analyze_calibs
import glob
from sys import argv
from os import path
import numpy as np
from matplotlib import dates, pyplot

def analyze_thermal_stability():
    """Function to analyze the thermal stability of a camera"""
    
    params = parse_args_thermal()
    
    current_ts = []
    currents = []
    for night_dir in params['dir_list']:
        print 'Analyzing data from '+ path.basename(night_dir)
        params['data_dir'] = path.join( night_dir, 'raw' )
        frames = analyze_calibs.FrameSet(params)
        frames.make_frame_listings()
        
        if len(frames.darks) > 0:
            (night_ts, night_currents) = frames.measure_dark_current()
            current_ts = current_ts + night_ts
            currents = currents + night_currents
        
        if frames.nframes > 0:
            (temp_ts,ccdatemp, ccdstemp) = frames.get_temperatures()
    
    plot_dark_current(params, current_ts, currents)
    plot_temperature(params,temp_ts,ccdatemp,ccdstemp)
    
def plot_dark_current(params, ts, currents):
    """Function to create a plot of a set of dark current measurements and 
    corresponding timestamps"""
    
    xplot = np.array(ts)
    yplot = np.array(currents)
    
    hfmt = dates.DateFormatter('%Y-%m-%d\n%H:%M')
    
    fig = pyplot.figure(1)
    pyplot.rcParams['font.size'] = 10.0
    ax = pyplot.subplot(111)
    pyplot.subplots_adjust(bottom=0.2)
    pyplot.plot(xplot,yplot,'b.')
    ax.xaxis.set_major_formatter(hfmt)
    pyplot.xticks(rotation=60.0)
    pyplot.xlabel('Date/time [UTC]')
    pyplot.ylabel('Median e-/pixel/s')
    pyplot.title('Dark current as a function of time')
    pyplot.axis([xplot.min(),xplot.max(),yplot.min()*0.98,yplot.max()*1.2])
    plotfile = path.join(params['out_dir'],'darkcurrent.png')
    pyplot.savefig(plotfile)
    pyplot.close(1)

def plot_temperature(params,temp_ts,ccdatemp,ccdstemp):
    """Function to create a plot of a camera's measured temperature as a function
    of time (CCDATEMP) compared with its set-point temperature (CCDSTEMP)"""
    
    xplot = np.array(temp_ts)
    yplot = np.array(ccdatemp)
    lplot = np.array(ccdstemp)
    set_point = np.median(lplot)
    
    fig = pyplot.figure(1)
    pyplot.rcParams['font.size'] = 10.0
    ax = pyplot.subplot(111)
    pyplot.subplots_adjust(bottom=0.15)
    pyplot.plot(xplot,yplot,'b.')
    pyplot.plot(xplot,np.array([ set_point ] * len(xplot)),'r-',label='Setpoint temp')
    ydiff = 0.25
    if (yplot.max()-yplot.min()) > ydiff: ydiff = (yplot.max()-yplot.min())/2.0
    temp_upper_limit = lplot.mean() + ydiff
    temp_lower_limit = lplot.mean() - ydiff
    pyplot.plot(xplot,np.array([ temp_upper_limit ] * len(xplot)),'r-.')
    pyplot.plot(xplot,np.array([ temp_lower_limit ] * len(xplot)),'r-.',label='Warning threshold')
    pyplot.gcf().autofmt_xdate()
    pyplot.xlabel('Date/time [UTC]')
    pyplot.ylabel('CCDATEMP [degC]')
    pyplot.title('Temperature as a function of time')
    pyplot.legend(loc='best',frameon=False)
    plotfile = path.join(params['out_dir'],'ccd_temp.png')
    pyplot.savefig(plotfile)
    pyplot.close(1)
    
    
def parse_args_thermal():
    """Function to harvest the parameters required for thermal properties 
    analysis from the commandline or prompts"""
    
    params = {}
    if len(argv) != 4:
        params['top_data_dir'] = raw_input('Please enter the path to the data directory: ')
        params['out_dir'] = raw_input('Please enter the path to the output directory: ')
        params['date_search_string'] = raw_input('Please enter the date search string [e.g. yyyymmd?]: ')
    else:
        params['top_data_dir'] = argv[1]
        params['out_dir'] = argv[2]
        params['date_search_string'] = argv[3]
    
    params['dir_list'] = glob.glob(path.join(params['top_data_dir'], params['date_search_string']))
    
    return params
    
if __name__ == '__main__':
    analyze_thermal_stability()