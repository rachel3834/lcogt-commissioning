# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 09:29:28 2016

@author: rstreet
"""

from os import path, remove
from sys import argv, exit
import glob
from astropy.io import fits
from datetime import datetime
import archive_access
import archive_access
import prepraw3d
import numpy as np
from shutil import copy
import noise_analysis

class FrameSet:
    
    def __init__(self, params):
        self.biases = []
        self.darks = []
        self.flats = {}
        self.exposures = []
        self.naxis1 = None
        self.naxis2 = None
        self.master_stats = {}
        self.masterbias = None
        self.masterbias_file = None
        self.masterdark = None
        self.masterdark_file = None
        self.masterflats = {}
        self.masterflat_files = {}
        self.exclude_file = None
        self.exclude_list = []
        
        if path.isdir(params['data_dir']) == True:
            self.data_dir = params['data_dir']
            self.out_dir = params['out_dir']
            self.exclude_file = path.join(self.out_dir,'exclude.frames')
        else:
            self.data_dir = None
            print('ERROR: Cannot find data directory ' + params['data_dir'])
            exit()
            
    def make_frame_listings(self):
        
        # Check to see if there is an exclude list:
        if path.isfile(self.exclude_file) == True:
            lines = open(self.exclude_file,'r').readlines()
            for line in lines:
                self.exclude_list.append(path.basename(line.replace('\n','')))
        
        frames = glob.glob( path.join(self.data_dir, '*.fits.fz') )
        for frame in frames:
            if path.basename(frame) not in self.exclude_list:
            
                # Always skip the first frame of the night, because un-flushed
                # charge means it will be close to saturation
                f_number = int(path.basename(frame).split('-')[3])
                if f_number > 1:
                
                    # Copy frame to working directory and uncompress - a 
                    # step required due to limited access to the data
                    uframe = archive_access.fetch_frame(frame,self.out_dir)
                    
                    hdr = fits.getheader(uframe)
                    if self.naxis1 == None:
                        datasec = int(hdr['TRIMSEC'].split(',')[0].split(':')[-1])
                        self.naxis1 = datasec
                        self.naxis2 = datasec
                    
                    if hdr['OBSTYPE'] == 'BIAS':
                        self.biases.append(uframe)
                    elif hdr['OBSTYPE'] == 'DARK':
                        self.darks.append(uframe)
                    elif hdr['OBSTYPE'] == 'SKYFLAT':
                        if hdr['FILTER'] not in self.flats.keys():
                            self.flats[hdr['FILTER']] = []
                        self.flats[hdr['FILTER']].append(uframe)
                    elif hdr['OBSTYPE'] == 'EXPOSE':
                        self.exposures.append(uframe)
        self.nframes = len(frames)
        
    def make_frame_listings_from_file(self,frames_file):
        
        frames = []
        if '.fits' in frames_file:
            frames = [ frames_file ]
        elif path.isfile( frames_file ) == True:
            frames = open( frames_file,'r').readlines()
        else:
            print('ERROR: Cannot find input frame or file list')
            exit()
        
        for frame in frames:
            
            # Copy frame to working directory and uncompress - a 
            # step required due to limited access to the data
            uframe = archive_access.fetch_frame(frame,self.out_dir)
            
            hdr = fits.getheader(uframe)
            if self.naxis1 == None:
                datasec = int(hdr['TRIMSEC'].split(',')[0].split(':')[-1])
                self.naxis1 = datasec
                self.naxis2 = datasec
            
            if hdr['OBSTYPE'] == 'BIAS':
                self.biases.append(uframe)
            elif hdr['OBSTYPE'] == 'DARK':
                self.darks.append(uframe)
            elif hdr['OBSTYPE'] == 'SKYFLAT':
                if hdr['FILTER'] not in self.flats.keys():
                    self.flats[hdr['FILTER']] = []
                self.flats[hdr['FILTER']].append(uframe)
        self.nframes = len(frames)
    
    def make_master(self,master_type,bandpass=None):

        if master_type == 'BIAS':
            frame_list = self.biases
            file_path = path.join(self.out_dir, 'masterbias.fits')
        elif master_type == 'DARK':
            frame_list = self.darks
            file_path = path.join(self.out_dir, 'masterdark.fits')
        elif master_type == 'FLAT':
            frame_list = self.flats[bandpass]
            print bandpass, frame_list
            file_path = path.join(self.out_dir, 'masterflat_'+bandpass+'.fits')
            
        if len(frame_list) <= 2:
            print('Warning: insufficient ('+str(len(frame_list))+\
                ') frames to build a master'+master_type.lower())
        
        else:
            (image_data, exp_times, master_header) = \
                read_frame_set(frame_list,self.naxis1,self.naxis2)
            
            if len(exp_times) > 0:
                if master_type in ['DARK','FLAT'] and self.masterbias is not None:
                    image_data = subtract_calib(image_data,self.masterbias)
                elif master_type in ['DARK','FLAT'] and self.masterbias is None:
                    print('Warning: No masterbias data available to subtract')
                    
                if master_type == 'FLAT' and self.masterdark is not None:
                    image_data = subtract_calib(image_data,self.masterdark)
                elif master_type == 'FLAT' and self.masterdark is None:
                    print('Warning: No masterdark data available to subtract')
                
                master_data = np.median(image_data,axis=0)
            
                if master_type == 'DARK':
                    master_data = master_data / np.median(exp_times)
            
                if master_type == 'FLAT':
                    flat_median = np.median(master_data[55:4065,100:4050])
                    master_data = master_data / flat_median
                    
                if master_type not in self.master_stats.keys():
                    self.master_stats[master_type] = {}
                if master_type != 'FLAT':
                    bandpass = 'None'
                
                self.master_stats[master_type][bandpass] = {}
                self.master_stats[master_type][bandpass]['mean'] = master_data.mean()
                self.master_stats[master_type][bandpass]['median'] = np.median(master_data)
                self.master_stats[master_type][bandpass]['stddev'] = master_data.std()
                
                output_frame(master_header, master_data, file_path)
    
                if master_type == 'BIAS':
                    self.masterbias = master_data
                    self.masterbias_file = file_path
                    print('Built masterbias')
                elif master_type == 'DARK':
                    self.masterdark = master_data
                    self.masterdark_file = file_path
                    print('Built masterdark')
                elif master_type == 'FLAT':
                    self.masterflats[bandpass] = master_data
                    self.masterflat_files[bandpass] = file_path
                    print('Built masterflat for '+ bandpass+' filter')
            else:
                print('ERROR: Insufficient valid data to build master'+master_type.lower())
                
                
    def fft_frame(self,frame_type,bandpass=None, data=None, file_path=None):
        
        if frame_type == 'MASTERBIAS' and self.masterbias_file != None:
            params = {  'image_data': self.masterbias, \
                        'image_path': self.masterbias_file, \
                        'out_dir': self.out_dir }
                        
        elif frame_type == 'BIAS' and file_path != None:
            params = {  'image_data': data, \
                        'image_path': file_path, \
                        'out_dir': self.out_dir }
            
        elif frame_type == 'MASTERDARK' and self.masterdark_file != None:
            params = {  'image_data': self.masterdark, \
                        'image_path': self.masterdark_file, \
                        'out_dir': self.out_dir }
                        
        elif frame_type == 'MASTERFLAT' and self.masterflat_files[bandpass] != None:
            params = {  'image_data': self.masterflats[bandpass], \
                        'image_path': self.masterflat_files[bandpass], \
                        'out_dir': self.out_dir }
        noise_analysis.plot_quadrant_ffts(params)
    
    def hist_frame(self,frame='MASTERBIAS',data=None,file_path=None,logy=True,\
                    xrange=None):
        if frame == 'MASTERBIAS' and xrange==None:
            params = {  'image_data': self.masterbias, \
                    'image_path': self.masterbias_file, \
                    'out_dir': self.out_dir,
                    'nbins': 22 }
        elif frame == 'MASTERBIAS' and xrange!= None:
            params = {  'image_data': self.masterbias, \
                    'image_path': self.masterbias_file, \
                    'out_dir': self.out_dir,
                    'nbins': 200 }
        else:
            params = {  'image_data': data, \
                    'image_path': file_path, \
                    'out_dir': self.out_dir }
        noise_analysis.plot_quadrant_hist(params,logy=logy,xrange=xrange)
    
    def stats_summary(self):
        
        for master_type in self.master_stats.keys():
            if master_type in ['BIAS', 'DARK']:
                bandpass = 'None'
                output = master_type
                for key,value in self.master_stats[master_type][bandpass].items():
                    output = output + ' ' + bandpass + ' ' + key + '=' + str(value)+'\n'
                print(output)
            else:
                for bandpass in self.master_stats['FLAT'].keys():
                    output = master_type
                    for key,value in self.master_stats[master_type][bandpass].items():
                        output = output + ' ' + bandpass + ' ' + key + '=' + str(value)+'\n'
                    print(output)
    
    def measure_frame_stats(self,image_data,params):
        stats_log = file(path.join(params['out_dir'],params['log_name']),'w')
        stats_log.write('# Frame Quadrant     Mean    Median   Stddev [ADU]\n')
    
        stats = np.zeros([len(params['frame_list']),12])
        for i in range(0,len(params['frame_list']),1):
            self.hist_frame(frame='SINGLE',data=image_data[i],\
                    file_path=params['frame_list'][i])
            self.fft_frame('BIAS', data=image_data[i],\
                    file_path=params['frame_list'][i])
            
            image_params = {'image_data': image_data[i], \
                            'image_path': params['frame_list'][i]}
            image_stats = noise_analysis.quadrant_stats(image_params)
            
            for q in range(0,4,1):
                stats[i,(3*q)] = image_stats[q]['mean']
                stats[i,(3*q+1)] = image_stats[q]['median']
                stats[i,(3*q+2)] = image_stats[q]['std']
                stats_log.write(path.basename(params['frame_list'][i])+' '+\
                        str(q+1)+' '+\
                        str(image_stats[q]['mean'])+'  '+\
                        str(image_stats[q]['median'])+'  '+\
                        str(image_stats[q]['std'])+'\n')
        
        stats_log.write('\n# Average values:\n')
        stats_log.write('# Quadrant     Mean    Median   Mean_stddev std_stddev [ADU]\n')
        for q in range(0,4,1):
            mean = stats[:,(3*q)].mean()
            median = np.median(stats[:,(3*q+1)])
            meanstd = stats[:,(3*q+2)].mean()
            std = stats[:,(3*q+2)].std()
            stats_log.write(str(q+1)+'     '+str(mean)+'  '+str(median)+\
                '  '+str(meanstd)+'  '+str(std)+'\n')
        stats_log.close()

    def measure_dark_current(self):
        
        (image_data, exp_times, master_header) = \
                read_frame_set(self.darks,self.naxis1,self.naxis2)
        
        frame_ts = []
        dark_current = []
        for i in range(0,len(self.darks),1):
            hdr = fits.getheader(self.darks[i])
            frame_ts.append( datetime.strptime(hdr['DATE-OBS'],"%Y-%m-%dT%H:%M:%S.%f") )
            dc = np.median(image_data[i, 70:4050, 65:4025] )
            dc = dc / exp_times[i]
            dark_current.append( dc )
        
        return frame_ts,dark_current
    
    def get_temperatures(self):
        ccdatemp = [] 
        ccdstemp = []
        frame_ts = []
        frames = self.biases + self.darks + self.exposures
        for bandpass, flats in self.flats.items():
            frames = frames + flats
            
        for frame in frames:
            hdr = fits.getheader(frame)
            frame_ts.append( datetime.strptime(hdr['DATE-OBS'],"%Y-%m-%dT%H:%M:%S.%f") )
            ccdatemp.append( float(hdr['CCDATEMP']) )
            ccdstemp.append( float(hdr['CCDSTEMP']) )
        return frame_ts,ccdatemp, ccdstemp
    
def analyze_night_calibs():
    """Driver function to analyze the calibration frames taken in a single
    night from a single instrument.  The expected data format is raw Sinistro
    frames in fpack-compressed FITS.
    
    The function will build a set of master calibration frames (bias, dark, 
    and flat fields for all available filters) wherever possible. 
    """
    
    params = parse_args_night()
    
    frames = FrameSet(params)
    frames.make_frame_listings()
    
    frames.make_master('BIAS')
    frames.fft_frame('MASTERBIAS')
    frames.hist_frame(frame='MASTERBIAS')
    frames.hist_frame(frame='MASTERBIAS',xrange=[-10.0,10.0])
    frames.make_master('DARK')
    frames.fft_frame('MASTERDARK')
    for bandpass in frames.flats.keys():
        frames.make_master('FLAT',bandpass)

    frames.stats_summary()

def analyze_bias_frames():
    """Function to perform a statistical analyze of a single bias frame, 
    including plotting a histogram of the pixel values per quadrant and
    taking an FFT"""
    
    params = parse_args_biases()
    frames = FrameSet(params)
    frames.make_frame_listings_from_file(params['frames_file'])
        
    (image_data, exp_times, master_header) = \
                read_frame_set(frames.biases,frames.naxis1,frames.naxis2)
    
    params['frame_list'] = frames.biases
    params['log_name'] = 'biases_statistics.data'
    frames.measure_frame_stats(image_data,params)
    
    
def parse_args_night():
    """Parse the commandline arguments and harvest the data directory 
    location.  Prompt the user if none are given."""
    params = {}
    if len(argv) != 3:
        params['data_dir'] = raw_input('Please enter the path to the data directory: ')
        params['out_dir'] = raw_input('Please enter the path to the output directory: ')
    else:
        params['data_dir'] = argv[2]
        params['out_dir'] = argv[3]
        
    return params
    
def parse_args_biases():
    """Parse the commandline arguments and harvest the data directory 
    location.  Prompt the user if none are given."""
    params = {}
    if len(argv) != 5:
        params['data_dir'] = raw_input('Please enter the path to the data directory: ')
        params['out_dir'] = raw_input('Please enter the path to the output directory: ')
        params['frames_file'] = raw_input('Please enter the path to the bias frame or list of frames: ')
    else:
        params['data_dir'] = argv[2]
        params['out_dir'] = argv[3]
        params['frames_file'] = argv[4]
    
    return params

def read_frame_set(frame_list,naxis1,naxis2):
    """Function to read into memory the image data from a set of frames.
    Returns a 3D numpy array of 2D images."""
    
    exp_times = []
    master_header = None
    max_frames = 50
    nframes = min(max_frames,len(frame_list))
    if nframes < len(frame_list):
        print('Warning: number of frames ('+str(len(frame_list))+\
            ') exceeds maximum (' + str(max_frames)+'), dataset will be limited')
    image_data = np.zeros([nframes,naxis2,naxis1])
    
    for i in range(0,nframes,1):
        (imstats, image, hdr) = prepraw3d.prepraw3d(frame_list[i])
        
        # Check for zero-array frames:
        izeros = (image == 0)
        if len(izeros.flatten()) > (0.1*len(image.flatten())):
            print('Warning: High number of zero-array entries in '+frame_list[i]+\
                    ' - frame skipped')
        else:
            image_data[i,:,:] = image
            if master_header == None:
                master_header = hdr
            exp_times.append( float(hdr['EXPTIME']) )
        
        print('Read in '+str(len(exp_times))+' frame(s)')
        
    return image_data, exp_times, master_header


def output_frame(header, data, file_path):
    """Function to output a new frame"""
    
    hdu = fits.PrimaryHDU(data,header=header)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(file_path, clobber=True)
    
def subtract_calib(image_data,master_frame):
    """Function to subtract a pre-existing master calibration frame from 
    a cube of image data."""
    
    for i in range(0,image_data.shape[0],1):
        image_data[i,:,:] = image_data[i,:,:] - master_frame
    return image_data

def command_menu():
    """Function to handle the menu of available options and harvesting of
    commandline arguments"""
    
    main_menu = """
    Analyze a set of calibration frames from one night.........N
    Analyze a one or more bias frames..........................B
    
    Exit.......................................................X
    """
    
    if len(argv) == 1:
        print main_menu
        opt = raw_input('Option: ')
    else:
        opt = argv[1]
    opt = str(opt).upper()
    
    if opt == 'N':
        analyze_night_calibs()
    elif opt == 'B':
        analyze_bias_frames()
        
        
if __name__ == '__main__':
    
    command_menu()
    