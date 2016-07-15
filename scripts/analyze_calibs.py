# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 09:29:28 2016

@author: rstreet
"""

from os import path, remove
from sys import argv, exit
import glob
from astropy.io import fits
import compression_handler
import prepraw3d
import numpy as np
from shutil import copy

class CalibFrameSet:
    
    def __init__(self, params):
        self.biases = []
        self.darks = []
        self.flats = {}
        self.naxis1 = None
        self.naxis2 = None
        self.master_stats = {}
        self.masterbias = None
        self.masterdark = None
        self.masterflats = {}
        
        if path.isdir(params['data_dir']) == True:
            self.data_dir = params['data_dir']
            self.out_dir = params['out_dir']
            self.make_frame_listings()
        else:
            self.data_dir = None
            print('ERROR: Cannot find data directory ' + params['data_dir'])
            exit()
            
    def make_frame_listings(self):
        
        frames = glob.glob( path.join(self.data_dir, '*.fits.fz') )
        
        for frame in frames:
            
            # Copy frame to working directory and uncompress - a 
            # step required due to limited access to the data
            wframe = path.join(self.out_dir, path.basename(frame))
            if path.isfile(wframe) == False:            
                copy(frame, wframe)
            uframe = wframe.replace('.fz','')
            if path.isfile(uframe) == False:
                compression_handler.uncompress( wframe )
            remove(wframe)
            
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
    
    
    def make_master(self,master_type,bandpass=None):

        if master_type == 'BIAS':
            frame_list = self.biases
            file_path = path.join(self.out_dir, 'masterbias.fits')
        elif master_type == 'DARK':
            frame_list = self.darks
            file_path = path.join(self.out_dir, 'masterdark.fits')
        elif master_type == 'FLAT':
            frame_list = self.flats[bandpass]
            file_path = path.join(self.out_dir, 'masterflat_'+bandpass+'.fits')
            
        if len(frame_list) <= 2:
            print('Warning: insufficient ('+str(len(frame_list))+\
                ') frames to build a master'+master_type.lower())
        
        else:
            (image_data, exp_times, master_header) = \
                read_frame_set(frame_list,self.naxis1,self.naxis2)
            
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
                print('Built masterbias')
            elif master_type == 'DARK':
                self.masterdark = master_data
                print('Built masterdark')
            elif master_type == 'FLAT':
                self.masterflats[bandpass] = master_data
                print('Built masterflat for '+ bandpass+' filter')
    
    def stats_summary(self):
        
        for master_type in [ 'BIAS', 'DARK', 'FLAT' ]:
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
                    
def analyze_night_calibs():
    """Driver function to analyze the calibration frames taken in a single
    night from a single instrument.  The expected data format is raw Sinistro
    frames in fpack-compressed FITS.
    
    The function will build a set of master calibration frames (bias, dark, 
    and flat fields for all available filters) wherever possible. 
    """
    
    params = parse_args()
    
    frame_set = CalibFrameSet(params)
    
    frame_set.make_master('BIAS')
    frame_set.make_master('DARK')
    for bandpass in frame_set.flats.keys():
        frame_set.make_master('FLAT',bandpass)

    frame_set.stats_summary()

def parse_args():
    """Parse the commandline arguments and harvest the data directory 
    location.  Prompt the user if none are given."""
    params = {}
    if len(argv) == 1:
        params['data_dir'] = raw_input('Please enter the path to the data directory: ')
        params['out_dir'] = raw_input('Please enter the path to the output directory: ')
    else:
        params['data_dir'] = argv[1]
        params['out_dir'] = argv[2]
        
    return params

def read_frame_set(frame_list,naxis1,naxis2):
    """Function to read into memory the image data from a set of frames.
    Returns a 3D numpy array."""
    
    image_data = np.zeros([len(frame_list),naxis2,naxis1])
    exp_times = []
    master_header = None
    for i in range(0,len(frame_list),1):
        (imstats, image, hdr) = prepraw3d.prepraw3d(frame_list[i])
        image_data[i,:,:] = image
        if master_header == None:
            master_header = hdr
        exp_times.append( float(hdr['EXPTIME']) )
    
    return image_data, exp_times, master_header

def output_frame(header, data, file_path):
    """Function to output a new frame"""
    
    hdu = fits.PrimaryHDU(data,header=header)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(file_path)
    
def subtract_calib(image_data,master_frame):
    """Function to subtract a pre-existing master calibration frame from 
    a cube of image data."""
    
    for i in range(0,image_data.shape[0],1):
        image_data[i,:,:] = image_data[i,:,:] - master_frame
    return image_data
    
    
if __name__ == '__main__':
    analyze_night_calibs()
    