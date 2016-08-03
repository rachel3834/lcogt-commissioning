# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 20:09:04 2016

@author: rstreet
"""

import archive_access
from os import path
import glob

class DataCounter:
    def __init__(self):
        self.night_dir = None
        self.nbiases = 0
        self.ndarks = 0
        self.nflats_per_filter = {}
        self.nflats = 0
        self.nscience_raw = 0
        self.nscience_quicklook = 0
        self.nscience_reduced = 0
        self.nraw = 0
    
    def summary(self):
        output = self.night_dir+' N biases='+str(self.nbiases)+' N darks='+str(self.ndarks)+\
                ' N flats='+str(self.nflats)+' N science (raw)='+str(self.nscience_raw)+\
                ' N science (quicklook)='+str(self.nscience_quicklook)+\
                ' N science (reduced)='+str(self.nscience_reduced)
        return output
    def flats_summary(self):
        if len(self.nflats_per_filter) == 0:
            output = self.night_dir+' no flat fields'
        else:
            output = self.night_dir+' flat fields: '
            for f,fcount in self.nflats_per_filter.items():
                output = output + ' N flats('+f+')='+str(fcount)
                return output
        
    def calc_nightly_data_totals(self,night_dir,out_dir):
        """Mothod to calculate the number of frames of different types 
        available within a single nights data directory
        """
        
        raw_path = path.join( night_dir, 'raw' )
        quicklook_path = path.join( night_dir, 'preview' )
        reduced_path = path.join( night_dir, 'processed' )
        
        raw_frames = glob.glob( path.join( raw_path, '*fits.fz' ) )
        quicklook_frames = glob.glob( path.join( quicklook_path, '-e11.fits.fz'))
        reduced_frames = glob.glob( path.join( reduced_path, '-e91.fits.fz'))
        self.nraw = len(raw_frames)
        
        key_list = ['FILTER']
        params = {'out_dir':out_dir}
        
        for frame in raw_frames:
            if '-b00' in frame:
                self.nbiases = self.nbiases + 1
            elif '-d00' in frame:
                self.ndarks = self.ndarks + 1
            elif '-f00' in frame:
                self.nflats = self.nflats + 1
                keywords= archive_access.get_frame_params(params, frame, key_list)
                if keywords['FILTER'] in self.nflats_per_filter.keys():
                    self.nflats_per_filter = self.nflats_per_filter + 1
                else:
                    self.nflats_per_filter[keywords['FILTER']] = 1
            elif '-e00' in frame:
                self.nscience_raw = self.nscience_raw + 1
        
        for frame in quicklook_frames:
            if '-e11' in frame:
                self.nscience_quicklook = self.nscience_quicklook + 1
        
        for frame in reduced_frames:
            if '-e91' in frame:
                self.nscience_reduced = self.nscience_reduced + 1
    
    def sum_nightly_data_totals(self,night_counter):
        """Method to add the totals from a second counter to the 
        corresponding attributes of the current counter"""
        
        attr_list = [ 'nraw', 'nbiases','ndarks','nflats','nscience_raw',\
                    'nscience_quicklook', 'nscience_reduced' ]
        for attr in attr_list:
            fcount = getattr(night_counter,attr)
            setattr(self,attr,(getattr(self,attr) + fcount))
        
        for f, fcount in night_counter.nflats_per_filter.items():
            if f in self.nflats_per_filter.keys():
                self.nflats_per_filter[f] = self.nflats_per_filter[f] + fcount
            else:
                self.nflats_per_filter[f] = fcount
                
