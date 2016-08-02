# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 20:09:04 2016

@author: rstreet
"""

class DataCounter:
    def __init__(self):
        self.nbiases = 0
        self.ndarks = 0
        self.nflats = 0
        self.nscience = 0
    
    def summary(self):
        output = 'N biases='+str(self.nbiases)+' N darks='+str(self.ndarks)+\
                ' N flats='+str(self.nflats)+' N science='+str(self.nscience)
        return output
