#!/usr/bin/env python
##############################################################################
#     	      	      	  CROSSTALK
#
# Software to analyse and remove the crosstalk signature
#
# Future improvements:
# - Currently the code subtracts a sigmaClipped mean to remove the sky background
#   (in addition to the overscan) from each frame before determining the crosstalk 
#   parameters.  Owing to the high degree of structure in the flat field, this is
#   not a good approximation, and yet dividing by the flat field before the determination
#   scales the flux.  
# - An alternative approach might be to use limited regions 
# - Is an iterative approach necessary?  Or minimal improvement?
##############################################################################

###################################
# IMPORTED FUNCTIONS
from astropy.io import fits
import numpy as np
from scipy import optimize
from sys import argv, exit
from os import path, mkdir
from shutil import copy, move
import glob
import warnings
from datetime import datetime, timedelta
import statistics
from matplotlib import use as useBackend
useBackend('Agg')
from matplotlib import pyplot
import archive_access

###################################
# DECLARATIONS

status_code = { 0: 'OK', 
  	      -1: 'Error: Cannot find frame ',\
		-2: 'Error: Missing data location',\
		-3: 'Error: More processed frames than raw - corrupted data location?',\
		-4: 'Warning: frame not in 3D data format',\
		-5: 'Warning: frame zero-length (will be ignored)',\
		-6: 'Error: No coefficients for camera and binning combination'
      	      }
verbose = True

############################################
# IMAGE CLASS
class CrossImage:
    
    def __init__(self,image_file):
    
        fitsobj = fits.open(image_file)
        self.header = fitsobj[0].header
        self.datacube = fitsobj[0].data.astype('float')
        fitsobj.close()

        if int(self.header['NAXIS']) != 3:
            self.istat = -4
            if verbose==True: print 'ERROR: Input image is not in the required raw 3D data format'
            exit()

        if '1 1' in self.header['CCDSUM']:
            self.raw_frame_dims = ( 2048,2080 )  # Ymax,Xmax
            self.bias_sec = ( 1, 2048, 2055, 2080 )   # Windowed; actually starts at col 2049
            self.image_sec = ( 30, 2048, 1, 2048 )
        elif '2 2' in self.header['CCDSUM']:
            self.raw_frame_dims = ( 1024,1040 )
            self.bias_sec = ( 1, 1024, 1025, 1040 )
            self.image_sec = ( 16, 1024, 1, 1024 )
        self.naxis1 = self.image_sec[1]*2
        self.naxis2 = self.image_sec[3]*2
	
        gainhdr = str(self.header['GAIN'])
        gainlist = gainhdr.replace('[','').replace(']','').split(',')
        self.gain = []
        if '[' in gainhdr or ',' in gainhdr or len(gainlist) > 1:
            for value in gainlist:
                self.gain.append(float(value))
        else:
            self.gain.append(float(gainhdr))

        # Sanity check: There's either one gain for the whole frame, or one
        # value per quadrant.  No other option makes sense here:
        if len(self.gain) != 1 and (len(self.gain)) != 4:
            self.istat = -7
            if verbose==True: 
                print 'ERROR: Wrong number of GAIN value(s) in the FITS header keyword'
            
        self.quad_flux = {
                            1: [],
                            2: [],
                            3: [],
                            4: []
                            }
        self.quad_background = {
                            1: 0.0,
                            2: 0.0,
                            3: 0.0,
                            4: 0.0
                            }
        self.flux_min = 2000.0
        self.flux_max = 50000.0       	# Default: 75000.0
        self.model = 'linear' 	     # One of { linear, polynomial, broken_power_law }
	
    def overscan_statistics(self,iquad):
        self.overscan_median = np.median(self.datacube[iquad,self.bias_sec[0]:self.bias_sec[1],\
                                        self.bias_sec[2]:self.bias_sec[3]])
	
        regiondata = self.datacube[iquad,self.bias_sec[0]:self.bias_sec[1],\
                                        self.bias_sec[2]:self.bias_sec[3]]
        idx = np.where(regiondata < 1e9)
        for it in range(1,3,1):
            mean = regiondata[idx].mean()
            std = regiondata[idx].std()
            idx1 = np.where(regiondata >= (mean-3.0*std))
            idx2 = np.where(regiondata <= (mean+3.0*std))
            idx = np.intersect1d(idx1[0],idx2[0])
        self.overscan_mean = regiondata[idx].mean()
        self.overscan_std = regiondata[idx].std()
        print 'Overscan stats: ', self.overscan_mean, \
                                self.overscan_median,self.overscan_std
	
    def debias_quadrants(self):
        for iquad in range(0,4,1):
            self.overscan_statistics(iquad)
            self.datacube[iquad,:,:] = self.datacube[iquad,:,:] - self.overscan_mean
            print 'Debiased all quadrants'
	
    def subtract_bkgd(self):
        for iquad in range(0,4,1):
            imagedata = self.datacube[iquad,self.image_sec[0]:self.image_sec[1],self.image_sec[2]:self.image_sec[3]]
            (mean,stddev) = statistics.calcRMSclip2D(imagedata,3.0,3)
            self.datacube[iquad,:,:] = self.datacube[iquad,:,:] - mean
            self.quad_background[iquad] = mean
            print 'Subtracted '+str(mean)+' background from Q'+(str(iquad+1))
	    
    def normalize_gain(self):
        for iquad in range(0,4,1):
            if len(self.gain) == 1: 
                self.datacube[iquad,:,:] = ( self.quadimage * self.gain[0] )
            elif len(self.gain) == 4: 
                self.datacube[iquad,:,:] = ( self.datacube[iquad,:,:] * self.gain[iquad] )
        self.header['GAIN'] = 1.0
    
    def create_mask(self,Quadrant):
        print 'Datacube min max: ',self.datacube[Quadrant-1,:,:].min(), self.datacube[Quadrant-1,:,:].max()
        self.region = self.datacube[Quadrant-1,self.image_sec[0]:self.image_sec[1],self.image_sec[2]:self.image_sec[3]]
        self.maskidx = statistics.select_pixels_in_flux_range(self.region,2000.0,45000.0)
        self.region = self.region[self.maskidx]
        if len(self.region) == 0:
            print 'Warning: All pixels masked from quadrant!'
	    
    def correlate_flux(self):
        for iquad in range(0,4,1):
            print 'Correlating flux for quadrant '+str(iquad+1)
            qregion = self.datacube[iquad,self.image_sec[0]:self.image_sec[1], \
                        self.image_sec[2]:self.image_sec[3]]
            xdata = self.region
            ydata = qregion[self.maskidx]
            ydata = ydata.flatten()
            self.quad_flux[iquad+1] = [xdata,ydata]
	    
    def apply_correction(self,Camera,Binning,Coefficients,debugname=None):
        if debugname != None:
            (nz,nx,ny) = self.datacube.shape
            diffcube = np.zeros([nz,nx,ny])
        
        for pQ in range(0,4,1):
            pquadrant = self.datacube[pQ,self.image_sec[0]:self.image_sec[1], \
                        self.image_sec[2]:self.image_sec[3]]
            # Loop over the other quadrants, and correct their image data for the
            # crosstalk from the primary quadrant:
            for iquad in range(0,4,1):
                if pQ != iquad:
                    print 'Correcting quad '+str(iquad+1)+' for pixels from Q'+str(pQ+1)+\
                        ', coefficient='+str(Coefficients[Camera][Binning][pQ,iquad])
                    qregion = self.datacube[iquad,self.image_sec[0]:self.image_sec[1],\
                        self.image_sec[2]:self.image_sec[3]]
                    
                    dqregion = Coefficients[Camera][Binning][pQ,iquad] * pquadrant
                    
                    qregion = qregion - dqregion
                    
                    self.datacube[iquad,self.image_sec[0]:self.image_sec[1], \
                        self.image_sec[2]:self.image_sec[3]] = qregion
                
                    if debugname != None:
                        diffcube[iquad,self.image_sec[0]:self.image_sec[1],\
                            self.image_sec[2]:self.image_sec[3]] = dqregion
        
        if debugname != None:
            iexec = output3Dimage(self.header,diffcube,debugname)
		    
def read_coefficients(ConfigFile):
    """Function to read the configuration file containing the measured 
    crosstalk coefficients"""
    
    def initcoeffs(Coefficients,camera):
        if camera not in Coefficients.keys():
            Coefficients[camera] = { 1: np.zeros([4,4]), 2: np.zeros([4,4]) }
        return Coefficients
    
    if path.isfile(ConfigFile) == False:
        print 'ERROR: Cannot read co-efficients configuration file'
        print 'Looking for '+ConfigFile
        exit()
    
    Coefficients = {}
    fileobj = open(ConfigFile,'r')
    linelist = fileobj.readlines()
    fileobj.close()
    
    for line in linelist:
        if line[0:1] != '#' and len(line.replace(' ','').replace('\n','')) > 0:
            (camera, binning, pQ, c1,c2,c3,c4) = line.split()
            binning = int(float(binning))
            pQ = int(float(pQ)) - 1
            c1 = float(c1)
            c2 = float(c2)
            c3 = float(c3)
            c4 = float(c4)
            if camera not in Coefficients.keys(): 
                Coefficients = initcoeffs(Coefficients,camera)
            Coefficients[camera][binning][pQ,0] = c1
            Coefficients[camera][binning][pQ,1] = c2
            Coefficients[camera][binning][pQ,2] = c3
            Coefficients[camera][binning][pQ,3] = c4
    
    return Coefficients

def outputimage(Header,ImageData,Filename):
    newhdu = fits.PrimaryHDU(ImageData)
    newhdulist = fits.HDUList([newhdu])
    for key,value in Header.items():
        if key not in [ 'NAXIS1', 'NAXIS2' ]:
            newhdulist[0].header[key] = value
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        newhdulist.writeto(Filename,clobber=True,output_verify='ignore')
    newhdulist.close()
    return 0

def output3Dimage(Header,ImageData,Filename):
    newhdu = fits.PrimaryHDU(ImageData)
    newhdulist = fits.HDUList([newhdu])
    for key,value in Header.items():
        newhdulist[0].header[key] = value
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        newhdulist.writeto(Filename,clobber=True,output_verify='ignore')
    newhdulist.close()
    return 0

def fit_gradient(xdata,ydata,pinit):
    """Function to fit a straight line function with an intercept of zero on 
    the y-axis to the given datasets."""
    
    fitfunc = lambda p, x: p[1]*x
    errfunc = lambda p, x, y: fitfunc(p,x) - y
    
    try:
        (p1,istat) = optimize.leastsq(errfunc, pinit, args=(xdata,ydata))
    except TypeError:
        print xdata
        print ydata
        exit()
        
    y = fitfunc(p1,xdata)
    rms = np.sqrt( (( ydata - y )*( ydata - y )).sum() / float( len(y) ) )
    
    return p1,fitfunc,errfunc, rms
    
def fit_polynomial_zero(xdata,ydata,pinit):
    """Function to fit a 2nd order polynomial function with an intercept of 
    zero on the y-axis to the given datasets."""
    
    fitfunc = lambda p, x: p[1]*x + p[2]*x*x
    errfunc = lambda p, x, y: fitfunc(p,x) - y
    
    try:
        (p1,istat) = optimize.leastsq(errfunc, pinit, args=(xdata,ydata))
    except TypeError:
        print xdata
        print ydata
        exit()
    
    return p1,fitfunc,errfunc

def fit_broken_power_law(xdata,ydata,pinit):
    """Function to fit a 2nd order polynomial function with an intercept of 
    zero on the y-axis to the given datasets."""
    
    def fitfunc(p,x):
        y = zeros(len(x))
        idx = np.where(x < 40000.0)
        y[idx] = p[1]*x
        idx = np.where(x > 40000.0)
        y[idx] = p[2]*x[idx]
        return y
    
    errfunc = lambda p, x, y: fitfunc(p,x) - y
    try:
        (p1,istat) = optimize.leastsq(errfunc, pinit, args=(xdata,ydata))
    except TypeError:
        print xdata
        print ydata
        exit()
    
    return p1,fitfunc,errfunc

def iterative_model_fit(xdata,ydata,pinit,fit_function,sigclip=3.0):
    """Fit a function iteratively with sigma clipping"""
    
    a1 = 9e36
    afit = []
    for i in range(0,len(pinit)-1,1):
        afit.append(0.0)
    i = 0
    while (abs(a1-afit[1]) > 1.0):
        i = i + 1
        a1 = afit[1]
        (afit,fitfunc, errfunc, rms) = fit_function(xdata,ydata,pinit)
        yfit = fitfunc(afit,xdata)
        resids = ydata - yfit
        stddev = resids.std()
        #idx = np.where( resids <= sigclip*rms)
        jdx = np.where( resids <= 200.0 )
        kdx = np.where( resids >= -200.0 )
        idx = np.intersect1d(jdx,kdx)
        resids = resids[idx]
        xdata = xdata[idx]
        ydata = ydata[idx]
        stddev = resids.std()
        print ydata
        print i,a1,afit
    
    return afit,fitfunc, errfunc, stddev

def crossanalysis1(ImageFile,Quadrant,verbose=False):
    """Function to analyse a single raw Sinistro frame"""
    
    if verbose==True:
        print 'Processing SINISTRO frame '+ImageFile
        print 'Primary quadrant is '+str(Quadrant)
	
    if path.isfile(ImageFile) == False:
        print 'Error: Cannot find frame '+ImageFile
        exit()
    
    imageobj = CrossImage(ImageFile)
    imageobj.debias_quadrants()
    imageobj.subtract_bkgd()
    imageobj.create_mask(Quadrant)
    imageobj.correlate_flux()
    
    return status_code[0], imageobj

def multicrossanalysis(data_dir,out_dir,ImageList,Quadrant,PlotFile,verbose=False):
    """Function to analyse a set of increasing exposures of a single pointing with a bright
    star in one quadrant."""
    
    if path.isfile(ImageList) == False:
        print 'ERROR: Cannot find imagelist '+ImageList
        exit()
    fileobj = open(ImageList,'r')
    FrameList = fileobj.readlines()
    fileobj.close()
    
    xplot = {
        1: 0.0,
        2: 0.0,
        3: 0.0,
        4: 0.0
        }
    yplot = {
        1: 0.0,
        2: 0.0,
        3: 0.0,
        4: 0.0
        }
    for i,imagefile in enumerate(FrameList):
        uframe = archive_access.fetch_frame(path.join(data_dir,imagefile),\
                                            out_dir)
        (status, imageobj) = crossanalysis1(uframe,Quadrant,verbose=True)
        for iquad in range(0,4,1):
            (xdata,ydata) = imageobj.quad_flux[iquad+1]
            if i == 0:
                xplot[iquad+1] = xdata
                yplot[iquad+1] = ydata
            else:
                xplot[iquad+1] = np.concatenate((xplot[iquad+1],xdata))
                yplot[iquad+1] = np.concatenate((yplot[iquad+1],ydata))

    print 'Completed analysis, plotting...'

    fmt = [ 'r', 'b', 'm', 'g' ]
    fig = pyplot.figure(1)
    pyplot.rcParams['font.size'] = 10.0
    plotord = [ 2, 3, 1, 4 ]
    if imageobj.model == 'linear': 
        pinit = [ 0.0, 0.0 ]
    elif imageobj.model == 'polynomial' or imageobj.model == 'broken_power_law':
        pinit = [ 0.0, 0.0, 0.0 ]
    coeffs = {}
    for q,iquad in enumerate(plotord):
        xdata = xplot[iquad]
        ydata = yplot[iquad]
        coeffs[iquad] = []
    
        for i in range(0,1,1):
            ax = pyplot.subplot(2,2,q+1)
            pyplot.subplots_adjust(left=0.125, bottom=0.1, right=0.9, \
                        top=0.9,wspace=0.3,hspace=0.35)
            pyplot.plot(xdata,ydata,fmt[iquad-1]+'.')

            fileobj= open(path.join(out_dir, 'crosstalk_data_Q'+str(iquad)+'.txt'),'w')
            for k in range(0,len(xdata),1):
                fileobj.write(str(xdata[k])+'  '+str(ydata[k])+'\n')
            fileobj.close()
            
            if iquad != Quadrant:
                idx = statistics.select_entries_within_bound(xdata,\
                            imageobj.flux_min,imageobj.flux_max)
            
                if imageobj.model == 'linear':
                    (afit,fitfunc,errfunc,stddev) = iterative_model_fit(xdata[idx],\
                                ydata[idx],pinit,fit_gradient)
                    label = 'p[1]='+str(round(afit[1],8))+'\nsig='+str(round(stddev,8))
                
                elif imageobj.model == 'polynomial':
                    (afit,fitfunc,errfunc,stddev) = iterative_model_fit(xdata[idx],\
                            ydata[idx],pinit,fit_polynomial_zero)
                    label = 'p[1]='+str(round(afit[1],10)) + '\np[2]='+str(round(afit[2],10))
                
                elif imageobj.model == 'broken_power_law':
                    (afit,fitfunc,errfunc,stddev) = iterative_model_fit(xdata[idx],\
                            ydata[idx],pinit,fit_broken_power_law)
                    label = 'p[1]='+str(round(afit[1],10)) + '\np[2]='+str(round(afit[2],10))
                    
                if afit[1] > 0.0: 
                    coeffs[iquad].append(afit[1])
                else: 
                    coeffs[iquad].append(0.0)
                
                xmodel = np.arange(0,xdata[idx].max(),100)
                pyplot.plot(xmodel,fitfunc(afit,xmodel),'k-',label=label)
                ymodel = fitfunc(afit,xdata)
                ydata = ydata - ymodel
            
        if iquad in [ 1, 4 ]: 
            pyplot.xlabel('Quadrant '+str(Quadrant)+' pixel value [ADU]')
        pyplot.ylabel('Pixel value [ADU]')
        pyplot.xticks(rotation=45)
        (xmin,xmax,ymin,ymax) = pyplot.axis()
        if iquad != Quadrant: 
            pyplot.axis([xmin,xmax,-600.0,800.0])
        pyplot.title('Quadrant '+str(iquad))
        if iquad != Quadrant: pyplot.legend(loc='best')
    pyplot.savefig(PlotFile)
    pyplot.close(1)
    
    fig = pyplot.figure(2)
    pyplot.rcParams['font.size'] = 10.0
    plotord = [ 2, 3, 1, 4 ]
    if imageobj.model == 'linear':
        pinit = [ 0.0, 0.0 ]
    elif imageobj.model == 'polynomial' or imageobj.model == 'broken_power_law':
        pinit = [ 0.0, 0.0, 0.0 ]
    for q,iquad in enumerate(plotord):
        pyplot.subplot(2,2,q+1)
        pyplot.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, \
                    wspace=0.3,hspace=0.35)
        xdata = xplot[iquad]
        ydata = yplot[iquad]
        pyplot.scatter(xdata,ydata,c=fmt[iquad-1], marker='o',s=0.02)
    
        if iquad != Quadrant:
            idx = statistics.select_entries_within_bound(xdata,\
                    imageobj.flux_min,imageobj.flux_max)
            if imageobj.model == 'linear':
                (afit,fitfunc,errfunc,stddev) = iterative_model_fit(xdata[idx],\
                ydata[idx],pinit,fit_gradient)
                label = 'p[1]='+str(round(afit[1],8))+'\nRMS='+str(round(stddev,8))
            
            elif imageobj.model == 'polynomial':
                (afit,fitfunc,errfunc,stddev) = iterative_model_fit(xdata[idx],\
                            ydata[idx],pinit,fit_polynomial_zero)
                label = 'p[1]='+str(round(afit[1],5)) + '\np[2]='+str(round(afit[2],10))
            
            elif imageobj.model == 'broken_power_law':
                (afit,fitfunc,errfunc,stddev) = iterative_model_fit(xdata[idx],\
                            ydata[idx],pinit,fit_broken_power_law)
                label = 'p[1]='+str(round(afit[1],5)) + '\np[2]='+str(round(afit[2],10))
            
            xmodel = np.arange(0,xdata[idx].max(),100)
            pyplot.plot(xmodel,fitfunc(afit,xmodel),'k-',label=label)
            print 'Primary Quadrant='+str(Quadrant)+' quad='+str(iquad)+' parameters='+label
            
        if iquad in [ 1, 4 ]: 
            pyplot.xlabel('Quadrant '+str(Quadrant)+' pixel value [ADU]')
        pyplot.ylabel('Pixel value [ADU]')
        pyplot.xticks(rotation=45)
        (xmin,xmax,ymin,ymax) = pyplot.axis()
        if iquad != Quadrant: 
            pyplot.axis([xmax-20000,xmax,-200.0,200.0])
        else: 
            pyplot.axis([0.0,xmax,0.0,ymax])
        pyplot.title('Quadrant '+str(iquad))
        if iquad != Quadrant: pyplot.legend(loc='best')
    pyplot.savefig(PlotFile.replace('.png','_zoom.png'))
    pyplot.close(2)
    
    return status_code[0]

def correct_crosstalk(data_dir,out_dir,ImageList,ConfigFile):
    if path.isfile(ImageList) == False:
        print 'ERROR: Cannot find imagelist '+ImageList
        exit()
    fileobj = open(ImageList,'r')
    FrameList = fileobj.readlines()
    fileobj.close()
    
    Coefficients = read_coefficients(ConfigFile)
    
    for i,frame in enumerate(FrameList):
        uframe = archive_access.fetch_frame(path.join(data_dir,frame),\
                                            out_dir)
        imageobj = CrossImage(uframe)
        
        camera = imageobj.header['INSTRUME']
        binning = int(float(str(imageobj.header['CCDSUM']).split(' ')[0]))
        if camera not in Coefficients.keys():
            print 'ERROR: No calibration available for camera '+camera
            return status_code[-6]
        if binning not in Coefficients[camera].keys():
            print 'ERROR: No calibration available for camera '+camera+\
                    ' and binning '+str(binning)+'x'+str(binning)
            return status_code[-6]
        
        else:
            print 'Applying calibration for '+camera+' and binning '+\
                str(binning)+'x'+str(binning)+' to frame '+path.basename(uframe)
            imageobj.apply_correction(camera,binning,Coefficients,\
                    debugname=uframe.replace('.fits','_diff.fits'))
            filename = uframe.replace('e00.fits','e01c.fits')
            iexec = output3Dimage(imageobj.header,imageobj.datacube,filename)
    
    return status_code[0]

if __name__ == '__main__':
    debug = False
    HelpText = """                       CROSSTALK
	
	This script is designed to analyse raw Sinistro frames to measure the crosstalk
	between quadrants.
	
	Call sequence:
	python crosstalk.py -option [arguments]
	
	where:
	    -measure : to measure the crosstalk coefficients for one primary quadrant
	      	Arguments: DataLoc ImageList Quadrant PlotFile
	      	  data_dir is the full path to the input data directory
                   out_dir is the full path to the output data directory
	      	  ImageList is an ASCII list of frames to be analysed together
	      	  Quadrant indicates which quadrant {1-4} should be used as a reference
	      	  PlotFile is the full path to the output plot
	    -read_config : [DEBUG] read the configuration
	      	Arguments: ConfigFile
		  ConfigFile is the full path to the configuration file
      	    -correct : to apply the correction for crosstalk to a set of images
	      	Arguments: DataLoc ImageList ConfigFile
	      	  DataLoc is the full path to the data directory
	      	  ImageList is an ASCII list of frames to be analysed together
		  ConfigFile is the full path to the configuration file
    """
    if len(argv) < 2:
        print HelpText
        exit()
    
    elif argv[1].lower() == '-help' or argv[1].lower() not in ['-measure', '-correct', '-read_config']:
        print HelpText
    
    elif argv[1].lower() == '-measure':
        data_dir = argv[2]
        out_dir = argv[3]
        ImageList = argv[4]
        Quadrant = int(argv[5])
        PlotFile = argv[6]
        status = multicrossanalysis(data_dir,out_dir,ImageList,Quadrant,PlotFile)
        print status
    
    elif argv[1].lower() == '-read_config':
        ConfigFile = argv[2]
        Coefficients = read_coefficients(ConfigFile)
        print Coefficients
    
    elif argv[1].lower() == '-correct':
        data_dir = argv[2]
        out_dir = argv[3]
        ImageList = argv[4]
        ConfigFile = argv[5]
        status = correct_crosstalk(data_dir, out_dir,ImageList,ConfigFile)
        print status