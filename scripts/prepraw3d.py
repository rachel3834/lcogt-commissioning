##############################################################################
#     	      	      	  PREPARE RAW 3D SINISTRO FRAMES
#
# Software to prepare raw 3D Sinistro datacubes by computing the overscan
# from the edge of each cubeslice (i.e. quadrant), subtracting this
# from that quadrant, trimming it and patching the quadrants together into 
# the normal configuration
##############################################################################

###################################
# IMPORTED FUNCTIONS
from astropy.io import fits
from numpy import array, median, zeros, where, arange, fft,sqrt, flipud,fliplr,intersect1d
from scipy import optimize
#from matplotlib import use as useBackend
#useBackend('Agg')
#from matplotlib import pyplot
from sys import argv, exit
from os import path
import glob
import warnings
from shutil import move

import statistics

###################################
# DECLARATIONS

# Maximum x, y pixel values
# For fq04: RawFrameDims = ( 4296, 4096 )
# For fq02:
RawFrameDims = ( 2080, 2048 )

# Right-hand edge bias section:
# Pixel range: xmin, xmax, ymin, ymax
# (neglects bias section at bottom of frame)
BiasSec = ( 1, 2048, 2049, 2080 )

# Image data section in array (numbers from zero):
# Pixel range: xmin, xmax, ymin, ymax
ImageSec = ( 1, 2048, 1, 2048 )

# Quadrant definitions:
NAXIS1 = (ImageSec[3] - ImageSec[2] + 1)*2
NAXIS2 = (ImageSec[1] - ImageSec[0] + 1)*2
Regions = {}
# Transformations to switch the quadrants so that the image sections line up
Regions['Quad_transformations'] = {
        1: [],
        2: ['fliplr'],
        3: ['fliplr','flipud'],
        4: ['flipud']
        }
# Image coordinates (numbers from 1,1) such that the whole frame is aligned
# N up and E left. 
#Regions['Quad_dimensions'] = {
#      	  1: [ImageSec[1],NAXIS1,ImageSec[3],NAXIS2],
#	  2: [ImageSec[0]-1,ImageSec[1],ImageSec[3],NAXIS2],
#	  3: [ImageSec[0]-1,ImageSec[1],ImageSec[2]-1,ImageSec[3]],
#	  4: [ImageSec[1],NAXIS1,ImageSec[2]-1,ImageSec[3]]
#	  }
# Image coordinates (numbers from 1,1) such that the whole frame is aligned
# N down and E right, as expected by ORAC-DR. 
Regions['Quad_dimensions'] = {
        1: [ImageSec[0]-1,ImageSec[1],ImageSec[2]-1,ImageSec[3]],
        2: [ImageSec[1],NAXIS1,ImageSec[2]-1,ImageSec[3]],
        3: [ImageSec[1],NAXIS1,ImageSec[3],NAXIS2],
        4: [ImageSec[0]-1,ImageSec[1],ImageSec[3],NAXIS2]
        }

class FrameStats:
    def __init__(self):
        self.imagename = None
        self.overscan_mean = []
        self.overscan_median = []
        self.overscan_stddev = []

    def summary(self):
        outstr = self.imagename+'  '
        for qid in range(0,4,1):
            outstr = outstr + str(self.overscan_mean[qid])+' '
        for qid in range(0,4,1):
            outstr = outstr + str(self.overscan_median[qid])+' '
        for qid in range(0,4,1):
            outstr = outstr + str(self.overscan_stddev[qid])+' '
        return outstr

def outputimage(header,image_data,filename):
    newhdu = fits.PrimaryHDU(image_data)
    newhdulist = fits.HDUList([newhdu])
    for key,value in header.items():
        if key not in [ 'NAXIS1', 'NAXIS2' ]:
            try:
                newhdulist[0].header.update( (key,value) )
            except:
                newhdulist[0].header[key] = value
                
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        newhdulist.writeto(filename,clobber=True,output_verify='ignore')
    newhdulist.close()
    return 0
    
def prepraw3d(frame_path,preserve_overscan=False,dbg=False):
    """Function to prepare raw 3D Sinistro datacubes by computing the overscan
    from the edge of each cubeslice (i.e. quadrant), subtracting this
    from that quadrant, trimming it and patching the quadrants together into 
    the normal configuration."""
    
    if path.isfile(frame_path) == False:
        print 'Error: Cannot find frame '+frame_path
        exit()

    imageobj = fits.open(frame_path)
    image_header = imageobj[0].header
    data_cube = imageobj[0].data
    data_cube = data_cube.astype('float')
    imageobj.close()
    
    image_stats = FrameStats()
    image_stats.imagename = path.basename(frame_path)
    
    image_data = zeros([NAXIS2,NAXIS1])
    for iquad in range(0,4,1):
        quadimage = data_cube[iquad,:,:]
        if preserve_overscan == True:
            scanimage = quadimage[BiasSec[0]:BiasSec[1],BiasSec[2]:BiasSec[3]]
            fname = path.splitext(frame_path)[0]+'_overscan_Q'+str(iquad+1)+'.fits'
            iexec = outputimage(image_header,scanimage,fname)
        
        overscan_median = median(quadimage[BiasSec[0]:BiasSec[1],BiasSec[2]:BiasSec[3]])
        (overscan_mean,stddev) = statistics.calcRMSclip2D(quadimage[BiasSec[0]:BiasSec[1],BiasSec[2]:BiasSec[3]],3.0,3)
        image_stats.overscan_mean.append(overscan_mean)
        image_stats.overscan_median.append(overscan_median)
        image_stats.overscan_stddev.append(stddev)
        
        if dbg==True: print '--> Quadrant '+str(iquad+1)+\
            ' overscan sigClip mean value = '+str(overscan_mean)+\
            ' overscan median value = '+str(overscan_median)
            
        quadimage = quadimage - overscan_mean
        quadimage = quadimage[ImageSec[2]-1:ImageSec[3],ImageSec[0]-1:ImageSec[1]]
        
        # Apply transformation to orient the image the correct way up. 
        transform_list = Regions['Quad_transformations'][iquad+1]
        for transform in transform_list:
            if transform == 'flipud': quadimage = flipud(quadimage)
            if transform == 'fliplr': quadimage = fliplr(quadimage)
            
        # Add this quadrant into the appropriate section of the combined image:
        (xmin,xmax,ymin,ymax) = Regions['Quad_dimensions'][iquad+1]
        if dbg==True: 
            print 'Array dims: ',xmin,xmax,'(',(xmax-xmin),')'\
                                ,ymin,ymax,'(',(ymax-ymin),\
                ') quad shape: ',quadimage.shape, \
                ' final image dimensions: ',NAXIS2,NAXIS1
        image_data[ymin:ymax,xmin:xmax] = quadimage
            
    return image_stats,image_data, image_header

def output_2d(frame_path,preserve_overscan=False,dbg=False):
    
    (image_stats, image_data,image_header) = prepraw3d(frame_path,\
                            preserve_overscan=preserve_overscan,dbg=dbg)
    fname = path.splitext(frame_path)[0]+'.fits'
    move(frame_path,path.splitext(frame_path)[0]+'_3d.fits')
    iexec = outputimage(image_header,image_data,fname)
    
    return image_stats
    
def parseset(frame_path,preserve_overscan=False,debug=False):
    """Function to evaluate whether a single frame or a directory of frames 
    should be processed"""
    
    if path.isdir(frame_path) == True:
        framelist = glob.glob(path.join(frame_path,'*fits'))
        print framelist
        for frame in framelist:
            print 'Converting frame: '+path.basename(frame)
            (imagestats,imagedata, hdr) = prepraw3d(frame,\
                            preserve_overscan=preserve_overscan,dbg=debug)
            print imagestats.summary()
    
    else:
        print 'Converting frame: '+path.basename(frame_path)
        imagestats = output_2d(frame_path,\
                            preserve_overscan=preserve_overscan,dbg=debug)
    return 0

if __name__ == '__main__':
    debug=True
    keep_over_scan=False
    
    if len(argv) < 2:
        print 'Call sequence: python prepraw3d.py [FramePath]'
        print 'where FramePath can be either a single FITS datacube or a directory'
        print 'If a directory is given, script will attempt to convert all FITS files within that directory'
        print 'Optional argument:'
        print '-keep-overscan: cut the side overscan into a separate frame'
        exit()
    else:
        frame_path = argv[1]
        if len(argv) > 2 and '-keep-overscan' in argv[2]:
            keep_over_scan = True
    iexec = parseset(frame_path,preserve_overscan=keep_over_scan,debug=debug)
