########################################################################################################
#     	      	      	      	      	STATISTICS FUNCTIONS
########################################################################################################

# IMPORT FUNCTIONS
from numpy import array, where, nan, isnan, logical_not
from numpy import median, zeros, arange, intersect1d, fft
from matplotlib import pyplot
from sys import argv, exit
from scipy import meshgrid, stats
from astropy.io import fits

###################################
# CALCRMSCLIP
def calcRMSclip(RegionData,sigClip,niter):
    '''Function to calculate the mean and standard deviation of a set of pixel values
    given an image region data array, using an iterative, sigma-clipping function
    Function name is a bit of a misnomer - returns STDDEV now'''
    
    RegionData = RegionData[logical_not(isnan(RegionData))]
    
    idx = where(RegionData < 1e9)
    #print RegionData[idx].mean(dtype='float64')
    for it in range(0,niter,1):
        mean = RegionData[idx].mean(dtype='float64')
    std = RegionData[idx].std(dtype='float64')
    idx1 = where(RegionData >= (mean-sigClip*std))
    idx2 = where(RegionData <= (mean+sigClip*std))
    idx = intersect1d(idx1[0],idx2[0])
    #print it,mean,std,(mean-sigClip*std),(mean+sigClip*std),len(idx)
    mean = RegionData[idx].mean(dtype='float64')
    std = RegionData[idx].std(dtype='float64')
    
    return mean,std

###################################
# CALCRMSCLIP
def calcRMSclip2Didx(RegionData,sigClip,niter):
    '''Function to calculate the mean and standard deviation of a set of pixel values
    given an image region data array, using an iterative, sigma-clipping function
    Function name is a bit of a misnomer - returns STDDEV now'''
    
    RegionData = RegionData.flatten()
    RegionData = RegionData[logical_not(isnan(RegionData))]
    
    idx = where(RegionData < 1e9)
    #print RegionData[idx].mean(dtype='float64')
    for it in range(0,niter,1):
        mean = RegionData[idx].mean(dtype='float64')
    std = RegionData[idx].std(dtype='float64')
    idx1 = where(RegionData >= (mean-sigClip*std))
    idx2 = where(RegionData <= (mean+sigClip*std))
    idx = intersect1d(idx1[0],idx2[0])
    #print it,mean,std,(mean-sigClip*std),(mean+sigClip*std),len(idx)
    mean = RegionData[idx].mean(dtype='float64')
    std = RegionData[idx].std(dtype='float64')
    
    return mean,std,idx
    
###################################
# CALCRMSCLIP
def calcRMSclip2D(RegionData,sigClip,niter):
    '''Function to calculate the mean and standard deviation of a set of pixel values
    given an image region data array, using an iterative, sigma-clipping function
    Function name is a bit of a misnomer - returns STDDEV now
    '''
    
    RegionData = RegionData.flatten()
    (mean,std) = calcRMSclip(RegionData,sigClip,niter)
    
    return mean,std

###################################
# CALCRMSCLIP
def calcMAD2D(RegionData):
    '''Function to calculate the MAD of a 2D set of pixel values
    given an image region data array, using an iterative, sigma-clipping function'''
    
    RegionData = RegionData[logical_not(isnan(RegionData))]
    
    idx = where(RegionData < 1e9)
    region_median = median(RegionData[idx])
    residuals = RegionData - region_median
    MAD = median( abs(residuals[idx]) )
    
    return MAD
    
###################################
# CALCALPHACLIP2D
def calcAlphaclip2D(RegionData,alpha):
    '''Function to calculate the mean, median and standard deviation of a set of pixel values
    given an image region data array, using an iterative, sigma-clipping function'''
    
    # Exclude NaNs and flatten the array for easier processing:
    RegionData = RegionData[logical_not(isnan(RegionData))]
    RegionData = RegionData.flatten()
    
    # Exclude extreme values and compute the initial statistics:
    idx = where(RegionData < 1e9)
    immedian = median(RegionData[idx])
    immean = RegionData[idx].mean(dtype='float64')
    imstd = RegionData[idx].std(dtype='float64')
    
    # Exclude the extreme wings of the distribution:
    jdx = RegionData[idx].argsort()
    ihi = int(len(jdx)* alpha)
    ilo = int((1-alpha)*len(jdx))
    clipHi = RegionData[idx][jdx[ihi]]
    clipLo = RegionData[idx][jdx[ilo]]
    RegionData = RegionData[idx][jdx[ilo:ihi]]
    
    # Recalculate the statistics:
    immedian = median(RegionData)
    immean = RegionData.mean(dtype='float64')
    imstd = RegionData.std(dtype='float64')
    #print immedian, immean, imstd, len(RegionData),clipHi,clipLo
    
    return immean,immedian,imstd,clipHi,clipLo

###################################
# MASKEDMEDIAN3D
def maskedmedian3D(DataCube,alpha):
    '''Function to compute the median 2D image of a datacube of 2D images, applying masking
    to exclude pixel values outside the clipLimits'''
    
    # Compute clip limits, using alpha=fraction of pixels at extremes of ranges to exclude:
    (immean,immedian,imstddev,clipHi,clipLo) = calcAlphaclip2D(DataCube,alpha)
    
    # Mask the 3D array using the clip limits:
    amask = stats.mstats.trim(DataCube,(clipLo,clipHi))
    
    # Calculate the median:
    median_frame = median(amask,axis=0)
    print (median_frame)
    exit()
    
    return median_frame
    
###################################
# HISTPIX
def histpix(RegionData, plotfile,xlimit_min=None,xlimit_max=None,mean=None,stddev=None):
    '''Function to plot a histogram of the image pixel data within a given region'''
    
    # Plot histogram of pixel data itself:
    fig = pyplot.figure(1)
    vector = RegionData.flatten()
    fullvector = vector[logical_not(isnan(vector))]
    vector = fullvector.copy()
    nbins = max(10,((vector.max() - vector.min())/0.05))
    #print 'Nbins = ',nbins
    pyplot.hist(vector,bins=nbins)
    (xmin,xmax,ymin,ymax) = pyplot.axis()
    if xlimit_min != None: xmin = xlimit_min
    if xlimit_max != None: xmax = xlimit_max
    pyplot.axis([xmin,xmax,ymin,ymax])
    
    if mean != None and stddev != None:
        # Indicate the mean value with a vertical line:
        xdata = array( [mean-stddev, mean+stddev] )
        ydata = zeros(2)
        ydata.fill( (ymin + (ymax-ymin)/2.0) )
        pyplot.plot(xdata,ydata,'r-')
    
        # Indicate the stddev with a horizontal line:
        xdata = zeros(2)
        xdata.fill(mean)
        ydata = array( [ymin,ymax] )
        pyplot.plot(xdata,ydata,'r-')
    
    pyplot.xlabel('Pixel value [ADU]')
    pyplot.ylabel('Frequency')
    pyplot.title('Histogram of image pixel values')
    pyplot.savefig(plotfile)
    pyplot.close(1)

    return 0

############################################
# CALCDELTASTDDEV
def calcdeltastddev(DataArray,sigClip,niter):
    '''Function to calculate the stddev of an a array of values both before and after
    sigma clipping it'''
    
    mean = DataArray.mean(dtype='float64')
    stddev = DataArray.std(dtype='float64')
    init_stats = (mean,stddev)
    
    idx = where(DataArray < 1e9)
    for it in range(0,niter,1):
        mean = DataArray[idx].mean(dtype='float64')
    stddev = DataArray[idx].std(dtype='float64')
    idx1 = where(DataArray >= (mean-sigClip*stddev))
    idx2 = where(DataArray <= (mean+sigClip*stddev))
    idx = intersect1d(idx1[0],idx2[0])
    print ('stats: ', it,mean,stddev,(mean-sigClip*stddev),(mean+sigClip*stddev),len(idx))
    mean = DataArray[idx].mean(dtype='float64')
    stddev = DataArray[idx].std(dtype='float64')  
    final_stats = (mean,stddev)
    
    return init_stats, final_stats

############################################
# EXCLUDEPOINTS
def excludepoints(DataArray,npts):
    '''Function to calculate the stddev of an a array of values both before and after
    excluding the top npts outlyers'''
    
    mean = DataArray.mean(dtype='float64')
    stddev = DataArray.std(dtype='float64')
    init_stats = (mean,stddev)
    
    darray = abs(DataArray - mean)
    idx = darray.argsort()
    
    mean = DataArray[idx[:-2]].mean()
    stddev = DataArray[idx[:-2]].std()
    
    final_stats = (mean,stddev)
    
    return init_stats, final_stats

############################################
# SELECT PIXELS IN FLUX RANGE
def select_pixels_in_flux_range(imagedata,flux_min,flux_max):
    '''Function to select those pixels from a 2D image with values between the fluxes given.
    Returns an np.where-like array.'''


    idx = where( (imagedata > flux_min) & (imagedata < flux_max))
    return idx
    
############################################
# SELECT ENTRIES WITH BOUNDS
def select_entries_within_bound(vector,entry_min,entry_max):
    '''Function to select those pixels from a 1D vector with values between the fluxes given.
    Returns an np.where-like array.'''
    
    # Identify all entries with values that are too low to include.  Set these 
    # pixels to values where they will be excluded by the next selection cut:
    #mask = vector.copy()
    #idx = where(vector < entry_min)
    #mask[idx] = entry_max + 1000.0
    
    # Now select everything below the maximum threshold, which excludes both too high and 
    # too low entries:
    idx = where( (vector < entry_max) & (vector > entry_min))
   # print 'Number of selected pixels=',len(idx[0])
    return idx
    
############################################
# COMMANDLINE RUN SECTION
if __name__ == '__main__':

    if len(argv) < 3:
        print ('Call sequence: python statistics [image path] sigclip niter')
        exit()
    else:
        ImageFile = argv[1]
        sigClip = float(argv[2])
        niter = int(float(argv[3]))

        # Read the input image:
        imageheader = fits.getheader()
        imagedata = fits.getdata()

        #(mean,stddev) = calcRMSclip2D(imagedata,sigClip,niter)
        #print 'Mean = ',mean
        #print 'Std. deviation = ',stddev
        
        select_pixels_in_flux_range(imagedata,500.0,40000.0)
