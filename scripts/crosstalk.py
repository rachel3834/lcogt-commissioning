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


import argparse
import logging
import os
import re
import sys
from os import path
from sys import exit

# useBackend('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy import optimize

import statistics
from Image import Image

_logger = logging.getLogger(__name__)

###################################
# DECLARATIONS

status_code = {0: 'OK',
               -1: 'Error: Cannot find frame ', \
               -2: 'Error: Missing data location', \
               -3: 'Error: More processed frames than raw - corrupted data location?', \
               -4: 'Warning: frame not in 3D data format', \
               -5: 'Warning: frame zero-length (will be ignored)', \
               -6: 'Error: No coefficients for camera and binning combination'
               }
verbose = True


#
#     def apply_correction(self,Camera,Binning,Coefficients,debugname=None):
#         if debugname != None:
#             (nz,nx,ny) = self.datacube.shape
#             diffcube = np.zeros([nz,nx,ny])
#
#         for pQ in range(0,4,1):
#             pquadrant = self.datacube[pQ,self.image_sec[0]:self.image_sec[1], \
#                         self.image_sec[2]:self.image_sec[3]]
#             # Loop over the other quadrants, and correct their image data for the
#             # crosstalk from the primary quadrant:
#             for iquad in range(0,4,1):
#                 if pQ != iquad:
#                     print 'Correcting quad '+str(iquad+1)+' for pixels from Q'+str(pQ+1)+\
#                         ', coefficient='+str(Coefficients[Camera][Binning][pQ,iquad])
#                     qregion = self.datacube[iquad,self.image_sec[0]:self.image_sec[1],\
#                         self.image_sec[2]:self.image_sec[3]]
#
#                     dqregion = Coefficients[Camera][Binning][pQ,iquad] * pquadrant
#
#                     qregion = qregion - dqregion
#
#                     self.datacube[iquad,self.image_sec[0]:self.image_sec[1], \
#                         self.image_sec[2]:self.image_sec[3]] = qregion
#
#                     if debugname != None:
#                         diffcube[iquad,self.image_sec[0]:self.image_sec[1],\
#                             self.image_sec[2]:self.image_sec[3]] = dqregion
#
#         if debugname != None:
#             iexec = output3Dimage(self.header,diffcube,debugname)

def read_coefficients(ConfigFile):
    """Function to read the configuration file containing the measured 
    crosstalk coefficients"""

    def initcoeffs(Coefficients, camera):
        if camera not in Coefficients.keys():
            Coefficients[camera] = {1: np.zeros([4, 4]), 2: np.zeros([4, 4])}
        return Coefficients

    if path.isfile(ConfigFile) == False:
        print ('ERROR: Cannot read co-efficients configuration file')
        print ('Looking for ' + ConfigFile)
        exit()

    Coefficients = {}
    fileobj = open(ConfigFile, 'r')
    linelist = fileobj.readlines()
    fileobj.close()

    for line in linelist:
        if line[0:1] != '#' and len(line.replace(' ', '').replace('\n', '')) > 0:
            (camera, binning, pQ, c1, c2, c3, c4) = line.split()
            binning = int(float(binning))
            pQ = int(float(pQ)) - 1
            c1 = float(c1)
            c2 = float(c2)
            c3 = float(c3)
            c4 = float(c4)
            if camera not in Coefficients.keys():
                Coefficients = initcoeffs(Coefficients, camera)
            Coefficients[camera][binning][pQ, 0] = c1
            Coefficients[camera][binning][pQ, 1] = c2
            Coefficients[camera][binning][pQ, 2] = c3
            Coefficients[camera][binning][pQ, 3] = c4

    return Coefficients


def outputimage(Header, ImageData, Filename):
    newhdu = fits.PrimaryHDU(ImageData)
    newhdulist = fits.HDUList([newhdu])
    for key, value in Header.items():
        if key not in ['NAXIS1', 'NAXIS2']:
            newhdulist[0].header[key] = value
    try:
        newhdulist.writeto(Filename, clobber=True, output_verify='ignore')
    except Exception as  e:
        _logger.error ("While saving output: " + e)

    newhdulist.close()
    return 0


def output3Dimage(Header, ImageData, Filename):
    newhdu = fits.PrimaryHDU(ImageData)
    newhdulist = fits.HDUList([newhdu])
    for key, value in Header.items():
        newhdulist[0].header[key] = value
        try:
            newhdulist.writeto(Filename, clobber=True, output_verify='ignore')
        except Exception as e:
            _logger.error ("while writing " + e )
    newhdulist.close()
    return 0


def fit_gradient(xdata, ydata, pinit):
    """Function to fit a straight line function with an intercept of zero on 
    the y-axis to the given datasets."""

    fitfunc = lambda p, x: p[1] * x
    errfunc = lambda p, x, y: fitfunc(p, x) - y

    try:
        (p1, istat) = optimize.leastsq(errfunc, pinit, args=(xdata, ydata))
    except TypeError as te:
        print ("There was an error!", te)
        print (xdata)
        print (ydata)
        exit()

    y = fitfunc(p1, xdata)
    rms = np.sqrt(((ydata - y) ** 2).sum() / float(len(y)))

    return p1, fitfunc, errfunc, rms


def fit_polynomial_zero(xdata, ydata, pinit):
    """Function to fit a 2nd order polynomial function with an intercept of 
    zero on the y-axis to the given datasets."""

    fitfunc = lambda p, x: p[1] * x + p[2] * x * x
    errfunc = lambda p, x, y: fitfunc(p, x) - y

    try:
        (p1, istat) = optimize.leastsq(errfunc, pinit, args=(xdata, ydata))
    except TypeError:
        print (xdata)
        print (ydata)
        exit()

    return p1, fitfunc, errfunc


def fit_broken_power_law(xdata, ydata, pinit):
    """Function to fit a 2nd order polynomial function with an intercept of 
    zero on the y-axis to the given datasets."""

    def fitfunc(p, x):
        y = np.zeros(len(x))
        idx = np.where(x < 40000.0)
        y[idx] = p[1] * x
        idx = np.where(x > 40000.0)
        y[idx] = p[2] * x[idx]
        return y

    errfunc = lambda p, x, y: fitfunc(p, x) - y
    try:
        (p1, istat) = optimize.leastsq(errfunc, pinit, args=(xdata, ydata))
    except TypeError:
        print (xdata)
        print (ydata)
        exit()

    return p1, fitfunc, errfunc


def iterative_model_fit(xdata, ydata, pinit, fit_function, sigclip=3.0):
    """Fit a function iteratively with sigma clipping"""

    def calc_resids(afit, ydata, fitfunc,nsig):
        yfit = fitfunc(afit, xdata)
        resids = ydata - yfit
        rms = np.std (resids)
        idx = np.where( np.abs(resids) <= nsig * rms)

        return idx, resids

    a1 = 9e36
    afit = [0.0, 0.001]
    (afit, fitfunc, errfunc, rms) = fit_function(xdata, ydata, pinit)
    (idx, resids) = calc_resids(afit, ydata, fitfunc, sigclip)
    i = 0
    cont = True
    while cont == True:
        i = i + 1
        a1 = afit[1]
        (afit, fitfunc, errfunc, rms) = fit_function(xdata[idx], ydata[idx], pinit)
        (idx, resids) = calc_resids(afit, ydata, fitfunc, sigclip)
        stddev = resids.std()

        if (abs(a1 - afit[1]) > 1e-5) or i == 10:
            cont = False

    return afit, fitfunc, errfunc, stddev, idx


def create_mask(image, Quadrant, fluxmin, fluxmax):
    """
    CValcualte a mask that contains all pixels in a quadrant that are within flux_min and flux_max level
    :param Quadrant:
    :param flux_min:
    :param flux_max:
    :return:
    """

    region = image.getccddata(Quadrant)
    _logger.debug("Input datacube min max %f %f " % (region.min(), region.max()))
    maskidx = statistics.select_pixels_in_flux_range(region, fluxmin, fluxmax)
    if len(maskidx) == 0:
        _logger.warn("No pixels in selected flux range!")
    return maskidx


def correlate_flux(image, sourcequadrant, maskidx):
    quadfluxes = []

    for iquad in range(4):
        xdata = image.getccddata(sourcequadrant)[maskidx].flatten()
        ydata = image.getccddata(iquad)[maskidx].flatten()

        # now make sure that we only see crostalk, and not additional stellar stuff in the fied.
        selectidx = ydata < 1e-3 * xdata # TODO: paramterize maximum allowable crosstalk.
        quadfluxes.append([xdata[selectidx], ydata[selectidx]])

    return quadfluxes


def crossanalysis1(_image, args):
    """Function to analyse a single raw Sinistro frame"""

    maskidx = create_mask(_image, args.opt_quadrant, args.fluxmin, args.fluxmax)
    quadfluxes = correlate_flux(_image, args.opt_quadrant, maskidx)

    return quadfluxes


def multicrossanalysis(args):
    """Function to analyse a set of increasing exposures of a single pointing with a bright
    star in one quadrant."""

    xplot = {}
    yplot = {}

    # build up statistics for all files we have given to us
    for ii, imagefile in enumerate(args.fitsfile):

        _image = Image(imagefile, gaincorrect=False, overscancorrect=True, skycorrect=True)

        quadfluxes = crossanalysis1(_image, args)

        for iquad in range(4):
            (xdata, ydata) = quadfluxes[iquad]
            if ii == 0:
                xplot[iquad + 1] = xdata
                yplot[iquad + 1] = ydata
            else:
                xplot[iquad + 1] = np.concatenate((xplot[iquad + 1], xdata))
                yplot[iquad + 1] = np.concatenate((yplot[iquad + 1], ydata))

    _logger.debug('Completed data fetching, plotting and fitting next...')

    fmt = ['r', 'b', 'm', 'g']
    fig = plt.figure(1)
    plt.rcParams['font.size'] = 10.0
    plotord = [2, 3, 1, 4]


    for zoomlevel in [1,2]:
        if args.linear:
            pinit = [0.0, 0.0]
        elif args.poly:
            pinit = [0.0, 0.0, 0.0]
        coeffs = {}

        for q, iquad in enumerate(plotord):
            xdata = xplot[iquad]
            ydata = yplot[iquad]

            coeffs[iquad] = []

            for ii in range(0, 1, 1):
                ax = plt.subplot(2, 2, q + 1)
                plt.subplots_adjust(left=0.125, bottom=0.15, right=0.9, top=0.9, wspace=0.3, hspace=0.35)
                plt.plot(xdata, ydata, 'k,')

                if iquad != args.opt_quadrant + 1:
                    idx = statistics.select_entries_within_bound(xdata, args.fluxmin, args.fluxmax)

                    if args.linear:
                        (afit, fitfunc, errfunc, stddev, kdx) = iterative_model_fit(xdata[idx], ydata[idx], pinit,
                                                                                fit_gradient, sigclip=5)
                        label = 'p[1]=' + str(round(afit[1], 10)) + '\nsig=' + str(round(stddev, 2))
                        if zoomlevel == 1:
                            print ('archon.header.CRSTLK%d%d = %9f' % ( (args.opt_quadrant + 1),  iquad, (round(afit[1],5))))
                            #print 'Primary Quadrant=' + str(args.opt_quadrant + 1) + ' quad=' + str(iquad) + ' parameters=' + label
                    elif args.poly:
                        (afit, fitfunc, errfunc, stddev, kdx) = iterative_model_fit(xdata[idx], ydata[idx], pinit,
                                                                                fit_polynomial_zero)
                        label = 'p[1]=' + str(round(afit[1], 10)) + '\np[2]=' + str(round(afit[2], 10))

                    if afit[1] > 0.0:
                        coeffs[iquad].append(afit[1])
                    else:
                        coeffs[iquad].append(0.0)

                    plt.plot(xdata[kdx], ydata[kdx], 'r,')
                    xmodel = np.arange(0, xdata[idx].max(), 100)
                    plt.plot(xmodel, fitfunc(afit, xmodel), 'k-', label=label)
                    ymodel = fitfunc(afit, xdata)
                    ydata = ydata - ymodel

            if iquad in [1, 4]:
                plt.xlabel('Quadrant ' + str(args.opt_quadrant + 1) + ' pixel value [ADU]')
            plt.ylabel('Pixel value [ADU]')
            plt.xticks(rotation=15)
            (xmin, xmax, ymin, ymax) = plt.axis()


            if iquad != args.opt_quadrant + 1:
                if zoomlevel == 1:
                    plt.axis([0, 65000, -60.0, 60.0])
                    plotfile = args.plotfile
                if zoomlevel == 2:
                    plt.axis([xmax - 10000, xmax + 1000, -100.0, 100.0])
                    plotfile = args.plotfile.replace('.png', '_zoom.png')
            else:
                plt.axis([0, 65000, 0, 65000])
            plt.title('Quadrant ' + str(iquad))
            if iquad != args.opt_quadrant + 1:
                plt.legend(loc='best')

        #plt.show()
        plt.savefig(plotfile)
        plt.close(1)


    return status_code[0]


def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='Measure cross talk of multi-amplifier FITS images.')

    parser.add_argument('fitsfile', type=str, nargs='+',
                        help='Fits files for cross talk measurement')

    parser.add_argument('--imagepath', dest='opt_imagepath', type=str, default=None)

    parser.add_argument('--log_level', dest='log_level', default='INFO', choices=['DEBUG', 'INFO'],
                        help='Set the debug level')

    parser.add_argument('--measure', dest='mode_measure', type=bool, default=True,
                        help='Measure crosstalk')

    parser.add_argument('--linear', action='store_true')
    parser.add_argument('--poly', action='store_true')

    parser.add_argument('--useheader', action="store_true", help="use OBJECT line to determine quadrant. ")

    parser.add_argument('--outputdir', dest='outdir', type=str, default='.',
                        help='Fie containing input or output xtalk measurements')

    parser.add_argument('--plotfile', dest='plotfile', type=str, default='xtalk.png',
                        help='File for crosstalk graph')

    parser.add_argument('--quadrant', dest='opt_quadrant', type=int,
                        help='Quadrant containing bright star. Start counting at 0')

    parser.add_argument('--minflux', dest='fluxmin', type=float, default='10000',
                        help='Minimum contaminationg flux')

    parser.add_argument('--maxflux', dest='fluxmax', type=float, default='50000',
                        help='Maximum contaminationg flux')

    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')

    for ii in range(len(args.fitsfile)):
        if args.opt_imagepath is not None:
            args.fitsfile[ii] = "%s/%s" % (args.opt_imagepath, args.fitsfile[ii])
        if not os.path.isfile(args.fitsfile[ii]):
            _logger.error("Fatal: file %s does not exists" % (args.fitsfile[ii]))
            sys.exit(0)

    return args


if __name__ == '__main__':

    args = parseCommandLine()

    if args.mode_measure:

        if args.opt_quadrant is not None:


            status = multicrossanalysis(args)

        else:
            print ("Trying to identofy location of contaminating source based on OBJECT name")

            QtoFitsdict = dict()
            fitsfiles = args.fitsfile

            for file in fitsfiles:
                hdul = fits.open (file)
                obj = hdul[0].header['OBJECT']
                hdul.close()
                m = re.search ("x talk q (\\d)", obj)
                if m:

                    quadrant = m.group(1)
                    print (obj, " -> ",    quadrant)

                    if quadrant in QtoFitsdict:
                        QtoFitsdict[quadrant].append (file)
                    else:
                        QtoFitsdict[quadrant] = [file]

            originalplotfile = args.plotfile
            for quadrant in QtoFitsdict:
                args.opt_quadrant = int (quadrant)
                args.fitsfile = QtoFitsdict[quadrant]
                args.plotfile = originalplotfile.replace(".png", "_%d.png" % int(quadrant))
                multicrossanalysis(args)

