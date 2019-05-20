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
import sys
from sys import exit
# useBackend('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
import statistics
from lcocommissioning.Image import Image

log = logging.getLogger(__name__)


def fit_linear_zero(xdata, ydata, pinit):
    """Function to fit a straight line function with an intercept of zero on 
    the y-axis to the given datasets."""

    fitfunc = lambda p, x: p[1] * x
    errfunc = lambda p, x, y: fitfunc(p, x) - y

    try:
        (p1, istat) = optimize.leastsq(errfunc, pinit, args=(xdata, ydata))
    except TypeError as te:
        log.error("While doing linear fit: {}".format(te))
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
    except TypeError as te:
        log.error("While doing polynomial fit: {}".format(te))
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
    except TypeError as te:
        log.error("While doing powerlaw fit: {}".format(te))
        exit()

    return p1, fitfunc, errfunc


def iterative_model_fit(xdata, ydata, pinit, fit_function, sigclip=3.0):
    """Fit a function iteratively with sigma clipping"""

    def calc_resids(afit, ydata, fitfunc, nsig):
        yfit = fitfunc(afit, xdata)
        resids = ydata - yfit
        rms = np.std(resids)
        idx = np.where(np.abs(resids) <= nsig * rms)

        return idx, resids

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
    log.debug("Input datacube min max %f %f " % (region.min(), region.max()))
    maskidx = statistics.select_pixels_in_flux_range(region, fluxmin, fluxmax)
    if len(maskidx) == 0:
        log.warning("No pixels in selected flux range!")
    return maskidx


def correlate_flux(image, sourcequadrant, maskidx):
    quadfluxes = []
    for iquad in range(4):
        xdata = image.getccddata(sourcequadrant)[maskidx].flatten()
        ydata = image.getccddata(iquad)[maskidx].flatten()

        # now make sure that we only see crostalk, and not additional stellar stuff in the field.
        selectidx = ydata < 1e-3 * xdata  # TODO: paramterize maximum allowable crosstalk.
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

    log.debug('Completed data fetching, plotting and fitting next...')

    plt.figure(1)
    plt.rcParams['font.size'] = 10.0
    plotord = [2, 3, 1, 4]

    pinit = [0.0, 0.0]
    coeffs = {}

    for q, iquad in enumerate(plotord):
        xdata = xplot[iquad]
        ydata = yplot[iquad]
        coeffs[iquad] = []

        plt.subplot(2, 2, q + 1)
        plt.subplots_adjust(left=0.125, bottom=0.15, right=0.9, top=0.9, wspace=0.3, hspace=0.35)

        plt.plot(xdata, ydata, 'k,')

        if iquad != args.opt_quadrant + 1:
            idx = statistics.select_entries_within_bound(xdata, args.fluxmin, args.fluxmax)

            (afit, fitfunc, errfunc, stddev, kdx) = iterative_model_fit(xdata[idx], ydata[idx], pinit,
                                                                        fit_linear_zero, sigclip=3)
            label = 'p[1]=' + str(round(afit[1], 7)) + '\nsig=' + str(round(stddev, 2))

            print('archon.header.CRSTLK%d%d = %9f' % (
                (args.opt_quadrant + 1), iquad, (round(afit[1], 6))))

            if afit[1] > 0.0:
                coeffs[iquad].append(afit[1])
            else:
                coeffs[iquad].append(0.0)

            plt.plot(xdata[kdx], ydata[kdx], 'r,')
            xmodel = np.arange(0, xdata[idx].max(), 100)
            plt.plot(xmodel, fitfunc(afit, xmodel), 'k-', label=label)
            # ymodel = fitfunc(afit, xdata)
            # ydata = ydata - ymodel

        if iquad in [1, 4]:
            plt.xlabel('Source ' + str(args.opt_quadrant + 1) + ' pixel value [ADU]')
        plt.ylabel('Victim Pixel value [ADU]')
        plt.xticks(rotation=15)

        if iquad != args.opt_quadrant + 1:

            plt.axis([0, 65000, -50.0, 50.0])
            plotfile = args.plotfile

        else:
            plt.axis([0, 65000, 0, 65000])
        plt.title('Quadrant ' + str(iquad))
        if iquad != args.opt_quadrant + 1:
            plt.legend(loc='best')

    # plt.show()
    plotfile = "{}_{}.png".format(args.plotfile, args.opt_quadrant + 1)
    plt.savefig(plotfile)
    plt.close(1)

    return


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

    parser.add_argument('--useheader', action="store_true", help="use OBJECT line to determine quadrant. ")

    parser.add_argument('--outputdir', dest='outdir', type=str, default='.',
                        help='Fie containing input or output xtalk measurements')

    parser.add_argument('--plotfile', dest='plotfile', type=str, default='xtalk',
                        help='Root file name for crosstalk graph')

    parser.add_argument('--quadrant', dest='opt_quadrant', type=int,
                        help='FITS extension containing contaminating source. Start counting at 0')

    parser.add_argument('--minflux', dest='fluxmin', type=float, default='5000',
                        help='Minimum contaminating flux')

    parser.add_argument('--maxflux', dest='fluxmax', type=float, default='40000',
                        help='Maximum contaminating flux')

    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')

    for ii in range(len(args.fitsfile)):
        if args.opt_imagepath is not None:
            args.fitsfile[ii] = "%s/%s" % (args.opt_imagepath, args.fitsfile[ii])
        if not os.path.isfile(args.fitsfile[ii]):
            log.error("Fatal: file %s does not exists" % (args.fitsfile[ii]))
            sys.exit(0)

    return args


def main():
    args = parseCommandLine()
    if args.mode_measure:
        if args.opt_quadrant is not None:
            multicrossanalysis(args)
        else:
            for quad in range(4):
                args.opt_quadrant = quad
                multicrossanalysis(args)

if __name__ == '__main__':
    main()
