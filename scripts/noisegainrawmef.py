"""
Program to calculate the noise and gain for each extension of a given mef file
"""

import sys
import os
import numpy as np
import math
import argparse
from Image import Image
import matplotlib.pyplot as plt

import logging

_logger = logging.getLogger(__name__)


def noisegainextension(flat1, flat2, bias1, bias2, minx=None, maxx=None, miny=None, maxy=None, showImages=False):
    """
    Measure the noise and gain from a pair of flat field , bias images. By default, the central
      2/8th square of the detector is used for measuring noise and levels.

    Both flat fields have to been observing in the same filter and must have the same exposure level.

    TODO: return actual values, not just print them out
    TODO: gross outlier rejection, e.g., cosmic ray hits

    :param flat1: flat field 1
    :param flat2: flat field 2
    :param bias1: bias 1
    :param bias2: bias 2
    :param minx:
    :param maxx:
    :param miny:
    :param maxy:
    :return:  (gain, readnoise) in e-/ADU, e-
    """

    if minx is None:
        minx = flat1.shape[1] * 8 / 16
    if maxx is None:
        maxx = flat1.shape[1] * 10 / 16
    if miny is None:
        miny = flat1.shape[0] * 8 / 16
    if maxy is None:
        maxy = flat1.shape[0] * 10 / 16

    flat1lvl = np.median(flat1[miny:maxy, minx:maxx])
    flat2lvl = np.median(flat2[miny:maxy, minx:maxx])
    bias1lvl = np.median(bias1[miny:maxy, minx:maxx])
    bias2lvl = np.median(bias2[miny:maxy, minx:maxx])
    _logger.debug(" Flat levels: % 8f % 8f" % (flat1lvl, flat2lvl))
    _logger.debug(" Bias levels: % 8f % 8f" % (bias1lvl, bias2lvl))

    # Measure noise of flat and biad differential iamges
    deltaflat = (flat1 - flat2)[miny:maxy, minx:maxx]
    deltabias = (bias1 - bias2)[miny:maxy, minx:maxx]
    biasnoise = np.std(deltabias)
    biasnoise = np.std(deltabias[np.abs(deltabias - np.median(deltabias)) < 10 * biasnoise])
    flatnoise = np.std(deltaflat)
    flatnoise = np.std(deltaflat[np.abs(deltaflat - np.median(deltaflat)) < 10 * flatnoise])

    flatlevel = (flat1lvl + flat2lvl) / 2 - (bias1lvl + bias2lvl) / 2
    gain = 2 * flatlevel / (flatnoise ** 2 - biasnoise ** 2)
    readnoise = gain * biasnoise / math.sqrt(2)

    if showImages:
        plt.imshow(deltaflat - np.median (deltaflat), clim=(-5 * flatnoise,5 * flatnoise))
        plt.colorbar()
        plt.title ("Delta flat")
        plt.show()
        plt.imshow(deltabias, clim=(-3 * biasnoise,3 * biasnoise))
        plt.colorbar()
        plt.title ("Delta Bias")
        plt.show()
    return (gain, readnoise)


def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='General purpose noise and gain measurement from a set of two flat fields and two biases.')

    parser.add_argument('fitsfile', type=str, nargs=4,
                        help='four fits files:  bias_1 bias_2 flat_1 flat_2')

    parser.add_argument('--imagepath', dest='opt_imagepath', type=str, default=None, help="pathname to prepend to fits file names.")
    parser.add_argument('--log_level', dest='log_level', default='WARN', choices=['DEBUG', 'INFO', 'WARN'],
                        help='Set the debug level')

    parser.add_argument ('--showimages', action='store_true', help="Show difference flat and bias images.")
    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')

    for ii in range(len(args.fitsfile)):
        if args.opt_imagepath is not None:
            args.fitsfile[ii] = "%s/%s" % (args.opt_imagepath, args.fitsfile[ii])
        if not os.path.isfile(args.fitsfile[ii]):
            _logger.error("Fatal: file %s does not exist. Giving up." % (args.fitsfile[ii]))
            sys.exit(0)

    return args


if __name__ == '__main__':

    args = parseCommandLine()

    bias1 = Image(args.fitsfile[0], overscancorrect=True)
    bias2 = Image(args.fitsfile[1], overscancorrect=True)
    flat1 = Image(args.fitsfile[2], overscancorrect=True)
    flat2 = Image(args.fitsfile[3], overscancorrect=True)

    for ii in range(len(flat1.data)):
        (gain, noise) = noisegainextension(flat1.data[ii], flat2.data[ii], bias1.data[ii], bias2.data[ii], showImages=args.showimages)
        print ("Extension %1d  Gain % 5.3f e-/ADU  Noise % 5.2f e-" % (ii, gain, noise))
