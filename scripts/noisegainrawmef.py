'''
Program to calculate the noise and gain for each extension of a given mef file
'''

from astropy.io import fits
import os
import numpy as np
import math
from sys import argv, exit
import matplotlib.pyplot as plt
import argparse
from Image import Image





def noisegainExtension (flat1, flat2, bias1, bias2, minx=None, maxx=None, miny = None, maxy=None):
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
    :return:
    """

    if minx is None:
        minx =flat1.shape[1] * 3/ 8
    if maxx is None:
        maxx = flat1.shape[1] * 5/8
    if miny is None:
        miny =flat1.shape[0] * 3/8
    if maxy is None:
        maxy = flat1.shape[0] * 5/8

    flat1lvl = np.mean (flat1[miny:maxy,minx:maxx])
    flat2lvl = np.mean (flat2[miny:maxy,minx:maxx])
    bias1lvl = np.mean (bias1[miny:maxy,minx:maxx])
    bias2lvl = np.mean (bias2[miny:maxy,minx:maxx])
    print (" Flat levels: % 8f % 8f" % (flat1lvl, flat2lvl))
    print (" Bias levels: % 8f % 8f" % (bias1lvl, bias2lvl))

    # plt.imshow(bias1[miny:maxy,minx:maxx], clim=(1000,1100))
    # plt.colorbar()
    # plt.show()
    deltaflat = flat1 - flat2
    deltabias = bias1 - bias2
    biasnoise = np.std (deltabias[miny:maxy,minx:maxx])
    flatnoise = np.std (deltaflat[miny:maxy,minx:maxx])

    flatlevel = (flat1lvl + flat2lvl) / 2 - (bias1lvl + bias2lvl) / 2

    gain = 2 * flatlevel / ( flatnoise**2 - biasnoise**2)
    readnoise = gain * biasnoise / math.sqrt(2)

    print "  Gain  [e-/ADU] : % 5.3f " % (gain)
    print "  Noise [e-]     : % 5.3f " % (readnoise)



def parseCommandLine():

    parser = argparse.ArgumentParser(
        description='General purpose noise and gain measurement from a set of two flat fields and two biases.')

    parser.add_argument('fitsfile', type=str, nargs=4,
                        help='four fits files:  bias_1 bias_2 flat_1 flat_2')

    parser.add_argument('--imagepath', dest = 'opt_imagepath' , type = str, default = None)
    args = parser.parse_args()

    for ii in range (len (args.fitsfile)):
        if args.opt_imagepath is not None:
            args.fitsfile[ii] = "%s/%s" % (args.opt_imagepath, args.fitsfile[ii])
        if not os.path.isfile( args.fitsfile[ii]):
            print "Fatal: file %s does not exists" % (args.fitsfile[ii])
            os.exit(0)

    return args


if __name__ == '__main__':


    args = parseCommandLine()

    bias1 =Image (args.fitsfile[0])
    bias2 =Image (args.fitsfile[1])
    flat1 =Image (args.fitsfile[2])
    flat2 =Image (args.fitsfile[3])

    for ii in range (len (flat1.data)):
        noisegainExtension(flat1.data[ii],flat2.data[ii],bias1.data[ii], bias2.data[ii])