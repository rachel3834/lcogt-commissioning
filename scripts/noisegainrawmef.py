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
from astropy.io import fits

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
        minx = (int)  (flat1.shape[1] * 3 / 8 )
    if maxx is None:
        maxx = (int) (flat1.shape[1] * 5 / 8)
    if miny is None:
        miny = (int) (flat1.shape[0] * 3 / 8 )
    if maxy is None:
        maxy = (int) (flat1.shape[0] * 5 / 8)

    flat1lvl = np.median(flat1[miny:maxy, minx:maxx])
    flat2lvl = np.median(flat2[miny:maxy, minx:maxx])
    bias1lvl = np.median(bias1[miny:maxy, minx:maxx])
    bias2lvl = np.median(bias2[miny:maxy, minx:maxx])


    # Measure noise of flat and biad differential iamges
    deltaflat = (flat1 - flat2)[miny:maxy, minx:maxx]
    deltabias = (bias1 - bias2)[miny:maxy, minx:maxx]
    biasnoise = np.std(deltabias)
    biasnoise = np.std(deltabias[np.abs(deltabias - np.median(deltabias)) < 7 * biasnoise])
    flatnoise = np.std(deltaflat)
    flatnoise = np.std(deltaflat[np.abs(deltaflat - np.median(deltaflat)) < 7 * flatnoise])
    _logger.debug(" Levels (flat,flat,bias,bias), and noise (flat, bias): [%d:%d, %d:%d]: % 7.2f, % 7.2f | % 5.2f, % 5.2f | % 7.2f % 5.2f" % (minx,maxx,miny,maxy,flat1lvl, flat2lvl,bias1lvl, bias2lvl, flatnoise,biasnoise))


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
    return (gain, readnoise, flatlevel)


def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='General purpose noise and gain measurement from a set of two flat fields and two biases.')

    parser.add_argument('fitsfile', type=str, nargs='+',
                        help='four fits files:  bias_1 bias_2 flat_1 flat_2')


    parser.add_argument('--minx', type=int, default=None)
    parser.add_argument('--maxx', type=int, default=None)
    parser.add_argument('--miny', type=int, default=None)
    parser.add_argument('--maxy', type=int, default=None)

    parser.add_argument('--imagepath', dest='opt_imagepath', type=str, default=None, help="pathname to prepend to fits file names.")
    parser.add_argument('--log_level', dest='log_level', default='INFO', choices=['DEBUG', 'INFO', 'WARN'],
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



def sortinputfitsfiles (listoffiles, sortby='exptime'):
    """ go through the list of input and sort files by bias and flat type. Find pairs of flats that are of the same exposure time / closest in illumination level."""

    sortedlistofFiles= {}

    filemetrics = {}

    for ii in listoffiles:
        hdu = fits.open (ii)
        if sortby == 'exptime':
            exptime = -1
            # f&*& fpack!
            if 'EXPTIME' in hdu[0].header:
                exptime = hdu[0].header['EXPTIME']
            if 'EXPTIME' in hdu[1].header:
                exptime = hdu[1].header['EXPTIME']
            if exptime > -1:
                filemetrics [ii] = str(exptime)

    for ii in filemetrics:
        _logger.debug ("%s -> %s" % (ii, filemetrics[ii]))


    # find the biases
    for filename in listoffiles:
        if 'b00' in filename:
            # identified a bias exposure
            if 'bias' not in sortedlistofFiles:
                sortedlistofFiles['bias'] = []
            if  len(sortedlistofFiles['bias']) < 2:
                sortedlistofFiles['bias'].append (filename)

    # pair the flat fields
    if sortby == 'exptime':
        unique = set(filemetrics.values())
        #print ("Unique exptimes: %s " % unique)
        for et in unique:
            if float(et) > 0.001:
                sortedlistofFiles[str(et)] = []

        for filename in filemetrics.keys():
            if 'x00' in filename:
                # identified a bias exposure
                #print ("%s, %s" % (filename, filemetrics[filename]))
                if  len( sortedlistofFiles[filemetrics[filename]]) < 2:
                    sortedlistofFiles[filemetrics[filename]].append (filename)

    return sortedlistofFiles


def dosingleLevelGain(fbias1,fbias2, fflat1, fflat2, overscancorrect = False):
    bias1 = Image(fbias1, overscancorrect=overscancorrect)
    bias2 = Image(fbias2, overscancorrect=overscancorrect)
    flat1 = Image(fflat1, overscancorrect=overscancorrect)
    flat2 = Image(fflat2, overscancorrect=overscancorrect)

    print ("Based on:\n Bias 1 %s\n Bias 2 %s\n Flat 1 %s\n Flat 2 %s\n" % (fbias1,fbias2,fflat1,fflat2))

    gains = []
    levels= []
    noises = []

    for ii in range(len(flat1.data)):
        (gain, noise, level) = noisegainextension(flat1.data[ii], flat2.data[ii], bias1.data[ii], bias2.data[ii], showImages=args.showimages,
                                              minx=args.minx, maxx=args.maxx, miny=args.miny, maxy=args.maxy, )


        print ("Extension %1d  Level: % 7.1f  Gain % 5.3f e-/ADU  Noise % 5.2f e-" % (ii, level, gain, noise))

        gains.append( gain)
        levels.append (level)
        noises.append (noise)

    # sanity check on gain and levels:
    retval = (gains,levels, noises)

    gains = gains / gains[0]
    levels = levels/levels[0]
    print ("Sanity checks of relative gain and levels above bias:")
    print ("Relative gains:  ", end="")
    for ii in range(len(gains)):
        print (" % 4.2f" % gains[ii], end="")
    print ("\nRelative levels: ", end="")
    for ii in range(len(levels)):
        print (" % 4.2f" % levels[ii], end="")
    print()
    return retval

def graphresults (alllevels, allgains, allnoises):


    _logger.debug ("Plotting gain vs level")
    plt.figure()
    for ext in alllevels:
        plt.plot ( alllevels[ext], allgains[ext], 'o', label = "extension %s" % ext)

    plt.legend()
    plt.xlabel(("Exposure level [ADU]"))
    plt.ylabel ("Gain [e-/ADU]")
    plt.savefig("levelgain.png")
    plt.close()

    _logger.debug ("Plotting ptc")
    plt.figure()
    print (alllevels)
    for ext in alllevels:
        plt.loglog (alllevels[ext], allnoises[ext], 'o', label = "extension %s" %ext )
    plt.legend()
    plt.xlim([1,64000])
    plt.ylim([1,300])
    plt.xlabel ("Exposure Level [ADU]")
    plt.ylabel ("Noise [ADU]")
    plt.savefig ("ptc.png")
    plt.close()

if __name__ == '__main__':

    args = parseCommandLine()

    sortedinputlist = sortinputfitsfiles(args.fitsfile)

    alllevels = {}
    allgains = {}
    allnoises = {}

    for ii in sortedinputlist:
        if 'bias' not in ii:
            if len (sortedinputlist[ii]) == 2:
                print ("\nNoise / Gain measuremnt based on metric %s" % ii)
                print ("===========================================")
                gains, levels, noises = dosingleLevelGain(sortedinputlist['bias'][0], sortedinputlist['bias'][1],sortedinputlist[ii][0],sortedinputlist[ii][1])

                for ii in range (len(levels)):
                    if ii not in alllevels:
                        alllevels[ii] = []
                        allgains[ii] = []
                        allnoises[ii] = []

                    alllevels[ii].append (levels[ii])
                    allgains[ii].append (gains[ii])
                    allnoises[ii].append(noises[ii])

    graphresults (alllevels, allgains, allnoises)