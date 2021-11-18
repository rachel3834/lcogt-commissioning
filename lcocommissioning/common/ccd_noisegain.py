import math
import matplotlib
import numpy as np
from astropy.io.fits import HDUList
from astropy.stats import sigma_clipped_stats
from matplotlib import pyplot as plt
from lcocommissioning.common.Image import Image
import logging

_logger = logging.getLogger(__name__)

# Utils to calculate CCD noise and gain from two flats and two biases.
# This module is concerened only with astropy.fits objects, and not their generation.


def noisegainextension(flat1, flat2, bias1, bias2, minx=None, maxx=None, miny=None, maxy=None, showImages=False):
    """
    Measure the noise and gain from a pair of flat field , bias images. By default, the central
      2/8th square of the detector is used for measuring noise and levels.

    Both flat fields have to been observing in the same filter and must have the same exposure level.

    TODO: return actual values, not just print them out
    TODO: gross outlier rejection, e.g., cosmic ray hits

    :param flat1: flat field 1  as numpy array
    :param flat2: flat field 2 as numpy array
    :param bias1: bias 1 as numpy array
    :param bias2: bias 2 as numpy array
    :param minx:
    :param maxx:
    :param miny:
    :param maxy:
    :return:  (gain, readnoise) in e-/ADU, e-
    """

    if minx is None:
        minx = (int)(flat1.shape[1] * 3 / 8)
    if maxx is None:
        maxx = (int)(flat1.shape[1] * 5 / 8)
    if miny is None:
        miny = (int)(flat1.shape[0] * 3 / 8)
    if maxy is None:
        maxy = (int)(flat1.shape[0] * 5 / 8)

    _,flat1lvl,_ = sigma_clipped_stats(flat1[miny:maxy, minx:maxx], sigma=5)
    _,flat2lvl,_ = sigma_clipped_stats(flat2[miny:maxy, minx:maxx], sigma=5)
    _,bias1lvl,_ = sigma_clipped_stats(bias1[miny:maxy, minx:maxx], sigma=3)
    _,bias2lvl,_ = sigma_clipped_stats(bias2[miny:maxy, minx:maxx], sigma=3)

    avgbiaslevel = 0.5 * (bias2lvl + bias2lvl)
    leveldifference = abs(flat1lvl - flat2lvl)
    avglevel = (flat1lvl - avgbiaslevel + flat2lvl - avgbiaslevel) / 2
    if (leveldifference > avglevel * 0.1):
        _logger.warning("flat level difference % 8f is large compared to level % 8f. Result will be questionable" % (
            leveldifference, flat1lvl - bias1lvl))
    # Measure noise of flat and bias differential images
    deltaflat = (flat1 - flat2)[miny:maxy, minx:maxx]
    deltabias = (bias1 - bias2)[miny:maxy, minx:maxx]
    #biasnoise = np.std(deltabias)
    #biasnoise = np.std(deltabias[np.abs(deltabias - np.median(deltabias)) < 10 * biasnoise])
    _,_,biasnoise = sigma_clipped_stats(deltabias, sigma=5)
    #flatnoise = np.std(deltaflat)
    #flatnoise = np.std(deltaflat[np.abs(deltaflat - np.median(deltaflat)) < 10 * flatnoise])
    _,_,flatnoise = sigma_clipped_stats(deltaflat, sigma=10)
    _logger.debug(
        " Levels (flat,flat,bias,bias), and noise (flat, bias): [%d:%d, %d:%d]: % 7.2f, % 7.2f | % 5.2f, % 5.2f | % 7.2f % 5.2f" % (
            minx, maxx, miny, maxy, flat1lvl, flat2lvl, bias1lvl, bias2lvl, flatnoise, biasnoise))

    flatlevel = (flat1lvl + flat2lvl) / 2 - (bias1lvl + bias2lvl) / 2
    gain = 2 * flatlevel / (flatnoise ** 2 - biasnoise ** 2)
    readnoise = gain * biasnoise / math.sqrt(2)

    if showImages:
        plt.switch_backend ('TkAgg')
        print ("Showing images", plt.get_backend())
        plt.imshow(deltaflat - np.median(deltaflat), clim=(-5 * flatnoise, 5 * flatnoise))
        plt.colorbar()
        plt.title("Delta flat")
        plt.show()
        plt.imshow(deltabias, clim=(-3 * biasnoise, 3 * biasnoise))
        plt.colorbar()
        plt.title("Delta Bias")
        plt.show()

    return (gain, readnoise, flatlevel, flatnoise, (flat1lvl - avgbiaslevel), (flat2lvl - avgbiaslevel))


def dosingleLevelGain(fbias1: HDUList, fbias2: HDUList, fflat1: HDUList, fflat2: HDUList, args, overscancorrect=True):
    """
    Calculate for each extension the noise and gain and print the result to console.

    :param fbias1: astropy.fits object
    :param fbias2:
    :param fflat1:
    :param fflat2:
    :param args:
    :param overscancorrect:
    :param alreadyopenhdu: if True, treat the input files as FITS hdu. Otherwise, trweat them as filenames.
    :return:
    """

    bias1 = Image(fbias1, overscancorrect=overscancorrect, alreadyopenedhdu=True)
    bias2 = Image(fbias2, overscancorrect=overscancorrect, alreadyopenedhdu=True)
    flat1 = Image(fflat1, overscancorrect=overscancorrect, alreadyopenedhdu=True)
    flat2 = Image(fflat2, overscancorrect=overscancorrect, alreadyopenedhdu=True)

    gains = []
    levels = []
    noises = []
    shotnoises = []
    level1s = []
    level2s = []
    exptimes = []

    for ii in range(len(flat1.data)):
        (gain, noise, level, shotnoise, level1, level2) = noisegainextension(flat1.data[ii], flat2.data[ii],
                                                                             bias1.data[ii], bias2.data[ii],
                                                                             showImages=args.showimages,
                                                                             minx=args.minx, maxx=args.maxx,
                                                                             miny=args.miny, maxy=args.maxy, )

        print("Extension %1d  Level: % 7.1f  Gain % 5.3f e-/ADU  Noise % 5.2f e-" % (ii, level, gain, noise))

        gains.append(gain)
        levels.append(level)
        noises.append(noise)
        shotnoises.append(shotnoise)
        level1s.append(level1)
        level2s.append(level2)
        exptimes.append(flat1.primaryheader['EXPTIME'])

    # sanity check on gain and levels:
    retval = (gains, levels, noises, shotnoises, level1s, level2s, exptimes)

    # Do some pretty console output with useful statistics.
    gains = gains / gains[0]
    levels = levels / levels[0]
    print("Sanity checks of relative gain and levels above bias:")
    print("Relative gains:  ", end="")
    for ii in range(len(gains)):
        print(" % 4.2f" % gains[ii], end="")
    print("\nRelative levels: ", end="")
    for ii in range(len(levels)):
        print(" % 4.2f" % levels[ii], end="")
    print()
    return retval
