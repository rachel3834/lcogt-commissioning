"""
Program to calculate the noise and gain for each extension of a given mef file
"""
import logging
import sys
import os
import os.path
import numpy as np
import math
import argparse

from lcocommissioning.common import common
from lcocommissioning.common.noisegaindbinterface import noisegaindbinterface
from lcocommissioning.common.Image import Image
import matplotlib.pyplot as plt
from astropy.io import fits

_logger = logging.getLogger(__name__)
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)

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
        minx = (int)(flat1.shape[1] * 3 / 8)
    if maxx is None:
        maxx = (int)(flat1.shape[1] * 5 / 8)
    if miny is None:
        miny = (int)(flat1.shape[0] * 3 / 8)
    if maxy is None:
        maxy = (int)(flat1.shape[0] * 5 / 8)

    flat1lvl = np.median(flat1[miny:maxy, minx:maxx])
    flat2lvl = np.median(flat2[miny:maxy, minx:maxx])
    bias1lvl = np.median(bias1[miny:maxy, minx:maxx])
    bias2lvl = np.median(bias2[miny:maxy, minx:maxx])

    avgbiaslevel = 0.5 * (bias2lvl + bias2lvl)

    leveldifference = abs(flat1lvl - flat2lvl)
    avglevel = (flat1lvl - avgbiaslevel + flat2lvl - avgbiaslevel) / 2
    if (leveldifference > avglevel * 0.1):
        _logger.warning("flat level difference % 8f is large compared to level % 8f. Result will be questionable" % (
            leveldifference, flat1lvl - bias1lvl))

    # Measure noise of flat and biad differential iamges
    deltaflat = (flat1 - flat2)[miny:maxy, minx:maxx]
    deltabias = (bias1 - bias2)[miny:maxy, minx:maxx]
    biasnoise = np.std(deltabias)
    biasnoise = np.std(deltabias[np.abs(deltabias - np.median(deltabias)) < 10 * biasnoise])
    flatnoise = np.std(deltaflat)
    flatnoise = np.std(deltaflat[np.abs(deltaflat - np.median(deltaflat)) < 10 * flatnoise])

    _logger.debug(
        " Levels (flat,flat,bias,bias), and noise (flat, bias): [%d:%d, %d:%d]: % 7.2f, % 7.2f | % 5.2f, % 5.2f | % 7.2f % 5.2f" % (
            minx, maxx, miny, maxy, flat1lvl, flat2lvl, bias1lvl, bias2lvl, flatnoise, biasnoise))

    flatlevel = (flat1lvl + flat2lvl) / 2 - (bias1lvl + bias2lvl) / 2
    gain = 2 * flatlevel / (flatnoise ** 2 - biasnoise ** 2)
    readnoise = gain * biasnoise / math.sqrt(2)

    if showImages:
        plt.imshow(deltaflat - np.median(deltaflat), clim=(-5 * flatnoise, 5 * flatnoise))
        plt.colorbar()
        plt.title("Delta flat")
        plt.show()
        plt.imshow(deltabias, clim=(-3 * biasnoise, 3 * biasnoise))
        plt.colorbar()
        plt.title("Delta Bias")
        plt.show()
    return (gain, readnoise, flatlevel, flatnoise, (flat1lvl - avgbiaslevel), (flat2lvl - avgbiaslevel))


def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='General purpose CCD noise and gain measurement from pairs of flat fields and biases.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('fitsfile', type=str, nargs='+',
                        help='Input fits files, must include at least two bias and two flat field frames.')

    group = parser.add_argument_group(
        'Optionally, specify the location of statistics window. All units in pixels with an FITS image extension')
    group.add_argument('--minx', type=int, default=None, help="minimum x.")
    group.add_argument('--maxx', type=int, default=None, help="maximum x.")
    group.add_argument('--miny', type=int, default=None, help="miniumm y.")
    group.add_argument('--maxy', type=int, default=None, help="maximum y.")

    parser.add_argument('--imagepath', dest='opt_imagepath', type=str, default=None,
                        help="pathname to prepend to fits file names.")
    parser.add_argument('--database', default="noisegain.sqlite", help="sqlite database where to store results.")
    parser.add_argument('--sortby', type=str, default="exptime", choices=['exptime', 'filterlevel'],
                        help="Automatically group flat fiel;ds by exposure time (great if using dome flas, or lab flats)."
                             ", or by measured light level (great when using sky flats, but more computing intensive")
    parser.add_argument('--readmode', default="full_frame")
    parser.add_argument('--noreprocessing', action='store_true',
                        help="Do not reprocess if data are already in database")

    parser.add_argument('--loglevel', dest='log_level', default='INFO', choices=['DEBUG', 'INFO', 'WARN'],
                        help='Set the debug level')
    parser.add_argument('--showimages', action='store_true', help="Interactively show difference flat and bias images.")
    parser.add_argument('--makepng', action='store_true', help="Create a png output image of noise, gain, and ptc.")

    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')

    for ii in range(len(args.fitsfile)):
        if args.opt_imagepath is not None:
            args.fitsfile[ii] = "%s/%s" % (args.opt_imagepath, args.fitsfile[ii])
        if not os.path.isfile(args.fitsfile[ii]):
            _logger.error("File %s does not exist. Giving up." % (args.fitsfile[ii]))
            sys.exit(0)

    return args


def findkeywordinhdul(hdulist, keyword):
    # f&*& fpack!
    for ext in hdulist:
        val = ext.header.get(keyword)
        if val is not None:
            return val
    return None


def sortinputfitsfiles(listoffiles, sortby='exptime', selectedreadmode="full_frame"):
    """ go through the list of input and sort files by bias and flat type. Find pairs of flats that are of the same exposure time / closest in illumination level."""

    sortedlistofFiles = {}

    filemetrics = {}

    for filecandidate in listoffiles:

        hdu = fits.open(filecandidate)

        ccdstemp = findkeywordinhdul(hdu, 'CCDSTEMP')
        ccdatemp = findkeywordinhdul(hdu, 'CCDATEMP')
        readoutmode = findkeywordinhdul(hdu, 'CONFMODE')

        if readoutmode not in selectedreadmode:
            _logger.debug ("Rejecting file as it is not in the correct readout mode ({} != {}".format(readoutmode, selectedreadmode))
            hdu.close()
            continue

        tempdiff = 0
        if (ccdstemp is not None) & (ccdatemp is not None):
            tempdiff = float(ccdatemp) - float(ccdstemp)

        if (abs(tempdiff) > 1):
            hdu.close()
            _logger.warning(
                "rejecting file {}: CCD temp is not near set point, delta = {:5.2f}".format(filecandidate, tempdiff))
            continue

        if ('b00' in filecandidate):  # it is a bias

            if 'bias' not in sortedlistofFiles:
                sortedlistofFiles['bias'] = []
            if len(sortedlistofFiles['bias']) < 2:
                sortedlistofFiles['bias'].append(filecandidate)

        else:  # it is a flat

            if (sortby == 'exptime') & (abs(tempdiff) < 0.5):
                exptime = findkeywordinhdul(hdu, 'EXPTIME')
                if exptime is not None:
                    filemetrics[filecandidate] = str(exptime)

            if (sortby == 'filterlevel') & (abs(tempdiff) < 0.5):
                filter = findkeywordinhdul(hdu, "FILTER")
                if (filter is not None) and ('b00' not in filecandidate):
                    image = Image(filecandidate, overscancorrect=True)
                    if image.data is None:
                        level = -1
                    else:
                        level = np.mean(image.data[0])
                    _logger.debug("Input file metrics %s %s %s" % (filecandidate, filter, level))
                    filemetrics[filecandidate] = (filter, level)

        hdu.close()

    if 'bias' not in sortedlistofFiles:
        _logger.fatal("No suitable bias frames found in list!")
        exit(1)

    # pair the flat fields
    if sortby == 'exptime':
        # task is simple: just find flats with the same exposure time
        unique = set(filemetrics.values())

        for et in sorted(unique, key=float):
            if float(et) > 0.001:
                sortedlistofFiles[str(et)] = []

        for filename in filemetrics.keys():
            if ('x00' in filename) or ('f00' in filename):
                # identified a bias exposure
                # print ("%s, %s" % (filename, filemetrics[filename]))
                if len(sortedlistofFiles[filemetrics[filename]]) < 2:
                    sortedlistofFiles[filemetrics[filename]].append(filename)

    if sortby == 'filterlevel':
        tempsortedListofFiles = {}

        for filename in filemetrics.keys():
            (filter, level) = filemetrics[filename]
            if level < 10:
                _logger.debug ("rejecting image {}, level is to low".format (filename))

            if (filter not in tempsortedListofFiles.keys()):
                tempsortedListofFiles[filter] = {}
                ## a new level detected
                tempsortedListofFiles[filter][level] = [filename, ]
                _logger.debug("Starting first entry in a hopeful pair: %s %s %s" % (filename, filter, level))
                continue

            matchfound = False
            for knownfilter in tempsortedListofFiles.keys():
                if (filter == knownfilter):
                    for knownlevel in tempsortedListofFiles[knownfilter].keys():
                        if (0.95 < knownlevel / level) and (knownlevel / level < 1.05):

                            if len(tempsortedListofFiles[knownfilter][knownlevel]) < 2:
                                _logger.debug("adding to set: %s level of %f is within 5 percent of level %f" % (
                                    filename, level, knownlevel))
                                tempsortedListofFiles[knownfilter][knownlevel].append(filename)
                                matchfound = True
                                continue

            if not matchfound:
                tempsortedListofFiles[filter][level] = [filename, ]
                _logger.debug("Starting new level pair with file %s " % (filename))

        for knownfilter in tempsortedListofFiles.keys():
            for knownlevel in tempsortedListofFiles[knownfilter].keys():
                sortedlistofFiles["%s% 8f" % (knownfilter, knownlevel)] = tempsortedListofFiles[knownfilter][knownlevel]

    return sortedlistofFiles


def dosingleLevelGain(fbias1, fbias2, fflat1, fflat2, args, overscancorrect=True):
    bias1 = Image(fbias1, overscancorrect=overscancorrect)
    bias2 = Image(fbias2, overscancorrect=overscancorrect)
    flat1 = Image(fflat1, overscancorrect=overscancorrect)
    flat2 = Image(fflat2, overscancorrect=overscancorrect)

    print("Based on:\n Bias 1 %s\n Bias 2 %s\n Flat 1 %s\n Flat 2 %s\n" % (fbias1, fbias2, fflat1, fflat2))

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


def graphresults(alllevels, allgains, allnoises, allshotnoises, allexptimes):
    _logger.debug("Plotting gain vs level")
    plt.figure()
    for ext in alllevels:
        gains = np.asarray(allgains[ext])
        levels = np.asarray(alllevels[ext])
        statdata = gains[(levels > 5000) & (levels < 55000) & (gains < 20)]
        for iter in range(2):
            mediangain = np.median(statdata)
            stdgain = np.std(statdata)
            goodgains = (np.abs(statdata - mediangain) < 1 * stdgain)
            bestgain = np.mean(statdata[goodgains])

        plt.plot(alllevels[ext], allgains[ext], 'o', label="extension %s data" % (ext))
        plt.hlines(bestgain, 0, 64000, label="Ext %d gain: %5.2f e-/ADU" % (ext, bestgain))
        print("Best gain for ext %d: %5.2f" % (ext, bestgain))

    plt.ylim([2, 20])

    plt.legend()
    plt.xlabel(("Exposure level [ADU]"))
    plt.ylabel("Gain [e-/ADU]")
    plt.savefig("levelgain.png")
    plt.close()

    _logger.debug("Plotting ptc")
    plt.figure()
    # print (alllevels)
    for ext in alllevels:
        plt.loglog(alllevels[ext], allshotnoises[ext], '.', label="extension %s" % ext)
    plt.legend()
    plt.xlim([1, 65000])
    plt.ylim([5, 300])
    plt.xlabel("Exposure Level [ADU]")
    plt.ylabel("Measured Noise [ADU]")
    plt.savefig("ptc.png")
    plt.close()

    _logger.debug("Plotting levle vs exptime")
    plt.figure()
    # print (alllevels)
    for ext in alllevels:
        plt.plot(allexptimes[ext], alllevels[ext], '.', label="extension %s" % ext)
    plt.legend()

    plt.xlabel("Exposure time [s]")
    plt.ylabel("Exposure label [ADU]")
    plt.savefig("texplevel.png")
    plt.close()


def main():
    args = parseCommandLine()
    database = noisegaindbinterface(args.database) if args.database is not None else None

    if (database is not None) and args.noreprocessing:
        for inputname in args.fitsfile:
            if database.checkifalreadyused(os.path.basename(inputname)):
                _logger.info("File %s was already used in a noise measurement. Skipping this entire batch." % inputname)
                exit(0)

    do_noisegain_for_fileset(args.fitsfile, database, args)

    database.close()



def do_noisegain_for_fileset(inputlist, database, args):
    alllevels = {}
    allgains = {}
    allnoises = {}
    allshotnoises = {}
    alllevel1s = {}
    alllevel2s = {}
    allexptimes = {}

    sortedinputlist = sortinputfitsfiles(inputlist, sortby=args.sortby, selectedreadmode=args.readmode)


    for pair_ii in sortedinputlist:

        if 'bias' not in pair_ii:
            if len(sortedinputlist[pair_ii]) == 2:
                print("\nNoise / Gain measuremnt based on metric %s" % pair_ii)
                print("===========================================")

                gains, levels, noises, shotnoises, level1s, level2s, exptimes = dosingleLevelGain(
                    sortedinputlist['bias'][0],
                    sortedinputlist['bias'][1],
                    sortedinputlist[pair_ii][0],
                    sortedinputlist[pair_ii][1], args)

                hdu = fits.open(sortedinputlist[pair_ii][0])
                dateobs = None
                if 'DATE-OBS' in hdu[0].header:
                    dateobs = hdu[0].header['DATE-OBS']
                if 'DATE-OBS' in hdu[1].header:
                    dateobs = hdu[1].header['DATE-OBS']
                camera = None

                if 'INSTRUME' in hdu[0].header:
                    camera = hdu[0].header['INSTRUME']
                if 'INSTRUME' in hdu[1].header:
                    camera = hdu[1].header['INSTRUME']

                filer = None
                if 'FILTER' in hdu[0].header:
                    filter = hdu[0].header['FILTER']
                if 'FILTER' in hdu[1].header:
                    filter = hdu[1].header['FILTER']

                readmode = None
                if 'CONFMODE' in hdu[0].header:
                    filter = hdu[0].header['CONFMODE']
                if 'CONFMODE' in hdu[1].header:
                    filter = hdu[1].header['CONFMODE']

                hdu.close()

                for extension in range(len(levels)):
                    if extension not in alllevels:
                        alllevels[extension] = []
                        allgains[extension] = []
                        allnoises[extension] = []
                        allshotnoises[extension] = []
                        alllevel1s[extension] = []
                        alllevel2s[extension] = []
                        allexptimes[extension] = []

                    if database is not None:
                        database.addmeasurement("%s-%s-%s" % (
                            os.path.basename(sortedinputlist[pair_ii][0]),
                            os.path.basename(sortedinputlist[pair_ii][1]),
                            extension), dateobs, camera, filter,
                                                extension, gains[extension], noises[extension], levels[extension],
                                                shotnoises[extension], level1s[extension], level2s[extension],
                                                readmode=readmode)

                    alllevels[extension].append(levels[extension])
                    allgains[extension].append(gains[extension])
                    allnoises[extension].append(noises[extension])
                    allshotnoises[extension].append(shotnoises[extension])
                    alllevel1s[extension].append(level1s[extension])
                    alllevel2s[extension].append(level2s[extension])
                    allexptimes[extension].append(exptimes[extension])
    if args.makepng:
        graphresults(alllevels, allgains, allnoises, allshotnoises, allexptimes)


if __name__ == '__main__':
    main()
