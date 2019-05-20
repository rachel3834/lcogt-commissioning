import argparse
import logging
import os
import scipy.ndimage as ndimage
import sys
import numpy as np
from astropy.io import fits
from Image import Image
import matplotlib.pyplot as plt

_logger = logging.getLogger(__name__)

def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='Create a bad pixel mask out of biases, darks, and flat fields')

    parser.add_argument('fitsfiles', type=str, nargs='+',
                        help='Input FITS files for bpm creation.')

    parser.add_argument('--loglevel', dest='log_level', default='INFO', choices=['DEBUG', 'INFO'],
                        help='Set the debug level')

    parser.add_argument('--outputfilename', dest='outputfilename', type=str, default='bpm.fits',
                        help='Outputfilename')

    parser.add_argument ('--showbpmimages', action='store_true', help="Show bpm of each extension for inspection")
    parser.add_argument ('--showstackedinput', action='store_true', help="Show stacked bias, dark, flat images for inspection")
    parser.add_argument ('--biassigma', type=float, default=30, help="rejection limit in sigma above background for bias")

    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')

    for ii in range(len(args.fitsfiles)):

        if not os.path.isfile(args.fitsfiles[ii]):
            _logger.error("Fatal: file %s does not exists" % (args.fitsfiles[ii]))
            sys.exit(0)

    return args


def combineData (listoffiles, extension=0, scale=False):
    """ Stack a set of images by median averaging"""

    _logger.debug("Averaging %d files..." % len(listoffiles))
    ccdsec = listoffiles[0].ccdsec[extension]
    stack = np.zeros( (len(listoffiles),  ccdsec[3]-ccdsec[2]+1,ccdsec[1] - ccdsec[0]+1), dtype=np.float32)
    for ii in range (len(listoffiles)):
        stack[ii] = listoffiles[ii].data[extension]

    averageimage = np.median (stack, axis=0)
    return averageimage


def createbpmFromBiasextension (extdata, sigma=15):
    """
    Build a bad pixel mask based on a bias image extension. Ideally, the input image is already a stack

    The strategy is to find pixels with a significant high level above background.

    """

    # estimate the variance of the images; subtract a smoothed version of image to get rid of gradients
    _logger.debug("Generating median-filter subtracted bias")

    smoothedbias = extdata - ndimage.median_filter(extdata,7)

    # level / variance iteration 1:
    variance = np.std(smoothedbias)
    level = np.median (smoothedbias)
    # now level / variance w/ outlier rejection
    level = np.median (extdata[np.abs(smoothedbias-level) < 5*variance])
    variance = np.std(extdata[np.abs(extdata-level) < 20*variance])

    baddata = np.abs (smoothedbias) > variance * sigma
    extdata *= 0
    extdata[baddata] = 1
    _logger.info ("BPM from Bias: Variance in extension data: % 8.2f +/- %6.3f" % (level,variance))

    return extdata


def createbpmFromDarkextension (extdata, sigma=30):

    variance = np.std(extdata)
    level = np.median (extdata)

    baddata = np.abs (extdata - level) > variance*sigma
    extdata *=0
    extdata[baddata] = 2
    _logger.info ("BPM from Dark: Variance in extension data: % 8.2f +/- %6.3f" % (level,variance))
    return extdata


def createbpmFromFlatextension (extdata):
    level = np.median (extdata)
    variance = np.std(extdata)
    baddata = extdata < level*0.2
    extdata *=0
    extdata[baddata] = 4

    _logger.info ("BPM from Flat: Variance in extension data: % 8.2f +/- %6.3f" % (level,variance))
    return extdata



def showanextenstion(data,title=None, minx=None,maxx=None):
    plt.imshow(data)
    if title is not None:
        plt.title(title)
    plt.show()

if __name__ == '__main__':
    args = parseCommandLine()

    biasfiles = []
    darkfiles = []
    flatfiles = []

    for ii in args.fitsfiles:
        if "-b00" in ii:
           biasfiles.append (Image(ii, overscancorrect=True, trim=False))
        if "-d00" in ii:
            darkfiles.append (Image(ii, overscancorrect=True, trim=False))
        if "-f00" in ii:
            flatfiles.append (Image(ii, overscancorrect=True, trim=False))

    if len (biasfiles) < 3:
        _logger.fatal ("*** There are not enough bias files defined to continue. Giving up.")
        exit(1)

    # we refer all image dimensions etc to the first bias in the list
    referenceImage = biasfiles[0]
    numext = referenceImage.data.shape[0]
    outputdata = np.zeros(referenceImage.data.shape)

    for ext in range(numext):

        extver = referenceImage.extver[ext]
        _logger.info ("Process extension # %d -> extver %s" % (ext+1, extver))
        if ext+1 != extver:
            _logger.warning ("Extension number and extver are out of sync. be aware!")

        extensionaveraged = combineData(biasfiles, extension = ext)
        if args.showstackedinput:
            showanextenstion(extensionaveraged, title="Stacked bias ext #%d" %ext)
        biasbpm = createbpmFromBiasextension(extensionaveraged, sigma = args.biassigma)

        darkbpm = None
        if len (darkfiles) > 3:
            extensionaveraged = combineData(darkfiles, extension = ext)
            darkbpm = createbpmFromDarkextension(extensionaveraged)
        else:
            _logger.info ("Ommitting DARK bpm due to lack of sufficient files")

        flatbpm = None
        if len(flatfiles) > 3:
            extensionaveraged = combineData(flatfiles, extension = ext)
            flatbpm = createbpmFromFlatextension(extensionaveraged)
        else:
            _logger.info ("Ommitting FLAT bpm due to lack of sufficient files")

        # add up all BPM fields.
        outputdata[ext] = (
                biasbpm +
                (flatbpm if flatbpm is not None else 0) +
                (darkbpm if darkbpm is not None else 0)
        ).astype(np.uint8)

        if args.showbpmimages:
            showanextenstion(outputdata[ext], title="Final BPM for ext #%d" %ext)

        nbad = sum (i > 0 for i in outputdata[ext].flatten())
        _logger.info ("Extension %d  has %d identified bad pixels." % (ext,nbad))


    hdul = fits.HDUList ()
    phdu = fits.PrimaryHDU(header=referenceImage.extension_headers[0])
    phdu.header['OBSTYPE'] = 'BPM'
    hdul.append (phdu)

    for ii in range (len (outputdata)):
        _logger.info ("extension %d  detsec %s shape %s" % (ii,referenceImage.extension_headers[ii]['DETSEC'], outputdata[ii].shape ))
        hdu = fits.ImageHDU (outputdata[ii].astype(np.uint8), header=referenceImage.extension_headers[ii])
        hdu.header['EXTNAME'] = 'BPM'
        hdu.header['OBSTYPE'] = 'BPM'
        hdu.header['EXTVER'] = '%d' % (ii+1)
        hdul.append (hdu)

    # BANZAI legacy from pre-fzcompress
    keywordd_to_ext1 = ['SITEID','INSTRUME','CCDSUM','DAY-OBS', 'DATE-OBS']
    for key in keywordd_to_ext1:
        _logger.debug ("Copy keyword %s to extention 1 for BANZAI" % key)
        #hdul[1].header[key] = referenceImage.header[key]

    _logger.info ("Writing output file %s to disk" % args.outputfilename)
    hdul.writeto(args.outputfilename, overwrite=True, output_verify='warn')