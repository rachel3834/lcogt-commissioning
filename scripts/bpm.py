import argparse
import logging
import os
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

    parser.add_argument('--log_level', dest='log_level', default='INFO', choices=['DEBUG', 'INFO'],
                        help='Set the debug level')

    parser.add_argument('--outputfilename', dest='outputfilename', type=str, default='bpm.fits',
                        help='Outputfilename')

    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')

    for ii in range(len(args.fitsfiles)):

        if not os.path.isfile(args.fitsfiles[ii]):
            _logger.error("Fatal: file %s does not exists" % (args.fitsfiles[ii]))
            sys.exit(0)

    return args


def combineData (listoffiles, extension=0, scale=False):

    _logger.debug("Averaging %d files..." % len(listoffiles))
    ccdsec = listoffiles[0].ccdsec[extension]
    stack = np.zeros( (len(listoffiles),  ccdsec[3]-ccdsec[2]+1,ccdsec[1] - ccdsec[0]+1), dtype=np.float32)
    for ii in range (len(listoffiles)):
        stack[ii] = listoffiles[ii].data[extension]
    averageimage = np.median (stack, axis=0)

    return averageimage


def createbpmFromBiasextension (extdata):

    variance = np.std(extdata)
    level = np.median (extdata)

    baddata = np.abs (extdata - level) > variance*20
    extdata *=0
    extdata[baddata] = 1
    _logger.info ("BPM from Bias: Variance in extension data: % 8.2f +/- %6.3f" % (level,variance))
    #plt.imshow (extdata, clim=(0,1))
    #plt.show();

    return extdata


def createbpmFromFlatextension (extdata):
    level = np.median (extdata)
    variance = np.std(extdata)
    baddata = extdata < level*0.5
    extdata *=0
    extdata[baddata] = 1
    #plt.imshow (extdata, clim=(0,1))
    #plt.show();
    _logger.info ("BPM from Flat: Variance in extension data: % 8.2f +/- %6.3f" % (level,variance))

    return extdata


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
        _logger.error ("There are not enough bias files defined to continue. Giving up.")
        exit(1)

    if len (flatfiles) < 3:
        _logger.error ("There are not enough flat files defined to continue. Giving up.")
        exit(1)


    numext = biasfiles[0].data.shape[0]
    outputdata = np.zeros(biasfiles[0].data.shape)

    for ext in range(numext):


        _logger.info("Process extension %d of %d" % (ext, numext))

        extensionaveraged = combineData(biasfiles, extension = ext)
        biasbpm = createbpmFromBiasextension(extensionaveraged)

        extensionaveraged = combineData(flatfiles, extension = ext)
        flatbpm = createbpmFromFlatextension(extensionaveraged)


        darkbpm = None
        if len (darkfiles) > 3:
            extensionaveraged = combineData(darkfiles, extension = ext)
            darkbpm = createbpmFromBiasextension(extensionaveraged)
        else:
            _logger.info ("Ommitting dark bpm due to lack of sufficient files")

        outputdata[ext] = biasbpm + flatbpm + (darkbpm if darkbpm is not None else 0)

        #plt.imshow (outputdata[ext], clim=(0,1))
        #plt.show();


    hdu = fits.PrimaryHDU(outputdata)
    hdul = fits.HDUList ([hdu])
    hdu.writeto(args.outputfilename)