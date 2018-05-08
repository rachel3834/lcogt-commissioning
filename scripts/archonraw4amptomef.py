import sys
import os
import numpy as np
import math
import argparse
from Image import Image
from astropy.io import fits

import matplotlib.pyplot as plt

import logging

_logger = logging.getLogger(__name__)

def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='Cut archon 4amp output into mef file.')

    parser.add_argument('inputfiles', type=str, nargs="+",
                        help='input files to convert')

    parser.add_argument('--log_level', dest='log_level', default='DEBUG', choices=['DEBUG', 'INFO', 'WARN'],
                        help='Set the debug level')

    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')


    return args


if __name__ == '__main__':

    args = parseCommandLine()
    for imagename in args.inputfiles:

        _logger.info ("Cutting apart: %s " % (imagename))

        inhdu = fits.open(imagename, 'readonly')

        header =  inhdu[1].header
        indata = inhdu[1].data
        _logger.info(indata.shape)
        naxis1 = indata.shape[1]
        naxis1half = int (naxis1 / 2)
        naxis2 = indata.shape[0]
        naxis2half = int (naxis2 / 2)

        # this is in fits corodia te conventions [x1:x2,y1:y2], but pixel count start at zero

        quadrants = [ [         0, naxis1half,           0, naxis2half,0,0],
                      [         0, naxis1half,  naxis2half, naxis2    ,0,1],
                      [naxis1half, naxis1    ,  naxis2half, naxis2    ,1,1],
                      [naxis1half, naxis1    ,            0,naxis2half,1,0],
                      ]
        outhdu = fits.HDUList()
        primary = fits.PrimaryHDU()

        primary.header = header
        outhdu.append(primary)

        for quadrant in quadrants:

            data = indata[ quadrant[2] : quadrant[3], quadrant[0]: quadrant[1]]
            if quadrant[4] == 1:
                    data = np.flip (data,1)
            if quadrant[5] == 1:
                data = np.flip (data,0)


            _logger.debug ("%s %s" % (quadrant, data.shape))
            imext = fits.ImageHDU(data)
          #  imext.header['DATASEC'] = "[%i:%i,%i,%i]" % (18,2048,1,2048)
          #  imext.header['BIASSEC'] = "[%i:%i,%i,%i]" % (data.shape[1]+1-10,data.shape[1]+1, 1, data.shape[0]+1)
            outhdu.append (imext)

        outname = os.path.basename(imagename).replace (".fits.fz", "_mef.fits")
        print ("Writing out to %s " % outname )
        outhdu.writeto (outname, clobber=True)





