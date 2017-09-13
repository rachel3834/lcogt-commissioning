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

class Image(object):

    filename = None


    def __init__(self, filename):
        """
        Load an image from a FITS file

        copypaste from banzai

        Parameters
        ----------
        filename: str
              Full path of the file to open

        Returns
        -------
        data: numpy array
            image data; will have 3 dimensions if the file was either multi-extension or
             a datacube
        header: astropy.io.fits.Header
              Header from the primary extension
        bpm: numpy array
            Array of bad pixel mask values if the BPM extension exists. None otherwise.
            extension_headers: list of astropy.io.fits.Header
                           List of headers from other SCI extensions that are not the
                           primary extension

        Notes
        -----
        The file can be either compressed or not. If there are multiple extensions,
        e.g. Sinistros, the extensions should be (SCI, 1), (SCI, 2), ...
        Sinsitro frames that were taken as datacubes will be munged later so that the
        output images are consistent
        """
        hdulist =  hdulist = fits.open(filename, 'readonly')

        # Get the main header
        self.header = hdulist[0].header

        # Check for multi-extension fits
        self.extension_headers = []
        sci_extensions = self.get_extensions_by_name(hdulist, 'SCI')
        if len(sci_extensions) >= 1:
            self. data = np.zeros((len(sci_extensions), sci_extensions[0].data.shape[0],
                         sci_extensions[0].data.shape[1]), dtype=np.float32)
            for i, hdu in enumerate(sci_extensions):
                self.data[i, :, :] = hdu.data[:, :]
                self.extension_headers.append(hdu.header)

        else:

            self.data = hdulist[0].data.astype(np.float32)

        try:
            self.bpm = hdulist['BPM'].data.astype(np.uint8)
        except KeyError:
            self.bpm = None




    def get_extensions_by_name(self, fits_hdulist, name):
        """
        Get a list of the science extensions from a multi-extension fits file (HDU list)

        Parameters
        ----------
        fits_hdulist: HDUList
                  input fits HDUList to search for SCI extensions

        name: str
          Extension name to collect, e.g. SCI

        Returns
        -------
        HDUList: an HDUList object with only the SCI extensions
        """
        # The following of using False is just an awful convention and will probably be
        # deprecated at some point
        extension_info = fits_hdulist.info(False)
        return fits.HDUList([fits_hdulist[ext[0]] for ext in extension_info if ext[1] == name])





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

def noisegainExtension (flat1, flat2, bias1, bias2, minx=None, maxx=None, miny = None, maxy=None):


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


if __name__ == '__main__':


    args = parseCommandLine()

    bias1 =Image (args.fitsfile[0])
    bias2 =Image (args.fitsfile[1])
    flat1 =Image (args.fitsfile[2])
    flat2 =Image (args.fitsfile[3])


    for ii in range (len (flat1.data)):
        noisegainExtension(flat1.data[ii],flat2.data[ii],bias1.data[ii], bias2.data[ii])