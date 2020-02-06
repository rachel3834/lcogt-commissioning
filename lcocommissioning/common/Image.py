import logging
import re

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

_log = logging.getLogger(__name__)

class Image(object):
    """ Generic class to read in all SCI extensions from a fits file, be it fz compressed or not.

    Code is taken from LCO Banzai pipeline.

    """

    filename = None

    def __init__(self, filename, alreadyopenedhdu=False , overscancorrect=False, gaincorrect=False, skycorrect=False, trim=True, minx=None, maxx=None,miny=None,maxy=None):
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
        if alreadyopenedhdu:
            hdulist = filename
        else:
            hdulist = fits.open(filename, 'readonly')

        # Get the main header
        self.primaryheader = hdulist[0].header
        if alreadyopenedhdu or filename.endswith(".fz") :
            for card in hdulist[1].header:
                self.primaryheader.append (card)

        # Check for multi-extension fits
        self.extension_headers = []
        self.ccdsec = []
        self.ccdsum = []
        self.extver = []

        sci_extensions = self.get_extensions_by_name(hdulist, ['SCI', 'COMPRESSED_IMAGE', 'SPECTRUM'])

        if len (sci_extensions) == 0:
            _log.warning ("No SCI extenstion found in image %s. Aborting." % filename)
            self.data = None
            return None

        # Find out whow big dat are. Warning: assumption is that all extensions have same dimensions.
        datasec = sci_extensions[0].header.get('DATASEC')
        _log.debug("DATASEC: {}".format(datasec ))

        if ( (datasec is None) or not trim):
            _log.info ("Not trimming")
            cs = [1,sci_extensions[0].header['NAXIS1'], 1,sci_extensions[0].header['NAXIS2']]
        else:
            cs = self.fitssection_to_slice(datasec)

        if len(sci_extensions) >0:
            # Generate internal storage array for pre-processed data
            self.data = np.zeros( ( len(sci_extensions), cs[3]-cs[2] + 1 ,cs[1] - cs[0] + 1), dtype=np.float32)

            for i, hdu in enumerate(sci_extensions):

                gain = 1.
                overscan = 0.
                extver = hdu.header.get('EXTVER', str(i+1))
                self.extver.append (extver)

                self.ccdsec.append(cs)

                if overscancorrect:
                    overscan = self.get_overscan_from_hdu(hdu)

                if gaincorrect:
                    gain = float(hdu.header.get('GAIN', "1.0"))
                    hdu.header['GAIN'] = "1.0"

                hdu.header['OVLEVEL'] = overscan

                self.data[i, :, :] = (hdu.data[cs[2]-1:cs[3],cs[0]-1:cs[1]] - overscan) * gain

                skylevel = 0
                if skycorrect:
                    imagepixels = self.data[
                                  i, 20:-20,20:-20]
                    skylevel = np.median(imagepixels)
                    std = np.std(imagepixels - skylevel)
                    skylevel = np.median(imagepixels[np.abs(imagepixels - skylevel) < 3 * std])
                    hdu.header['SKYLEVEL'] = skylevel
                    self.data[i, :, :] = self.data[i, :,:] - skylevel

                _log.debug(
                    "Correcting image extension #%d with gain / overscan / sky: % 5.3f % 8.1f  % 8.2f" % (
                        i, gain, overscan, skylevel))

                self.extension_headers.append(hdu.header)

        else:

            self.data = hdulist[0].data.astype(np.float32)

        try:
            self.bpm = hdulist['BPM'].data.astype(np.uint8)
        except KeyError:
            self.bpm = None

        if not alreadyopenedhdu:
            hdulist.close()

    def getccddata (self, extension):
        """
        Returns the data as define by DATASEC in ehader
        :param extension:
        :return:
        """
        retdata = self.data[extension]#, sect[2]-1:sect[3]-1, sect[0]-1:sect[1]-1]
        return retdata


    def get_overscan_from_hdu (self, hdu, sig_rej=2, biassecheader='BIASSEC'):
        """ Calculate the overscan level of an image extension.
            Calculation is based on slice defined by header keyword.
        """

        overkeyword = hdu.header.get(biassecheader)
        if overkeyword is None:
            return 0

        biassecslice = self.fitssection_to_slice(overkeyword)
        ovpixels = hdu.data[
                   biassecslice[2]+1:biassecslice[3]-1, biassecslice[0]+1: biassecslice[1]-1]
        overscanlevel = np.median(ovpixels)
        std = np.std(ovpixels)
        overscanlevel = np.mean(ovpixels[np.abs(ovpixels - overscanlevel) < sig_rej * std])
        return overscanlevel

    def fitssection_to_slice (self, keyword):
        integers =  [int(n) for n in re.split(',|:', keyword[1:-1])]
        return integers

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
        return fits.HDUList([fits_hdulist[ext[0]] for ext in extension_info if ext[1] in name])



