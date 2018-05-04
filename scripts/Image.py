from astropy.io import fits
import numpy as np
import re

import logging

_logger = logging.getLogger(__name__)


class Image(object):
    """ Generic class to read in all SCI extensions from a fits file, be it fz compressed or not.

    Code is taken from LCO Banzai pipeline.

    """

    filename = None

    def __init__(self, filename, overscancorrect=False, gaincorrect=False, skycorrect=False):
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
        hdulist = fits.open(filename, 'readonly')
        # Get the main header
        self.header = hdulist[0].header

        # Check for multi-extension fits
        self.extension_headers = []
        self.biassec = []
        self.ccdsec = []
        sci_extensions = self.get_extensions_by_name(hdulist, ['SCI', 'COMPRESSED_IMAGE'])
        if len(sci_extensions) >= 1:
            self.data = np.zeros((len(sci_extensions), sci_extensions[0].data.shape[0],
                                  sci_extensions[0].data.shape[1]), dtype=np.float32)
            for i, hdu in enumerate(sci_extensions):

                self.data[i, :, :] = hdu.data.astype('float32')[:, :]

                gain = 1.
                overscan = 0.
                if (hdu.header['BIASSEC'] == 'UNKNOWN'):
                    bs = None
                else:
                    bs = [int(n) for n in re.split(',|:', hdu.header['BIASSEC'][1:-1])]
                print ("CCDSEC: %s " % hdu.header['DATASEC'])
                if (hdu.header['DATASEC'] is None):
                    cs = [1,hdu.header['NAXIS1'], 1,hdu.header['NAXIS2']]
                else:
                    cs = [int(n) for n in re.split(',|:', hdu.header['DATASEC'][1:-1])]
                self.biassec.append(bs)
                self.ccdsec.append(cs)

                if overscancorrect & (bs is not None):
                    ovpixels = self.data[
                               i, bs[2]:bs[3], bs[0]: bs[1]]
                    overscan = np.median(ovpixels)
                    std = np.std(ovpixels)
                    overscan = np.mean(ovpixels[np.abs(ovpixels - overscan) < 3 * std])

                if gaincorrect:
                    gain = float(hdu.header['GAIN'])
                    hdu.header['GAIN'] = "1.0"

                hdu.header['OVLEVEL'] = overscan
                self.data[i, :, :] = (self.data[i, :, :] - overscan) * gain
                _logger.debug(
                    "Correcting image extension corrected #%d with gain / overscan: % 5.3f % 8.1f" % (
                    i, gain, overscan))

                if skycorrect:
                    imagepixels = self.data[
                                  i, cs[2]: cs[3], cs[0]: cs[1]]
                    skylevel = np.median(imagepixels)
                    std = np.std(imagepixels - skylevel)
                    skylevel = np.mean(imagepixels[np.abs(imagepixels - skylevel) < 5 * std])
                    hdu.header['SKYLEVEL'] = skylevel
                    self.data[i, :, :] = self.data[i, :,:] - skylevel
                    _logger.debug("Sky correct extension #%d with % 8.2f \pm  % 8.2f" % (i, skylevel, std))

                self.extension_headers.append(hdu.header)

        else:

            self.data = hdulist[0].data.astype(np.float32)

        try:
            self.bpm = hdulist['BPM'].data.astype(np.uint8)
        except KeyError:
            self.bpm = None

        hdulist.close()

    def getccddata (self, extension):
        """
        Returns the data as define by DATASEC in ehader
        :param extension:
        :return:
        """
        # TODO: Add safety inspection: is extension in range, is image sub area in image?
        datasec = self.ccdsec[extension]
        sect = self.ccdsec[extension]

        retdata = self.data[extension, sect[2]-1:sect[3]-1, sect[0]-1:sect[1]-1]
        return retdata


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


