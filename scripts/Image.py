from astropy.io import fits
import numpy as np


class Image(object):
    """ Generic class to read in all SCI extensions from a fits file, be it fz compressed or not.

    Code is taken from LCO Banzai pipeline.

    """

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

