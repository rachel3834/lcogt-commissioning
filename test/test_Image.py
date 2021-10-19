from astropy.io import fits

from lcocommissioning.common.Image import Image
import logging
logging.basicConfig(level=getattr(logging, 'DEBUG'),
                    format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')

def test_texp():
    image = Image ('testdata/coj1m003-fa19-20211018-0010-b00.fits.fz', overscancorrect=True)
    assert image.primaryheader['EXPTIME'] ==0
    image = Image ('testdata/coj1m003-fa19-20211018-0047-f00.fits.fz', overscancorrect=True)
    assert image.primaryheader['EXPTIME'] == 21.137
    image = Image ('testdata/coj0m401-sq01-20211018-0163-b00.fits.fz', overscancorrect=True)
    assert image.primaryheader['EXPTIME'] ==0
    image = Image ('testdata/coj0m401-sq01-20211018-0167-x00.fits.fz', overscancorrect=True)
    assert image.primaryheader['EXPTIME'] ==0.500



    hdu = fits.open ('testdata/coj1m003-fa19-20211018-0010-b00.fits.fz')
    image = Image (hdu, overscancorrect=True, alreadyopenedhdu=True)
    assert image.primaryheader['EXPTIME'] ==0

    hdu = fits.open ('testdata/coj0m401-sq01-20211018-0163-b00.fits.fz')
    image = Image (hdu, overscancorrect=True, alreadyopenedhdu=True)
    assert image.primaryheader['EXPTIME'] ==0

test_texp()
