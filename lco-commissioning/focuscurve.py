import sys
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from astropy.io.fits import ImageHDU
from scipy import optimize

from SourceCatalogProvider import SEPSourceCatalogProvider, getImageFWHM

L1FWHM = "L1FWHM"
FOCDMD = "FOCDMD"


# def getImageFWHM(imagename):
#     hdul = fits.open(imagename, 'readonly', ignore_missing_end=True)
#
#     deltaFocus = None
#     for ii in range(len(hdul)):
#         if FOCDMD in hdul[ii].header:
#             deltaFocus = hdul[ii].header[FOCDMD]
#             continue
#
#
#     catalog = SEPSourceCatalogProvider(refineWCSViaLCO=False)
#     fwhmcat = np.asarray([])
#     for ii in range(len(hdul)):
#         if 'EXTNAME' in hdul[ii].header:
#             if 'SCI' in hdul[ii].header['EXTNAME']:
#                 cat, wcs = catalog.get_source_catalog(imagename, ext=ii )
#                 fwhmcat = np.append(fwhmcat, cat['fwhm'])
#     hdul.close()
#     fwhm = np.median(fwhmcat)
#     print(imagename, deltaFocus, fwhm)
#     return deltaFocus, fwhm


def overplotfit(xdata, ydata, func, pinit, label=""):


    for iter in range(2):
        (p1, istat) = optimize.curve_fit(func, xdata, ydata)
        base = np.arange(-3.6, 3.6, 0.1)
        if len (p1) == 3:
            y = func(base, p1[0], p1[1], p1[2])
            fit = func(xdata, p1[0], p1[1], p1[2])


        if len (p1) == 2:
            y = func(base, p1[0], p1[1])
            fit = func(xdata, p1[0], p1[1])

        delta = ydata - fit
        s = np.std (delta)
        print ("Fit rms: {:5.3f}".format (s))
        good = delta < 2 * s
        xdata = xdata[good]
        ydata = ydata[good]



    if 'poly' in label:
        label = "{} curv={:4.2f}".format (label, p1[2])
    if 'sqrt' in label and len(p1)==3:
        label = '{} slope: {:4.2f}'.format (label, p1[2])

    plt.plot(base, y, "--", label=label)
    return p1


polyfit = lambda x, p0, p1, p2: p0 + p1 * x + p2 * x ** 2
polyinit = [2, 0, 1]

sqrtfit = lambda x, p0, p2, p1 : (p0 ** 2 + (p1 * (x - p2)) ** 2)** 0.5
#sqrtfit = lambda x, p0, p2: np.sqrt(p0 ** 2 + (2.5 * (x - p2)) ** 2)
sqrtinit = [2,0]
if __name__ == '__main__':

    focuslist = []
    fwhmlist = []
    for image in sys.argv[1:]:
        focus, fwhm = getImageFWHM(image)

        if fwhm < 15:
            focuslist.append(focus)
            fwhmlist.append(fwhm)

    of = focuslist
    os = fwhmlist
    for jackstart in range(1):
        focuslist = np.asarray(of)[jackstart:]
        fwhmlist = np.asarray(os)[jackstart:]

        plt.figure()
        plt.plot(focuslist, fwhmlist, 'o')
        plt.xlabel("FOCUS Demand [mm foc plane]")
        plt.ylabel("FWHM (Pixels")
        plt.xlim([-3.6, 3.6])
        plt.ylim([0, 14])
        p = overplotfit(focuslist, fwhmlist, polyfit, polyinit, label="poly(2)")
        print("Best poly focus at {:5.2f}  is {:5.2f}; curvature is {:5.2f}".format(-p[1] / (2 * p[2]), p[0] - p[1] ** 2 / 4 / p[2], p[2]))

        p = overplotfit(focuslist, fwhmlist, sqrtfit, sqrtinit, label="sqrt(a+df^2)")
        print("Best sqrt focus at {:5.2f} is {:5.2f}".format(p[1], p[0]))
        plt.axvline(x=p[1], label="best focus sqrt")
        plt.legend()
        plt.title("Sqrt fixed slope best focus found at {:5.2f}".format(p[1]))
        plt.savefig("focus_{}.png".format(jackstart))
