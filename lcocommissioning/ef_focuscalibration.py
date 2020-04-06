import logging
import math
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.io.fits import ImageHDU, CompImageHDU
from scipy import optimize

from lcocommissioning.common.SourceCatalogProvider import SEPSourceCatalogProvider

from matplotlib import rc
rc('text', usetex=True)

_log = logging.getLogger(__name__)

L1FWHM = "L1FWHM"
FOCDMD = "FOCDMD"

_LIMIT_EXPONENT_U = 0.7
_MIN_NUMBER_OF_POINTS = 5

# This describes our model for a focus curve: Seeing and defocus add in quadrature.
sqrtfit = lambda x, seeing, bestfocus, slope, tweak: (seeing ** 2 + (slope * (x - bestfocus)) ** 2) ** tweak
fermifit = lambda x, thetazero, x0, T, a: a / (1 + np.e ** ((x - x0) / T)) + thetazero


def focus_curve_fit(xdata, ydata, func=sqrtfit):
    """
    Generic iterative fit with sigma rejection.

    :param xdata:
    :param ydata:
    :param func:
    :param plot:
    :return:
    """

    # TODO: Boundaries and initial guess externally driven by context.
    initial_guess = None
    bound = None
    label = None

    if func == sqrtfit:
        initial_guess = [2, 0, 1, 0.6]
        bounds = [[0, -3, 0, 0.5], [5, 3, 5, _LIMIT_EXPONENT_U]]

    if func == fermifit:
        initial_guess = [0, 0, 100, 0]
        ydata = ydata * np.pi / 180.
        bounds = [[-np.pi, -np.inf, 0, -np.pi], [+np.pi, np.inf, np.inf, np.pi]]

    for iter in range(2):
        try:
            (paramset, istat) = optimize.curve_fit(func, xdata, ydata, p0=initial_guess, bounds=bounds)
        except:
            paramset = None
            istat = None
            _log.debug("fit error")

        if paramset is None:
            continue
        # sigma rejection and rms estimate:
        fit = func(xdata, *paramset)
        delta = ydata - fit
        s = np.std(delta)
        good = (delta < 2 * s)
        xdata = xdata[good]
        ydata = ydata[good]

    paramerrors = np.sqrt(np.diag(istat)) if istat is not None else None
    if (func == fermifit) & (paramset is not None):
        paramset[3] *= 180./np.pi
        paramset[0] *= 180./np.pi

    return paramset, paramerrors


def overplot_fit(func, paramset, color=None):
    if paramset is None:
        _log.info ("Not plotting function since plot parameters are None")
        return
    base = np.arange(-3.6, 3.6, 0.1)
    y = func(base, *paramset)
    if color is None:
        color = 'blue' if func == sqrtfit else 'grey'
    plt.plot(base, y, "--", color=color,
             label="sqrt {:5.2f}".format(paramset[3]) if func == sqrtfit else "parabola")


def getImageData(imagename, minarea=20, deblend=0.5):
    """ Measure the FWHM of an image, tuned to get a reasonable FWHM also for defocussed images.
    """

    hdul = fits.open(imagename, 'readonly', ignore_missing_end=True)

    deltaFocus = None
    pixelscale = None
    for ii in range(len(hdul)):
        if FOCDMD in hdul[ii].header:
            deltaFocus = hdul[ii].header[FOCDMD]
        if 'PIXSCALE' in hdul[ii].header:
            pixelscale = hdul[ii].header['PIXSCALE']

    catalog = SEPSourceCatalogProvider(refineWCSViaLCO=False)
    fwhmcat = np.asarray([])
    thetacat = np.asarray([])
    ellcat = np.asarray([])
    for ii in range(len(hdul)):
        if isinstance(hdul[ii], ImageHDU) or isinstance(hdul[ii], CompImageHDU):
            cat, wcs = catalog.get_source_catalog(imagename, ext=ii, minarea=minarea, deblend=deblend)
            fwhmcat = np.append(fwhmcat, cat['fwhm'])
            thetacat = np.append(thetacat, cat['theta'])
            ellcat = np.append(thetacat, cat['ellipticity'])
    hdul.close()

    # comprehension of the object catalog....
    good = fwhmcat > 0
    goodtheta = thetacat > -720
    _log.debug("Theta", np.min(thetacat), np.max(thetacat), np.mean(thetacat))
    goodell = (ellcat > 0) & (ellcat <= 1)
    meanfwhm = np.mean(fwhmcat[good])

    for iter in range(3):
        medianfwhm = np.median(fwhmcat[good])
        mediantheta = np.median(thetacat[goodtheta])
        medianell = np.median(ellcat[goodell])
        fwhmstd = np.std(fwhmcat[good])
        thetastd = np.std(thetacat[goodtheta])
        ellstd = np.std(ellcat[goodell])

        good = abs(fwhmcat - medianfwhm) < 2 * fwhmstd
        goodtheta = abs(thetacat - mediantheta) < 3 * thetastd
        goodell = abs(ellcat - medianell) < 2 * ellstd

        if np.sum(good) > 10:
            medianfwhm = np.median(fwhmcat[good])
        if np.sum(goodtheta) > 10:
            mediantheta = np.median(thetacat[good])
        if np.sum(goodell) > 10:
            medianell = np.median(ellcat[goodell])


    _log.debug ("{}  FOCCMD {: 5.3f} FWHM (mean med) ({: 5.2f} {: 5.2f}) \pm {:5.2f} pixel  {:5.2f}".format(imagename, deltaFocus,
                                                                                             meanfwhm, medianfwhm,
                                                                                             fwhmstd, mediantheta))

    if pixelscale is not None:
        medianfwhm *= pixelscale
    else:
        _log.warning("Pixel scale was not defined!")
    return deltaFocus, medianfwhm, mediantheta, medianell


def main():
    logging.basicConfig(level=getattr(logging, 'INFO'),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')

    focuslist = []
    fwhmlist = []
    thetalist = []
    elllist = []

    efimages = [i for i in sys.argv[1:] if ("ef" in os.path.basename(i))]
    faimages = [i for i in sys.argv[1:] if ("fa" in os.path.basename(i))]

    _log.debug("EF images:", efimages)
    _log.debug("FA images:", faimages)


    for image in efimages:
        focus, fwhm, theta, ell = getImageData(image)
        if np.isfinite(fwhm):
            focuslist.append(focus)
            fwhmlist.append(fwhm)
            thetalist.append(theta)
            elllist.append(ell)

    fafwhmlist = []
    fafocuslist = []

    for image in faimages:
        focus, fwhm, theta, ell = getImageData(image)
        fafwhmlist.append(fwhm)
        fafocuslist.append(focus)

    focuslist = np.asarray(focuslist)
    fwhmlist = np.asarray(fwhmlist)
    thetalist = np.asarray(thetalist)
    elllist = np.asarray(elllist)

    fafwhmlist = np.asarray(fafwhmlist)
    fafocuslist = np.asarray(fafocuslist)

    exponential_p, exponential_rms = focus_curve_fit(focuslist, fwhmlist, sqrtfit)
    fermi_p, fermi_rms = focus_curve_fit(focuslist, thetalist, fermifit)

    faexponential_p, faexponential_rms = focus_curve_fit(fafocuslist, fafwhmlist, sqrtfit)

    # we will need this a few times - meaningful references here
    if exponential_p is not None:
        bestfocus = exponential_p[1]
        bestfocus_error = exponential_rms[1]

        if not math.isfinite(bestfocus_error):
            error_string = "fit did not converge"
        if bestfocus_error > 0.25:
            error_string = "focus fit is too noisy"
        if abs(exponential_p[1]) > 2.5:
            error_string = "Focus offset too large to be credible."


    bestfocusfa = faexponential_p[1]
    bestfocusfa_error = faexponential_rms[1]

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    if math.isfinite(bestfocus_error):
        plt.axes().axvspan(bestfocus - bestfocus_error, bestfocus + bestfocus_error, alpha=0.1, color='grey')

    plt.xlabel("FOCUS Demand [mm foc plane]")
    plt.ylabel("FWHM ['']")
    plt.xlim([-2.6, 2.6])
    plt.ylim([0.5, 6])
    overplot_fit(sqrtfit, exponential_p)
    plt1, = plt.plot(focuslist, fwhmlist, 'o', color="blue", label="EF FWWM")
    overplot_fit(sqrtfit, faexponential_p, color="lightblue")
    plt1fa, = plt.plot(fafocuslist, fafwhmlist, '.', color="lightblue", label="FA FWHM")

    # Plot the position angle and present a nice fermi fit
    ax2 = ax1.twinx()
    plt2, = ax2.plot(focuslist, thetalist, '.', color='grey', label='EF theta')
    # fermi_p = [  1.41135 ,  0 , 0.06772695 ,-1.34700651 ]
    overplot_fit(fermifit, fermi_p)

    _log.info(f"fermi fit: {fermi_p}")
    ax2.set_ylabel("$\Theta$ [$^\circ$]")
    ax2.set_xlim([-2.6, 2.6])

    # Plot the ellipticity.
    ax3 = ax1.twinx()
    ax3.spines['right'].set_position(('outward', 60))
    plt3, = ax3.plot(focuslist, elllist, "o", color="lightgreen", label='EF ellipticity')

    ax3.set_ylabel("Ellipticity")
    ax3.set_xlim([-2.6, 2.6])
    ax3.set_ylim([0, 1])

    ax1.legend(handles=[plt1, plt1fa, plt2, plt3], loc="lower right", bbox_to_anchor=(1, -0.1),
               bbox_transform=fig.transFigure, ncol=4)

    ax1.yaxis.label.set_color(plt1.get_color())
    ax2.yaxis.label.set_color(plt2.get_color())
    ax3.yaxis.label.set_color(plt3.get_color())
    plt.title("Best EF focus at {:5.2f} +/- {:5.2f} mm\nBest FA focus at {:5.2f} +/- {:5.2f} mm".format(bestfocus,
                                                                                                        bestfocus_error,
                                                                                                        bestfocusfa,
                                                                                                        bestfocusfa_error) if math.isfinite(
        bestfocus_error) else "Fit failed")

    plt.savefig("{}".format("ef_focus.png"), bbox_inches="tight")
    _log.info ("All done")

if __name__ == '__main__':
    main()
