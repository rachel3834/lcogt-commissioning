import json
import sys
import math
import matplotlib.pyplot as plt
import numpy as np
import scipy
from astropy.io import fits
from astropy.io.fits import ImageHDU, CompImageHDU
from scipy import optimize

from lcocommissioning.common.SourceCatalogProvider import getImageFWHM, SEPSourceCatalogProvider

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

    # TODO: Verification that we have enough points.

    # TODO: Boundaries and intial guess externally driven by context.
    initial_guess = None
    bound = None
    label = None

    if func == sqrtfit:
        initial_guess = [2, 0, 1, 0.6]
        bounds = [[0, -3, 0, 0.5], [5, 3, 5, _LIMIT_EXPONENT_U]]

    if func == fermifit:
        initial_guess = [0, -0.2, 1, -1]
        bounds = [[-np.pi, -np.inf, 0, -2], [+np.pi, np.inf, np.inf, 2]]

    for iter in range(2):
        try:
            (paramset, istat) = optimize.curve_fit(func, xdata, ydata, p0=initial_guess, bounds=bounds)
        except:
            paramset = None
            istat = None
            print("fit error")

        if paramset is None:
            continue
        # sigma rejection and rms estimate:
        fit = func(xdata, *paramset)
        delta = ydata - fit
        s = np.std(delta)
        good = (delta < 1 * s)
        xdata = xdata[good]
        ydata = ydata[good]

    paramerrors = np.sqrt(np.diag(istat)) if istat is not None else None

    return paramset, paramerrors


def overplot_fit(func, paramset):
    if paramset is None:
        return
    base = np.arange(-3.6, 3.6, 0.1)
    y = func(base, *paramset)
    plt.plot(base, y, "--", color='blue' if func == sqrtfit else 'grey',
             label="sqrt {:5.2f}".format(paramset[3]) if func == sqrtfit else "parabola")


def getImageData(imagename, minarea=20, deblend=0.5):
    """ Measure the FWHM of an image, tuned to get a reasonable FWHM also for defocussed images.
    """

    hdul = fits.open(imagename, 'readonly', ignore_missing_end=True)

    deltaFocus = None
    for ii in range(len(hdul)):
        if FOCDMD in hdul[ii].header:
            deltaFocus = hdul[ii].header[FOCDMD]
            continue
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
    goodtheta = thetacat > -10
    goodell = (ellcat > 0 ) & (ellcat <=1)
    meanfwhm = np.mean(fwhmcat[good])

    for iter in range(3):
        medianfwhm = np.median(fwhmcat[good])
        mediantheta = np.median(thetacat[goodtheta])
        medianell = np.median(ellcat[goodell])
        fwhmstd = np.std(fwhmcat[good])
        thetastd = np.std(thetacat[goodtheta])
        ellstd = np.std(ellcat[goodell])

        good = abs(fwhmcat - medianfwhm) < 2 * fwhmstd
        goodtheta = abs(thetacat - mediantheta) < 2 * thetastd
        goodell = abs(ellcat - medianell) < 2 * ellstd

        if np.sum(good) > 10:
            medianfwhm = np.median(fwhmcat[good])
        if np.sum(goodtheta) > 10:
            mediantheta = np.median(thetacat[good])
        if np.sum(goodell) > 10:
            medianell = np.median(ellcat[goodell])

    print("{}  FOCCMD {: 5.3f} FWHM (mean med) ({: 5.2f} {: 5.2f}) \pm {:5.2f} pixel".format(imagename, deltaFocus,
                                                                                             meanfwhm, medianfwhm,
                                                                                             fwhmstd))
    return deltaFocus, medianfwhm, mediantheta, medianell


def main():
    error_string = None
    focuslist = []
    fwhmlist = []
    thetalist = []
    elllist = []

    for image in sys.argv[1:]:
        focus, fwhm, theta, ell = getImageData(image)
        if np.isfinite(fwhm):
            focuslist.append(focus)
            fwhmlist.append(fwhm)
            thetalist.append(theta)
            elllist.append(ell)

    of = focuslist
    os = fwhmlist
    focuslist = np.asarray(of)
    fwhmlist = np.asarray(os)
    thetalist = np.asarray(thetalist)
    elllist = np.asarray(elllist)
    # thetalist = np.mean(thetalist) - thetalist

    print("{}\n{}".format(np.round(focuslist, 2), np.round(fwhmlist, 2)))
    exponential_p, exponential_rms = focus_curve_fit(focuslist, fwhmlist, sqrtfit)
    fermi_p, fermi_rms = focus_curve_fit(focuslist, thetalist, fermifit)

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

        return_package = {'fitok': True if error_string is None else False,
                          'fit_seeing': round(exponential_p[0], 2),
                          'fit_focus': round(bestfocus, 2),
                          'fit_slope': round(exponential_p[2], 2),
                          'fit_exponent': round(exponential_p[3], 2),
                          'fit_rms': round(bestfocus_error, 2),
                          'errormsg': error_string}
    else:
        return_package = None
    print(json.dumps(return_package))

    fig= plt.figure()
    ax1 = fig.add_subplot (111)
    if math.isfinite(bestfocus_error):

        plt.axes().axvspan(bestfocus - bestfocus_error, bestfocus + bestfocus_error, alpha=0.1, color='grey')

    plt.xlabel("FOCUS Demand [mm foc plane]")
    plt.ylabel("FWHM [Pixels]")
    plt.xlim([-3.6, 3.6])
    plt.ylim([0, 30])
    overplot_fit(sqrtfit, exponential_p)
    plt1,=plt.plot(focuslist, fwhmlist, 'o', color="blue")

    ax2 = ax1.twinx()
    plt2,=ax2.plot(focuslist, thetalist, '.', color='grey', label='theta')
    # fermi_p = [  1.41135 ,  0 , 0.06772695 ,-1.34700651 ]
    overplot_fit(fermifit, fermi_p)

    print(fermi_p)
    ax2.set_ylabel("theta [radians]")
    ax2.set_xlim([-3.6, 3.6])

    ax3 = ax1.twinx()
    ax3.spines['right'].set_position(('outward', 60))
    plt3,=ax3.plot(focuslist, elllist, "o", color="lightgreen", label='ellipticity')
    print ("Ellipticity", elllist)
    ax3.set_ylabel ("Ellipticity")
    ax3.set_xlim([-3.6, 3.6])
    ax3.set_ylim([0,1])
    plt.legend()

    ax1.yaxis.label.set_color (plt1.get_color())
    ax2.yaxis.label.set_color (plt2.get_color())
    ax3.yaxis.label.set_color (plt3.get_color())
    plt.title("Sqrt best focus found at {:5.2f} +/- {:5.2f}".format(bestfocus, bestfocus_error) if math.isfinite(
        bestfocus_error) else "Fit failed")
    fig.tight_layout()
    plt.savefig("{}".format("ef_focus.png"))


if __name__ == '__main__':
    main()
