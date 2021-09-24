import json
import logging

import sys
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

from lcocommissioning.common.SourceCatalogProvider import getImageFWHM

L1FWHM = "L1FWHM"
FOCDMD = "FOCDMD"

_LIMIT_EXPONENT_U = 0.7
_MIN_NUMBER_OF_POINTS = 5

# This describes our model for a focus curve: Seeing and defocus add in quadrature.
sqrtfit = lambda x, seeing, bestfocus, slope, tweak: (seeing ** 2 + (slope * (x - bestfocus)) ** 2) ** tweak
polyfit = lambda x, seeing, bestfocus, slope: slope * (x - bestfocus) ** 2 + seeing


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
    initial_guess = [2, 0, 1, 0.6] if func == sqrtfit else None
    bounds = [[0, -3, 0, 0.5], [5, 3, 5, _LIMIT_EXPONENT_U]] if func == sqrtfit else [-math.inf, math.inf]

    for iter in range(2):
        try:
            (paramset, istat) = optimize.curve_fit(func, xdata, ydata, p0=initial_guess, bounds=bounds)
        except:
            paramset = None
            istat = None

        if paramset is None:
            continue
        # sigma rejection and rms estimate:
        fit = func(xdata, *paramset)
        delta = ydata - fit
        s = np.std(delta)
        good = (delta < 3 * s)
        xdata = xdata[good]
        ydata = ydata[good]

    paramerrors = np.sqrt(np.diag(istat)) if istat is not None else None

    return paramset, paramerrors


def overplot_fit(func, paramset):
    if paramset is None:
        return
    base = np.arange(-3.6, 3.6, 0.1)
    y = func(base, *paramset)
    plt.plot(base, y, "--", color='orange' if func == sqrtfit else 'grey',
             label="sqrt {:5.2f}".format(paramset[3]) if func == sqrtfit else "parabola")


def main():
    logging.basicConfig(level=getattr(logging, 'INFO'),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')

    error_string = None
    focuslist = []
    fwhmlist = []
    for image in sys.argv[1:]:
        focus, fwhm = getImageFWHM(image, minarea=5)
        if np.isfinite(fwhm):
            focuslist.append(focus)
            fwhmlist.append(fwhm)

    of = focuslist
    os = fwhmlist
    focuslist = np.asarray(of)
    fwhmlist = np.asarray(os)

    print ("{}\n{}".format (np.round(focuslist,2), np.round (fwhmlist,2)))

    parabola_p, parabola_e = focus_curve_fit(focuslist, fwhmlist, polyfit)
    exponential_p, exponential_rms = focus_curve_fit(focuslist, fwhmlist, sqrtfit)

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
    print (json.dumps(return_package))

    plt.figure()
    if math.isfinite(bestfocus_error):
        plt.axvline(x=bestfocus, color='orange', label="best focus sqrt")
        plt.axes().axvspan(bestfocus - bestfocus_error, bestfocus + bestfocus_error, alpha=0.1, color='grey')

    plt.xlabel("FOCUS Demand [mm foc plane]")
    plt.ylabel("FWHM (arcsec)")

    overplot_fit(polyfit, parabola_p)
    overplot_fit(sqrtfit, exponential_p)
    plt.plot(focuslist, fwhmlist, 'o')
    plt.legend()
    plt.xlim([-3.6, 3.6])
    plt.ylim([0, 6])
    plt.title("Sqrt best focus found at {:5.2f} +/- {:5.2f}".format(bestfocus, bestfocus_error) if math.isfinite(
        bestfocus_error) else "Fit failed")
    plt.savefig("{}".format("focus_0.png"), bbox_inches='tight')


if __name__ == '__main__':
    main()
