import sys
import os
import numpy as np
import math
import argparse
from Image import Image
import matplotlib.pyplot as plt

import logging
_logger = logging.getLogger(__name__)



def parseCommandLine():
    parser=argparse.ArgumentParser()
    parser.add_argument('fitsfile', type=str, nargs='+')
    parser.add_argument('--ron', default=9)
    parser.add_argument('--log_level', dest='log_level', default='DEBUG', choices=['DEBUG', 'INFO', 'WARN'],
                        help='Set the debug level')
    args = parser.parse_args()



    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')
    return args



def linfactorImagePair (image1, image2):
    pass




def nresBackIllumLinearity (args):
    images = []
    print (args.fitsfile)
    for ii in args.fitsfile:
        _logger.debug ("Loading image %s " % (ii))
        images.append (Image(ii, overscancorrect=True))



    select = (images[0].data[0] < 50000)  & (images[1].data[0] > 1000)

    ratio = images[0].data[0] / images[1].data[0]

    #select = select & (np.abs(ratio) < 10) & ( np.isfinite(ratio))
    expected = np.median (ratio[select])
    select = select & ( np.abs(ratio - expected) < 1)
    # fit

    x_g = images[0].data[0][select].flatten()
    y_g = ratio[select].flatten()

    print ("%f %f %f %f" % (x_g.min(), x_g.max(), y_g.min(), y_g.max()))


    A = np.vstack([x_g**2, x_g,np.ones(len(x_g))]).T
    fit = np.linalg.lstsq (A,y_g)[0]

    #fit = np.polyfit(x_g-x_g.mean(), y_g - y_g.mean(), 1)
    _logger.info ("fit: %s" % (fit))
    fit_fn = np.poly1d(fit)

    expectederror = 3 * np.sqrt (  (images[0].data[0]+9) / images[1].data[0]**2 + images[0].data[0]**2 / images[1].data[0]**3   )[select].flatten()

    plt.figure()
    #plt.plot (images[0].data[0].flatten(), ratio.flatten() - fit_fn(), ',', color='grey')
    plt.plot (x_g, y_g - fit_fn(x_g), ',', color='blue')

    plt.ylim ([-0.5,+0.5])
    print (expectederror)
    _x = range (0,60000,5000)
    plt.plot (x_g, expectederror , ",", color='yellow')

    plt.savefig('nreslinearity.png')



if __name__ == '__main__':

    args = parseCommandLine()

    nresBackIllumLinearity(args)