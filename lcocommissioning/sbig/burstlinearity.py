import logging

import astropy
import math
import sys
import os
import os.path
import numpy as np
import argparse

from scipy import optimize

from lcocommissioning.common import lco_archive_utilities
from lcocommissioning.common.ccd_noisegain import dosingleLevelGain
from lcocommissioning.common.noisegaindb_orm import NoiseGainMeasurement, noisegaindb
from lcocommissioning.common.Image import Image
import matplotlib.pyplot as plt
from astropy.io import fits

_logger = logging.getLogger(__name__)
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)



sqrtfunc = lambda x, gain, z , exp:  ((gain*x)**exp + z**exp)**(1/exp) - z


def linearityfit (exptime, signaladu, func = sqrtfunc ):

    bounds = [ 0,math.inf]
    (paramset, istat) = optimize.curve_fit (func, exptime, signaladu, bounds = bounds)
    paramerrors = np.sqrt(np.diag(istat)) if istat is not None else None
    _logger.info (f" Paramterset: {paramset} with errors: {paramerrors}")
    return paramset, paramerrors


def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='General purpose CCD noise and gain measurement from pairs of flat fields and biases.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('fitsfiles', type=str, nargs='+',
                        help='Input fits files, must include at least two bias and two flat field frames.')
    parser.add_argument('--camera', default=None)

    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, 'DEBUG'),
                    format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')


    return args


def getfiles(filenames):
    darks = []
    flats = []
    for file in filenames:
        hdu = fits.open(file)
        image = Image(hdu, overscancorrect=True, alreadyopenedhdu=True)
        exptime = float(image.primaryheader['OBJECT'].split()[3])
        if exptime == 0:
            darks.append (image)
        else:
            flats.append (image)
    return darks, flats

def combinedarks (listofdarks):

    shape3d = [len(listofdarks)] + list(listofdarks[0].data[0].shape)

    a = np.zeros(shape3d, dtype=listofdarks[0].data[0].dtype)


    for i, data in enumerate(listofdarks):
        #dark,_,_ = astropy.stats.sigma_clipped_stats(data.data[0][:,-50:-1])
        a[i, :, :] = data.data[0][:, :] #- dark
        #_logger.info (f"dark oversan: {dark}")
    stacked_data = a.sum(axis=0) / len(listofdarks)
    std = np.std(stacked_data)
    m = np.median(stacked_data)
    plt.figure()
    plt.imshow (stacked_data, vmin=m-1*std, vmax = m+1*std)
    plt.savefig('masterdark.png')
    _logger.info(f"Background level of 0 exposure image is {m}")
    return stacked_data


def getlevelforimage (image, masterdark):


    #imageov,_,_ = astropy.stats.sigma_clipped_stats(image[:,-50:-10])
    #_logger.info (f"Oversscan: {imageov}")
    imageov = 0
    zerocorrected = (image-imageov) - (masterdark)

    levell, mean, std = astropy.stats.sigma_clipped_stats (zerocorrected[1100:1500,500:800])
    levelr, mean, std = astropy.stats.sigma_clipped_stats (zerocorrected[1100:1500,-500:-200])
    return levell,levelr

def do_linearity_for_fileset (fitsfiles, args):

    darks,flats = getfiles(fitsfiles)
    masterdark = combinedarks(darks)

    exptimes = []
    levels = []
    levelsr = []
    for flat in flats:
        exptime = float(flat.primaryheader['OBJECT'].split()[3])
        levell,levelr = getlevelforimage(flat.data[0], masterdark)

        exptimes.append (exptime)
        levels.append(levell)

    exptimes = np.asarray(exptimes)
    levels = np.asarray (levels)
    print (exptimes, levels)

    plt.figure()
    plt.plot (exptimes, levels, '.')
    plt.title (args.camera)
    plt.xlabel ('illumination cycles')
    plt.ylabel ('Exposure level [ADU]')


    good = np.asarray(exptimes)>1000
    z = np.polyfit (exptimes[good], levels[good], 1)
    p = np.poly1d(z)
    plt.plot (exptimes, p(exptimes), label = p)
    plt.legend()
    title = args.camera.replace(" ","_") if args.camera is not None else None
    plt.savefig ("exptimelevel.png")
    plt.savefig (f'{title}_exptimelevel.png')
    plt.figure()


    plt.axhline (0, color='grey')
    plt.plot (levels, (levels - p(exptimes)), '.', label="linear fit")


    paramset, paramerror = linearityfit(exptimes, levels, sqrtfunc)

    fitted  = sqrtfunc(exptimes, paramset[0], paramset[1], paramset[2])

    plt.plot (levels, (levels - fitted) , '.' , label = f"(x**k+z**k)**(1/k)-z\nz={paramset[1]:5.1f} k={paramset[2]:5.3f}")

    plt.title (args.camera)
    plt.xlabel("level [ADU]")
    plt.ylabel("(Level - fit) [ADU]")
    plt.ylim([-10,60])
    plt.legend ()


    plt.savefig (f'{title}_exptimelevelresidual.png')





def main():
    args = parseCommandLine()
    do_linearity_for_fileset(args.fitsfiles, args)



if __name__ == '__main__':
    main()