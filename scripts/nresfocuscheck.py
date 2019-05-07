import numpy as np
import sep
from astropy.io import fits
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

import scipy.linalg

from SourceCatalogProvider import SEPSourceCatalogProvider


def readFile(file):
    sourcecatalogProvider = SEPSourceCatalogProvider(refineWCSViaLCO=False)

    cat, wcs = sourcecatalogProvider.get_source_catalog(
        file)
    return cat


inputs = ['/archive/engineering/lsc/nres01/20181226/raw/lscnrs01-fa09-20181226-0011-a00.fits.fz',
          '/archive/engineering/elp/nres02/20190315/raw/elpnrs02-fa17-20190315-0012-a00.fits.fz',
          '/archive/engineering/cpt/nres03/20190315/raw/cptnrs03-fa13-20190315-0021-a00.fits.fz',
          '/archive/engineering/tlv/nres04/20190313/raw/tlvnrs04-fa18-20190313-0008-a00.fits.fz',]

index = 0
for input in inputs[:]:

    index += 1

    cat1 = readFile(input)
    cat1 = cat1[cat1['x'] > 500]

    plt.hist(cat1['fwhm'], range=(0, 10), histtype='step', bins=40, label='after')
    # plt.hist(cat2['fwhm'], range=(0,10),histtype='step', bins=40, label='after')
    plt.legend()
    plt.savefig("nres_fwhm_histogram{}.png".format(index))
    plt.close()

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    X, Y = np.meshgrid(np.arange(0, 4000, 50), np.arange(0, 4000, 50))
    XX = X.flatten()
    YY = Y.flatten()

    good = (cat1['fwhm'] < 10) & (cat1['fwhm'] > 1)

    for ii in range(5):
        cat1 = cat1[good]

        A = np.c_[cat1['x'], cat1['y'], np.ones(len(cat1))]
        C, _, _, _ = scipy.linalg.lstsq(A, cat1['fwhm'])  # coefficients
        # evaluate it on grid
        Z = C[0] * X + C[1] * Y + C[2]
        print("{} -> {:8.6f} {:8.6f} {:8.6f}".format(len(cat1), C[0], C[1], C[2]))
        fitted = C[0] * cat1['x'] + C[1] * cat1['y'] + C[2]

        good = np.abs(cat1['fwhm'] - fitted) < 1

    ax.scatter(cat1['x'], cat1['y'], cat1['fwhm'], color='c', marker='o', alpha=0.5)

    plt.xlabel('x')
    plt.ylabel('y')
    ax.set_zlabel('FWHM')
    ax.axis('equal')
    ax.axis('tight')
    ax.plot_surface(X, Y, Z, alpha=0.5)
    #plt.show()
    plt.savefig('nresfocus3d_{}.png'.format(index))
    plt.close()

    plt.figure()
    plt.subplot(2, 1, 1)
    plt.plot(cat1['x'], cat1['fwhm'], '.')
    plt.xlabel('X')
    plt.ylim([1,7])
    plt.subplot(2, 1, 2)
    plt.plot(cat1['y'], cat1['fwhm'], '.')
    plt.ylim([1,7])
    plt.xlabel ("Y")
    plt.savefig('nresfocusxy_{}.png'.format(index))
    plt.close()
