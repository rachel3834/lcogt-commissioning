
import json
import sys
import math
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy import optimize

from lcocommissioning.common.SourceCatalogProvider import SEPSourceCatalogProvider

''' make  1D and 2D plots of focus / ellipticity / orientaion distriobitions

'''


def plotdistribution(imagecatalog, name=None):


    medianfwhm = np.median (imagecatalog['fwhm'])
    std = np.std (imagecatalog['fwhm'])
    medianfwhm = np.median (imagecatalog['fwhm'][ np.abs (imagecatalog['fwhm'] - medianfwhm) < std*2])
    good =  np.abs (imagecatalog['fwhm'] - medianfwhm) < std*2


    fig, ax = plt.subplots()
    ax.set_box_aspect(1)
    plt.scatter (x=imagecatalog['x'][good], y=imagecatalog['y'][good], c=imagecatalog['ellipticity'][good], vmin=0,vmax=0.15)
    plt.colorbar()
    plt.xlabel ("x center")
    plt.ylabel ("y center")
    plt.savefig ('ellipticity2d.png')

    fig = plt.figure()

    plt.subplot (2,1,2)
    plt.plot (imagecatalog['x'][good], imagecatalog['ellipticity'][good],'.')
    plt.xlabel ("x position")
    plt.ylabel ("ellipticity")


    plt.ylim([0,0.2])
    plt.subplot (2,1,1)
    plt.plot (imagecatalog['y'][good], imagecatalog['ellipticity'][good], '.')
    plt.xlabel ("y position")
    plt.ylabel ("ellipticity")
    plt.ylim([0,0.2])
    plt.savefig ('ellipticity1d.png')



def main():
    catalogmaker = SEPSourceCatalogProvider()

    for image in sys.argv[1:]:

        catalog, wcs = catalogmaker.get_source_catalog(image, )

        print (catalog['x'])
        plotdistribution(catalog)



if __name__ == '__main__':
    main()
