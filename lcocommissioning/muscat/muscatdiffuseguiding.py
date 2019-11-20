import logging
import sys
import matplotlib.pyplot as plt

from lcocommissioning.common.SourceCatalogProvider import SEPSourceCatalogProvider
import numpy as np

reference_x = [-1.210,
               -2.730
               - 2.940,
               -3.040,
               -2.640,
               -2.640,
               -2.330,
               -2.720,
               -3.080, ]

reference_y = [-0.710,
               -0.500,
               -0.140,
               -0.310,
               -0.470,
               -0.290,
               -1.230,
               -0.160,
               -0.340, ]


def dist(x1, y1, x2, y2, i):
    return np.sqrt((x1[i] - x2) ** 2 + (y1[i] - y2) ** 2)


def match(x1, y1, x2, y2, tolerance=10):
    assert (len(x1) == len(y1))
    assert (len(x2) == len(y2))

    xm = x1 * 0
    ym = y1 * 0

    for ii in range(len(x1)):

        d = dist(x1, y1, x2, y2, ii)
        nearest = np.argmin(d)
        if d[nearest] < tolerance:
            xm[ii] = x2[nearest]
            ym[ii] = y2[nearest]

    return xm, ym


if __name__ == '__main__':
    # TODO: Make this a test code
    logging.basicConfig(level=getattr(logging, 'INFO'),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')
    sourcecatalogProvider = SEPSourceCatalogProvider(refineWCSViaLCO=False)

    if len(sys.argv) < 2:
        exit(0)

    print(sys.argv[1])
    refcatalog, wcs = sourcecatalogProvider.get_source_catalog(sys.argv[1], ext=0)
    ii = 0
    for image in sys.argv[2:]:
        catalog, wcs = sourcecatalogProvider.get_source_catalog(image, ext=0)
        xm, ym = match(refcatalog['x'], refcatalog['y'], catalog['x'], catalog['y'])
        plt.plot([ii] * len(xm), ym - refcatalog['y'], 'o')
        print(ii, image)
        if ii < len(reference_x):
            plt.plot(ii, -reference_y[ii], 'x', color='black')
        ii += 1
    plt.ylim([-10, 10])
    plt.xlabel("Muscat file number")
    plt.ylabel("individual y-matched star pixel shift x")
    plt.savefig('muscatguide-y.png')
