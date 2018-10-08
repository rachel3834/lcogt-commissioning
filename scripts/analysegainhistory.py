import argparse
import logging
import numpy as np
from noisegainrawmef import noisegaindbinterface
import matplotlib.pyplot as plt
_logger = logging.getLogger(__name__)
import matplotlib.dates as mdates
import datetime


starttime = datetime.datetime(2016, 1, 1)

endtime = datetime.datetime.utcnow().replace(day=28) + datetime.timedelta(days=31+4)
endtime.replace(day =1)

def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='Analyse long term gain behaviour in LCO cameras')


    parser.add_argument('--log_level', dest='log_level', default='INFO', choices=['DEBUG', 'INFO', 'WARN'],
                        help='Set the debug level')


    parser.add_argument ('--database', default="noisegain.sql")
    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')


    return args


def dateformat (starttime=None,endtime=None):
    """ Utility to prettify a plot with dates.
    """

    #plt.xlim([starttime, endtime])
    plt.gcf().autofmt_xdate()
    years = mdates.YearLocator()   # every year
    months = mdates.MonthLocator(bymonth=[4, 7, 10])  # every month
    yearsFmt = mdates.DateFormatter('%Y %b')
    monthformat = mdates.DateFormatter('%b')
    plt.gca().xaxis.set_major_locator(years)
    plt.gca().xaxis.set_major_formatter(yearsFmt)
    plt.gca().xaxis.set_minor_locator(months)
    plt.gca().xaxis.set_minor_formatter(monthformat)
    plt.setp(plt.gca().xaxis.get_minorticklabels(), rotation=45)
    plt.setp(plt.gca().xaxis.get_majorticklabels(), rotation=45)
    plt.gca().grid(which='minor')


if __name__ == '__main__':

    args = parseCommandLine()
    plt.style.use('ggplot')
    database = None
    database = noisegaindbinterface(args.database)

    for camera in database.getcameras():

        dataset = database.readmeasurements(camera)
        extensions = sorted(set(dataset['extension']))

        plt.figure()

        for ext in extensions:
            _logger.info ("ptc for extension %s" %ext)
            l = dataset['level'][dataset['extension'] == ext]
            n = (dataset['diffnoise'][dataset['extension'] == ext])
            plt.loglog (l,n, '.', label="ext %s" % ext)

        plt.title ("Photon transfer curve %s" % camera)
        plt.xlabel('Level [ADU]')
        plt.ylabel('delta flat Noise [ADU]')
        plt.ylim([5,500])
        plt.xlim([1,64000])
        plt.legend()
        plt.savefig ("ptchist-%s.png" % camera)
        plt.cla()
        plt.close()


        ######################################
        plt.figure()
        for ext in extensions:
            d = dataset['dateobs'][dataset['extension'] == ext]
            g = dataset['gain'][dataset['extension'] == ext]
            plt.plot (d,g, '.', label="ext %s" % ext)

        if 'fa' in camera:
            plt.ylim([2,4])
        if 'fl' in camera:
            plt.ylim([1,3])

        plt.xlabel('Date Obs')
        plt.ylabel('Gain [e-/ADU]')
        plt.title ('Gain history for %s' % camera)
        dateformat(starttime, endtime)
        plt.legend()
        plt.savefig ("gainhist-%s.png" % camera)
        plt.cla()
        plt.close()

        ##################################3

    database.close()




