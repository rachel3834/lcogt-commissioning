import errno
import multiprocessing
import os
from concurrent.futures import wait
from concurrent.futures.thread import ThreadPoolExecutor

import matplotlib
matplotlib.use('Agg')
import argparse
import logging
from lcocommissioning.common.noisegaindbinterface import noisegaindbinterface
import matplotlib.dates as mdates
from multiprocessing import Pool

import matplotlib.pyplot as plt
_logger = logging.getLogger(__name__)
import datetime


starttimeall = datetime.datetime(2016, 1, 1)
starttimefa = datetime.datetime(2018, 7, 1)

endtime = datetime.datetime.utcnow().replace(day=28) + datetime.timedelta(days=31+4)
endtime.replace(day =1)

def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='Analyse long term gain behaviour in LCO cameras')


    parser.add_argument('--loglevel', dest='log_level', default='INFO', choices=['DEBUG', 'INFO', 'WARN'],
                        help='Set the debug level')

    parser.add_argument ('--outputdir', default='gainhistory', help="directory for output graphs")

    parser.add_argument ('--database', default="noisegain.sqlite")
    parser.add_argument ('--ncpu', default=1, type=int)
    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')



    if not os.path.exists(args.outputdir):
        _logger.info ("Creating output directory [%s]" % args.outputdir)
        try:
            os.makedirs(args.outputdir)
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise

    return args


def dateformat (starttime=None,endtime=None):
    """ Utility to prettify a plot with dates.
    """

    plt.xlim([starttime, endtime])
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



def renderHTMLPage (args, cameras):
    _logger.info ("Now rendering output html page")

    outputfile = "%s/index.html" % (args.outputdir)

    message = """<html>
<head></head>
<body><title>LCO Gain History Plots</title>
"""
    message += "<p/>Figures updated %s UTC <p/>\n"  % (datetime.datetime.utcnow())
    message += """
<h1> Details by Camera: </h1>
"""

    for camera in cameras:
        message = message + " <h2> %s </h2>\n" % (camera)

        historyname = "gainhist-%s.png" % camera
        ptcname = "ptchist-%s.png" % camera
        lvlname = "levelgain-%s.png" % camera
        noisename = "noise-%s.png" % camera
        line = f'<a href="{historyname}"><img src="{historyname}" height="450"/></a>  <a href="{ptcname}"><img src="{ptcname}" height="450"/> </a> <a href="{lvlname}"><img src="{lvlname}" height="450"/> </a>  <a href="{noisename}"><img src="{noisename}" height="450"/> </a>'
        message = message + line

    message = message + "</body></html>"

    with open (outputfile, 'w+') as f:
        f.write (message)
        f.close()


def make_plots_for_camera(camera,  args):
    database = noisegaindbinterface(args.database)
    starttime = starttimeall
    if 'fa' in camera:
        starttime = starttimefa

    dataset = database.readmeasurements(camera, levelratio = 0.02, filters=goodfilters)
    if dataset is None:
        return
    extensions = sorted(set(dataset['extension']))

    plt.figure()

    for ext in extensions:

        l = dataset['level'][dataset['extension'] == ext]
        n = (dataset['diffnoise'][dataset['extension'] == ext])
        plt.loglog (l,n, '.', label="ext %s" % ext, markersize=1)

    plt.title ("Photon transfer curve %s" % camera)
    plt.xlabel('Level [ADU]')
    plt.ylabel('delta flat Noise [ADU]')
    plt.ylim([5,1000])
    plt.xlim([1,70000])
    plt.legend()
    plt.savefig ("%s/ptchist-%s.png" % (args.outputdir,camera))
    plt.cla()
    plt.close()


    ######################################
    plt.figure()
    for ext in extensions:
        d = dataset['dateobs'][dataset['extension'] == ext]
        g = dataset['gain'][dataset['extension'] == ext]
        plt.plot (d,g, '.', label="ext %s" % ext, markersize=1)

    if 'fa' in camera:
        plt.ylim([2.5,4])
    if 'fl' in camera:
        plt.ylim([1,3])
    if 'fs' in camera:
        plt.ylim([6,8])
    if 'kb' in camera:
        plt.ylim([1,3])


    plt.xlabel('Date Obs')
    plt.ylabel('Gain [e-/ADU]')
    plt.title ('Gain history for %s' % camera)
    dateformat(starttime, endtime)
    plt.legend()
    plt.savefig ("%s/gainhist-%s.png" % (args.outputdir,camera))
    plt.cla()
    plt.close()

    ##################################3

    plt.figure()
    for ext in extensions:
        d = dataset['level'][dataset['extension'] == ext]
        g = dataset['gain'][dataset['extension'] == ext]
        plt.plot (d,g, '.', label="ext %s" % ext, markersize=1)

    if 'fa' in camera:
        plt.ylim([2.5,4])
    if 'fl' in camera:
        plt.ylim([1,3])
    if 'fs' in camera:
        plt.ylim([6,8])
    if 'kb' in camera:
        plt.ylim([1,3])

    plt.xlim([0,70000])
    plt.xlabel('Avg. Level [ADU]')
    plt.ylabel('Gain [e-/ADU]')
    plt.title ('Level vs Gain for %s' % camera)

    plt.legend()
    plt.savefig ("%s/levelgain-%s.png" % (args.outputdir,camera))
    plt.cla()
    plt.close()


    ######################################
    plt.figure()
    for ext in extensions:
        d = dataset['dateobs'][dataset['extension'] == ext]
        g = dataset['readnoise'][dataset['extension'] == ext]
        plt.plot (d,g, '.', label="ext %s" % ext, markersize=1)

    if 'fa' in camera:
        plt.ylim([5,10])
    if 'fl' in camera:
        plt.ylim([5,15])
    if 'fs' in camera:
        plt.ylim([5,15])
    if 'kb' in camera:
        plt.ylim([7,20])


    dateformat(starttime, endtime)
    plt.xlabel('Date')
    plt.ylabel('Readnoise [e-]')
    plt.title ('Readnoise vs time for %s' % camera)

    plt.legend()
    plt.savefig ("%s/noise-%s.png" % (args.outputdir,camera))
    plt.cla()
    plt.close()

    database.close()


goodfilters = ['up','gp','rp','ip','zp','U','B','V','R','I']


def main():

    args = parseCommandLine()
    plt.style.use('ggplot')
    matplotlib.rcParams['savefig.dpi'] = 200

    database = noisegaindbinterface(args.database)
    cameras = database.getcameras()
    database.close()

    for camera in cameras:
        make_plots_for_camera( camera,  args)


    renderHTMLPage(args, sorted(database.getcameras()))


if __name__ == '__main__':
    main()


