import errno
import os
import matplotlib
from lcocommissioning.common.common import dateformat

matplotlib.use('Agg')
import argparse
import logging
from lcocommissioning.common.noisegaindbinterface import noisegaindbinterface, darkcurrentdbinterface
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.FATAL)
import matplotlib.pyplot as plt

_logger = logging.getLogger(__name__)
import datetime


def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='Analyse long term gain behaviour in LCO cameras')

    parser.add_argument('--loglevel', dest='log_level', default='WARN', choices=['DEBUG', 'INFO', 'WARN'],
                        help='Set the debug level')

    parser.add_argument('--outputdir', default='darkhistory', help="directory for output graphs")

    parser.add_argument('--database', default="darkcurrent.sqlite")
    parser.add_argument('--ncpu', default=1, type=int)
    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')

    if not os.path.exists(args.outputdir):
        _logger.info("Creating output directory [%s]" % args.outputdir)
        try:
            os.makedirs(args.outputdir)
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise

    return args


def make_plots_for_camera(camera, args):
    database = darkcurrentdbinterface(args.database)
    dataset = database.readmeasurements(camera)

    plt.figure()

    plt.title ('Dark Current for {}'.format (camera))
    plt.xlabel ("DATE-OBS")
    plt.ylabel ("dark current [e-/s]")
    plt.ylim ([0,0.01])
    plt.plot (dataset['dateobs'], dataset['darkcurrent'], '.')
    dateformat()
    plt.savefig ("{}/darkcurrenthist-{}.png".format (args.outputdir, camera))
    plt.close()

def renderHTMLPage(args, cameras):
    _logger.info("Now rendering output html page")

    outputfile = "%s/index.html" % (args.outputdir)

    message = """<html>
<head></head>
<body><title>LCO Dark Current History Plots</title>
"""
    message += "<p/>Figures updated %s UTC <p/>\n" % (datetime.datetime.utcnow())
    message += """
<h1> Details by Camera: </h1>
"""

    for camera in cameras:
        if 'fl' in camera:
            continue
        message = message + " <h2> %s </h2>\n" % (camera)

        historyname = "darkcurrenthist-%s.png" % camera

        line = f'<a href="{historyname}"><img src="{historyname}" height="450"/></a>  '
        message = message + line

    message = message + "</body></html>"

    with open(outputfile, 'w+') as f:
        f.write(message)
        f.close()




def main():
    args = parseCommandLine()
    plt.style.use('ggplot')
    matplotlib.rcParams['savefig.dpi'] = 200

    database = darkcurrentdbinterface(args.database)
    cameras = database.getcameras()
    _logger.debug("Cameras: {}".format(cameras))
    for camera in cameras:
        make_plots_for_camera(camera, args)

    renderHTMLPage(args, sorted(database.getcameras()))
    database.close()


if __name__ == '__main__':
    main()
