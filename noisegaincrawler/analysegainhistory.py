import errno
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
import logging
import datetime
import boto3
import io
from lcocommissioning.common.common import dateformat
from lcocommissioning.common.noisegaindb_orm import noisegaindb

_logger = logging.getLogger(__name__)
logging.getLogger('matplotlib').setLevel(logging.FATAL)

starttimeall = datetime.datetime(2016, 1, 1)
starttimefa = datetime.datetime(2018, 7, 1)
endtime = datetime.datetime.utcnow().replace(day=28) + datetime.timedelta(days=31 + 4)
endtime.replace(day=1)

fareadmodes = [['full_frame', None], ['central_2k_2x2', ]]


def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='Analyse long term gain behaviour in LCO cameras')

    parser.add_argument('--loglevel', dest='log_level', default='INFO', choices=['DEBUG', 'INFO', 'WARN'],
                        help='Set the debug level')

    parser.add_argument('--outputdir', default=None,  help="directory for output graphs")

    parser.add_argument('--cameras', default=None, nargs="*")
    parser.add_argument('--database', default="sqlite:///noisegain.sqlite")
    parser.add_argument('--ncpu', default=1, type=int)
    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')

    if not aws_enabled():
        if (args.outputdir is not None) and (not os.path.exists(args.outputdir)):
            _logger.info("Creating output directory [%s]" % args.outputdir)
            try:
                os.makedirs(args.outputdir)
            except OSError as exc:  # Guard against race condition
                if exc.errno != errno.EEXIST:
                    raise
    else:
        _logger.info ("Working in AWS environment, not using file system storage backend.")

    return args


def aws_enabled():
    '''Return True if AWS support is configured'''
    access_key = os.environ.get('AWS_ACCESS_KEY_ID', None)
    secret_key = os.environ.get('AWS_SECRET_ACCESS_KEY', None)
    s3_bucket = os.environ.get('AWS_S3_BUCKET', None)
    region = os.environ.get('AWS_DEFAULT_REGION', None)

    return access_key and secret_key and s3_bucket and region


def write_to_storage_backend(directory, filename, data, binary=True):
    _logger.info (f"writing object {filename}")
    if aws_enabled():
        # AWS S3 Bucket upload
        client = boto3.client('s3')
        bucket = os.environ.get('AWS_S3_BUCKET', None)
        with io.BytesIO(data) as fileobj:
            _logger.debug(f'Write data to AWS S3: {bucket}/{filename}')
            response = client.upload_fileobj(fileobj, bucket, filename)
            return response
    else:
        fullpath = os.path.join(directory, filename)
        _logger.info (f'writing to file system {fullpath}')
        with open(fullpath, 'wb' if binary else 'w') as fileobj:
            fileobj.write(data)
            return True


def renderHTMLPage(args, cameras, filenames):
    _logger.info("Now rendering output html page")

    outputfile = "%s/index.html" % (args.outputdir)

    message = """<html>
<head></head>
<body><title>LCO Gain History Plots</title>
"""
    message += "<p/>Figures updated %s UTC <p/>\n" % (datetime.datetime.utcnow())
    message += """
<h1> Details by Camera: </h1>
"""

    for camera in cameras:

        readmodes = [None, ]
        if 'fl' in camera:
            continue
        if 'fa' in camera:
            readmodes = fareadmodes

        message = message + " <h2> %s </h2>\n" % (camera)
        for readmode in readmodes:
            if readmode is not None:
                readmode = [x if x is not None else 'None' for x in readmode]

            readmode = "".join(readmode) if readmode is not None else ""

            historyname = "gainhist-%s%s.png" % (camera, readmode)
            ptcname = "ptchist-%s%s.png" % (camera, readmode)
            lvlname = "levelgain-%s%s.png" % (camera, readmode)
            noisename = "noise-%s%s.png" % (camera, readmode)
            flatlevelname = "flatlevel-%s%s.png" % (camera, readmode)
            line = f'<a href="{historyname}"><img src="{historyname}" height="450"/></a>  <a href="{ptcname}"><img src="{ptcname}" height="450"/> </a> <a href="{lvlname}"><img src="{lvlname}" height="450"/> </a>  <a href="{noisename}"><img src="{noisename}" height="450"/> </a> <a href="{flatlevelname}"><img src="{flatlevelname}" height="450"/><br/> </a>'
            message = message + line

    message = message + "</body></html>"

    with io.StringIO() as fileobj:
        fileobj.write(message)
        filename = 'index.html'
        write_to_storage_backend(args.outputdir, filename, fileobj.getvalue(), binary=False)



def make_plots_for_camera(camera, args, database):
    ''' Crate a list of plots for all cameras matching pattern and retuirn a list of filenames of the resulting plots.
    '''

    filenames = []
    outputdir = args.outputdir
    _logger.info("readmodes: {} {}".format(camera, database.getReadmodesFroCamera(camera)))
    readmodes = [None, ]

    starttime = starttimeall
    if 'fa' in camera:
        starttime = starttimefa
        readmodes = fareadmodes

    for readmode in readmodes:
        dataset = database.getMeasurementsForCamera(camera, levelratio=0.02, filters=goodfilters, readmode=readmode)
        if dataset is None:
            return
        extensions = sorted(set(dataset['extension']))
        if readmode is not None:
            readmode = [x if x is not None else 'None' for x in readmode]

        filenames.append(plot_ptc(camera, dataset, extensions, outputdir, readmode))
        filenames.append(plot_gainhist(camera, dataset, extensions, outputdir, starttime, readmode))
        filenames.append(plot_levelgain(camera, dataset, extensions, outputdir, readmode))
        filenames.append(plotnoisehist(camera, dataset, extensions, outputdir, starttime, readmode))
        filenames.append(plot_flatlevelhist(camera, dataset, extensions, outputdir, starttime, readmode))

    return filenames


def plot_flatlevelhist(camera, dataset, extensions, outputdir, starttime, readmode=None):
    plt.figure()
    for ext in extensions:
        d = dataset['dateobs'][dataset['extension'] == ext]
        l = dataset['level'][dataset['extension'] == ext]
        plt.plot(d, l, '.', label="ext %s" % ext, markersize=1)
    plt.ylim([1, 100000])
    dateformat(starttime, endtime)
    plt.xlabel('Date')
    plt.ylabel('Flat Level [ADU]')
    plt.title('Flat level vs time for %s' % camera)
    plt.legend()
    readmode = "".join(readmode) if readmode is not None else ""

    with io.BytesIO() as fileobj:
        plt.savefig(fileobj, format='png', bbox_inches='tight')
        plt.close()
        filename = "flatlevel-{}{}.png".format(camera, readmode)
        write_to_storage_backend(outputdir, filename, fileobj.getvalue())
    plt.cla()
    return filename


def plotnoisehist(camera, dataset, extensions, outputdir, starttime, readmode=None):
    readmode = "".join(readmode) if readmode is not None else ""
    plt.figure()
    for ext in extensions:
        d = dataset['dateobs'][dataset['extension'] == ext]
        g = dataset['readnoise'][dataset['extension'] == ext]
        plt.plot(d, g, '.', label="ext %s" % ext, markersize=1)
    if 'fa' in camera:
        plt.ylim([5, 10])
    if 'fl' in camera:
        plt.ylim([5, 15])
    if 'fs' in camera:
        plt.ylim([5, 15])
    if 'kb' in camera:
        plt.ylim([7, 20])
    dateformat(starttime, endtime)
    plt.xlabel('Date')
    plt.ylabel('Readnoise [e-]')
    plt.title('Readnoise vs time for %s %s' % (camera, readmode))
    plt.legend()

    with io.BytesIO() as fileobj:
        plt.savefig(fileobj, format='png', bbox_inches='tight')
        plt.close()
        filename = "noise-{}{}.png".format(camera, readmode)
        write_to_storage_backend(outputdir, filename, fileobj.getvalue())
    plt.cla()
    return filename


def plot_levelgain(camera, dataset, extensions, outputdir, readmode=None):
    readmode = "".join(readmode) if readmode is not None else ""
    plt.figure()
    for ext in extensions:
        d = dataset['level'][dataset['extension'] == ext]
        g = dataset['gain'][dataset['extension'] == ext]
        plt.plot(d, g, '.', label="ext %s" % ext, markersize=1)
    if 'fa' in camera:
        plt.ylim([2.5, 4])
        if '2x2' in readmode:
            plt.ylim([5.5, 7])
    if 'fl' in camera:
        plt.ylim([1, 3])
    if 'fs' in camera:
        plt.ylim([6, 8])
    if 'kb' in camera:
        plt.ylim([1, 3])
    plt.xlim([0, 70000])
    plt.xlabel('Avg. Level [ADU]')
    plt.ylabel('Gain [e-/ADU]')
    plt.title('Level vs Gain for %s %s' % (camera, readmode))
    plt.legend()

    with io.BytesIO() as fileobj:
        plt.savefig(fileobj, format='png', bbox_inches='tight')
        plt.close()
        filename = "levelgain-{}{}.png".format(camera, readmode)
        write_to_storage_backend(outputdir, filename, fileobj.getvalue())
    plt.cla()
    return filename


def plot_gainhist(camera, dataset, extensions, outputdir, starttime, readmode=None):
    readmode = "".join(readmode) if readmode is not None else ""
    plt.figure()
    for ext in extensions:
        d = dataset['dateobs'][dataset['extension'] == ext]
        g = dataset['gain'][dataset['extension'] == ext]
        plt.plot(d, g, '.', label="ext %s" % ext, markersize=1)
    if 'fa' in camera:
        plt.ylim([2.5, 4])
        if '2x2' in readmode:
            plt.ylim([5.5, 7])
    if 'fl' in camera:
        plt.ylim([1, 3])
    if 'fs' in camera:
        plt.ylim([6, 8])
    if 'kb' in camera:
        plt.ylim([1, 3])
    plt.xlabel('Date Obs')
    plt.ylabel('Gain [e-/ADU]')
    plt.title('Gain history for %s %s ' % (camera, readmode))
    dateformat(starttime, endtime)
    plt.legend()

    with io.BytesIO() as fileobj:
        plt.savefig(fileobj, format='png', bbox_inches='tight')
        plt.close()
        filename = "gainhist-{}{}.png".format(camera, readmode)
        write_to_storage_backend(outputdir, filename, fileobj.getvalue())

    plt.cla()
    return filename


def plot_ptc(camera, dataset, extensions, outputdir, readmode=None):
    readmode = "".join(readmode) if readmode is not None else ""

    plt.figure()
    for ext in extensions:
        l = dataset['level'][dataset['extension'] == ext]
        n = (dataset['diffnoise'][dataset['extension'] == ext])
        plt.loglog(l, n, '.', label="ext %s" % ext, markersize=1)
    plt.title("Photon transfer curve %s %s " % (camera, readmode))
    plt.xlabel('Level [ADU]')
    plt.ylabel('delta flat Noise [ADU]')
    plt.ylim([5, 1000])
    plt.xlim([1, 70000])
    plt.legend()
    with io.BytesIO() as fileobj:
        plt.savefig(fileobj, format='png', bbox_inches='tight')
        plt.close()
        filename = "ptchist-{}{}.png".format(camera, readmode)
        write_to_storage_backend(outputdir, filename, fileobj.getvalue())
    plt.cla()
    return filename


goodfilters = ['up', 'gp', 'rp', 'ip', 'zp', 'U', 'B', 'V', 'R', 'I']


def main():
    args = parseCommandLine()
    plt.style.use('ggplot')
    matplotlib.rcParams['savefig.dpi'] = 400

    database = noisegaindb(args.database)
    cameras = args.cameras if args.cameras is not None else database.getCameras()
    #we are not intersted in old fl data that may be included
    cameras = [ c for c in cameras if not c.startswith('fl')]
    _logger.info("Cameras: {}".format(cameras))

    filenames = []
    for camera in cameras:
        filenames += make_plots_for_camera(camera, args, database)

    renderHTMLPage(args, sorted(cameras), filenames)
    database.close()


if __name__ == '__main__':
    main()
