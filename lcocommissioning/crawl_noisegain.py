import argparse
import logging
import multiprocessing as mp
import os
from concurrent.futures.process import ProcessPoolExecutor

from lcocommissioning.noisegainrawmef import do_noisegain_for_fileset
from lcocommissioning.common.noisegaindbinterface import noisegaindbinterface
from lcocommissioning.common.lcoarchivecrawler import ArchiveCrawler, get_frames_for_noisegainanalysis, \
    filename_to_archivepath

log = logging.getLogger(__name__)


def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='Crawl LCO archive tyo measure noise, gain from paitrs of biases and darks',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--cameratype', type=str, nargs='+', default=['fa', 'fs', 'kb'],
                        help='Type of cameras to parse')

    parser.add_argument('--instrument', type=str,  default=None,
                        help='Individual camera to request')
    group = parser.add_argument_group(
        'Optionally, specify the location of statistics window. All units in pixels with an FITS image extension')
    group.add_argument('--minx', type=int, default=None, help="minimum x.")
    group.add_argument('--maxx', type=int, default=None, help="maximum x.")
    group.add_argument('--miny', type=int, default=None, help="miniumm y.")
    group.add_argument('--maxy', type=int, default=None, help="maximum y.")

    parser.add_argument('--ndays', default=3, type=int, help="How many days to look into the past")
    parser.add_argument('--database', default="noisegain.sqlite", help="sqlite database where to store results.")
    parser.add_argument('--readmode', default="full_frame")

    parser.add_argument('--ncpu', default=2, type=int)
    parser.add_argument('--noreprocessing', action='store_true',
                        help="Do not reprocess if datra are already in database")

    parser.add_argument('--makepng', action='store_true', help="Create a png output image of noise, gain, and ptc.")
    parser.add_argument('--showimages', action='store_true', help="Interactively show difference flat and bias images.")
    parser.add_argument('--ignoretemp', action='store_true',
                        help="ignore if actual temperature differs from set point temperature. Reject by default.")
    parser.add_argument('--loglevel', dest='log_level', default='DEBUG', choices=['DEBUG', 'INFO', 'WARN'],
                        help='Set the debug level')

    args = parser.parse_args()
    args.sortby = "filterlevel"
    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')
    return args


def findfilesanddonoisegain(date, args, camera=None, cameratype=None, useElasticsearch=True):


    if useElasticsearch:
        files = get_frames_for_noisegainanalysis(date, camera=camera, cameratype=cameratype[0] if cameratype else None)
        filedict =  filename_to_archivepath(files)
    else:
        if camera is None:
            raise ("Camera is not specified")
        filedict = {camera:  ArchiveCrawler.findfiles_for_camera_dates(camera, date, 'raw', "*[bf]00.fits*")}

    for camera in filedict:
        files = filedict[camera]
        if len(files) > 3:
            database = noisegaindbinterface(args.database) if args.database is not None else None
            if (database is not None) and args.noreprocessing:
                for inputname in files:
                    if database.checkifalreadyused(os.path.basename(inputname)):
                        log.info("File %s was already used in a noise measurement. Skipping this entire batch." % inputname)
                        if database is not None:
                            database.close()
                        return

            log.debug("{} {} # files: {}".format(camera, date, len(files)))
            try:
                do_noisegain_for_fileset(files, database, args)
            except Exception as e:
                log.error ('While launching task:',e)
            if database is not None:
                database.close()


if __name__ == '__main__':
    args = parseCommandLine()

    exec = ProcessPoolExecutor (max_workers = args.ncpu)

    dates = ArchiveCrawler.get_last_N_days(args.ndays)
    for date in dates:
        exec.submit (findfilesanddonoisegain, date, args, args.instrument, args.cameratype)

    log.info("Waiting for all threads to complete.")
    exec.shutdown()
    log.info ("All Done")