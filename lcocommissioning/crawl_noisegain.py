import argparse
import logging
import multiprocessing as mp
import os
from concurrent.futures.process import ProcessPoolExecutor

from lcocommissioning.noisegainrawmef import do_noisegain_for_fileset
from lcocommissioning.common.noisegaindbinterface import noisegaindbinterface
from lcocommissioning.common.lcoarchivecrawler import ArchiveCrawler

log = logging.getLogger(__name__)


def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='Crawl LCO archive tyo measure noise, gain from paitrs of biases and darks',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--cameratype', type=str, nargs='+', default=['fa??', 'fs??', 'kb??'],
                        help='Type of cameras to parse')

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

    parser.add_argument('--loglevel', dest='log_level', default='DEBUG', choices=['DEBUG', 'INFO', 'WARN'],
                        help='Set the debug level')

    args = parser.parse_args()
    args.sortby = "filterlevel"
    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')
    return args


def findfilesanddonoisegain(camera, date, args):
    files = ArchiveCrawler.findfiles_for_camera_dates(camera, date, 'raw', "*[x00|f00].fits*")

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
    c = ArchiveCrawler()
    dates = c.get_last_N_days(args.ndays)
    cameras = c.find_cameras(cameras=args.cameratype)

    for camera in cameras:
        for date in dates:
            exec.submit (findfilesanddonoisegain, camera, date, args)

    log.info("Waiting for all threads to complete.")
    exec.shutdown()
    log.info ("All Done")