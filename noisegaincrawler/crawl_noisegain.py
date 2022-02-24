import argparse
import logging
import os
from lcocommissioning.common.lco_archive_utilities import ArchiveDiskCrawler, get_frames_for_noisegainanalysis, \
    filename_to_archivepath_dict
from lcocommissioning.common.noisegaindb_orm import noisegaindb
from lcocommissioning.noisegainrawmef import do_noisegain_for_fileset

log = logging.getLogger(__name__)

def find_files_and_invoke_noisegain(date, args, camera=None, cameratype=None):
    """ Identify viable files for noisegain measurement from OpenSearch  and execute the noise gain calculation on those.
    """
    files = get_frames_for_noisegainanalysis(date, camera=camera, cameratype=cameratype if cameratype else None,
                                             readmode=args.readmode)
    filedict = filename_to_archivepath_dict(files)
    for camera in filedict:
        files = filedict[camera]
        if len(files) >= 4:  # Chances we have two flats and two biases..... if less, no need to go on.
            database = noisegaindb(args.database) if args.database is not None else None
            if (database is not None) and args.noreprocessing:
                # Remove duplicates - do not touch data we already analyzed and have a database record of.
                for inputname in files['FILENAME']:
                    if database.checkifalreadyused(os.path.basename(inputname)):
                        log.info(
                            "File %s was already used in a noise measurement. Skipping this entire batch." % inputname)
                        if database is not None:
                            database.close()
                        return
            log.info("{} {} # files: {}".format(camera, date, len(files)))
            try:
                do_noisegain_for_fileset(files, database, args, frameidtranslationtable=files)
            except Exception as e:
                log.error(f'While doing noisegain for file set: {e}')
            if database is not None:
                database.close()


def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='Crawl LCO archive tyo measure noise, gain from paitrs of biases and darks',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--cameratype', type=str, nargs='+', default=['fa', 'fs', 'kb', 'ep', 'sq'],
                        help='Type of cameras to parse')
    parser.add_argument('--instrument', type=str, default=None,
                        help='Individual camera to request')

    group = parser.add_argument_group(
        'Optionally, specify the location of statistics window. All units in pixels with an FITS image extension')
    group.add_argument('--minx', type=int, default=None, help="minimum x.")
    group.add_argument('--maxx', type=int, default=None, help="maximum x.")
    group.add_argument('--miny', type=int, default=None, help="miniumm y.")
    group.add_argument('--maxy', type=int, default=None, help="maximum y.")

    mutex = parser.add_mutually_exclusive_group()
    mutex.add_argument('--date', dest='date', nargs='+', help='Specific date to process.')
    mutex.add_argument('--ndays', type=int, help='Look back for the last 3 days, including current one, as in DAY-OBS')

    parser.add_argument('--useaws', action='store_true')
    parser.add_argument('--database', default="sqlite:///noisegain.sqlite", help="sqlite database where to store results.")

    parser.add_argument('--readmode', default="full_frame",
                        help="CCD readmode, typically full_frame, default, or central_2k_2x2")
    parser.add_argument('--ncpu', default=2, type=int, help='How many threads to use (One thread per day to process)')
    parser.add_argument('--noreprocessing', action='store_true',
                        help="Do not reprocess if data are show up in database")

    parser.add_argument('--makepng', action='store_true', help="Create a png output image of noise, gain, and ptc.")
    parser.add_argument('--showimages', action='store_true', help="Interactively show difference flat and bias images.")
    parser.add_argument('--ignoretemp', action='store_true',
                        help="ignore if actual temperature differs from set point temperature. Reject by default.")
    parser.add_argument('--loglevel', dest='log_level', default='INFO', choices=['DEBUG', 'INFO', 'WARN'],
                        help='Set the debug level')

    args = parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')

    args.sortby = "filterlevel"
    if args.ndays is not None:
        log.info (f"Determining dates from ndays argument {args.ndays}")
        args.dates = ArchiveDiskCrawler.get_last_n_days(args.ndays)
    else:
        args.dates = args.date



    return args

def main():
    args = parseCommandLine()

    if args.cameratype is None:
        args.cameratype = args.camera[0:2]
    log.info (f"These are the dates to process: {args.dates}")
    for date in args.dates:
        for ct in args.cameratype:
            log.debug ("Processing from date, cameratype, camera {} {} {}".format (date, args.cameratype, args.instrument))
            find_files_and_invoke_noisegain(date, args, camera=args.instrument, cameratype=ct)
    log.info("All Done")

if __name__ == '__main__':
    main()
