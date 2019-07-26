import argparse

from astropy.io import fits
import numpy as np
import logging

from common.noisegaindbinterface import darkcurrentdbinterface
from lcocommissioning.common.lcoarchivecrawler import ArchiveCrawler

log = logging.getLogger(__name__)

def findkeywordinhdul(hdulist, keyword):
    # f&*& fpack!
    for ext in hdulist:
        val = ext.header.get(keyword)
        if val is not None:
            return val
    return None

def get_dark_current_from_master(file):
    hdul = fits.open (file)
    data = hdul[1].data[50:-50,:]
    readmode = findkeywordinhdul(hdul, 'CONFMODE')
    dateobs = findkeywordinhdul(hdul, 'DATE-OBS')
    instrument = findkeywordinhdul(hdul, 'INSTRUME')
    hdul.close()
    sigma = np.std (data)
    dark_current = np.average ( data [ np.abs ( data - np.nanmedian(data)) < 5 * sigma ] )
    log.info ("dark current for image {} is {} ".format (file, dark_current))
    return dark_current, dateobs, readmode, instrument



def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='Crawl LCO archive tyo measure dark current from master darks',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--cameratype', type=str, nargs='+', default=['fa??', 'fs??', 'kb??'],
                        help='Type of cameras to parse')

    parser.add_argument('--ndays', default=3, type=int, help="How many days to look into the past")
    parser.add_argument('--database', default="darkcurrent.sqlite", help="sqlite database where to store results.")
    parser.add_argument('--readmode', default="full_frame")

    parser.add_argument('--ncpu', default=2, type=int)
    parser.add_argument('--noreprocessing', action='store_true',
                        help="Do not reprocess if data are already in database")

    parser.add_argument('--loglevel', dest='log_level', default='INFO', choices=['DEBUG', 'INFO', 'WARN'],
                        help='Set the debug level')

    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')
    return args


def find_darks_and_process (camera, dates, args):

    files = []
    for date in dates:
        datefiles = ArchiveCrawler.findfiles_for_camera_dates(camera, date, 'processed', "*dark*.fits*")
        files.extend(datefiles)

    if len(files) > 0:
        database = darkcurrentdbinterface(args.database) if args.database is not None else None
        if (database is not None):

            for file in files:
                if not database.checkifalreadyused(file):
                    darkcurrent, dateobs, readmode, instrument = get_dark_current_from_master(file)
                    database.addmeasurement (file, dateobs, instrument, darkcurrent, readmode)
                else:
                    log.debug ("file already measured, skipping")

            database.close()


def main ():
    pass




if __name__ == '__main__':
    args = parseCommandLine()
    c = ArchiveCrawler()
    dates = c.get_last_N_days(args.ndays)
    cameras = c.find_cameras(cameras=args.cameratype)

    for camera in cameras:

            find_darks_and_process(camera, dates, args)