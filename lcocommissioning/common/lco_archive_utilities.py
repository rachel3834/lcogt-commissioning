import datetime
import glob
import io
import logging
import os

import numpy as np
import requests
from astropy.io import fits
from astropy.table import Table

from lcocommissioning.common.common import lco_site_lonlat
from opensearchpy import OpenSearch
from opensearch_dsl import Search

log = logging.getLogger(__name__)
logging.getLogger('opensearch').setLevel(logging.WARNING)
logging.getLogger('connectionpool').setLevel(logging.WARNING)

ARCHIVE_ROOT = "/archive"
ARCHIVE_API_TOKEN = os.getenv('ARCHIVE_API_TOKEN', '')


class ArchiveDiskCrawler:
    ''' Legacy code from the good old times (2019) when the /archive mount was accessible, and everybody was happy in
    the file system land.

    Now everybody moved on to databases, object stores, and other witchcraft.

    '''

    archive_root = None

    def __init__(self, rootdirectory=ARCHIVE_ROOT):
        self.archive_root = rootdirectory

    def find_cameras(self, sites=lco_site_lonlat, cameras=["fa??", "fs??", "kb??"]):
        sitecameras = []
        for site in sites:
            for camera in cameras:
                dir = "{}/{}/{}".format(self.archive_root, site, camera)
                candiates = glob.glob(dir)
                sitecameras.extend(candiates)
        return sitecameras

    @staticmethod
    def get_last_n_days(lastNdays):
        ''' Utility method to return the last N days as string YYYYMMDD, nicely arranged in an array.'''

        date = []
        today = datetime.datetime.utcnow()
        for ii in range(lastNdays):
            day = today - datetime.timedelta(days=ii)
            date.append(day.strftime("%Y%m%d"))
        return date[::-1]

    @staticmethod
    def findfiles_for_camera_dates(sitecamera, date, raworprocessed, filetempalte, prefix=""):
        dir = "{}{}/{}/{}/{}".format(prefix, sitecamera, date, raworprocessed, filetempalte)
        files = glob.glob(dir)
        if (files is not None) and (len(files) > 0):
            myfiles = np.asarray ([[f, "-1"] for f in files])
            return Table(myfiles, names=['FILENAME', 'FRAMEID'])
        return None



def make_opensearch(index, filters, queries=None, exclusion_filters=None, range_filters=None, prefix_filters=None,
                    terms_filters=None,
                    opensearch_url='https://opensearch.lco.global'):
    """
    Make an ElasticSearch query

    Parameters
    ----------
    index : str
            Name of index to search
    filters : list of dicts
              Each dict has a criterion for an OpenSearch "filter"
    queries : list of dicts
              Each dict has a "type" and "query" entry. The 'query' entry is a dict that has a criterion for an
              ElasticSearch "query"
    exclusion_filters : list of dicts
                        Each dict has a criterion for an OpenSearch "exclude"
    range_filters: list of dicts
                   Each dict has a criterion an openSearch "range filter"
    opensearch_url : str
             URL of the openSearch host

    Returns
    -------
    search : opensearch_dsl.Search
             The OpenSearch object
    """
    if queries is None:
        queries = []
    if exclusion_filters is None:
        exclusion_filters = []
    if range_filters is None:
        range_filters = []
    if terms_filters is None:
        terms_filters = []
    if prefix_filters is None:
        prefix_filters = []
    opensearch = OpenSearch(opensearch_url)
    s = Search(using=opensearch, index=index)
    for f in filters:
        s = s.filter('term', **f)
    for f in terms_filters:
        s = s.filter('terms', **f)
    for f in range_filters:
        s = s.filter('range', **f)
    for f in prefix_filters:
        s = s.filter('prefix', **f)
    for f in exclusion_filters:
        s = s.exclude('term', **f)
    for q in queries:
        s = s.query(q['type'], **q['query'])
    return s


def get_frames_for_noisegainanalysis(dayobs, site=None, cameratype=None, camera=None, readmode='full_frame',
                                     obstype=['BIAS', 'SKYFLAT'], opensearch_url='https://opensearch.lco.global'):
    """ Queries for a list of processed LCO images that are viable to get a photometric zeropoint in the griz bands measured.

        Selection criteria are by DAY-OBS, site, by camera type (fs,fa,kb), what filters to use, and minimum exposure time.
        Only day-obs is a mandatory fields, we do not want to query the entire archive at once.
     """
    log.debug ("Starting opensearch query")
    query_filters = [{'DAY-OBS': dayobs}, {'RLEVEL': 0}, {'CONFMODE': readmode}]
    range_filters = []
    terms_filters = [{'OBSTYPE': obstype}]
    prefix_filters = []

    if site is not None:
        query_filters.append({'SITEID': site})
    if camera is not None:
        query_filters.append({'INSTRUME': camera})
    if cameratype is not None:
        prefix_filters.append({'INSTRUME': cameratype})

    queries = []
    records = make_opensearch('fitsheaders', query_filters, queries, exclusion_filters=None, opensearch_url=opensearch_url,
                              range_filters=range_filters, prefix_filters=prefix_filters,
                              terms_filters=terms_filters).scan()
    records_sanitized = [[record['filename'], record['SITEID'], record['INSTRUME'], record['RLEVEL'], record['DAY-OBS'],
                          record['frameid']] for record in records]

    t = Table(np.asarray(records_sanitized), names=('FILENAME', 'SITEID', 'INSTRUME', 'RLEVEL', 'DAY-OBS', 'frameid'))

    return t


def filename_to_archivepath_dict(filenametable, rootpath=ARCHIVE_ROOT):
    ''' Return a dictionary with camera -> list of FileIO-able path of imagers from an elastic search result.
        We are still married to /archive file names here - because reasons. Long term we should go away from that.
    '''

    cameras = set(filenametable['INSTRUME'])
    returndict = {}
    for camera in cameras:
        returndict[camera] = [['{}/{}/{}/{}/{}/{}'.format(rootpath, record['SITEID'], record['INSTRUME'],
                                                          record['DAY-OBS'], 'raw', record['FILENAME']),
                               record['frameid']] for record in filenametable[filenametable['INSTRUME'] == camera]]

        returndict[camera] = Table(np.asarray(returndict[camera]), names=['FILENAME', 'frameid'])
    return returndict


def download_from_archive(frameid):
    """
    Download a file from the LCO archive by frame id.
    param frameid: Archive API frame ID
    return: Astropy HDUList
    """
    url = f'https://archive-api.lco.global/frames/{frameid}'
    log.info("Downloading image frameid {} from URL: {}".format(frameid, url))
    headers = {'Authorization': 'Token {}'.format(ARCHIVE_API_TOKEN)}
    response = requests.get(url, headers=headers)
    response.raise_for_status()
    response_dict = response.json()
    if response_dict == {}:
        log.warning("No file url was returned from id query")
        raise Exception('Could not find file remotely.')
    frame_url = response_dict['url']
    log.debug(frame_url)
    file_response = requests.get(frame_url)
    file_response.raise_for_status()
    f = fits.open(io.BytesIO(file_response.content))
    return f


if __name__ == '__main__':

    camera = 'fa15'
    dates = ArchiveDiskCrawler.get_last_n_days(20)

    for dayobs in dates:
        listofframes = get_frames_for_noisegainanalysis(dayobs, cameratype='fa')
        filelist = filename_to_archivepath_dict(listofframes)
        print("{} {} ".format(dayobs, filelist.keys()))

    # c = ArchiveCrawler()
    # for dayobs in dates:
    #     listofframes = c.findfiles_for_camera_dates("/archive/engineering/lsc/{}".format(camera), dayobs, 'raw', "*[xbf]00.fits*")
    #     #print ("{} {} ".format (dayobs, listofframes))
