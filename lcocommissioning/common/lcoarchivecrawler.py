import datetime
import glob
import logging
import numpy as np
from astropy.table import Table

from lcocommissioning.common.common import lco_site_lonlat
from elasticsearch import Elasticsearch
from elasticsearch_dsl import Search


log = logging.getLogger(__name__)

ARCHIVE_ROOT="/archive/engineering"



class ArchiveCrawler:

    archive_root = None

    def __init__(self, rootdirectory = ARCHIVE_ROOT):
        self.archive_root = rootdirectory



    def find_cameras (self, sites=lco_site_lonlat, cameras=["fa??", "fs??", "kb??"]):
        sitecameras = []

        for site in sites:
            for camera in cameras:
                dir = "{}/{}/{}".format (self.archive_root, site, camera)
                candiates=glob.glob (dir)
                sitecameras.extend (candiates)
        return sitecameras


    @staticmethod
    def get_last_N_days (lastNdays):
        date=[]
        today = datetime.datetime.utcnow()
        for ii in range (lastNdays):
            day = today - datetime.timedelta(days=ii)
            date.append (day.strftime("%Y%m%d"))
        return date[::-1]

    @staticmethod
    def findfiles_for_camera_dates (sitecamera, date, raworprocessed, filetempalte):
        dir="{}/{}/{}/{}".format ( sitecamera, date, raworprocessed, filetempalte)
        files = glob.glob (dir)
        return (files)



def make_elasticsearch(index, filters, queries=None, exclusion_filters=None, range_filters=None, prefix_filters=None,
                       terms_filters=None,
                       es_url='http://elasticsearch.lco.gtn:9200'):
    """
    Make an ElasticSearch query

    Parameters
    ----------
    index : str
            Name of index to search
    filters : list of dicts
              Each dict has a criterion for an ElasticSearch "filter"
    queries : list of dicts
              Each dict has a "type" and "query" entry. The 'query' entry is a dict that has a criterion for an
              ElasticSearch "query"
    exclusion_filters : list of dicts
                        Each dict has a criterion for an ElasticSearch "exclude"
    range_filters: list of dicts
                   Each dict has a criterion an ElasticSearch "range filter"
    es_url : str
             URL of the ElasticSearch host

    Returns
    -------
    search : elasticsearch_dsl.Search
             The ElasticSearch object
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
    es = Elasticsearch(es_url)
    s = Search(using=es, index=index)
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


def get_frames_for_noisegainanalysis(dayobs, site=None, cameratype=None, camera=None, mintexp=30,
                                     filterlist=['gp', 'rp', 'ip', 'zp'], obstype=['BIAS','SKYFLAT'],  es_url='http://elasticsearch.lco.gtn:9200'):

    """ Queries for a list of processed LCO images that are viable to get a photometric zeropoint in the griz bands measured.

        Selection criteria are by DAY-OBS, site, by camaera type (fs,fa,kb), what filters to use, and minimum exposure time.
        Only day-obs is a mandatory fields, we do not want to query the entire archive at once.
     """

    query_filters = [{'DAY-OBS': dayobs}, {'RLEVEL': 0}, ]
    range_filters = []
    terms_filters = [ {'OBSTYPE': obstype}]
    prefix_filters = []

    if site is not None:
        query_filters.append({'SITEID': site})
    if camera is not None:
        query_filters.append({'INSTRUME': camera})
    if cameratype is not None:
        prefix_filters.append({'INSTRUME': cameratype})

    queries = []
    records = make_elasticsearch('fitsheaders', query_filters, queries, exclusion_filters=None, es_url=es_url,
                                 range_filters=range_filters, prefix_filters=prefix_filters,
                                 terms_filters=terms_filters).scan()
    records_sanitized = [[record['filename'], record['SITEID'], record['INSTRUME'], record['RLEVEL'], record['DAY-OBS']]
                         for record in records if (record['FILTER'] in filterlist) or (record['OBSTYPE']=='BIAS')]

    t = Table (np.asarray(records_sanitized), names=('FILENAME','SITEID','INSTRUME','RLEVEL','DAY-OBS'))
    return t


def filename_to_archivepath (filenametable, rootpath='/archive/engineering'):
    '''return a dictionary with camera -> list of FileIO-able path of iamgers
    '''
    cameras = set (filenametable['INSTRUME'])
    returndict = {}
    for camera in cameras:
        returndict[camera] = ['{}/{}/{}/{}/{}/{}'.format(rootpath, record['SITEID'],record['INSTRUME'],record['DAY-OBS'],'raw',record['FILENAME']) for record in filenametable[filenametable['INSTRUME']== camera] ]
    return returndict


if __name__ == '__main__':

    camera = 'fa15'
    dates = ArchiveCrawler.get_last_N_days(20)

    for dayobs in dates:
        listofframes = get_frames_for_noisegainanalysis(dayobs, cameratype='fa')
        filelist = filename_to_archivepath(listofframes)
        print ("{} {} ".format (dayobs, filelist.keys()))

    # c = ArchiveCrawler()
    # for dayobs in dates:
    #     listofframes = c.findfiles_for_camera_dates("/archive/engineering/lsc/{}".format(camera), dayobs, 'raw', "*[xbf]00.fits*")
    #     #print ("{} {} ".format (dayobs, listofframes))




