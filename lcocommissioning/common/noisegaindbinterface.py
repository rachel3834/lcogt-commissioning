import logging
import sqlite3
import numpy as np
from astropy.table import Table
import astropy.time as astt

_logger = logging.getLogger(__name__)


class noisegaindbinterface:
    ''' Storage model for flat field based noise gain measurements'''

    createstatement = "CREATE TABLE IF NOT EXISTS noisegain (" \
                      "name TEXT PRIMARY KEY, " \
                      "dateobs text," \
                      " camera text," \
                      " filter text," \
                      " readmode text,"\
                      " extension integer," \
                      " gain real," \
                      " readnoise real," \
                      " level real," \
                      " differencenoise real," \
                      " level1 real," \
                      " level2 real);"


    def __init__(self, fname):
        _logger.debug("Open data base file %s" % (fname))
        self.sqlite_file = fname
        self.conn = sqlite3.connect(self.sqlite_file)
        self.conn.execute(self.createstatement)
        self.conn.execute("PRAGMA journal_mode=WAL;")
        self.conn.commit()

    def addmeasurement(self, identifier, dateobs, camera, filter, extension, gain, readnoise, level, diffnoise, level1,
                       level2, readmode='default', commit=True):
        with self.conn:
            _logger.debug("Inserting: %s\n %s %s %s %s %s %s %s %s %s %s" % (
                identifier, dateobs, camera, filter, extension, gain, readnoise, level, diffnoise, level1, level2, readmode))
            self.conn.execute("insert or replace into noisegain values (?,?,?,?,?,?,?,?,?,?,?)",
                              (identifier, dateobs, camera, filter, extension, gain, readnoise, level, float(diffnoise),
                               float(level1), float(level2), readmode))

        if (commit):
            self.conn.commit()

    def checkifalreadyused(self, flat1):

        """ Check if a noisegain measurement based on two flat field expsoures already exists or not"""

        with self.conn:
            query = 'select name from noisegain where (name like ?)'
            cursor = self.conn.execute(query, ("%{}%".format(flat1),))
            allmatch = cursor.fetchall()
            if len(allmatch) > 0:
                _logger.debug("match found for %s" % (flat1))
                return True

        _logger.debug("no match found for %s" % (flat1))
        return False

    def getcameras(self):
        query = "select distinct camera from noisegain"

        cursor = self.conn.execute(query)
        retarray = []

        allcameras = (cursor.fetchall())
        if len(allcameras) == 0:
            _logger.warning("Zero results returned from query")

        for c in allcameras:
            retarray.append(c[0])
        _logger.debug("Distinct cameras: %s" % retarray)
        return retarray


    def get_readmodes_for_cameras(self, camera):
        query = "select distinct readmode from noisegain where camera like (?)"
        cursor = self.conn.execute(query,[camera,])
        retarray = []
        allreadmodes = (cursor.fetchall())
        if len(allreadmodes) == 0:
            _logger.warning("Zero results returned from query")
        for c in allreadmodes:
            retarray.append(c[0])
        _logger.debug("Distinct readmodes for camera %s: %s" % (camera, retarray))
        return retarray

    def readmeasurements(self, camera=None, filters=None, readmode=None, levelratio=None):
        """

        :param camera: name of the camera to read
        :param filters: array of filters to use. None if no filter selection
        :param levelratio: maximum ratio how much the two flat field levels may vary
        :return: astropy.table with query results. May be None if no results are returned.
        """

        query = "select name,dateobs,camera,filter,extension,gain,readnoise,level,differencenoise,level1,level2,readmode from noisegain " \
                "where (camera like ?) __FILTER__ __READMODE__ ORDER BY dateobs"

        queryargs = [camera if camera is not None else '%', ]

        if filters is not None:
            filtercondition = 'AND (filter in (%s))' % ','.join('?' * len(filters))
            query = query.replace("__FILTER__", filtercondition)
            queryargs.extend(filters)
        else:
            query = query.replace("__FILTER__", "")

        if readmode is not None:
            readmodecondition = 'AND (readmode like (?))'
            query = query.replace("__READMODE__", readmodecondition)
            queryargs.extend(readmode)
        else:
            query = query.replace("__READMODE__", "")

        _logger.debug("Read from database with query\n  %s\n  %s\n" % (query, queryargs))

        cursor = self.conn.execute(query, queryargs)

        allrows = np.asarray(cursor.fetchall())
        if len(allrows) == 0:
            _logger.warning("Zero results returned from query")
            return None

        t = Table(allrows,
                  names=['identifier', 'dateobs', 'camera', 'filter', 'extension', 'gain', 'readnoise', 'level',
                         'diffnoise', 'level1', 'level2','readmode'])
        t['dateobs'] = t['dateobs'].astype(str)
        t['dateobs'] = astt.Time(t['dateobs'], scale='utc', format=None).to_datetime()
        t['gain'] = t['gain'].astype(float)
        t['level'] = t['level'].astype(float)
        t['diffnoise'] = t['diffnoise'].astype(float)
        t['readnoise'] = t['readnoise'].astype(float)
        t['level1'] = t['level1'].astype(float)
        t['level2'] = t['level2'].astype(float)

        if levelratio is not None:
            t = t[abs((t['level1'] - t['level2']) / t['level']) < levelratio]

        return t

    def close(self):
        _logger.debug("Closing data base file %s " % (self.sqlite_file))
        self.conn.commit()
        self.conn.close()
