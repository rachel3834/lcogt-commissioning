import logging
import sqlite3
import numpy as np
from astropy.table import Table
import astropy.time as astt

_logger = logging.getLogger(__name__)


class lcoimagerpropertydatabase:
    createstatement = None
    indexstatements = None
    tablename = None

    def __init__(self, fname):
        _logger.debug("Open data base file %s" % (fname))
        self.sqlite_file = fname
        self.conn = sqlite3.connect(self.sqlite_file)

        if self.createstatement is not None:
            self.conn.execute(self.createstatement)
            self.conn.execute("PRAGMA journal_mode=WAL;")
            self.conn.commit()

        if self.indexstatements is not None:
            for statement in self.indexstatements:
                self.conn.execute(statement)

            self.conn.commit()

    def close(self):
        _logger.debug("Closing data base file %s " % (self.sqlite_file))
        self.conn.commit()
        self.conn.close()

    def checkifalreadyused(self, fitsimage, tablename):

        """ Check if a noisegain measurement based on two flat field expsoures already exists or not"""

        with self.conn:
            query = 'select name from {} where (name like %?%)'.format(tablename)
            cursor = self.conn.execute(query, ("%{}%".format(fitsimage),))
            allmatch = cursor.fetchall()
            if len(allmatch) > 0:
                return True
        return False

    def getcameras(self):
        if self.tablename is None:
            return []
        query = "select distinct camera from {}".format(self.tablename)

        cursor = self.conn.execute(query)
        retarray = []

        allcameras = (cursor.fetchall())
        if len(allcameras) == 0:
            _logger.warning("Zero results returned from query")

        for c in allcameras:
            retarray.append(c[0])
        _logger.debug("Distinct cameras: %s" % retarray)
        return retarray


class darkcurrentdbinterface(lcoimagerpropertydatabase):
    tablename = 'darkcurrent'
    createstatement = "CREATE TABLE IF NOT EXISTS darkcurrent (" \
                      "name TEXT PRIMARY KEY, " \
                      "dateobs text," \
                      " camera text," \
                      " readmode text," \
                      " darkcurrent real" \
                      " );"
    indexstatements = ["CREATE INDEX IF NOT EXISTS camera_idx on darkcurrent(camera);",]

    def checkifalreadyused(self, fitsimage):
        return super().checkifalreadyused(fitsimage, 'darkcurrent')

    def addmeasurement(self, identifier, dateobs, camera, darkcurrent, readmode='default', commit=True):
        with self.conn:
            _logger.debug("Inserting: %s\n %s %s %s %s" % (
                identifier, dateobs, camera, darkcurrent, readmode))
            self.conn.execute("insert or replace into darkcurrent values (?,?,?,?,?)",
                              (identifier, dateobs, camera,  readmode, float(darkcurrent)))
        if (commit):
            self.conn.commit()


    def readmeasurements(self, camera=None, readmode=None):
        """

        :param camera: name of the camera to read
        :param filters: array of filters to use. None if no filter selection
        :param levelratio: maximum ratio how much the two flat field levels may vary
        :return: astropy.table with query results. May be None if no results are returned.
        """

        query = "select name,dateobs,camera,darkcurrent,readmode from darkcurrent " \
                "where (camera like ?)  __READMODE__ ORDER BY dateobs"

        queryargs = [camera if camera is not None else '%', ]
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
              names=['identifier', 'dateobs', 'camera', 'darkcurrent', 'readmode'])
        t['dateobs'] = t['dateobs'].astype(str)
        t['dateobs'] = astt.Time(t['dateobs'], scale='utc', format=None).to_datetime()
        t['darkcurrent'] = t['darkcurrent'].astype(float)
        return t


