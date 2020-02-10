import logging
import sys
import datetime
import numpy as np
from astropy.table import Table
import astropy.time as astt
import math
from sqlalchemy.orm import sessionmaker
from sqlalchemy_utils import database_exists, create_database
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String, Float, create_engine

assert sys.version_info >= (3, 5)
_logger = logging.getLogger(__name__)

Base = declarative_base()

class NoiseGainMeasurement(Base):
    __tablename__ = 'noisegain'

    def __init__(self, rec):
        self.name=rec.name
        self.dateobs = rec.dateobs
        self.camera = rec.camera
        self.filter = rec. filter
        self.extension = rec.extension
        self.gain = rec.gain
        self.radnoise = rec.readnoise
        self.level = rec.level
        self.differencenoise = rec.differencenoise
        self.level1 = rec.level1
        self.level2 = rec.level2
        self.readmode = rec.readmode

    name = Column(String, primary_key=True)
    dateobs = Column(String)
    camera = Column(String, index=True)
    filter = Column(String, index=True)
    extension = Column(Integer)
    gain = Column(Float)
    readnoise = Column(Float)
    level = Column(Float)
    differencenoise = Column(Float)
    level1 = Column(Float)
    level2 = Column(Float)
    readmode = Column(String, index=True)

    def __repr__(self):
        return f'{self.name} {self.filter} {self.gain} {self.readnoise}'


class noisegaindb():
    def __init__(self, fname):
        _logger.debug("Open data base file %s" % (fname))
        self.engine = create_engine(f'{fname}', echo=False)
        if not database_exists(self.engine.url):
            create_database(self.engine.url)
        NoiseGainMeasurement.__table__.create(bind=self.engine, checkfirst=True)
        self.session = sessionmaker(bind=self.engine)()

    def close(self):
        """ Close the database safely"""
        _logger.debug("Closing data base session")
        self.session.close()

    def exists(self, name):
        return self.session.query(NoiseGainMeasurement).filter_by(name=name).first()

    def getCameras(self):
        q = self.session.query(NoiseGainMeasurement.camera).distinct()
        allrows = np.asarray([e.camera for e in q.all()])
        return allrows

    def getReadmodesFroCamera(self, camera):
        q = self.session.query(NoiseGainMeasurement.readmode).filter(NoiseGainMeasurement.camera == camera).distinct()
        allrows = np.asarray([e.readmode for e in q.all()])
        return allrows

    def addMeasurement(self, m, commit=True):
        existingEntry = self.exists(m.name)
        if (existingEntry):
            existingEntry.dateobs = m.dateobs
            # TODO: etc...
        else:
            self.session.add(m)
        if commit:
            self.session.commit()

    def getMeasurementsForCamera(self, camera=None, readmode=None, filters=None, extension=None, levelratio=None):

        q = self.session.query(NoiseGainMeasurement).filter(NoiseGainMeasurement.camera == camera)
        if readmode is not None:
            q = q.filter(NoiseGainMeasurement.readmode.in_(readmode))
        if filters is not None:
            q = q.filter(NoiseGainMeasurement.filter.in_(filters))
        if extension is not None:
            q = q.filter(NoiseGainMeasurement.extension == extension)

        allrows = [[e.name, e.dateobs, e.camera, e.filter, e.extension, e.gain, e.readnoise, e.level, e.differencenoise,
                    e.level1, e.level2, e.readmode] for e in q.all()]

        t = Table(np.asarray(allrows),
                  names=['identifier', 'dateobs', 'camera', 'filter', 'extension', 'gain', 'readnoise', 'level',
                         'diffnoise', 'level1', 'level2', 'readmode'])

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

    def checkifalreadyused(self, filename):
        q = self.query(NoiseGainMeasurement.name).filter (NoiseGainMeasurement.name.like(f'%{filename}%'))
        return q.all().count()


if __name__ == '__main__':
    # TODO: Move this stuff into a test routine
    c = noisegaindb('noisegain.sqlite')
    print(c.getCameras())
    print(c.getReadmodesFroCamera('fa03'))
    print(c.getMeasurementsForCamera('fa03', readmode= 'full_frame', filter=['gp','rp']))
    c.close()
