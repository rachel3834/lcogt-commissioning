import datetime
import glob

from lcocommissioning.common.common import lco_site_lonlat

ARCHIVE_ROOT="/archive/engineering"


cameras="fa?? fs?? kb??"               # Camera types to crawl.
sites="elp cpt lsc ogg coj tfn"        # Which sites to crawl
dates="2018121?"

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


    def get_last_N_days (self, lastNdays):
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


if __name__ == '__main__':
    c = ArchiveCrawler()

    dates = c.get_last_N_days(2)
    cameras = c.find_cameras(['elp',])
    for camera in cameras:
        for date in dates:
            files = c.findfiles_for_camera_dates(camera, date, 'raw', "*[xbf]00.fits*")
            print (files)




