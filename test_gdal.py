from osgeo import gdal

from ekskursja import track as tr
from ekskursja import lightr as lr

import logging
FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
logging.basicConfig(level=logging.DEBUG, format=FORMAT)

#initiate logger
logger = logging.getLogger(__name__)


version_num = int(gdal.VersionInfo('VERSION_NUM'))
logger.info('gdal.version: %s'%version_num)

t = tr.Track()
t.readGpx('data/170705-Miedzywodzie.gpx')

lr.applyMileage('data/170705-photos.txt', t)