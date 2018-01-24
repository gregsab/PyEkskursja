'''
Created on 6 lip 2015

@author: sabakgrz
'''

from osgeo import ogr
from osgeo import osr
from osgeo import gdal

from datetime import datetime
from dateutil import tz
from math import floor

import os.path

import tools

import logging
FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
logging.basicConfig(level=logging.DEBUG, format=FORMAT)

#initiate logger
logger = logging.getLogger(__name__)


#initiate timezone global variables
tz_utc = tz.tzutc()
tz_wawa = tz.tzlocal()


track_colors = ["#4385f5", "#d65777", "#fe6500", "#8465ab", "#eeb4d1", "#51bb85"]
track_colors_karkonosze = ["#1b57e2", "#2acadf", "#e31a1c", "#e9e60c", "#1b9d5a", "#0dd538", "#6767d8"]

DEFAULT_LINE_WIDTH = 6
DEFAULT_POINT_SIZE = 6
DEFAULT_START_STOP_SIZE = 8

'''
Converts geometry from WGS84 projection to Poland 92 


Args:
    geom (ogr.Geometry) : Point in  WGS84 projection
    
Returns:
    arg.Geomtery: Point in EPSG2180 projection

'''


def _joinLineStrings(geom):
    logger.info('No of geometries to join: %d' % geom.GetGeometryCount())
    
    lstring = ogr.Geometry(ogr.wkbLineString)
    
    for i in range(0, geom.GetGeometryCount()):
        g = geom.GetGeometryRef(i)
        
        for j in range(0, g.GetPointCount()):
            pt = g.GetPoint(j)
            lstring.AddPoint(pt[0], pt[1])
                
    logger.info('%d points joined to one LineString' % lstring.GetPointCount())
                
    return lstring    

def _saveTrack(filename, geom, points, popup, distance, color="#4385f5"):
    driverName = "GeoJSON"
    driver = ogr.GetDriverByName(driverName)

    if driver is None:
        logger.error("%s driver not available." % driverName)
        return        
    # Create the output GeoJSON
    outDataSource = driver.CreateDataSource(filename)
    
    if outDataSource is None:
        logger.error('Could not create file %s' % filename)
        return
    
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)      
    
    outLayer = outDataSource.CreateLayer('track', srs, ogr.wkbLineString, ["COORDINATE_PRECISION=6","WRITE_BBOX=YES"] )

    _createLayerFields(outLayer)

    _saveSimpleTrack(outLayer, geom, srs, popup, color)
    _saveMileage(outLayer, points, srs, distance, color)

    outDataSource.Destroy()

def _createLayerFields(outLayer):
    field_color = ogr.FieldDefn("color", ogr.OFTString)
    field_color.SetWidth(12)
    outLayer.CreateField(field_color)
    outLayer.CreateField(ogr.FieldDefn("opacity", ogr.OFTReal))
    outLayer.CreateField(ogr.FieldDefn("radius", ogr.OFTReal))
    field_fillColor = ogr.FieldDefn("fillColor", ogr.OFTString)
    field_fillColor.SetWidth(12)
    outLayer.CreateField(field_fillColor)
    outLayer.CreateField(ogr.FieldDefn("weight", ogr.OFTReal))
    outLayer.CreateField(ogr.FieldDefn("opacity", ogr.OFTReal))
    outLayer.CreateField(ogr.FieldDefn("fillOpacity", ogr.OFTReal))
    field_popup = ogr.FieldDefn("popup", ogr.OFTString)
    field_popup.SetWidth(12)
    outLayer.CreateField(field_popup)


def _saveSimpleTrack(outLayer, geom, srs, popup, color="#4385f5"):
    
    outFeature = ogr.Feature(outLayer.GetLayerDefn())

    outFeature.SetField("color", color)
    outFeature.SetField("weight", DEFAULT_LINE_WIDTH)
    outFeature.SetField("opacity", 0.9)
    outFeature.SetField("popup", popup)
        
    # Set new geometry
    outFeature.SetGeometry(tools.from2180to4326(geom))
    
    # Add new feature to output Layer
    outLayer.CreateFeature(outFeature)
    outFeature.Destroy()
    
#     logger.info('Track saved to file %s.' % filename)


def _saveMileage(outLayer, points, srs, distance, color="#4385f5"):
    color_first = "#ffe168"
    color_last = "#f48159"

    for i, pt in enumerate(points):
        rr = DEFAULT_LINE_WIDTH-2
        cc = color
        fc = color_first
         
        pp = "%dkm" % (i*distance/1000)
         
        if i == 0:
            cc = fc = color_first
            rr = DEFAULT_START_STOP_SIZE
            pp = "start"
        elif i == len(points) - 1:
            cc = fc = color_last
            rr = DEFAULT_START_STOP_SIZE
            pp = "koniec"
        
        outFeature = ogr.Feature(outLayer.GetLayerDefn())
    
        outFeature.SetField("radius", rr)
        outFeature.SetField("color", cc)
        outFeature.SetField("fillColor", fc)
        outFeature.SetField("weight", 2)
        outFeature.SetField("opacity", 0.9)
        outFeature.SetField("fillOpacity", 0.9)

        outFeature.SetField("popup", pp)

            
        # Set new geometry
        outFeature.SetGeometry(tools.from2180to4326(pt))
        
        # Add new feature to output Layer
        outLayer.CreateFeature(outFeature)
        outFeature.Destroy

    logger.info('%d mileage points written.' % len(points))


def _pointAtDistance(lstring, distance):
    if (distance < 0):
        return None

    from_start=0

    p0 = ogr.Geometry(ogr.wkbPoint)
    pp = lstring.GetPoint(0)
    
    p0.AddPoint(pp[0], pp[1])

    for i in range(1, lstring.GetPointCount()-1):
        p1 = ogr.Geometry(ogr.wkbPoint)
        pq = lstring.GetPoint(i)
        p1.AddPoint(pq[0], pq[1])
        
        if (distance - from_start < 0.01):
            return p0.Clone()
        
        d = p0.Distance(p1)
        
        if (from_start <= distance and from_start+d >= distance):
            r = (distance - from_start) / d
            
            xr = p0.GetX() + r * (p1.GetX() - p0.GetX())
            yr = p0.GetY() + r * (p1.GetY() - p0.GetY())
            
            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(xr, yr)
            
            return point
        
        from_start = from_start + d
        
        p0 = p1
    
    return None
    
    
def _mileage(lstring, distance=1000):
    points = []

    d = 0
    
    p = _pointAtDistance(lstring, d)

    while (p is not None):
        points.append(p)
        d = d + distance
        p = _pointAtDistance(lstring, d)
            
    p0 = lstring.GetPoint(lstring.GetPointCount()-1)
    p = ogr.Geometry(ogr.wkbPoint)
    p.AddPoint(p0[0], p0[1])
    points.append(p)
            
    return points

def _openGPX(filename):
    pass


def processGarminBaseCamp(filename, popup = "", colors="#4385f5", distance=10000):
    version_num = int(gdal.VersionInfo('VERSION_NUM'))

    logger.info("GDAL version: %s" % version_num)

    driverName = "GPX"
    driver = ogr.GetDriverByName(driverName)

    if driver is None:
        logger.error("%s driver not available." % driverName)
        return

    dataSource = driver.Open(filename, 0)

    if dataSource is None:
        logger.error('Could not open %s.' % filename)
        return
    else:
        logger.info('File opened %s.' % filename)
        layer = dataSource.GetLayerByName("tracks")

    i = 0

    for feature in layer:
        geom = feature.GetGeometryRef()

        geom2 = _joinLineStrings(geom)

        geom2180 = tools.from4326to2180(geom2)

        logger.info('%d points before simplification.' % geom2180.GetPointCount())
        geom2180s = geom2180.Simplify(10)
        logger.info('%d points after simplification.' % geom2180s.GetPointCount())

        mileage = _mileage(geom2180s, distance)

        fname = "%s-%d.geojson" % (os.path.splitext(filename)[0], i)

        _saveTrack(fname, geom2180s, mileage, popup, distance, colors[i])
        i = i+1

def process4mileage(filename, popup = "", distance = 1000, color="#4385f5"):
    version_num = int(gdal.VersionInfo('VERSION_NUM'))

    logger.info("GDAL version: %s" % version_num)

    driverName = "GPX"
    driver = ogr.GetDriverByName(driverName)

    if driver is None:
        logger.error("%s driver not available." % driverName)
        return

    dataSource = driver.Open(filename, 0)

    if dataSource is None:
        logger.error('Could not open %s.' % filename)
        return
    else:
        logger.info('File opened %s.' % filename)
        layer = dataSource.GetLayerByName("tracks")
    for feature in layer:
        geom = feature.GetGeometryRef()

        geom2 = _joinLineStrings(geom)

        geom2180 = tools.from4326to2180(geom2)

        logger.info('%d points before simplification.' % geom2180.GetPointCount())
        geom2180s = geom2180.Simplify(10)
        logger.info('%d points after simplification.' % geom2180s.GetPointCount())

        mileage = _mileage(geom2180s, distance)

        _saveTrack(os.path.splitext(filename)[0]+".geojson", geom2180s, mileage, popup, distance, color)



class TrackPoint(object):
    '''
    classdocs
    '''
    def distance(self, tpoint, base_csr=2180):
        if base_csr == 3035:
            return self.point3035.Distance(tpoint.point3035)
        else:
            return self.point2180.Distance(tpoint.point2180)
    
    def __init__(self, point, tstamp, height=0):
        self.point4326 = point
        self.point2180 = tools.from4326to2180(point)
        self.point3035 = tools.from4326to3035(point)
        self.height = height
                
        temp = datetime(tstamp[0], tstamp[1], tstamp[2], tstamp[3], tstamp[4], int(tstamp[5]), 0, tz_utc)
        self.timestamp = temp.astimezone(tz_wawa)
        
        pass
    
    def __str__(self):
        buf = '%s, %s, %s, %s' % (self.timestamp, self.point4326, self.point2180, self.height)
        return buf
    
    pass

class Track(object):
    '''
    classdocs
    '''
    def __init__(self):
        self.trackpoints = None
        self.tracks = None


    @staticmethod
    def saveSlices2GeoJson(filename, sliced_track):

        # Create the output Driver
        outDriver = ogr.GetDriverByName('GeoJSON')
        
        # Create the output GeoJSON
        outDataSource = outDriver.CreateDataSource(filename)
        
        if outDataSource is None:
            print 'Could not create file'
        
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)        
        
        outLayer = outDataSource.CreateLayer('track', srs, ogr.wkbLineString )
        
        field_color = ogr.FieldDefn("color", ogr.OFTString)
        field_color.SetWidth(12)
        outLayer.CreateField(field_color)
        outLayer.CreateField(ogr.FieldDefn("opacity", ogr.OFTReal))
        outLayer.CreateField(ogr.FieldDefn("section", ogr.OFTInteger))
        
        i=0;
        for line in sliced_track:
            # create a new feature
            outFeature = ogr.Feature(outLayer.GetLayerDefn())

            outFeature.SetField("opacity", 0.8)
            outFeature.SetField("section", i%2==0)
            
            if (i%2==0):
                outFeature.SetField("color", '#e51c55')
            else:
                outFeature.SetField("color", '#800f31')
                
            # Set new geometry
            outFeature.SetGeometry(tools.from2180to4326(line))
            
            # Add new feature to output Layer
            outLayer.CreateFeature(outFeature)
            outFeature.Destroy
            
            i += 1
        
        # Close DataSources
        outDataSource.Destroy()


    @staticmethod
    def saveMileage2GeoJson(filename, mileage):
        driverName = "GeoJSON"
        driver = ogr.GetDriverByName(driverName)

        if driver is None:
            logger.error("%s driver not available." % driverName)
            return        
        # Create the output GeoJSON
        outDataSource = driver.CreateDataSource(filename)
        
        if outDataSource is None:
            logger.error('Could not create file %s' % filename)
            return
        
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)        
        
        outLayer = outDataSource.CreateLayer('mileage', srs, ogr.wkbLineString )
        
        field_color = ogr.FieldDefn("color", ogr.OFTString)
        field_color.SetWidth(12)
        outLayer.CreateField(field_color)
        outLayer.CreateField(ogr.FieldDefn("size", ogr.OFTReal))
        
        i=0;
        for point in mileage:
            # create a new feature
            outFeature = ogr.Feature(outLayer.GetLayerDefn())

            outFeature.SetField("color", '#e51c55')
            outFeature.SetField("size", 4)
                
            # Set new geometry
            outFeature.SetGeometry(tools.from2180to4326(point))
            
            # Add new feature to output Layer
            outLayer.CreateFeature(outFeature)
            outFeature.Destroy
            
            i += 1
        
        # Close DataSources
        outDataSource.Destroy()
        logger.info('%d mileage point written to %s' % (i, filename))


    def getFirstPoint(self):
        p1 = self.trackpoints[0]
        
        return p1
    
    def getLastPoint(self):
        return self.trackpoints[len(self.trackpoints)-1];

    def mileageByDist(self, distance=1000):
        multipoint = ogr.Geometry(ogr.wkbMultiPoint)
        
        multipoint.AddGeometry(self.getFirstPoint().point2180)
        
        multipoint.AddGeometry(self.getLastPoint().point2180)
        
        return multipoint


    def getPointAtDist(self, dist):
       
        if (dist < 0 or dist > self.length()):
            return None


        for i in range(len(self.trackpoints)-1):
            if (self.trackpoints[i].distFromStart <= dist and  self.trackpoints[i+1].distFromStart >= dist):
                
                delta = self.trackpoints[i+1].distFromStart - self.trackpoints[i].distFromStart
                
                if delta < 0.01:
                    return self.trackpoints[i].point2180.Clone()
                
                r = (dist - self.trackpoints[i].distFromStart) / delta
                tp1 = self.trackpoints[i].point2180
                tp2 = self.trackpoints[i+1].point2180
                xr = tp1.GetX() + r * (tp2.GetX() - tp1.GetX())
                yr = tp1.GetY() + r * (tp2.GetY() - tp1.GetY())
                
                point = ogr.Geometry(ogr.wkbPoint)
                point.AddPoint(xr, yr)
                
                return point
        
        return None

# Uwazac na timezone, zwraca punkt GPS
    def getPointAtTimestamp(self, tstamp):
        
        for i in range(len(self.trackpoints)-1):
            if (self.trackpoints[i].timestamp <= tstamp and  self.trackpoints[i+1].timestamp >= tstamp):
                
                delta = (self.trackpoints[i+1].timestamp - self.trackpoints[i].timestamp).total_seconds()
                
                if delta <= 1 :
                    return self.trackpoints[i].point2180.Clone()
                
                r = (tstamp - self.trackpoints[i].timestamp).total_seconds() / delta
                tp1 = self.trackpoints[i].point2180
                tp2 = self.trackpoints[i+1].point2180
                xr = tp1.GetX() + r * (tp2.GetX() - tp1.GetX())
                yr = tp1.GetY() + r * (tp2.GetY() - tp1.GetY())
                
                point = ogr.Geometry(ogr.wkbPoint)
                point.AddPoint(xr, yr)
                
                return tools.from2180to4326(point)
        
        return None    
    
    def saveProfile(self, filename=None, init_dist=0):
        f = open(filename,'w')
        for i in range(len(self.trackpoints)-3):
            tp0 = self.trackpoints[i]
            tp1 = self.trackpoints[i+2]
            
            dtime = (tp1.timestamp - tp0.timestamp).total_seconds()
            dist = tp0.point3035.Distance(tp1.point3035)
            
            v = dist * 3.6 / dtime
            
            f.write("%s, %f, %f, %f\n" % (self.trackpoints[i].timestamp.strftime("%Y-%m-%d %H:%M:%S"),
                                          (init_dist+self.trackpoints[i].distFromStart)/1000.0,
                                          self.trackpoints[i].height, v))
        f.close()
    
    def sliceByDist(self, slice_length):
        multiline = ogr.Geometry(ogr.wkbMultiLineString)
        
        for i in range(int(floor(self.length() / slice_length))+1):
            line = self._sliceByDist(i*slice_length, (i+1)*slice_length)
            
            line.Simplify(10)
            
            multiline.AddGeometry(line)
        
        return multiline    
    
    def length(self):
        if (self.trackpoints is None):
#             logger.warning('No trackpoints in track. Length = 0')
            return 0
        
        return self.trackpoints[len(self.trackpoints)-1].distFromStart
    
    def _sliceByDist(self, start, finish):
        line = ogr.Geometry(ogr.wkbLineString)
        
        p1 = self.getPointAtDist(start)
        line.AddPoint(p1.GetX(), p1.GetY())
        
        for tp in self.trackpoints:
            if (tp.distFromStart >= start and tp.distFromStart <= finish):
                line.AddPoint(tp.point2180.GetX(), tp.point2180.GetY())

        p2 = self.getPointAtDist(finish)
        
        if p2 is not None:
            line.AddPoint(p2.GetX(), p2.GetY())

                                
        return line

    def _layer_meta(self, layer):   
        capabilities = [
            ogr.OLCRandomRead,
            ogr.OLCSequentialWrite,
            ogr.OLCRandomWrite,
            ogr.OLCFastSpatialFilter,
            ogr.OLCFastFeatureCount,
            ogr.OLCFastGetExtent,
            ogr.OLCCreateField,
            ogr.OLCDeleteField,
            ogr.OLCReorderFields,
            ogr.OLCAlterFieldDefn,
            ogr.OLCTransactions,
            ogr.OLCDeleteFeature,
            ogr.OLCFastSetNextByIndex,
            ogr.OLCStringsAsUTF8,
            ogr.OLCIgnoreFields
        ]

        print("Layer Capabilities:")
        for cap in capabilities:
            print("  %s = %s" % (cap, layer.TestCapability(cap)))

        print
        
        lyrDefn = layer.GetLayerDefn()
        for i in range( lyrDefn.GetFieldCount() ):
            fieldName =  lyrDefn.GetFieldDefn(i).GetName()
            fieldTypeCode = lyrDefn.GetFieldDefn(i).GetType()
            fieldType = lyrDefn.GetFieldDefn(i).GetFieldTypeName(fieldTypeCode)
            fieldWidth = lyrDefn.GetFieldDefn(i).GetWidth()
            GetPrecision = lyrDefn.GetFieldDefn(i).GetPrecision()

            print fieldName + " - " + fieldType+ " " + str(fieldWidth) + " " + str(GetPrecision)
            
        print

    
    def readGpx(self, filename, base_csr=2180):
        driverName = "GPX"
        driver = ogr.GetDriverByName(driverName)

        if driver is None:
            logger.error("%s driver not available." % driverName)
            return

        dataSource = driver.Open(filename, 0)
        
        if dataSource is None:
            logger.error('Could not open %s.' % filename)
            return
        else:
            self.trackpoints = []
            layer = dataSource.GetLayerByName("track_points")

            distance = 0
            duration = 0
            lasttp = None

            i = 0
    
            for feature in layer:
                geom = feature.GetGeometryRef()
                
                idx = feature.GetFieldIndex("time") 
                tstamp = feature.GetFieldAsDateTime(idx)
                
                h = feature.GetFieldIndex("ele")
                height = feature.GetFieldAsDouble(h)
                            
                tp = TrackPoint(geom, tstamp, height)
                
                if (i == 0):
                    start = tp.timestamp
#                     logger.info("First location at: %s" % start)
                else:
                    d = lasttp.distance(tp, base_csr)
                    distance += d  
                    duration = tp.timestamp - start
                    
                tp.distFromStart = distance
                tp.durFromStart = duration
                
                self.trackpoints.append(tp)
                i += 1
                lasttp = tp
                
            logger.info('%d points read from ''%s''' % (i, os.path.basename(filename)))
            logger.info('Dist. from beginning [m]: %d ' % lasttp.distFromStart)

            
    def readGpxSegments(self, filename):
        driverName = "GPX"
        driver = ogr.GetDriverByName(driverName)

        if driver is None:
            logger.error("%s driver not available." % driverName)
            return

        dataSource = driver.Open(filename, 0)
        
        if dataSource is None:
            logger.error('Could not open %s.' % filename)
            return
        else:
            self.tracks = []
            self.track_names = []
            layer = dataSource.GetLayerByName("tracks")

            i = 0
    
            for feature in layer:
                geom = feature.GetGeometryRef()
                name = feature.GetFieldAsString("name")
                                
                self.tracks.append(geom.Clone())
                self.track_names.append(name)
                
                i += 1
                
            logger.info('%d tracks read from ''%s''' % (i, os.path.basename(filename)))
            
            
    def saveTracks(self, filename):    
        driverName = "GeoJSON"
        driver = ogr.GetDriverByName(driverName)

        if driver is None:
            logger.error("%s driver not available." % driverName)
            return        
        # Create the output GeoJSON
        outDataSource = driver.CreateDataSource(filename)
        
        if outDataSource is None:
            logger.error('Could not create file %s' % filename)
            return
        
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)        
        
        outLayer = outDataSource.CreateLayer('tracks', srs, ogr.wkbLineString, ["COORDINATE_PRECISION=6","WRITE_BBOX=YES"] )
        
        field_color = ogr.FieldDefn("color", ogr.OFTString)
        field_color.SetWidth(12)
        outLayer.CreateField(field_color)
        outLayer.CreateField(ogr.FieldDefn("opacity", ogr.OFTReal))
        outLayer.CreateField(ogr.FieldDefn("weight", ogr.OFTReal))
        outLayer.CreateField(ogr.FieldDefn("opacity", ogr.OFTReal))
        field_popup = ogr.FieldDefn("popup", ogr.OFTString)
        field_popup.SetWidth(128)
        outLayer.CreateField(field_popup)
                
        for i, t in enumerate(self.tracks):
            outFeature = ogr.Feature(outLayer.GetLayerDefn())

            outFeature.SetField("weight", DEFAULT_LINE_WIDTH)
            outFeature.SetField("opacity", 0.9)
            outFeature.SetField("color", track_colors_karkonosze[i])
            outFeature.SetField("popup", self.track_names[i])

            geom2180 = tools.from4326to2180(t)
            geom2180s = geom2180.Simplify(10)

            # Set new geometry
            outFeature.SetGeometry(tools.from2180to4326(geom2180s))
            
            # Add new feature to output Layer
            outLayer.CreateFeature(outFeature)
            outFeature.Destroy()

        logger.info('%d tracks saved to ''%s''' % (len(self.tracks), os.path.basename(filename)))   
        outDataSource.Destroy()

    def snap(self, x, y, acc = 100):
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(x, y)
        point2180 = tools.from4326to2180(point)

        s = -1
        d = float("inf")

        i = 0
        while i <= self.length():
            p = self.getPointAtDist(i)
            d1 = point2180.Distance(p)

            if d1 < d :
                d = d1
                s = i

            # logger.info('Distance at %d: %.1f' %(i, d1))

            i += acc

        return s