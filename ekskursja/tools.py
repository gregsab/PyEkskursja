from osgeo import osr

target3035 = osr.SpatialReference()
target3035.ImportFromEPSG(3035)


def from4326to2180(geom):
    source = osr.SpatialReference()
    source.ImportFromEPSG(4326)

    target = osr.SpatialReference()
    target.ImportFromEPSG(2180)

    transform = osr.CoordinateTransformation(source, target)

    gtemp = geom.Clone()

    gtemp.Transform(transform)

    return gtemp


def from4326to3035(geom):
    source = osr.SpatialReference()
    source.ImportFromEPSG(4326)

    transform = osr.CoordinateTransformation(source, target3035)

    gtemp = geom.Clone()
    gtemp.Transform(transform)

    return gtemp


def transform(geom, fromSRS=4326, toSRS=2180):
    source = osr.SpatialReference()
    source.ImportFromEPSG(fromSRS)

    target = osr.SpatialReference()
    target.ImportFromEPSG(toSRS)

    transform = osr.CoordinateTransformation(source, target)

    gtemp = geom.Clone()

    gtemp.Transform(transform)

    return gtemp


def from2180to4326(geom):
    source = osr.SpatialReference()
    source.ImportFromEPSG(2180)

    target = osr.SpatialReference()
    target.ImportFromEPSG(4326)

    transform = osr.CoordinateTransformation(source, target)

    gtemp = geom.Clone()

    gtemp.Transform(transform)

    return gtemp
