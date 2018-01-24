import os

import logging
FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
logging.basicConfig(level=logging.DEBUG, format=FORMAT)

#initiate logger
logger = logging.getLogger(__name__)

def applyMileage(fin, t):
    logger.info('applyMileage to file: %s', fin)
    fout = fin[:-4] + '-km.txt'
    logger.info('output file: %s', fout)

    # Open the file with read only permit
    f = open(fin, "r")
    fo = open(fout, "w")

    s = ";"

    # iterate over lines
    for line in f.readlines():
        ll = line.split(";")
        logger.info("%s", ll)
        x = float(ll[2])
        y = float(ll[3])
        km = t.snap(x,y)
        skm = "%.1f" % (km/1000.0)
        ll.append(ll[6])
        ll[6] = skm
        logger.info("%s", ll)
        fo.write(s.join(ll))

    # close the file after reading the lines.
    f.close()
    fo.close()