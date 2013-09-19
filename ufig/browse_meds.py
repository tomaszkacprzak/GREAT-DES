import os
import pyfits
import logging
import sys
import numpy
import sys
import math
import argparse
import yaml
import meds
import pylab

if __name__ == "__main__":

    description = 'quick and simple browse for meds files'

    # parse arguments
    parser = argparse.ArgumentParser(description=description, add_help=True)
    parser.add_argument('filepath_meds', type=str, help='input meds file')
    parser.add_argument('object_id', type=int, help='number of the object to view')
    parser.add_argument('-v', '--verbosity', type=int, action='store', default=2, choices=(0, 1, 2, 3 ), help='integer verbosity level: min=0, max=3 [default=2]')
    
    args = parser.parse_args()

    # Parse the integer verbosity level from the command line args into a logging_level string
    logging_levels = { 0: logging.CRITICAL, 
                       1: logging.WARNING,
                       2: logging.INFO,
                       3: logging.DEBUG }
    logging_level = logging_levels[args.verbosity]
    logging.basicConfig(format="%(message)s", level=logging_level, stream=sys.stdout)
    global logger 
    logger = logging.getLogger("add_weight") 
    config_logging_level = logging_levels[args.verbosity - 1]
    
    m = meds.MEDS(args.filepath_meds)
    logger.debug(m)
    logger.info('opened meds file with size %d' % (len(m._cat)))

    images = m.get_mosaic(args.object_id,type='image')
    weights = m.get_mosaic(args.object_id,type='weight')
    segs = m.get_mosaic(args.object_id,type='seg')

    pylab.subplot(1,3,1)
    pylab.imshow(images,interpolation='nearest')
    pylab.colorbar()
    pylab.subplot(1,3,2)
    pylab.imshow(weights,interpolation='nearest')
    pylab.colorbar()
    pylab.subplot(1,3,3)
    pylab.imshow(segs,interpolation='nearest')
    pylab.colorbar()
    pylab.show()


    image1 = m.get_cutout(args.object_id,0,type='image')
    image2 = m.get_cutout(args.object_id,1,type='image')
    pylab.subplot(1,3,1)
    pylab.imshow(image1,interpolation='nearest')
    pylab.colorbar()
    pylab.subplot(1,3,2)
    pylab.imshow(image2,interpolation='nearest')
    pylab.colorbar()
    pylab.show()
    