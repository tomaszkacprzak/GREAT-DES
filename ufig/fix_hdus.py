import os
import pyfits
import logging
import sys
import numpy
import sys
import math
import argparse
import yaml
import copy

# ufig
HDU_UFIG_IMAGE = 0
# coadd
COADD_HDU     = 2; # that is for image 
COADD_WT_HDU  = 3;
COADD_SEG_HDU = 1;
# single exp
SE_HDU        = COADD_HDU; # that is for image
SE_WT_HDU     = COADD_WT_HDU;
SE_BADPIX_HDU = 4;
# other
SKY_HDU       = 1; # in the bkg file
SEG_HDU       = 1; # in the seg file

# total
NHDU_OUT = 4

if __name__ == "__main__":

    description = ''

    # parse arguments
    parser = argparse.ArgumentParser(description=description, add_help=True)
    parser.add_argument('filepath_ufig', type=str, help='input ufig image')
    parser.add_argument('filepath_seg', type=str, help='input segmap image')
    parser.add_argument('filepath_out', type=str, help='fits file with added weight map')
    parser.add_argument('-v', '--verbosity', type=int, action='store', default=2, choices=(0, 1, 2, 3 ), help='integer verbosity level: min=0, max=3 [default=2]')
    parser.add_argument('--weight_value', type=float, action='store', default=0.01, help='weight map will be filled in with this number')
    
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
    
    hdus_ufig = pyfits.open(args.filepath_ufig)
    logger.debug(hdus_ufig)
    nx = hdus_ufig[HDU_UFIG_IMAGE].shape[0]
    ny = hdus_ufig[HDU_UFIG_IMAGE].shape[1]
    logger.info('ufig image size %d %d' % (nx,ny))

    hdus_seg = pyfits.open(args.filepath_seg)
    logger.info('opened seg file of size %s' % str(hdus_seg[SEG_HDU-1].shape))

    # get all the data we need
    image = hdus_ufig[HDU_UFIG_IMAGE].data
    weight = numpy.ones([nx,ny])*args.weight_value
    badpix = numpy.zeros([nx,nx],dtype=numpy.int32)
    # the convention for the hdu indices in settings iterates from 1, and in pyfits from 0, so subtract one
    seg = hdus_seg[SEG_HDU-1].data

    # get the output hdus, the convention for the hdu indices in settings iterates from 1, and in pyfits from 0, so subtract one
    out_list = [None]*NHDU_OUT
    out_list[0] = pyfits.PrimaryHDU(header=hdus_ufig[0].header)
    out_list[COADD_HDU-1]      = pyfits.ImageHDU(data=image,header=hdus_ufig[HDU_UFIG_IMAGE].header)
    out_list[COADD_SEG_HDU-1]  = pyfits.ImageHDU(data=seg,header=hdus_ufig[HDU_UFIG_IMAGE].header)
    out_list[COADD_WT_HDU-1]   = pyfits.ImageHDU(data=weight,header=hdus_ufig[HDU_UFIG_IMAGE].header)
    out_list[SE_BADPIX_HDU-1]  = pyfits.ImageHDU(data=badpix,header=hdus_ufig[HDU_UFIG_IMAGE].header)

    hdulist_out = pyfits.HDUList(out_list)
    logger.info('saving..')
    hdulist_out.writeto(args.filepath_out,clobber=True)
    logger.info('saved %s' % args.filepath_out)


