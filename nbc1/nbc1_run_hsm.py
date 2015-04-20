import numpy, galsim, sys, logging, yaml, argparse, time, pylab, copy, itertools, meds, pyfits
from nbc1_dtypes import *


ERR_FLAG = -99;

def image_array_to_galsim(array):


    nx,ny = array.shape
    img_gs = galsim.ImageD(ny,nx)
    for ix in range(nx):
        for iy in range(ny):
            img_gs.setValue(iy+1,ix+1,float(array[ix][iy]))

    return img_gs

def write_result(file_result,hsmresult=None,ig=None,header_only=False):

    if header_only:
        line= '#index,g1,g2,size,x0,y0'
    else:

        fmt_str = '%d\t' + '% f\t'*5 + '\n'

        if hsmresult.meas_type == 'g':
            hsmresult_corrected_g1 = hsmresult.corrected_g1
            hsmresult_corrected_g2 = hsmresult.corrected_g2
        elif hsmresult.meas_type == 'e':
            hsmresult_corrected_g1 = hsmresult.corrected_e1/2.
            hsmresult_corrected_g2 = hsmresult.corrected_e2/2.
        else:
            logger.error('meas type: %s' % hsmresult.meas_type)


        line = fmt_str % (ig,
                        hsmresult_corrected_g1,
                        hsmresult_corrected_g2,
                        hsmresult.moments_sigma,
                        hsmresult.moments_centroid.x,
                        hsmresult.moments_centroid.y)

    
    file_result.write(line)        

def write_empty_result(file_result,ig):

    fmt_str = '%d\t' + '%f\t'*5 + '\n'    
    vals = [ig] + [ERR_FLAG]*5
    line = fmt_str % tuple(vals)
    file_result.write(line)

def run():

    cat = numpy.loadtxt(args.filename_input,dtype=dtype_table_cat)
    logger.info('opened %s' % args.filename_input)

    method = 'REGAUSS'
    hsmparams = galsim.hsm.HSMParams()
    hsmparams.max_mom2_iter=1000 
 
    for filename_meds in cat['filename_meds']:

        filename_result = '%s.hsm.cat' % (filename_meds)
        file_result = open(filename_result,'w')
        write_result(header_only=True)

        m = meds.MEDS(filename_meds)
        logger.info('opened meds file %s with %d galaxies',filename_meds,m.size)

        filename_psf = '%s.psf' % filename_meds
        psf = pyfits.getdata(filename_psf)
        img_psf = image_array_to_galsim(psf)
        
        for ig in range(m.size):

            gal = m.get_cutout(ig,0) # first cutout

            img_gal = image_array_to_galsim(gal)

            try:
                hsmresult = galsim.hsm.EstimateShear(img_gal,img_psf,shear_est=method,strict=False,hsmparams=hsmparams)

                if hsmresult.error_message != "":                                   
                    logger.debug(('%3d hsm failed with msg: %s' % (ig,hsmresult.error_message)).strip())
                    write_empty_result(file_result,ig)
                else:
                    write_result(file_result,hsmresult,ig)

            except Exception, errmsg:
                logger.debug(('%3d hsm failed with msg: %s' % (ig,errmsg)).strip())
                write_empty_result(file_result,ig)

            if ig % 100 == 0 :
                logger.info('passing %5d' % ig)


        file_result.close()
        logger.info('closed %s' % filename_result)
   
    

def main():

    import numpy, galsim, sys, logging, yaml, argparse, time

    global logger , config , args

    description = 'Run on sha1 data with HSM'
    parser = argparse.ArgumentParser(description=description, add_help=True)
    parser.add_argument('-v', '--verbosity', type=int, action='store', default=2, choices=(0, 1, 2, 3 ), help='integer verbosity level: min=0, max=3 [default=2]')
    parser.add_argument('-i', '--filename_input', default='sha1-O1.cat',type=str, action='store', help='name of the output catalog')
    
    args = parser.parse_args()
    # Parse the integer verbosity level from the command line args into a logging_level string
    logging_levels = { 0: logging.CRITICAL, 
                       1: logging.WARNING,
                       2: logging.INFO,
                       3: logging.DEBUG }
    logging_level = logging_levels[args.verbosity]
    logging.basicConfig(format="%(message)s", level=logging_level, stream=sys.stdout)
    logger = logging.getLogger("run_hsm_sha1") 
    logger.setLevel(logging_level)

    logger.info(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))

    run()


    logger.info(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))


main()