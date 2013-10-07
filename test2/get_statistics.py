import numpy
import glob

#index filename_meds hlr snr ig g1 g2

dtype_stats =  { 'names'   : ['index', 'name_meds' , 'hlr' ,'snr', 'ig','g1','g2'],
                        'formats' : ['i8']*1 + ['a100'] + ['f4']*2 + ['i8']*1 + ['f4']*2 } 

dtype_im3shape =  { 'names'   : [
                                'id',
                                'ra',
                                'dec',
                                'e1',
                                'e2',
                                'rgp_rp',
                                'stats_signal_to_noise',
                                'likelihood',
                                'time_taken',
                                'ra_pix',
                                'dec_pix',
                                'sersic_parameter_radius',
                                'sersic_bulge_ampl',
                                'sersic_disc_ampl',
                                'sersic_bulge_flux',
                                'sersic_disc_flux',
                                'sersic_flux_ratio',
                                'stats_min_residuals',
                                'stats_max_residuals',
                                'stats_model_min',
                                'stats_model_max',
                                'levmar_number_of_likelihood_evaluations',
                                'levmar_number_of_iterations',
                                'levmar_reason_of_termination'
                                ],
                            'formats': 
                                ['i8']*1 + ['f4']*20 + ['i8']*3

                        }


dtype_im3shape_new =  { 'names'   : [
                                        'identifier',
                                        'ra',
                                        'dec',
                                        'ra_pix',
                                        'dec_pix',
                                        'e1',
                                        'e2',
                                        'radius',
                                        'radius_ratio',
                                        'bulge_A',
                                        'disc_A',
                                        'bulge_index',
                                        'disc_index',
                                        'delta_e_bulge',
                                        'delta_theta_bulge',
                                        'time',
                                        'bulge_flux',
                                        'disc_flux',
                                        'flux_ratio',
                                        'snr',
                                        'min_residuals',
                                        'max_residuals',
                                        'model_min',
                                        'model_max',
                                        'likelihood',
                                
                                                                ],
                            'formats': 
                                ['i8']*1 + ['f4']*24

                        }



def _getLineFit(x,y,sig):
        """
        @brief get linear least squares fit with uncertainity estimates
        y(X) = b*X + a
        see numerical recipies 15.2.9
        @param X    function arguments 
        @param y    function values
        @param sig  function values standard deviations   
        @return a - additive 
        @return b - multiplicative
        @return C - covariance matrix
        """
        
        import numpy

        invsig2 = sig**-2;
        
        S  = numpy.sum(invsig2)
        Sx = numpy.inner(x,invsig2)
        Sy = numpy.inner(y,invsig2)
        Sxx = numpy.inner(invsig2*x,x)
        Sxy = numpy.inner(invsig2*x,y)

        D = S*Sxx - Sx**2
        a = (Sxx*Sy - Sx*Sxy)/D
        b = (S*Sxy  - Sx*Sy)/D
        
        Cab = numpy.zeros((2,2))
        Cab[0,0] = Sxx/D
        Cab[1,1] = S/D
        Cab[1,0] = -Sx/D
        Cab[0,1] = Cab[1,0]
        
        return a,b,Cab

def select_usable_shears(e1,e2):

    select = numpy.abs(e1 + 1j*e2) < 1.
    e1s = e1[select]
    e2s = e2[select]
    n_gals_use = len(e1s)
    return e1s, e2s,  n_gals_use


def get_stats():

    if args.results_format == 'im3shape':
        dtype_results = dtype_im3shape
    elif args.results_format == 'im3shape_new':
        dtype_results = dtype_im3shape_new
    else:
        logger.error('so far only im3shape results are used')
     
    truth_cat = numpy.loadtxt(args.filename_input,dtype=dtype_stats)
    n_truth_cat = len(truth_cat)
    logger.info('opened truth catalog file %s with %d lines' % (args.filename_input,n_truth_cat))

    for itr,tr in enumerate(truth_cat):

        glob_search = '%s/%s.*.cat' % (args.dirname_results,tr['name_meds'])
        file_list = glob.glob(glob_search)
        n_file_list = len(file_list)

        logger.debug('test case %s has %d files' % (glob_search,n_file_list))

        e1_est = []
        e2_est = []

        if len(file_list) == 0:
            logger.error('no results for %s' % glob_search)
            continue

        for fl in file_list:
            results = numpy.loadtxt(fl,dtype=dtype_results,usecols=range(25))
            # results = numpy.loadtxt(fl,usecols=range(10))
            n_results = len(results)
            logger.debug('opened file %s with %d galaxies' % (fl,n_results))

            e1_est.extend(results['e1'])
            e2_est.extend(results['e2'])

        n_gals_total = len(e1_est)
        e1_est = numpy.array(e1_est)
        e2_est = numpy.array(e2_est)

        e1_est,e2_est,n_gals_use = select_usable_shears(e1_est,e2_est)
        # n_gals_use  =n_gals_total
        
        mean_g1 = numpy.mean(e1_est)
        mean_g2 = numpy.mean(e2_est)
        stdv_g1 = numpy.std(e1_est,ddof=1) 
        stdv_g2 = numpy.std(e2_est,ddof=1) 
        stdm_g1 = stdv_g1 / numpy.sqrt(n_gals_use)
        stdm_g2 = stdv_g2 / numpy.sqrt(n_gals_use)

        logger.info('%10d \tsnr=% 2.4f\thlr=% 2.4f\tg1=% 2.4f\tg2=% 2.4f\tg1e=% 2.4f\tg1e=% 2.4f\tstdm_g1e=% 2.4f\tstdm_g1e=% 2.4f\tstdv_g1e=% 2.4f\tstdv_g1e=% 2.4f \t %5d / %5d' % (itr,tr['snr'],tr['hlr'],tr['g1'],tr['g2'],mean_g1,mean_g2,stdm_g1,stdm_g2,stdv_g1,stdv_g2,n_gals_use,n_gals_total))

        if logger.level == logging.DEBUG:

            import pylab
            n_bins = 50
            title_line = 'g_true = (%.4f,%.4f) g_est = (%.4f,%.4f)' % (tr['g1'],tr['g2'],mean_g1,mean_g2)
            pylab.subplot(1,2,1)
            pylab.hist(e1_est,n_bins)
            pylab.subplot(1,2,2)
            pylab.hist(e2_est,n_bins)
            pylab.suptitle(title_line)

            filename_fig = 'figures/fig-hlr%2.2f-snr%02d-shear%02d.png' % (tr['hlr'],tr['snr'],tr['ig'])
            pylab.savefig(filename_fig)
            logger.debug('saved %s' % filename_fig)
            pylab.close()


        # logger.info('m1 = % 2.4f \t +/- % 2.4f' % ( m1, std_m1))
        # logger.info('m2 = % 2.4f \t +/- % 2.4f' % ( m2, std_m2))
        # logger.info('c1 = % 2.4f \t +/- % 2.4f' % ( c1, std_c1))
        # logger.info('c2 = % 2.4f \t +/- % 2.4f' % ( c2, std_c2))


        




def main():


    global logger , config , args

    description = 'produce m and c estimates from shear_test2 results'
    parser = argparse.ArgumentParser(description=description, add_help=True)
    parser.add_argument('-v', '--verbosity', type=int, action='store', default=2, choices=(0, 1, 2, 3 ), help='integer verbosity level: min=0, max=3 [default=2]')
    parser.add_argument('-o', '--filename_input', default='test2.cat',type=str, action='store', help='name of the input catalog')
    parser.add_argument('-r', '--dirname_results', default='results',type=str, action='store', help='name of the dir with all the results files')
    parser.add_argument('-f', '--results_format', default='im3shape',type=str, action='store', help='which format to use, so far implemented only [im3shape] ')
    # parser.add_argument('-c', '--filename_config', default='test2.yaml',type=str, action='store', help='name of the yaml config file')
    # parser.add_argument('-d', '--dry', default=False,  action='store_true', help='Dry run, dont generate data')

    args = parser.parse_args()
    # Parse the integer verbosity level from the command line args into a logging_level string
    logging_levels = { 0: logging.CRITICAL, 
                       1: logging.WARNING,
                       2: logging.INFO,
                       3: logging.DEBUG }
    logging_level = logging_levels[args.verbosity]
    logging.basicConfig(format="%(message)s", level=logging_level, stream=sys.stdout)
    logger = logging.getLogger("ell_noise") 
    logger.setLevel(logging_level)

    get_stats()

import yaml, argparse, sys, logging 
main()