import os
import pdb
import pyfits
import logging
import sys
import numpy
import sys
import math
import pylab
import argparse
import yaml
import galsim
import copy
import datetime  

dtype_results   = { 'names'  : ['g1','g2'],
                    'formats': ['f4']*2 }

global logger , config , args


figure_size_x = 8
figure_size_y = 6

lines_g1 = ['bx-','rx-','cx-']
lines_g2 = ['b+-','rx-','cx-']

def savePlots():

	n_angles = config['gal']['num']

	g1_true = config['eval_variables']['fg1']
	g2_true = config['eval_variables']['fg2']

	ellip = [(config['eval_variables']['ffirst_ellip'] + i*config['eval_variables']['fellip_step']) for i in range(config['eval_variables']['inellip'])]
	hlr   =	[(config['eval_variables']['ffirst_hlr']   + i*config['eval_variables']['fhlr_step'])   for i in range(config['eval_variables']['inhlr'])]
	profile = range(1,len(config['gal']['first']['items'])+1)  
	method = range(len(args.methods))

	n_ellip = len(ellip)
	n_hlr = len(hlr)
	n_profile = len(profile)
	n_method = len(args.methods)

	results_g1_mean = numpy.empty((n_method,n_profile,n_ellip,n_hlr)) 
	results_g2_mean = numpy.empty((n_method,n_profile,n_ellip,n_hlr)) 
	results_g1_stdm = numpy.empty((n_method,n_profile,n_ellip,n_hlr)) 
	results_g2_stdm = numpy.empty((n_method,n_profile,n_ellip,n_hlr)) 

	for im,m in enumerate(args.methods):
		for ip,p in enumerate(profile):
			for ie,e in enumerate(ellip):
				for ih,h in enumerate(hlr):
					filename_results = 'profile%02dellip%02dhlr%02d.meds.%s.cat' % (p, e*10, h*10, m)
					filepath_results = os.path.join(args.dirpath_results,filename_results)
					# print filename_results
					results = numpy.loadtxt(filepath_results,dtype=dtype_results)
					n_gals = results.shape[0]

					# the errorbars are taken from the 0 intrinsic ellipticity case, which doesn't contain shape noise
					# and then these errorbars are used for all ellipticities
					# this is a very crude workaround for the shape noise, and can easilly be misleading.
					if ie == 0:			
						results_g1_stdm[im,ip,ie,ih] = numpy.std(results['g1'],ddof=1) / numpy.sqrt(n_gals)
						results_g2_stdm[im,ip,ie,ih] = numpy.std(results['g2'],ddof=1) / numpy.sqrt(n_gals)
					else:
						results_g1_stdm[im,ip,ie,ih] = results_g1_stdm[im,ip,0,ih]
						results_g2_stdm[im,ip,ie,ih] = results_g2_stdm[im,ip,0,ih]

					results_g1_mean[im,ip,ie,ih] = numpy.mean(results['g1'])
					results_g2_mean[im,ip,ie,ih] = numpy.mean(results['g2'])

					logger.info('ngals=%d, method=%s, profile=%d, ellip=%2.2f, hlr=%2.2f, g1=%.4f +/- %2.4e, g2 = %.4f +/- %2.4e' % (
						n_gals , m , p , e , h , results_g1_mean[im,ip,ie,ih] , results_g1_stdm[im,ip,ie,ih] , results_g2_mean[im,ip,ie,ih] , results_g2_stdm[im,ip,ie,ih]))

	
	# there are 6x3 plots, one for each profile and ellipticity
	# on the x axis we vary the hlr

	for ie,e in enumerate(ellip):
		for ip,p in enumerate(profile):
			
			pylab.figure()
			pylab.clf

			for im in method:
				y1 = (results_g1_mean[im,ip,ie,:] - g1_true)  / g1_true
				y2 = (results_g2_mean[im,ip,ie,:] - g2_true)  / g2_true
				s1 = results_g1_stdm[im,ip,ie,:] / g1_true
				s2 = results_g2_stdm[im,ip,ie,:] / g2_true
				pylab.errorbar(hlr,y1,yerr=s1,fmt=lines_g1[im])
				pylab.errorbar(hlr,y2,yerr=s2,fmt=lines_g2[im])
			
			profile_name = config['gal']['first']['items'][ip]['type']
			pylab.xlabel('hlr')
			pylab.ylabel('dgi/gi')
			pylab.title('profile %s , ellipticity %2.2f' % (profile_name,e))
			legend = [val for val in args.methods for _ in (0, 1)]
			pylab.legend(legend)
			filename_fig = 'fig.profile%02dellip%02d.png' % (ip,e*10)
			pylab.savefig(filename_fig)
			logger.info('saved figure %s' % filename_fig)











if __name__ == "__main__":

    
    description = 'Compare reconvolved and directly created galaxies.'

    # parse arguments
    parser = argparse.ArgumentParser(description=description, add_help=True)
    parser.add_argument('--methods',  type=str, nargs='+',  action='store', help= 'use those methods to create plots')  
    parser.add_argument('--filepath_config',  type=str, default='shear_test1.meds.yaml',  action='store', help= 'use those methods to create plots')  
    parser.add_argument('--dirpath_results',  type=str, default='./results/',  action='store', help= 'directory containing results')     
    parser.add_argument('-v', '--verbosity', type=int, action='store', default=2, choices=(0, 1, 2, 3 ), help='integer verbosity level: min=0, max=3 [default=2]')
    args = parser.parse_args()
    # Parse the integer verbosity level from the command line args into a logging_level string
    logging_levels = { 0: logging.CRITICAL, 
                       1: logging.WARNING,
                       2: logging.INFO,
                       3: logging.DEBUG }
    logging_level = logging_levels[args.verbosity]
    logging.basicConfig(format="%(message)s", level=logging_level, stream=sys.stdout)
    logger = logging.getLogger("nmb_main") 

    # load the configuration file
    config = yaml.load(open(args.filepath_config,'r')) 

    savePlots()