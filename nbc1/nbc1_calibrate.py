import numpy, galsim, sys, logging, yaml, argparse, time, copy, itertools, tabletools; 
import numpy as np; import pylab as pl
from nbc1_dtypes import *



def get_results_sample(filelist_svclusters,filename_all_svclusters_results):

    files_svclusters = [line.strip() for line in open(filelist_svclusters)]
    
    total_gals= 0
    list_results = []

    for f in files_svclusters:

        results_svclusters = np.loadtxt(f, dtype=dtype_table_results_sv)
        total_gals += len(results_svclusters)
        print total_gals
        list_results.append(results_svclusters)

    all_results = np.concatenate(list_results)
    print len(all_results)
    tabletools.saveTable(filename_all_svclusters_results,all_results)


def create_histogram(ell,size,snr,filename_hist='hist.test.png',show=False):

    bins_ell = np.linspace(-1,1,30)
    bins_size = np.linspace(1.2,2.2,30)
    # bins_snr = [5.,7.5,10,12.5,15,17.5,20,25,30,40]
    bins_snr = np.linspace(0,100,50)
    samples=np.concatenate([ell[:,None],size[:,None],snr[:,None]],axis=1)

    h3d, edges = pl.histogramdd(samples, bins = [bins_ell,bins_size,bins_snr])
    print edges
    print h3d.shape

    # for e_slice in range(h3d.shape[0]):
    #     pl.imshow(h3d[e_slice,:,:],interpolation='nearest')
    #     pl.colorbar()
    #     pl.show()

    pl.figure(1)
    pl.hist(ell,bins=bins_ell,histtype='step',normed=True)
    pl.title('ell')
    pl.legend(filename_hist)
    filename_hist_current = filename_hist.replace('.png','.ell.png')
    pl.savefig(filename_hist_current)
    logger.info('saved %s' % filename_hist_current)

    pl.figure(2)
    pl.hist(size,bins=bins_size,histtype='step',normed=True)
    pl.title('size')
    filename_hist_current = filename_hist.replace('.png','.size.png')
    pl.savefig(filename_hist_current)
    pl.legend(filename_hist)
    logger.info('saved %s' % filename_hist_current)

    pl.figure(3)
    pl.hist(snr,bins=bins_snr,histtype='step',normed=True)
    pl.title('snr')
    filename_hist_current = filename_hist.replace('.png','.snr.png')
    pl.savefig(filename_hist_current)
    logger.info('saved %s' % filename_hist_current)

    if show:
        pl.show()

    return h3d, edges



def match_properties(filename_results_calib,filename_all_svclusters_results):

    results_calib = tabletools.loadTable(filename_results_calib)
    results_svclusters = tabletools.loadTable(filename_all_svclusters_results)

    filename_hist_cal = 'hist.cal.png'
    filename_hist_sv = 'hist.sv.png'
    create_histogram(ell=results_calib['g1'],size=results_calib['size'],snr=results_calib['snr'],filename_hist=filename_hist_cal)
    create_histogram(ell=results_svclusters['g1'],size=results_svclusters['size'],snr=results_svclusters['snr'],filename_hist=filename_hist_sv)

    pl.show()



def get_calibration_sample(filename_results_calib):

    list_results = []
    for index in cat['index']:

        filename_results = 'calib.v4.2013.12.20/cleaned_calib.v4.2013.12.20/nbc_%03d.fits.im3.cleaned.cat' % index
        current_results = tabletools.loadTable(filename_results,dtype=dtype_table_results_im3shape_cleaned,logger=logger)
        logger.info('loaded %s' % filename_results)
        n_gals = len(current_results)
        g1_true = np.ones(n_gals) * cat[index]['g1']
        g2_true = np.ones(n_gals) * cat[index]['g2']
        current_results = tabletools.appendColumn(rec=current_results,arr=g1_true,name='g1_true',dtype= 'f8')
        current_results = tabletools.appendColumn(rec=current_results,arr=g2_true,name='g2_true',dtype= 'f8')
        list_results.append(current_results)

    all_results = np.concatenate(list_results)
    tabletools.saveTable(filename_results_calib,all_results)


def main():

    import numpy, galsim, sys, logging, yaml, argparse, time

    global logger , config , args , cat

    description = 'Run on sha1 data with HSM'
    parser = argparse.ArgumentParser(description=description, add_help=True)
    parser.add_argument('-v', '--verbosity', type=int, action='store', default=2, choices=(0, 1, 2, 3 ), help='integer verbosity level: min=0, max=3 [default=2]')
    parser.add_argument('-i', '--filename_input', default='sha1-O1.cat',type=str, action='store', help='name of the output catalog')
    parser.add_argument('-c', '--filename_config', default='nbc1.tiled.yaml',type=str, action='store', help='name of the yaml config file')
    parser.add_argument('-m', '--method_id', default='hsm',type=str, action='store', help='name of the yaml config file')
    parser.add_argument('-t', '--tag', default='def',type=str, action='store', help='tag the stats (maybe for different cuts etc')
    
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

    filelist_svclusters = 'filelist_svclusters.txt'
    filename_all_svclusters_results = 'sv_clusters.all.fits'
    filename_results_calib = 'calib.v4.2013.12.20.all.fits'

    cat = tabletools.loadTable(args.filename_input,dtype=dtype_table_cat)
    config = yaml.load(open(args.filename_config))

    # get_results_sample(filelist_svclusters,filename_all_svclusters_results)
    # get_calibration_sample(filename_results_calib)
    match_properties(filename_results_calib,filename_all_svclusters_results)

    logger.info(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))


main()
