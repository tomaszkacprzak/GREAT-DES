import galsim, sys, logging, yaml, argparse, time, copy, itertools, tabletools, fitting, pickle
import numpy as np; import pylab as pl; 
from nbc1_dtypes import * ; 

bins_size = np.arange(1.2,2.2,0.1).tolist()
bins_snr = [2.5,7.5,12.5,17.5,22.5,40]

calib_struct = {}
calib_struct['bins_size'] = bins_size
calib_struct['bins_snr'] = bins_snr
calib_struct['bias_m'] =  np.zeros([ len(bins_size)-1 , len(bins_snr) - 1 ]).tolist()
calib_struct['bias_m_std'] = np.zeros([ len(bins_size)-1 , len(bins_snr) - 1 ]).tolist()

filelist_svclusters = 'filelist_svclusters.txt'
filename_results_calib = 'run-tiled-002/cleaned_calib.all.fits'
filename_calibration = 'nbc_bins.json'

"""
Calculate the uncertainty on multiplicative correction for selected results sample.
Results should have fields bin_id and nbc_m_std.
"""
def get_m_std_for_sample(results):

    # sum of variances for all bins
    var_all = 0.
    # number of results in this sample, which fall into bins
    n_all = 0
    # iterate over all bins
    for ibin in range(max(results['bin_id'])+1):
        # select galaxies falling in that bin
        results_bin = results[results['bin_id'] == ibin]
        # get number of galaxies which fell in that bin
        n_bin = len(results_bin)     
        # if no galaxis in current bin, continue
        if n_bin < 1: continue
        # calculate variance for that bin 
        var_bin = results_bin['nbc_m_std'][0]**2
        # add to the sum of variances, weighted by the number of galaxies in that bin
        var_all += var_bin*n_bin**2
        # update the total number of galaxies which fall into bins
        n_all += n_bin

    # square root the variance to get standard deviation, divide by n_all as we are taking a mean
    std_all = np.sqrt(var_all) / n_all
    return std_all

def analyse_population(results):

    calib_struct = pickle.load(open(filename_calibration))
    bins_size = calib_struct['bins_size']
    bins_snr = calib_struct['bins_snr']

    log.info('analysing population')
    for isize,vsize in enumerate(bins_size[:-1]):
        for isnr,vsnr in enumerate(bins_snr[:-1]):

            nbc_m=calib_struct['bias_m'][isize][isnr]
            bin_size_min = bins_size[isize]
            bin_size_max = bins_size[isize+1]
            bin_snr_min = bins_snr[isnr]
            bin_snr_max = bins_snr[isnr+1]
            select_size = (results['size'] < bin_size_max) * (results['size'] > bin_size_min)
            select_snr = (results['snr'] < bin_snr_max) * (results['snr'] > bin_snr_min)
            select = select_size * select_snr
            results_select=results[select]
            n_gals=len(results_select)
            log.info("size=[%5.2f,%5.2f]\tsnr=[%5.2f,%5.2f]\tn_gals=%d\tnbc_m=%5.2f" % (bin_size_min,bin_size_max,bin_snr_min,bin_snr_max,n_gals,nbc_m))

    g1_mean = np.mean(results['g1'])
    g1_bias_m = np.mean(results['nbc_m'])
    g1_bias_m_std = get_m_std_for_sample(results)
    g1_mean_calibrated = g1_mean / g1_bias_m

    log.info( 'g1_mean            %10.4f' % g1_mean  )
    log.info( 'g1_mean_calibrated %10.4f' % g1_mean_calibrated )
    log.info( 'g1_bias_m          %10.4f' % g1_bias_m )
    log.info( 'g1_bias_m_std      %10.4f' % g1_bias_m_std )


            # nbc_m = results_select['nbc_m'][0]
            # log.info("size=[%5.2f,%5.2f]\tsnr=[%5.2f,%5.2f]\tn_gals=%d\tnbc_m=%5.2f" % (bin_size_min,bin_size_max,bin_snr_min,bin_snr_max,n_gals,nbc_m))

def get_marg(key,results):

    calib_struct = pickle.load(open(filename_calibration))

    bins_key = 'bins_%s' % key 
    list_marg_m = []
    list_marg_m_std = []
    list_bin_centers = []
    bins = calib_struct[bins_key]
    for ib,vb in enumerate(bins[:-1]):

        bins_min = bins[ib]
        bins_max = bins[ib+1]

        select = (results[key] < bins_max) * (results[key] > bins_min)
        marg_m = np.mean(results['nbc_m'][select])
        marg_m_std = get_m_std_for_sample(results[select])
        list_marg_m.append(marg_m)
        list_marg_m_std.append(marg_m_std)
        list_bin_centers.append(bins_min+(bins_max-bins_min)/2.) 

    return list_marg_m, list_marg_m_std, list_bin_centers


def test_calibration_procedure():

    filelist = np.loadtxt(filelist_svclusters,dtype='a')
    calib_struct = pickle.load(open(filename_calibration))


    for ifile,filename_results in enumerate(filelist):

        filename_results_calibrated = filename_results.replace('.cat','.nbc.cat')
        results_sv = tabletools.loadTable(filename_results_calibrated,dtype=dtype_table_results_sv_calibrated,log=1)
        select = results_sv['flag'] == 0
        results_sv = results_sv[select]

        # analyse_population(results_sv)

        g1_mean = np.mean(results_sv['g1'])
        g1_bias_m = np.mean(results_sv['nbc_m'])
        g1_bias_m_std = get_m_std_for_sample(results_sv)
        g1_mean_calibrated = g1_mean / g1_bias_m

        log.info( 'g1_mean            %10.4f' % g1_mean  )
        log.info( 'g1_mean_calibrated %10.4f' % g1_mean_calibrated )
        log.info( 'g1_bias_m          %10.4f' % g1_bias_m )
        log.info( 'g1_bias_m_std      %10.4f' % g1_bias_m_std )

        marg_size_m, marg_size_m_std, marg_size_centers = get_marg('size',results_sv)
        marg_snr_m, marg_snr_m_std, marg_snr_centers = get_marg('snr',results_sv)

        pl.figure()
        pl.suptitle('%s\ntotal_calibration = %0.4f +/- %0.4f' % (filename_results,g1_bias_m,g1_bias_m_std))
        pl.subplot(1,2,1)
        hh,bh,_=pl.hist(results_sv['size'],bins=calib_struct['bins_size'],histtype='step')
        # hh,bh,_=pl.hist(results_sv['size'],bins=np.linspace(1.2,3,20),histtype='step')
        pl.xlabel('FWHM_RATIO')
        pl.ylabel('histogram')
        ax2=pl.gca().twinx()
        ax2.errorbar(marg_size_centers,marg_size_m,yerr=marg_size_m_std,fmt='md')
        ax2.set_ylabel('multiplicative bias m')

        pl.subplot(1,2,2)
        # pl.hist(results_sv['snr'],bins=calib_struct['bins_snr'],histtype='step')
        pl.hist(results_sv['snr'],bins=np.linspace(0,60,20),histtype='step')
        pl.xlabel('SNR')
        pl.ylabel('histogram')
        ax2=pl.gca().twinx()
        ax2.errorbar(marg_snr_centers,marg_snr_m,yerr=marg_snr_m_std,fmt='cd')
        ax2.set_ylabel('multiplicative bias m')

        left  = 0.125  # the left side of the subplots of the figure
        right = 0.9    # the right side of the subplots of the figure
        bottom = 0.1   # the bottom of the subplots of the figure
        top = 0.9      # the top of the subplots of the figure
        wspace = 0.7   # the amount of width reserved for blank space between subplots
        hspace = 0.5   # t
        pl.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

        filename_fig = ('figs/fig.marg.%s.png' % filename_results ).replace('.fits.im3.cleaned.wcs.cat','').replace('sv-clusters-shears/','')
        pl.savefig(filename_fig)
        pl.close()
        log.info('saved %s' % filename_fig)

        
        # pl.show()

        # n_gals = len(results_sv)

        # n_gals_select = int(n_gals*0.25)
        # select = np.random.permutation(n_gals)[0:n_gals_select]
        # n_selected = len(select)

        select = results_sv['snr'] < 15
        g1_mean = np.mean(results_sv[select]['g1'])
        g1_bias_m = np.mean(results_sv[select]['nbc_m'])
        g1_bias_m_std = get_m_std_for_sample(results_sv[select])
        g1_mean_calibrated = g1_mean / g1_bias_m

        print 'g1_mean' , g1_mean 
        print 'g1_mean_calibrated' , g1_mean_calibrated
        print 'g1_bias_m' , g1_bias_m
        print 'g1_bias_m_std' , g1_bias_m_std

        # import pdb; pdb.set_trace()


        # h,b,_=pl.hist(results_sv['nbc_m'][select],histtype='step')
        # pl.ylim([-10,max(h)])
        # pl.yscale('log')
        # pl.hist(results_sv['nbc_m'],histtype='step')
        # results_sv['nbc_m'].sort()
        # pl.plot(results_sv['nbc_m'])
        # pl.figure()
        # pl.scatter(results_sv['nbc_m'],results_sv['size'])
        # pl.plot()
        # pl.figure()
        # pl.scatter(results_sv['nbc_m'],results_sv['snr'])

        pl.show()




def get_calibration_bins():

    results_calib = tabletools.loadTable(filename_results_calib)  

    n_all = 0
    for isize,vsize in enumerate(bins_size[:-1]):
            for isnr,vsnr in enumerate(bins_snr[:-1]):

                bin_size_min = bins_size[isize]
                bin_size_max = bins_size[isize+1]
                bin_snr_min = bins_snr[isnr]
                bin_snr_max = bins_snr[isnr+1]
                select_size = (results_calib['size'] < bin_size_max) * (results_calib['size'] > bin_size_min)
                select_snr = (results_calib['snr'] < bin_snr_max) * (results_calib['snr'] > bin_snr_min)
                select = select_size * select_snr
                current_results = results_calib[select]
                # pl.subplot(1,3,1)
                # pl.hist(current_results['size'])
                # pl.subplot(1,3,2)
                # pl.hist(current_results['snr'])
                # pl.subplot(1,3,3)
                # pl.hist(current_results['isnr_true'])
                # pl.show()
                n_gals = len(current_results)
                n_all+=n_gals


                g1_tru = current_results['g1_true']
                g2_tru = current_results['g2_true']
                g1_est = current_results['g1']
                g2_est = current_results['g2']
                g1_err  = np.ones(n_gals)*np.std(current_results['g1'],ddof=1)
                g2_err  = np.ones(n_gals)*np.std(current_results['g2'],ddof=1)

                [c1,m1,C1cm] = fitting.get_line_fit(g1_tru,g1_est,g1_err)
                [c2,m2,C2cm] = fitting.get_line_fit(g2_tru,g2_est,g2_err)
                m1_std = np.sqrt(C1cm[1,1])
                m2_std = np.sqrt(C2cm[1,1])
                c1_std = np.sqrt(C1cm[0,0])
                c2_std = np.sqrt(C2cm[0,0])
                m_mean = (m1 + m2)/2.
                m_mean_std = np.sqrt((m1_std**2 + m2_std**2)/2.)

                calib_struct['bias_m'][isize][isnr] = m_mean
                calib_struct['bias_m_std'][isize][isnr] = m_mean_std

                del(select_size)
                del(select_snr)
                del(current_results)              

                log.info('size=%d [%5.2f,% 5.2f]\tsnr=%d [%5.2f,%5.2f]\tn_gals=%d\tm_mean=%2.4f\t(% 2.4f)' % (isize,bin_size_min,bin_size_max,isnr,bin_snr_min,bin_snr_max,n_gals,m_mean,m_mean_std))

    file_calibration = open(filename_calibration,'w')   
    pickle.dump(calib_struct,file_calibration)  
    file_calibration.close()  
    log.info('pickled %s' % filename_calibration)


def calibrate_results():

    filelist = np.loadtxt(filelist_svclusters,dtype='a')

    for ifile,vfile in enumerate(filelist):

        results_sv = tabletools.loadTable(vfile,dtype=dtype_table_results_sv,log=1)
        results_sv = add_calibration_columns(results_sv)

        filename_calibrated = vfile.replace('.cat','.nbc.cat')
        tabletools.saveTable(filename_calibrated,results_sv)
        log.info('calibrated file %3d %s' %(ifile,filename_calibrated)) 

def add_calibration_columns(results):

        calib_struct = pickle.load(open(filename_calibration))

        column_m = results['g1']*0 + 1.
        column_m_std = results['g1']*0 + 1e-9
        column_bin_id_other = results['g1']*0 - 1
        column_flag = 0
        results = tabletools.appendColumn(rec=results,arr=column_m,dtype='f8',name='nbc_m')
        results = tabletools.appendColumn(rec=results,arr=column_m_std,dtype='f8',name='nbc_m_std')
        results = tabletools.appendColumn(rec=results,arr=column_bin_id_other,dtype='i4',name='bin_id_size')
        results = tabletools.appendColumn(rec=results,arr=column_bin_id_other,dtype='i4',name='bin_id_snr')
        results = tabletools.appendColumn(rec=results,arr=column_bin_id_other,dtype='i4',name='bin_id')
        results = tabletools.appendColumn(rec=results,arr=column_flag,dtype='i4',name='flag')

        bins_size = calib_struct['bins_size']
        bins_snr = calib_struct['bins_snr']

        select = results['size'] < 1.2
        results['flag'][select] = 10

        iall=0
        for isize,vsize in enumerate(bins_size[:-1]):
            for isnr,vsnr in enumerate(bins_snr[:-1]):

                bin_size_min = bins_size[isize]
                bin_size_max = bins_size[isize+1]
                bin_snr_min = bins_snr[isnr]
                bin_snr_max = bins_snr[isnr+1]
                select_size = (results['size'] < bin_size_max) * (results['size'] > bin_size_min)
                select_snr = (results['snr'] < bin_snr_max) * (results['snr'] > bin_snr_min)
                select = select_size * select_snr

                nbc_m = calib_struct['bias_m'][isize][isnr]
                nbc_m_std = calib_struct['bias_m_std'][isize][isnr]
                results['nbc_m'][select] = nbc_m
                results['nbc_m_std'][select] = nbc_m_std
                results['bin_id_size'][select] = isize
                results['bin_id_snr'][select] = isnr
                results['bin_id'][select] = iall

                log.debug('size=%d [%5.2f,%5.2f] , snr=%d [%5.2f,%5.2f] n_gals=%d nbc_m=%2.4f (%2.4f)' % (isize,bin_size_min,bin_size_max,isnr,bin_snr_min,bin_snr_max,len(results[select]),nbc_m,nbc_m_std))
                iall+=1

        return results




    
    

def main():

    import numpy, galsim, sys, logging, yaml, argparse, time

    global log , config , args , cat

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
    
    log = logging.getLogger("nbc1_calibrate") 
    log.setLevel(logging_level)  
    log_formatter = logging.Formatter("%(asctime)s %(name)s %(levelname)s   %(message)s", "%Y-%m-%d %H:%M:%S")
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(log_formatter)
    log.addHandler(stream_handler)

    log.info(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))


    # get_calibration_bins()
    # calibrate_results()
    test_calibration_procedure()

    log.info(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))


main()
