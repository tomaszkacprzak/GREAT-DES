import numpy, galsim, sys, logging, yaml, argparse, time, copy, itertools, tabletools, nbc1_plots
import numpy as np; import pylab as pl
from nbc1_dtypes import * ; 



def get_results_sample(filelist_svclusters,filename_all_svclusters_results):

    files_svclusters = [line.strip() for line in open(filelist_svclusters)]
    
    total_gals= 0
    list_results = []

    for f in files_svclusters:

        results_svclusters = tabletools.loadTable(f, dtype=dtype_table_results_sv)
        total_gals += len(results_svclusters)
        log.info('opened %s with %d gals, total number %d' % (f, len(results_svclusters),total_gals) )
        list_results.append(results_svclusters)

    all_results = np.concatenate(list_results)
    print len(all_results)
    tabletools.saveTable(filename_all_svclusters_results,all_results)


def create_histogram(ell,size,snr,filename_hist='hist.test.png',show=False,close=True):


    bins_ell = np.linspace(0.,1.,10)
    bins_size = np.linspace(1.2,2.2,10)
    bins_snr = [2.5,7.5,12.5,17.5,22.5,40]
    # bins_snr = np.linspace(0,100,50)
    samples=np.concatenate([ell[:,None],size[:,None],snr[:,None]],axis=1)

    h3d, edges = pl.histogramdd(samples, bins = [bins_ell,bins_size,bins_snr],normed=True)
    
    pl.figure(1)
    pl.hist(ell,bins=bins_ell,histtype='step',normed=True)
    pl.title('ell')
    pl.legend(filename_hist)
    filename_hist_current = filename_hist.replace('.png','.ell.png')
    pl.savefig(filename_hist_current)
    log.info('saved %s' % filename_hist_current)

    pl.figure(2)
    pl.hist(size,bins=bins_size,histtype='step',normed=True)
    pl.title('size')
    filename_hist_current = filename_hist.replace('.png','.size.png')
    pl.savefig(filename_hist_current)
    pl.legend(filename_hist)
    log.info('saved %s' % filename_hist_current)

    pl.figure(3)
    pl.hist(snr,bins=bins_snr,histtype='step',normed=True)
    pl.title('snr')
    filename_hist_current = filename_hist.replace('.png','.snr.png')
    pl.savefig(filename_hist_current)
    log.info('saved %s' % filename_hist_current)

    h,b,_ = pl.hist(snr,bins=bins_snr,histtype='step',normed=True)
    # np.savetxt('hist_snr.txt',h)
    # np.savetxt('bins_snr.txt',b)
    # print 'saved hists'


    if show:
        pl.show()
    if close:
        pl.close()

    return h3d, edges



def match_properties(filename_results_calib,filename_all_svclusters_results):

    results_calib = tabletools.loadTable(filename_results_calib)
    results_svclusters = tabletools.loadTable(filename_all_svclusters_results)

    filename_hist_cal = 'hist.cal.png'
    filename_hist_sv = 'hist.sv.png'
    ell_calib = np.abs(results_calib['g1'] + 1j*results_calib['g2'])
    ell_sv = np.abs(results_svclusters['g1'] + 1j*results_svclusters['g2'])
    hist_cal, edges = create_histogram(ell=ell_calib,size=results_calib['size'],snr=results_calib['snr'],filename_hist=filename_hist_cal)
    hist_sv, edges = create_histogram(ell=ell_sv,size=results_svclusters['size'],snr=results_svclusters['snr'],filename_hist=filename_hist_sv)

    bins_ell , bins_size , bins_snr = edges

    # save plots
    for e_slice in range(hist_cal.shape[0]):
        
        pl.suptitle('%2.2f < |e| < %2.2f' % (edges[0][e_slice],edges[0][e_slice+1]))
        
        vmax=max([ max(hist_cal[e_slice,:,:].flatten()) , max(hist_sv[e_slice,:,:].flatten()) ] )

        pl.subplot(1,2,1)
        pl.imshow(hist_cal[e_slice,:,:],interpolation='nearest',origin='lower',extent=(min(bins_snr),max(bins_snr),min(bins_size),max(bins_size)),aspect='auto',vmin=0,vmax=vmax)
        pl.colorbar()
        pl.title('calibration data')

        pl.subplot(1,2,2)
        pl.imshow(hist_sv[e_slice,:,:],interpolation='nearest',origin='lower',extent=(min(bins_snr),max(bins_snr),min(bins_size),max(bins_size)),aspect='auto',vmin=0,vmax=vmax)
        pl.colorbar()
        pl.title('sv coadd v4')

        filename_fig = 'figs/hist2d.counts.ell%d.png' % e_slice
        pl.savefig(filename_fig)
        log.info('saved %s' % filename_fig)
        # pl.show()
        pl.close()

    cal_ell = np.abs(results_calib['g1'] + 1j*results_calib['g2'])
    cal_size = results_calib['size']
    cal_snr = results_calib['snr']

    hist_m = np.ones(hist_cal.shape)
    hist_m_std = np.ones(hist_cal.shape)
    hist_c = np.ones(hist_cal.shape)
    hist_c_std = np.ones(hist_cal.shape)
    hist_ngal = np.ones(hist_cal.shape)

    # now calculate m_j
    iall=0
    n_bins_total = len(bins_ell[:-1])*len(bins_size[:-1])*len(bins_snr[:-1])
    log.info('getting m and for %d bins' % n_bins_total)

    for isize,vsize in enumerate(bins_size[:-1]):
        for isnr,vsnr in enumerate(bins_snr[:-1]):

            # get the standard deviation of ellipcitiy estimator before binning on ellipticity
            select = (cal_size > bins_size[isize]) * (cal_size < bins_size[isize+1]) *(cal_snr > bins_snr[isnr]) * (cal_snr < bins_snr[isnr+1])
            current_sampels = results_calib[select]
            g_err = np.std(current_sampels['g1'],ddof=1)

            g1_true = current_sampels['g1_true']
            g1_est  = current_sampels['g1']
            g2_true = current_sampels['g2_true']
            g2_est  = current_sampels['g2']
            g_std = np.ones(current_sampels['g1'].shape)*g_err
            
            [c1,m1,C1cm] = nbc1_plots.get_line_fit(g1_true,g1_est,g_std)
            [c2,m2,C2cm] = nbc1_plots.get_line_fit(g2_true,g2_est,g_std)
            m1_std = numpy.sqrt(C1cm[1,1])
            m2_std = numpy.sqrt(C2cm[1,1])
            c1_std = numpy.sqrt(C1cm[0,0])
            c2_std = numpy.sqrt(C2cm[0,0])
            m = np.mean([m1,m2])
            m_std = m1_std/np.sqrt(2)
            c = np.mean([c1,c2])
            c_std = np.sqrt((c1_std**2 + c2_std**2)/4.)

            alle_m = m
            alle_m_std = m_std

            for iell,vell in enumerate(bins_ell[:-1]):

                select = (cal_ell > bins_ell[iell]) * (cal_ell < bins_ell[iell+1]) * (cal_size > bins_size[isize]) * (cal_size < bins_size[isize+1]) *(cal_snr > bins_snr[isnr]) * (cal_snr < bins_snr[isnr+1])
                current_sampels = results_calib[select]
                n_gal = len(current_sampels)
                if n_gal<10:
                    log.info('%4d ell=[%2.2f,%2.2f] size=[%2.2f,%2.2f] snr=[%2.2f,%2.2f] ngals=%d' % (
                                    iall,
                                    bins_ell[iell],bins_ell[iell+1],
                                    bins_size[isize],bins_size[isize+1],
                                    bins_snr[isnr],bins_snr[isnr+1], 
                                    n_gal
                                    ))
                    continue

                g1_true = current_sampels['g1_true']
                g1_est  = current_sampels['g1']
                g2_true = current_sampels['g2_true']
                g2_est  = current_sampels['g2']
                g_std = np.ones(current_sampels['g1'].shape)*g_err
                
                [c1,m1,C1cm] = nbc1_plots.get_line_fit(g1_true,g1_est,g_std)
                [c2,m2,C2cm] = nbc1_plots.get_line_fit(g2_true,g2_est,g_std)
                m1_std = numpy.sqrt(C1cm[1,1])
                m2_std = numpy.sqrt(C2cm[1,1])
                c1_std = numpy.sqrt(C1cm[0,0])
                c2_std = numpy.sqrt(C2cm[0,0])
                m = np.mean([m1,m2])
                m_std = m1_std/np.sqrt(2)
                c = np.mean([c1,c2])
                c_std = np.sqrt((c1_std**2 + c2_std**2)/4.)
                
                pl.figure()
                pl.errorbar(g1_true,g1_est,yerr=g_std,fmt='b.',label='m1=%2.3f +/- %2.3f' % (m1,m1_std))
                pl.errorbar(g2_true,g2_est,yerr=g_std,fmt='r.',label='m1=%2.3f +/- %2.3f' % (m2,m2_std))
                pl.plot(g1_true,g1_true*m + c,'b-',label='mean m=%2.3f +/- %2.3f' % (m,m_std))
                pl.title('snr=%2.2f size=%2.2f ell=%2.2f' % (vsnr,vsize,vell))
                pl.legend()
                filename_fig = 'figs/fig.bias.ell%02d.size%02d.snr%02d.png' % (iell,isize,isnr)
                pl.savefig(filename_fig)
                pl.close()
                log.debug('saved %s' % filename_fig)

                hist_m[iell,isize,isnr]     = m
                hist_m_std[iell,isize,isnr] = m_std
                hist_c[iell,isize,isnr]     = c
                hist_c_std[iell,isize,isnr] = c_std
                hist_ngal[iell,isize,isnr] = n_gal

                # log.info('size=%2.2f\tsnr=%2.2f\tm=% 2.2e\tc=% 2.2f\tm_err=%2.2e\tc_err=%2.2e\tg_std=%2.2f' % (vsize, vsnr,m1,c1,m1_std,c1_std,g_std))                

                log.info('%4d\tell=[%2.2f,%2.2f]\tsize=[%2.2f,%2.2f]\tsnr=[%2.2f,%2.2f]\tngals=%d\tm=% 2.4f (%2.4f)\tc=% 2.4f (%2.4f)\tg_std=%2.3f' % (
                                    iall,
                                    bins_ell[iell],bins_ell[iell+1],
                                    bins_size[isize],bins_size[isize+1],
                                    bins_snr[isnr],bins_snr[isnr+1], 
                                    n_gal,
                                    m,
                                    m_std,
                                    c,
                                    c_std,
                                    g_err
                                    ))


                iall+=1

            mean_ebin = sum( hist_m[:,isize,isnr] / hist_m_std[:,isize,isnr]**2 ) / sum(1./hist_m_std[:,isize,isnr]**2)
            # std_ebin = np.std(hist_m[:,isize,isnr],ddof=1)))
            log.info('m calculated for that size and snr bin=%2.4f (%2.4f)' % (alle_m,alle_m_std))
            log.info('mean from bins of e for that size and snr bin=%2.4f' % mean_ebin)

    tabletools.saveTable('hist_m.fits', hist_m,logger=log)
    tabletools.saveTable('hist_m_std.fits', hist_m_std,logger=log)
    tabletools.saveTable('hist_c.fits', hist_m,logger=log)
    tabletools.saveTable('hist_c_std.fits', hist_c_std,logger=log)
    tabletools.saveTable('hist_ngal.fits', hist_ngal,logger=log)

    for e_slice in range(hist_cal.shape[0]):
        
        pl.suptitle('%2.2f < |e| < %2.2f' % (edges[0][e_slice],edges[0][e_slice+1]))
        
        vmax_m= max( hist_m[e_slice,:,:].flatten() )
        vmax_m_std= max( hist_m_std[e_slice,:,:].flatten() )
        vmax_c= max( hist_c[e_slice,:,:].flatten() )
        vmax_c_std= max( hist_c_std[e_slice,:,:].flatten() )

        pl.subplot(2,2,1)
        pl.imshow(hist_m[e_slice,:,:],interpolation='nearest',origin='lower',extent=(min(bins_snr),max(bins_snr),min(bins_size),max(bins_size)),aspect='auto',vmin=0,vmax=vmax_m)
        pl.colorbar()
        pl.title('m calibration data')

        pl.subplot(2,2,2)
        pl.imshow(hist_m_std[e_slice,:,:],interpolation='nearest',origin='lower',extent=(min(bins_snr),max(bins_snr),min(bins_size),max(bins_size)),aspect='auto',vmin=0,vmax=vmax_m_std)
        pl.colorbar()
        pl.title('m_std sv coadd v4')

        pl.subplot(2,2,3)
        pl.imshow(hist_c[e_slice,:,:],interpolation='nearest',origin='lower',extent=(min(bins_snr),max(bins_snr),min(bins_size),max(bins_size)),aspect='auto',vmin=0,vmax=vmax_c_std)
        pl.colorbar()
        pl.title('c calibration data')

        pl.subplot(2,2,4)
        pl.imshow(hist_c_std[e_slice,:,:],interpolation='nearest',origin='lower',extent=(min(bins_snr),max(bins_snr),min(bins_size),max(bins_size)),aspect='auto',vmin=0,vmax=vmax_c_std)
        pl.colorbar()
        pl.title('c_std sv coadd v4')


        filename_fig = 'figs/hist2d.mc.ell%d.png' % e_slice
        pl.savefig(filename_fig)
        log.info('saved %s' % filename_fig)
        # pl.show()
        pl.close()

                # bias_result = {}
                # bias_result['index']  = iall
                # bias_result['iell'] = iell
                # bias_result['isize'] = isize
                # bias_result['isnr']  = isnr
                # bias_result['ell'] = vell
                # bias_result['size'] = vsize
                # bias_result['snr'] = vsnr
                # bias_result['m'] = m
                # bias_result['c'] = c
                # bias_result['m_std'] = m_std
                # bias_result['c_std'] = c_std
                # bias_result['g_std'] = g_std
                # # bias_result['n_fail'] = n_fail













def get_calibration_sample(filename_results_calib,dirpath_results):

    list_results = []
    for index in cat['index']:

        filename_results = '%s/nbc_%03d.fits.im3.cleaned.cat' % (dirpath_results,index)
        current_results = tabletools.loadTable(filename_results,dtype=dtype_table_results_im3shape_cleaned,logger=log)
        log.info('loaded %s' % filename_results)
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
    log_formatter = logging.Formatter("%(asctime)s %(name)s %(levelname)s:   %(message)s ","%Y-%m-%d %H:%M:%S")
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(log_formatter)
    log.addHandler(stream_handler)

    log.info(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))

    filelist_svclusters = 'filelist_svclusters.txt'
    filename_all_svclusters_results = 'sv_clusters.all.fits'
    filename_results_calib = 'calib.v4.2013.12.20.all.fits'
    dirpath_results = 'run-tiled-001/calib.v4.2013.12.20/cleaned_calib.v4.2013.12.20'

    cat = tabletools.loadTable(args.filename_input,dtype=dtype_table_cat)
    config = yaml.load(open(args.filename_config))

    # get_results_sample(filelist_svclusters,filename_all_svclusters_results)
    # get_calibration_sample(filename_results_calib,dirpath_results)
    match_properties(filename_results_calib,filename_all_svclusters_results)

    log.info(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))


main()
