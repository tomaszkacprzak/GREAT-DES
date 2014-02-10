import numpy as np
import pylab as pl
import galsim, sys, logging, yaml, argparse, time, copy, itertools, fitting
from nbc1_dtypes import *
sys.path.append('/home/tomek/code/tktools')
import tabletools
import plotstools

          

dtype_results = { 'hsm' : dtype_table_results , 'im3' : dtype_table_results_im3shape, 'im3_cleaned' : dtype_table_results_calib}

req1_dg = 0.003
req2_dg = 0.02
req3_dg = 0.1
req1_nfail = 0.01
req2_nfail = 0.05

DES_PIXEL_SCALE = 0.27

def get_line_fit(x,y,sig):
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
        
        invsig2 = sig**-2;
        
        S  = np.sum(invsig2)
        Sx = np.inner(x,invsig2)
        Sy = np.inner(y,invsig2)
        Sxx = np.inner(invsig2*x,x)
        Sxy = np.inner(invsig2*x,y)

        D = S*Sxx - Sx**2
        a = (Sxx*Sy - Sx*Sxy)/D
        b = (S*Sxy  - Sx*Sy)/D
        
        Cab = np.zeros((2,2))
        Cab[0,0] = Sxx/D
        Cab[1,1] = S/D
        Cab[1,0] = -Sx/D
        Cab[0,1] = Cab[1,0]
        
        return a,b,Cab


def surf_interp(arg1,arg2,value):

    grid_x, grid_y = np.mgrid[min(arg1):max(arg1):100j,min(arg2):max(arg2):100j]
    array_hlr  = np.array(arg1,ndmin=2).T
    array_snr  = np.array(arg2,ndmin=2).T
    points = np.concatenate((array_hlr,array_snr),axis=1)
    import scipy.interpolate
    grid = scipy.interpolate.griddata(points, value, (grid_x, grid_y), method='cubic')
    return grid.T


def plot_add_req(mult=1.):

    corner = pl.xlim()[0]
    length = abs(pl.xlim()[1]) + abs(pl.xlim()[0])
    pl.gca().add_patch(pl.Rectangle(  (corner,1-mult*req3_dg), length , 2*mult*req3_dg , facecolor = '0.9' , edgecolor='k' ))
    pl.gca().add_patch(pl.Rectangle(  (corner,1-mult*req2_dg), length , 2*mult*req2_dg , facecolor = '0.8' , edgecolor='k' ))
    pl.gca().add_patch(pl.Rectangle(  (corner,1-mult*req1_dg), length , 2*mult*req1_dg , facecolor = '0.7' , edgecolor='k' ))

def plot_add_nfail_req():

    corner = pl.xlim()[0]
    length = abs(pl.xlim()[1]) + abs(pl.xlim()[0])
    pl.gca().add_patch(pl.Rectangle(  (corner,-req2_nfail), length , 2*req2_nfail , facecolor = '0.9' , edgecolor='k' ))
    pl.gca().add_patch(pl.Rectangle(  (corner,-req1_nfail), length , 2*req1_nfail , facecolor = '0.7' , edgecolor='k' ))

def write_stats(file_stats,stats):

    stats_line_fmt = '%d\t'*3 + '% .8f\t'*9 + '\n'
    line = stats_line_fmt % (
            stats['index'],
            stats['n_gals'],
            stats['n_fail'],
            stats['est_g1'],
            stats['est_g2'],
            stats['est_size'],
            stats['est_stdv_g1'],
            stats['est_stdv_g2'],
            stats['est_stdm_g1'],
            stats['est_stdm_g2'],
            stats['est_stdv_size'],
            stats['est_stdm_size']
        )
    file_stats.write(line)


def get_shear(filename_result,id_result):

    use_dtype = dtype_results[args.method_id]
    res = np.loadtxt(filename_result,dtype=use_dtype)
    n_res = len(res)
    select_successful = np.abs(res['g1'] + 1j*res['g2'])<1
    n_success = sum(select_successful)
    n_fail = n_res - n_success
    res = res[select_successful]

    stats = {}
    
    if n_fail == n_res:
        stats['est_g1']  = 0.01
        stats['est_g2'] = 0.01
        stats['est_size'] = 1.
        stats['est_stdv_g1'] = 0.1
        stats['est_stdv_g2'] = 0.1
        stats['est_stdm_g1'] = 0.1
        stats['est_stdm_g2'] = 0.1
        stats['est_stdv_size'] = 0.1
        stats['est_stdm_size'] = 0.1
        stats['n_fail'] = n_res 
        stats['n_gals'] = n_res
        log.error('file %s - all measurements are error' % filename_result)
    else:
        stats['est_g1'] = np.mean(res['g1'])
        stats['est_g2'] = np.mean(res['g2'])
        stats['est_size'] = np.mean(res['size'])
        stats['est_stdv_g1'] = np.std(res['g1'],ddof=1)
        stats['est_stdv_g2'] = np.std(res['g2'],ddof=1)
        stats['est_stdm_g1'] = np.std(res['g1'],ddof=1)/np.sqrt(n_success)
        stats['est_stdm_g2'] = np.std(res['g2'],ddof=1)/np.sqrt(n_success)
        stats['est_stdv_size'] = np.std(res['size'],ddof=1)
        stats['est_stdm_size'] = np.std(res['size'],ddof=1)/np.sqrt(n_success)
        stats['n_fail'] = n_fail
        stats['n_gals'] = n_res

    # save histograms
    filename_fig = 'figs/fig.hist.%03d.png' % id_result
    save_histograms(res,filename_result,filename_fig)
    
    return stats

def save_histograms(res,title,filename_fig):

    n_bins = 20
    bins_e = np.linspace(-1,1,n_bins)
    bins_size = np.linspace(1.,3.,100)


    pl.figure()
    pl.subplot(1,2,1)
    pl.suptitle(title)
    try:
        pl.hist(res['g1'],bins=bins_e,histtype='step',color='r',label='g1')
        pl.hist(res['g2'],bins=bins_e,histtype='step',color='b',label='g2')
    except Exception,errmsg:
        log.error('hist failed')
        pl.hist(np.random.rand(1000)*2-1,bins=bins_e,histtype='step',color='r',label='g1')
        pl.hist(np.random.rand(1000)*2-1,bins=bins_e,histtype='step',color='b',label='g2')

    pl.legend()
    pl.xlabel('ellip [g]')
    pl.ylabel('histogram')

    pl.subplot(1,2,2)
    try:
        pl.hist(res['size'],bins=bins_size,histtype='step',color='r',label='FWHM_GAL/FWHM_PSF')
    except Exception,errmsg:
        log.error('hist failed')
        pl.hist(np.random.rand(1000)*2-1,bins=bins_e,histtype='step',color='r',label='FWHM_GAL/FWHM_PSF')

    pl.legend()
    pl.xlabel('size [FWHM_GAL/FWHM_PSF]')
    pl.ylabel('histogram')

    pl.savefig(filename_fig)
    log.info('saved %s' % filename_fig)
    pl.close()


def get_stats():

    truth_cat = np.loadtxt(args.filename_input,dtype=dtype_table_cat)

    filename_stats = '%s.%s.stats.cat' % (args.filename_input,args.method_id)
    file_stats = open(filename_stats,'w')
    log.debug('opened file %s' % filename_stats)

    stats_header = '# index n_gals n_fail g1 g2 size stdv_g1 stdv_g2 stdm_g1 stdm_g2 stdv_size stdm_size\n'
    file_stats.write(stats_header)

    for it,vt in enumerate(truth_cat):

        filename_result = results_filename_fmt % (it)
        log.debug('getting stats for file %s' % filename_result)
        stats = get_shear(filename_result,it)
        stats['index'] = it
        write_stats(file_stats,stats)

    file_stats.close()
    log.info('saved file %s' % filename_stats)

def write_bias_results(file_bias,bias_result=None,header=False):

    if header:
        header = '#  index ipsf_fwhm isnr n_fail psf_fwhm snr m1 m2 c1 c2 m1_std m2_std c1_std c2_std g_std \n'
        file_bias.write(header)
    else:
        fmt = '%d\t'*4 + '% 2.8f\t'*11 + '\n'
        line = fmt % (
            bias_result['index'] ,
            bias_result['ipsf_fwhm'],
            bias_result['isnr'] ,
            bias_result['n_fail'],
            bias_result['psf_fwhm'],
            bias_result['snr'],
            bias_result['m1'],
            bias_result['m2'],
            bias_result['c1'],
            bias_result['c2'],
            bias_result['m1_std'],
            bias_result['m2_std'],
            bias_result['c1_std'],
            bias_result['c2_std'],
            bias_result['g_std']
        )

        file_bias.write(line)


def get_weights_with_histograms():

    truth_cat = np.loadtxt(args.filename_input,dtype=dtype_table_cat)
    n_gals_total = sum(truth_cat['n_gals'])
    log.info('opened %s with %d rows and %d galaxies total' % (args.filename_input,len(truth_cat),n_gals_total))

    iall = 0
    
        
    for ipsf_fwhm,vpsf_fwhm in enumerate(config['grid']['psf_fwhm']):
           for isnr,vsnr in enumerate(config['grid']['snr']):

            select = np.array(truth_cat['ipsf_fwhm'] == ipsf_fwhm) * np.array(truth_cat['isnr'] == isnr)
            current_truth_cat = truth_cat[select]

            list_results = []
            for itruth,vtruth in enumerate(current_truth_cat):
                filename_result = results_filename_fmt % (vtruth['index'])
                current_res = tabletools.loadTable(filename_result,dtype=dtype_table_results_calib)
                list_results.append(current_res)

            results_ca_current = np.concatenate(list_results)

            current_snr = current_truth_cat['snr'][0]
            current_psf_fwhm = current_truth_cat['psf_fwhm'][0]

            select1 = (results_sv_example['fwhm_psf']*DES_PIXEL_SCALE < current_psf_fwhm + 0.05) * (results_sv_example['fwhm_psf']*DES_PIXEL_SCALE > current_psf_fwhm - 0.05)
            select2 = (results_sv_example['snr'] < current_snr + 2.5) * (results_sv_example['snr'] > current_snr - 2.5)
            results_sv_current = results_sv_example[select1*select2]
            n_current_gals = len(results_sv_current)
            log.info('snr=%2.2f psf_fwhm=%2.2f n_gals in that bin=%d' % (current_truth_cat['snr'][0],current_truth_cat['psf_fwhm'][0],n_current_gals))

            n_bins = 20
            bins_size=np.linspace(1.2,3,n_bins)
            # log.info('size of size bin %2.2f' % (bins_size[1]-bins_size[0]))
            title = 'snr=%2.2f psf_fwhm=%2.2f' % (current_truth_cat['snr'][0],current_truth_cat['psf_fwhm'][0])
            if n_current_gals > 2000:
                
                pl.subplot(1,2,1)
                h1,b1,_=pl.hist(results_ca_current['size'],bins=bins_size,histtype='step',normed=True,color='r',label='cal')
                h2,b2,_=pl.hist(results_sv_current['size'],bins=bins_size,histtype='step',normed=True,color='b',label='sv')
                pl.ylim([0,1.1*max([max(h1),max(h2)])])
                pl.legend()
                pl.xlabel('measured fwhm_ratio')

                pl.subplot(1,2,2)
                # abs_g_ca = np.abs(results_ca_current['g1']+1j*results_ca_current['g2'])
                # abs_g_sv = np.abs(results_sv_current['g1']+1j*results_sv_current['g2'])
                h1,b1,_=pl.hist(results_ca_current['g1'],bins=np.linspace(-1,1,n_bins),histtype='step',normed=True,color='r',label='cal')
                h2,b2,_=pl.hist(results_sv_current['g1'],bins=np.linspace(-1,1,n_bins),histtype='step',normed=True,color='b',label='sv')
                pl.legend()
                pl.xlabel('measured ellipticity')
                pl.ylim([0,1.1*max([max(h1),max(h2)])])
                
                pl.suptitle(title)
                filename_fig = 'figs/fig.hists.snr%02d.psf%02d.png' % (isnr,ipsf_fwhm)
                pl.savefig(filename_fig)
                pl.close()
                log.info('saved %s' % filename_fig)
               

            else:
                log.info('not enough gals to bother')



            # title = 'snr=%2.2f psf_fwhm=%2.2f' % (current_truth_cat['snr'][0],current_truth_cat['psf_fwhm'][0])
            # filename_fig = 'figs/fig.hist.allg.snr%03d.psf%03d.png' % (current_truth_cat['isnr'][0] , current_truth_cat['ipsf_fwhm'][0])
            # save_histograms(results_sv_current,title=title,filename_fig=filename_fig)


def plot_mc_cuts():

    size_cut = 1.2
    list_snr = []
    list_m_mean = []
    list_m_mean_std = []

    for isnr,vsnr in enumerate(config['grid']['snr']):

        # select  = (results_cal_all['size'] > size_cut) * (results_cal_all['isnr_true'] == isnr)
        select  = results_cal_all['isnr_true'] == isnr
        current_results_cal = results_cal_all[select]
        n_gals=len(current_results_cal)
        if n_gals < 5000:
            log.warning('number of galaxies in that bin is low')
            continue

        g1_true = current_results_cal['g1_true']
        g2_true = current_results_cal['g2_true']
        g1_bias = (current_results_cal['g1'] - g1_true)
        g2_bias = (current_results_cal['g2'] - g2_true)
        g1_err  = np.ones(n_gals)*np.std(current_results_cal['g1'],ddof=1)
        g2_err  = np.ones(n_gals)*np.std(current_results_cal['g2'],ddof=1)


        [c1,m1,C1cm] = get_line_fit(g1_true,g1_bias,g1_err)
        [c2,m2,C2cm] = get_line_fit(g2_true,g2_bias,g2_err)
        m1_std = np.sqrt(C1cm[1,1])
        m2_std = np.sqrt(C2cm[1,1])
        c1_std = np.sqrt(C1cm[0,0])
        c2_std = np.sqrt(C2cm[0,0])
        m_mean = (m1 + m2)/2.
        m_mean_std = np.sqrt((m1_std**2 + m2_std**2)/2.)

        list_snr.append(vsnr)
        list_m_mean.append(m_mean)
        list_m_mean_std.append(m_mean_std)

        log.info('n_gals=% 5d snr=%2.2f\tm=% 2.2e\tc=% 2.2f\tm_err=%2.2e\tc_err=%2.2e\tg_std=%2.2f' % (n_gals,vsnr,m1,c1,m1_std,c1_std,0))                

        
    pl.errorbar(list_snr,list_m_mean,yerr=list_m_mean_std)
    # pl.legend(mode='expand',ncol=5,loc='lower center')
    plot_add_req()
    pl.show()


def plot_ellipticity_variance():

    size_bin_width=0.2
    bins_size = np.arange(1.0,2.2,size_bin_width)

    bins_ell = np.linspace(-1,1,100)

    colorscale = plotstools.get_colorscale(len(bins_size)-1)

    n_snr = len(config['grid']['snr'])
    n_size = len(bins_size[:-1])

    pl.figure(1)

    variance_matrix = np.zeros([n_snr,n_size])

    for ib,vb in enumerate(bins_size[:-1]):

        list_stde = []
        list_snr = []
               
        bin_size_min = bins_size[ib]
        bin_size_max = bins_size[ib+1]

        for isnr,vsnr in enumerate(config['grid']['snr']):

            select  = (results_cal_all['size'] > bin_size_min) * (results_cal_all['size'] < bin_size_max) * (results_cal_all['isnr_true'] == isnr)
            current_results_cal = results_cal_all[select]
            n_gals=len(current_results_cal)
            if n_gals < 5000:
                log.warning('number of galaxies in that bin is low')
                continue

            stde = np.std(current_results_cal['g1'],ddof=1)

            variance_matrix[isnr,ib] = stde

            list_snr.append(vsnr)
            list_stde.append(stde)

            pl.figure(2)
            pl.hist(current_results_cal['g1'],histtype='step',bins=bins_ell,color=colorscale[isnr],normed=True,label='SNR=%2.2f' % vsnr)

        pl.figure(1)
        pl.plot(list_snr,list_stde,'d-',color=colorscale[ib],label='FWHM_RATIO=[%2.1f,%2.1f]' % (bin_size_min,bin_size_max))

        pl.figure(2)
        pl.ylabel('normalised histogram')
        pl.xlabel('measured ellipticity')
        pl.title('measured ellipticity for FWHM_RATIO=[%2.1f,%2.1f]' % (bin_size_min,bin_size_max))
        pl.legend()
        filename_fig = 'figs/fig.hist.ell.size%d.png' % ib
        pl.savefig(filename_fig)
        pl.close()
        log.info('saved %s' % filename_fig)

    pl.figure(1)
    pl.ylabel('ellipticity standard deviation')
    pl.xlabel('SNR')
    pl.legend()
    plotstools.adjust_limits(0.1,0.05)
    filename_fig = 'figs/fig.stde.png'
    pl.savefig(filename_fig)
    pl.close()
    log.info('saved %s' % filename_fig)

    # pl.figure(3)
    # pl.imshow(variance_matrix,interpolation='nearest')
    # pl.show()


def plot_mc_in_fwhm_ratio_bins():

    # get the colorscale over full range for n_colors
    size_bin_width=0.2
    # bins_size = np.arange(1.2,2.2,size_bin_width)
    bins_size = np.arange(1.2,2.01,0.2).tolist()
    bins_snr = [2.5,7.5,12.5,17.5,22.5,40]


    import colorsys
    n_colors = len(bins_size)-1
    HSV_tuples = [(x*1.0/n_colors, 0.75, 0.75) for x in range(n_colors)]
    RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)

    filename_data = 'data.bias.size-snr.txt'
    file_data = open(filename_data,'w')
    header = '# fwhm_mid snr_mid m m_std \n'
    file_data.write(header)

    pl.figure(1)


    for ib,vb in enumerate(bins_size[:-1]):
        
        list_snr = []
        list_m_mean = []
        list_m_mean_std = []
        

        # for isnr,vsnr in enumerate(config['grid']['snr']):
        for isnr,vsnr in enumerate(bins_snr[:-1]):

            bin_size_min = bins_size[ib]
            bin_size_max = bins_size[ib+1]
            bin_snr_min = bins_snr[isnr]
            bin_snr_max = bins_snr[isnr+1]

            # select  = (results_cal_all['size'] > bin_size_min) * (results_cal_all['size'] < bin_size_max) * (results_cal_all['isnr_true'] == isnr)
            select  = (results_cal_all['size'] > bin_size_min) * (results_cal_all['size'] < bin_size_max) * (results_cal_all['snr'] < bin_snr_max) * (results_cal_all['snr'] > bin_snr_min)
            current_results_cal = results_cal_all[select]
            n_gals=len(current_results_cal)
            if n_gals < 5000:
                log.warning('number of galaxies in that bin is low')
                continue

            fwhm_mid = (bins_size[ib] + bins_size[ib+1])/2.
            snr_mid = (bins_snr[isnr] + bins_snr[isnr+1])/2.

            g1_tru = current_results_cal['g1_true']
            g2_tru = current_results_cal['g2_true']
            g1_est = current_results_cal['g1']
            g2_est = current_results_cal['g2']
            g1_err  = np.ones(n_gals)*np.std(current_results_cal['g1'],ddof=1)
            g2_err  = np.ones(n_gals)*np.std(current_results_cal['g2'],ddof=1)

            [c1,m1,C1cm] = fitting.get_line_fit(g1_tru,g1_est,g1_err)
            [c2,m2,C2cm] = fitting.get_line_fit(g2_tru,g2_est,g2_err)
            m1_std = np.sqrt(C1cm[1,1])
            m2_std = np.sqrt(C2cm[1,1])
            c1_std = np.sqrt(C1cm[0,0])
            c2_std = np.sqrt(C2cm[0,0])
            m_mean = (m1 + m2)/2.
            m_mean_std = np.sqrt((m1_std**2 + m2_std**2)/2.)
            list_snr.append( snr_mid )
            list_m_mean.append( m_mean )
            list_m_mean_std.append( m_mean_std )
            line = '%2.6f\t%2.6f\t%2.6f\t%2.6f\n' % (fwhm_mid,snr_mid,m_mean,m_mean_std)
            file_data.write(line)

            log.info('n_gals=% 5d snr=%2.2f\tfwhm=[%2.2f,%2.2f]\tm=% 2.2e\tc=% 2.2f\tm_err=%2.2e\tc_err=%2.2e\tg_std=%2.4f' % (n_gals,vsnr,bin_size_min,bin_size_max,m1,c1,m1_std,c1_std,0))                

            # pl.figure(2)
            # pl.hist(current_results_cal['g1'],bins=np.linspace(-1,1,50),histtype='step',label='g1')
            # pl.hist(current_results_cal['g2'],bins=np.linspace(-1,1,50),histtype='step',label='g2')
            # pl.legend()
            # pl.title('snr=%2.2f size=[%2.2f-%2.2f]' %(vsnr,bins_size[ib],bins_size[ib+1]))
            # pl.show()
    
        pl.figure(1)    
        pl.errorbar(list_snr,list_m_mean,yerr=list_m_mean_std,label='fwhm_ratio=[%2.1f,%2.1f]' % (bin_size_min,bin_size_max),color=RGB_tuples[ib])

    file_data.close()
    log.info('closed %s' % filename_data)

    pl.xlabel('SNR')
    pl.ylabel('multiplicative bias m')
    pl.title('GREAT-DES noise bias calibration')
    pl.legend(loc='lower right')
    plotstools.adjust_limits()    
    plot_add_req()
    filename_fig = 'figs/fig.bias.size-snr.png'
    pl.savefig(filename_fig)
    pl.close()
    log.info('saved %s' % filename_fig)





    

def get_mc():

    filename_stats = '%s.%s.stats.cat' % (args.filename_input,args.method_id)
    truth_cat = np.loadtxt(args.filename_input,dtype=dtype_table_cat)
    stats_cat = np.loadtxt(filename_stats,dtype=dtype_table_stats)
    n_gals_total = sum(truth_cat['n_gals'])
    log.info('opened %s with %d rows and %d galaxies total' % (args.filename_input,len(truth_cat),n_gals_total))
    log.info('opened %s with %d rows' % (filename_stats,len(stats_cat)))

    filename_bias = '%s.%s.bias.cat' % (args.filename_input,args.method_id)
    file_bias = open(filename_bias,'w')
    write_bias_results(file_bias,header=True)

    iall = 0
    
    for ipsf_fwhm,vpsf_fwhm in enumerate(config['grid']['psf_fwhm']):

        for isnr,vsnr in enumerate(config['grid']['snr']):

            select = np.array(truth_cat['ipsf_fwhm'] == ipsf_fwhm) * np.array(truth_cat['isnr'] == isnr)
            current_stats_cat = stats_cat[select]
            current_truth_cat = truth_cat[select]

            g1_true = current_truth_cat['g1']
            g2_true = current_truth_cat['g2']
            g1_bias = (current_stats_cat['g1'] - g1_true)
            g2_bias = (current_stats_cat['g2'] - g2_true)
            g1_err  = stats_cat[select]['stdm_g1']
            g2_err  = stats_cat[select]['stdm_g2']

            [c1,m1,C1cm] = get_line_fit(g1_true,g1_bias,g1_err)
            [c2,m2,C2cm] = get_line_fit(g2_true,g2_bias,g2_err)
            m1_std = np.sqrt(C1cm[1,1])
            m2_std = np.sqrt(C2cm[1,1])
            c1_std = np.sqrt(C1cm[0,0])
            c2_std = np.sqrt(C2cm[0,0])

            g_std = np.mean(current_stats_cat['stdv_g1'])
            n_fail = float(sum(current_stats_cat['n_fail'])) / float(sum(current_stats_cat['n_gals']))

            pl.figure()
            pl.errorbar(g1_true,g1_bias,yerr=g1_err,fmt='b.',label='m1=%2.3f +/- %2.3f' % (m1,m1_std))
            pl.plot(g1_true,g1_true*m1 + c1,'b-')
            pl.errorbar(g2_true,g2_bias,yerr=g2_err,fmt='r.',label='m1=%2.3f +/- %2.3f' % (m2,m2_std))
            pl.plot(g2_true,g2_true*m2 + c2,'r-')
            pl.title('snr=%2.2f psf_fwhm=%2.2f' % (vsnr,vpsf_fwhm))
            pl.legend()
            filename_fig = 'figs/fig.bias.snr%03d.psf%03d.png' % (isnr,ipsf_fwhm)
            pl.savefig(filename_fig)
            pl.close()
            log.info('saved %s' % filename_fig)

            bias_result = {}
            bias_result['index']  = iall
            bias_result['ipsf_fwhm'] = ipsf_fwhm
            bias_result['isnr']  = isnr
            bias_result['psf_fwhm'] = vpsf_fwhm
            bias_result['snr'] = vsnr
            bias_result['m1'] = m1
            bias_result['m2'] = m2
            bias_result['c1'] = c1
            bias_result['c2'] = c2
            bias_result['m1_std'] = m1_std
            bias_result['m2_std'] = m2_std
            bias_result['c1_std'] = c1_std
            bias_result['c2_std'] = c2_std
            bias_result['g_std'] = g_std
            bias_result['n_fail'] = n_fail
            log.info('psf_fwhm=%2.2f\tsnr=%2.2f\tm=% 2.2e\tc=% 2.2f\tm_err=%2.2e\tc_err=%2.2e\tg_std=%2.2f\tn_fail=%2.2f' % (vpsf_fwhm, vsnr,m1,c1,m1_std,c1_std,g_std,n_fail))                

            write_bias_results(file_bias, bias_result=bias_result)

            iall+=1


def plot_mc():

        filename_bias = '%s.%s.bias.cat' % (args.filename_input,args.method_id)
        results_mc = np.loadtxt(filename_bias,dtype=dtype_table_bias)  

        psf_fwhm_colors=['r','g','b']
        
        for ipsf_fwhm,vpsf_fwhm in enumerate(config['grid']['psf_fwhm']):
            # for isnr,vsnr in enumerate(config['grid']['snr']):

            select = (results_mc['ipsf_fwhm'] == ipsf_fwhm)
            results_select = results_mc[select]
            sort = np.argsort(results_select['snr'])
            results_select=results_select[sort] 

            pl.figure(1)
            pl.errorbar(results_select['snr'],results_select['m1'],yerr=results_select['m1_std'],fmt='r+--')
            pl.errorbar(results_select['snr'],results_select['m2'],yerr=results_select['m2_std'],fmt='rx:')
            xlim_add = (max(results_select['snr']) - min(results_select['snr']))*0.1
            pl.xlim([min(results_select['snr'])-xlim_add,max(results_select['snr'])+xlim_add])
            plot_add_req()
            pl.xlabel('SNR')
            pl.ylabel('multiplicative bias m')
            pl.xticks(results_select['snr'])
            title_str = 'PSF FWHM=%2.2f ' % (vpsf_fwhm)
            pl.title(title_str)
            pl.legend(['%s g1' % args.method_id,'%s g2' % args.method_id],loc='lower left',ncol=2,mode='expand')
            filename_fig = 'figs/fig.%s.psf_fwhm%d.png' % (filename_bias,ipsf_fwhm)
            pl.savefig(filename_fig)
            log.info('saved %s' % filename_fig)
            pl.close()

            m_mean = (results_select['m1'] + results_select['m2'])/2.
            m_mean_std = np.sqrt((results_select['m1_std']**2 + results_select['m2_std']**2)/2.)
            pl.figure(2)
            # pl.errorbar(results_select['snr'],results_select['m1'],yerr=results_select['m1_std'],fmt=psf_fwhm_colors[ipsf_fwhm]+'--',label='psf_fwhm=%2.2f m1' % vpsf_fwhm)
            # pl.errorbar(results_select['snr'],results_select['m2'],yerr=results_select['m2_std'],fmt=psf_fwhm_colors[ipsf_fwhm]+':',label='psf_fwhm=%2.2f m2' % vpsf_fwhm)
            pl.errorbar(results_select['snr'],m_mean,yerr=m_mean_std,fmt=psf_fwhm_colors[ipsf_fwhm]+'-',label='psf_fwhm=%2.2f m' % vpsf_fwhm)

        pl.figure(2)
        xlim_add = (max(results_select['snr']) - min(results_select['snr']))*0.1
        pl.xlim([min(results_select['snr'])-xlim_add,max(results_select['snr'])+xlim_add])
        plot_add_req()
        pl.xlabel('SNR')
        pl.ylabel('multiplicative bias m')
        pl.xticks(results_select['snr'])
        pl.legend(loc='lower left',ncol=2,mode='expand')
        filename_fig = 'figs/fig.%s.png' % (filename_bias)
        pl.savefig(filename_fig)
        log.info('saved %s' % filename_fig)
        pl.close()

        # pl.figure()
        # pl.plot(snr[sort],n_fail[sort],'ro-')
        # xlim_add = (max(snr) - min(snr))*0.1
        # pl.xlim([min(snr)-xlim_add,max(snr)+xlim_add])
        # plot_add_nfail_req()
        # pl.xlabel('SNR')
        # pl.ylabel('n_fail [%]')
        # pl.xticks(snr)
        # title_str = 'hlr=%2.2f' % vhlr
        # pl.title(title_str)
        # filename_fig = 'fig.%s.hlr%d.nfail.png' % (filename_stats,ihlr)
        # pl.savefig(filename_fig)
        # log.info('saved %s' % filename_fig)
        # pl.close()

    # # plot the surface

    # angrad = float(config['gal']['shear']['beta'].split()[0])/180.*np.pi
    # e_true = config['gal']['shear']['g'] * np.exp(1j*angrad*2.)   # angle of galaxy orientation, not complex shear angle

    # n_grid = 100
    # pl.figure()
    # bias = (stats_cat['g1'] - e_true.real)/e_true.real
    # grid=surf_interp(truth_cat['hlr'],truth_cat['snr'],bias)
    # pl.imshow(grid, origin='lower',extent=(min(truth_cat['hlr']),max(truth_cat['hlr']),min(truth_cat['snr']),max(truth_cat['snr'])), aspect='auto')
    # pl.ylabel('SNR')
    # pl.xlabel('hlr [arcsec]')
    # title_str = '%s' % args.method_id
    # pl.title(title_str)
    # filename_fig = 'fig.%s.hlr_snr.%s.png' % (filename_stats,args.method_id)
    # pl.colorbar()
    # pl.savefig(filename_fig)
    # log.info('saved %s' % filename_fig)
    # pl.close()


def main():

    global log , config , args

    description = 'Get statistics and plot results of noise bias calibration runs'
    parser = argparse.ArgumentParser(description=description, add_help=True)
    parser.add_argument('-v', '--verbosity', type=int, action='store', default=2, choices=(0, 1, 2, 3 ), help='integer verbosity level: min=0, max=3 [default=2]')
    parser.add_argument('-i', '--filename_input', default='sha1-O1.cat',type=str, action='store', help='name of the output catalog')
    parser.add_argument('-c', '--filename_config', default='sha1.yaml',type=str, action='store', help='name of the yaml config file')
    parser.add_argument('-m', '--method_id', default='hsm',type=str, action='store', help='name of the yaml config file')
    
    args = parser.parse_args()
    # Parse the integer verbosity level from the command line args into a logging_level string
    logging_levels = { 0: logging.CRITICAL, 
                       1: logging.WARNING,
                       2: logging.INFO,
                       3: logging.DEBUG }
    logging_level = logging_levels[args.verbosity]
    
    log = logging.getLogger("nbc1_plots") 
    log.setLevel(logging_level)  
    log_formatter = logging.Formatter("%(asctime)s  %(name)s  %(levelname)s  %(message)s ")
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(log_formatter)
    log.addHandler(stream_handler)    

    log.info(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))

    config = yaml.load(open(args.filename_config))

    filename_results_sv_example = '/home/tomek/projects/131118_nbc1/sv-clusters-shears/results_sv_clusters.all.fits'
    filename_cal_all = '/home/tomek/projects/131118_nbc1/run-tiled-002/cleaned_calib.all.fits'

    global results_filename_fmt
    # results_filename_fmt = 'calib.v4.2013.12.20/cleaned_calib.v4.2013.12.20/nbc_%03d.fits.im3.cleaned.cat'
    results_filename_fmt = 'calib.v4.2014.01.24/cleaned_calib.v4.2014.01.24/nbc_%03d.fits.im3.cleaned.cat'

    global results_sv_example
    results_sv_example = tabletools.loadTable(filename_results_sv_example)

    global results_cal_all
    results_cal_all = tabletools.loadTable(filename_cal_all)

    # get_stats()  
    # get_weights_with_histograms()
    plot_mc_in_fwhm_ratio_bins()
    plot_ellipticity_variance()
    # get_mc()
    # plot_mc()
    # plot_mc_cuts()

    log.info(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))


if __name__ == '__main__':
    main()