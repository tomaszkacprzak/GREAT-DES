import numpy as np
import pylab as pl
import galsim, sys, logging, yaml, argparse, time, copy, itertools, fitting
sys.path.append('/home/tomek/code/tktools')
import tabletools
import plotstools
from nbc2_dtypes import *

logging_level = logging.INFO
log = logging.getLogger("nbc2_analy") 
log.setLevel(logging_level)  
log_formatter = logging.Formatter("%(asctime)s %(name)s %(levelname)s   %(message)s", "%Y-%m-%d %H:%M:%S")
stream_handler = logging.StreamHandler(sys.stdout)
stream_handler.setFormatter(log_formatter)
if log.handlers == [] : log.addHandler(stream_handler)
log.propagate = False
       
req1_dg = 0.003
req2_dg = 0.02
req3_dg = 0.1
req1_nfail = 0.01
req2_nfail = 0.05

figsize = [10,7]
params = {'legend.fontsize': 10,
        'legend.linewidth': 2}
pl.rcParams.update(params)

def plot_add_requirements(mult=1.):

    req1_dg = 0.003
    req2_dg = 0.02
    req3_dg = 0.1
    req1_nfail = 0.01
    req2_nfail = 0.05

    corner = pl.xlim()[0]
    length = abs(pl.xlim()[1]) + abs(pl.xlim()[0])
    pl.gca().add_patch(pl.Rectangle(  (corner, 0 -mult*req3_dg), length , 2*mult*req3_dg , facecolor = '0.9' , edgecolor='k' ))
    pl.gca().add_patch(pl.Rectangle(  (corner, 0 -mult*req2_dg), length , 2*mult*req2_dg , facecolor = '0.8' , edgecolor='k' ))
    pl.gca().add_patch(pl.Rectangle(  (corner, 0 -mult*req1_dg), length , 2*mult*req1_dg , facecolor = '0.7' , edgecolor='k' ))

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


def get_shear_results(list_all_res, list_all_tru , use_ell_for_truth=False ):

    n_split = len(list_all_tru)
    list_shears = []

    for il in range(n_split):

        try:
            # stats=get_shear_estimator(list_all_res[il],list_all_tru[il])
            stats=get_shear_estimator(list_all_res[il])
        except Exception,errmsg:
            print errmsg
            log.error('not enough galaxies in sample %d' % len(list_all_res[il]))
            continue
        if use_ell_for_truth:
            true_e = (1-list_all_tru[il]['sf_q'])/(1+list_all_tru[il]['sf_q'])*np.exp(2.*1j*list_all_tru[il]['sf_phi'])
            stats['tru_g1'] = np.mean(true_e.real + list_all_tru[il]['g1_true'])
            stats['tru_g2'] = np.mean(true_e.imag + list_all_tru[il]['g2_true'])      
        else:
            stats['tru_g1'] = np.mean(list_all_tru[il]['g1_true'])
            stats['tru_g2'] = np.mean(list_all_tru[il]['g2_true'])

        if 'fake' in config['methods'][args.method]: 
            if config['methods'][args.method]['fake']:
                fake_e1 = 0.1*np.random.randn(len(cat_tru)) + stats['tru_g1']
                fake_e2 = 0.1*np.random.randn(len(cat_tru)) + stats['tru_g2']
                stats['est_g1'] = np.mean( fake_e1 )
                stats['est_g2'] = np.mean( fake_e2 )

        list_shears.append( stats )

    # import pdb; pdb.set_trace()
    # perm = np.random.permutation(len(all_tru))[0:10000]
    # pl.scatter(all_tru['snr'][perm],all_res['snr'][perm])
    # pl.show()

    shears = np.concatenate(list_shears)


    return shears




def get_mc(selection_string,tag):

    # these are the columns that we select for the results catalog (we are interested only in few, don't need all columns)
    col_g1 = config['methods'][args.method]['col_g1']
    col_g2 = config['methods'][args.method]['col_g2']
    col_size = config['methods'][args.method]['col_size']
    col_snr = config['methods'][args.method]['col_snr']
    cols_res = [col_g1, col_g2, col_size , col_snr]
    cols_tru = ['g1_true', 'g2_true' , 'snr','flux' , 'hsm_obs_sigma', 'sf_q' , 'sf_phi']
    
    split_all_res, split_all_tru = get_selection_split(selection_string, cols_res, cols_tru)
    all_tru = np.concatenate(split_all_tru)
    all_res = np.concatenate(split_all_res)

    stats = get_shear_results(split_all_res, split_all_tru )

    # for i in stats: print i['tru_g1'], i['est_g1'], i['est_g1'] - i['tru_g1'] 


    g1_tru = stats['tru_g1']
    g2_tru = stats['tru_g2']
    g1_est  = stats['est_g1'] 
    g2_est  = stats['est_g2'] 
    g1_err  = stats['est_stdm_g1']
    g2_err  = stats['est_stdm_g2']
    g1_bias = g1_est-g1_tru
    g2_bias = g2_est-g2_tru

    [c1,m1,C1cm] = get_line_fit(g1_tru,g1_bias,g1_err)
    [c2,m2,C2cm] = get_line_fit(g2_tru,g2_bias,g2_err)
    d1 = np.std(g1_bias - (g1_tru*m1+c1),ddof=1)
    d2 = np.std(g2_bias - (g2_tru*m2+c2),ddof=1)
    [c1,m1,C1cm] = get_line_fit(g1_tru,g1_bias,np.ones_like(g1_tru)*d1)
    [c2,m2,C2cm] = get_line_fit(g2_tru,g2_bias,np.ones_like(g2_tru)*d2)
    m1_std = np.sqrt(C1cm[1,1])
    m2_std = np.sqrt(C2cm[1,1])
    c1_std = np.sqrt(C1cm[0,0])
    c2_std = np.sqrt(C2cm[0,0])
    log.info('errorbar of each point sig=%2.4f , after jackknife sig=%2.4f , m_std=%2.4f' , g1_err[0], d1 , m1_std)


    g1_std = np.mean(stats['est_stdv_g1'])
    g2_std = np.mean(stats['est_stdv_g2'])
    n_fail = float(sum(stats['n_fail'])) / float(sum(stats['n_gals']))

    bias_result=np.empty(1,dtype=dtype_bias)
    bias_result['m1'] = m1
    bias_result['m2'] = m2
    bias_result['c1'] = c1
    bias_result['c2'] = c2
    bias_result['c1_std'] = c1_std
    bias_result['c2_std'] = c2_std
    bias_result['m1_std'] = m1_std
    bias_result['m2_std'] = m2_std
    bias_result['g1_std'] = g1_std
    bias_result['g2_std'] = g2_std
    bias_result['n_fail'] = n_fail

    global fig_counter
    pl.figure()
    pl.subplot(2,2,1)
    pl.errorbar(g1_tru,g1_bias,yerr=g1_err,fmt='b.',label='m1=%2.3f +/- %2.3f' % (m1,m1_std))
    pl.plot(g1_tru,g1_tru*m1 + c1,'b-')
    pl.errorbar(g2_tru,g2_bias,yerr=g2_err,fmt='r.',label='m2=%2.3f +/- %2.3f' % (m2,m2_std))
    pl.plot(g2_tru,g2_tru*m2 + c2,'r-')
    ylim=list(pl.ylim()); ylim[1]*=1.5; pl.ylim(ylim)
    pl.title(tag)
    pl.legend()
    # pl.show()

    pl.subplot(2,2,2)
    pl.hist(all_res[col_g1], bins=np.linspace(-1,1,100),histtype='step',normed=True , label='e1', color='b')
    pl.hist(all_res[col_g2], bins=np.linspace(-1,1,100),histtype='step',normed=True , label='e2', color='r')
    ylim=list(pl.ylim()); ylim[1]*=1.5; pl.ylim(ylim)
    pl.legend()

    pl.subplot(2,2,3)
    pl.hist(all_res[col_size], bins=np.linspace(0,4,100),histtype='step',normed=True , label='Rgpp/Rpp', color='r')
    pl.hist(all_tru['hsm_obs_sigma'], bins=np.linspace(2,4,100),histtype='step',normed=True , label='hsm_mom_sigma', color='b')
    ylim=list(pl.ylim()); ylim[1]*=1.5; pl.ylim(ylim)
    pl.legend()

    filename_fig = 'figs/fig.bias+hist.%s.png' % (tag)
    pl.savefig(filename_fig)
    pl.close()
    log.info('saved %s' % filename_fig)

    return bias_result, all_tru, all_res


def get_shear_estimator(res,tru=None):

    stats = np.zeros(1,dtype=dtype_stats)

    n_res = len(res)
    col_g1 = config['methods'][args.method]['col_g1']
    col_g2 = config['methods'][args.method]['col_g2']
    col_size = config['methods'][args.method]['col_size']
    flip_g1 =  config['methods'][args.method]['flip_g1']
    flip_g2 =  config['methods'][args.method]['flip_g2']
    select_successful = ((np.abs(res[col_g1] + 1j*res[col_g2]))<1) * (~np.isnan(res[col_size])) * (~np.isinf(res[col_size])) * (~np.isnan(res[col_g1])) * (~np.isinf(res[col_g1]))
    n_success = len(np.nonzero(select_successful)[0])
    n_fail = n_res - n_success
    res_success = res[select_successful]

    if (n_res==0) or (n_success<2): 
        raise Exception('no galaxies in sample')
   
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
        log.error('all measurements are error')
    else:
        res_g1 = flip_g1*res_success[col_g1]
        res_g2 = flip_g2*res_success[col_g2]
        if tru!=None:
            tru_success = tru[select_successful]
            true_e = (1-tru_success['sf_q'])/(1+tru_success['sf_q'])*np.exp(2.*1j*tru_success['sf_phi'])
            stats['est_g1'] = np.mean(true_e.real)
            stats['est_g2'] = np.mean(true_e.imag)      
            stats['est_size'] = np.mean(res_success[col_size])
            stats['est_stdv_g1'] = np.std(true_e.real,ddof=1)
            stats['est_stdv_g2'] = np.std(true_e.imag,ddof=1)
            stats['est_stdm_g1'] = np.std(true_e.real,ddof=1)/np.sqrt(n_success)
            stats['est_stdm_g2'] = np.std(true_e.imag,ddof=1)/np.sqrt(n_success)
            stats['est_stdv_size'] = np.std(res_success[col_size],ddof=1)
            stats['est_stdm_size']  =np.std(res_success[col_size],ddof=1)/np.sqrt(n_success)
            stats['n_fail'] = n_fail
            stats['n_gals'] = n_res

        else:
            stats['est_g1'] = np.mean(res_g1)
            stats['est_g2'] = np.mean(res_g2)
            stats['est_size'] = np.mean(res_success[col_size])
            stats['est_stdv_g1'] = np.std(res_g1,ddof=1)
            stats['est_stdv_g2'] = np.std(res_g2,ddof=1)
            stats['est_stdm_g1'] = np.std(res_g1,ddof=1)/np.sqrt(n_success)
            stats['est_stdm_g2'] = np.std(res_g2,ddof=1)/np.sqrt(n_success)
            stats['est_stdv_size'] = np.std(res_success[col_size],ddof=1)
            stats['est_stdm_size']  =np.std(res_success[col_size],ddof=1)/np.sqrt(n_success)
            stats['n_fail'] = n_fail
            stats['n_gals'] = n_res


    log.debug('got shear estimator for %d galaxies with est_stdm_g1=%2.4f' , n_success, stats['est_stdm_g1'])
 

    # check nans and infs
    for key in stats.dtype.names:
        if any(np.isinf(stats[key])) or any(np.isnan(stats[key])):
            import pdb; pdb.set_trace()


    return stats


def get_selection_DES(selection_string,cols,n_files=30):

    filelist_i = np.loadtxt('filelist_r.txt',dtype='a1024')
    list_results = []
    for filename_des in filelist_i[:n_files]:
        cat_res=tabletools.loadTable(filename_des,log=1)
        exec selection_string
        res_select = cat_res[select][cols]
        list_results.append(res_select)
    results = np.concatenate(list_results)

    log.info('selected DES galaxies %d' , len(results))

    return results

        

def get_selection(selection_string, cols_res, cols_tru):

    list_all_res, list_all_tru = get_selection_split(selection_string, cols_res, cols_tru)
    all_tru = np.concatenate(list_all_tru)
    all_res = np.concatenate(list_all_res)

    return all_res, all_tru



def get_selection_split(selection_string, cols_res, cols_tru,n_split=30):

    results_filename_fmt = config['methods'][args.method]['filename_results']  
    truth_filename_fmt = config['filename_truth']  
    list_shears = []
    list_all_res = []
    list_all_tru = []

    for ig,vg in enumerate(config['shear']):
        
        list_results = []

        id_first = args.first
        id_last = id_first + args.num

        for ip in range(id_first,id_last):
           
            filename_tru = truth_filename_fmt % (ip,ig)
            filename_res = results_filename_fmt % (ip,ig)
            try:
                cat_tru = tabletools.loadTable(filename_tru,log=1,remember=True)
            except:
                log.error('file %s not found' % filename_tru )
                continue

            for col in cols_tru:
                if col not in cat_tru.dtype.names:
                    raise Exception('column %s not found in truth catalog %s' % (col,filename_tru))

            try:
                cat_res = tabletools.loadTable(filename_res,log=1,remember=True)
            except:
                log.error('file %s not found' % filename_res )
                continue
 
            if 'patch_g2' in config['methods'][args.method]:
                if config['methods'][args.method]['patch_g2']:
                    log.info('patching g2--')
                    cat_res = tabletools.appendColumn(rec=cat_res, arr=cat_res['hsm_cor_g1'].copy(), name='hsm_cor_g2', dtype='f4')          
 

            for col in cols_res:
                if col not in cat_res.dtype.names:
                    raise Exception('column %s not found in results catalog %s' % (col,filename_res))

            exec selection_string
            if len(np.nonzero(select)[0]) < 1:
                log.debug('select didnt give any results %s' % selection_string)
                # import pdb;pdb.set_trace()
                # raise Exception('select didnt give any results %s' % selection_string)


            selected_res = cat_res[select][cols_res]
            selected_tru = cat_tru[select][cols_tru]


            # selected_res = selected_res.astype( ['f4']*len(cols_res) )
            selected_tru = selected_tru.astype( dtype={ 'formats':['f4']*len(cols_tru) , 'names': cols_tru })

            list_all_res.append(selected_res)
            list_all_tru.append(selected_tru)

    all_res = np.concatenate(list_all_res)
    all_tru = np.concatenate(list_all_tru)

    n_per_split = len(all_res)/n_split
    list_split_res = []
    list_split_tru = []
    for il in range(n_split):
        istart = il*n_per_split
        iend = (il+1)*n_per_split
        list_split_res.append(all_res[istart:iend])
        list_split_tru.append(all_tru[istart:iend])
   
    n_total = 0
    for res in list_all_res:
        n_total += len(res)

        
    log.info('selected %d parts with average %d, total %d' % (len(list_all_res) , float(n_total)/float(len(list_all_res)), n_total) )
    return list_split_res, list_split_tru




def get_plots():

    size_cut = 0.65
    snr_cut = 10

    list_plots = [22]

    # --------------------------------------------------------------------------------------------------
    n_total = 0
    if 1 in list_plots:
        list_bias_result = []
        list_psf = [0.8,0.9,1.0,1.1,1.2,1.3]
        for ipsf,vpsf in enumerate(list_psf):

            tag = 'psf_fwhm=%2.2f' % vpsf
            # selection_string = "select=  ( np.isclose(cat_tru['psf_fwhm'],%1.1f)) * ( (cat_tru['fwhm']*0.27/cat_tru['psf_fwhm']) > %f)" % (vpsf,size_cut)
            selection_string = "select=  ( np.isclose(cat_tru['psf_fwhm'],%1.1f,atol=0.01)) * (cat_tru['sf_hlr'] > %f) * (cat_tru['snr'] > %f) " % (vpsf,size_cut,snr_cut)
            log.info(selection_string)
            bias_result, all_res, all_tru = get_mc(selection_string,tag)
            n_total+=len(all_res)
            list_bias_result.append(bias_result)

        print 'total number in all bins n=' , n_total
        bias_results = np.concatenate(list_bias_result)

        m_mean = (bias_results['m1'] + bias_results['m2'])/2.
        m_mean_std = np.sqrt((bias_results['m1_std']**2 + bias_results['m2_std']**2))/np.sqrt(2.)

        pl.figure(figsize=figsize)
        pl.errorbar(list_psf,bias_results['m1'],yerr=bias_results['m1_std'],fmt=':r',label='m1')
        pl.errorbar(list_psf,bias_results['m2'],yerr=bias_results['m2_std'],fmt=':b',label='m2')
        pl.errorbar(list_psf,m_mean,m_mean_std,fmt='k',label='<m>')
        pl.legend()
        pl.xlabel('PSF_FWHM')
        pl.ylabel('multiplicative bias')
        plotstools.adjust_limits()
        plot_add_requirements()
        filename_fig = 'figs/bias_vs_psf.png' 
        pl.savefig(filename_fig)
        log.info('saved %s', filename_fig)

    # --------------------------------------------------------------------------------------------------
    # plot m vs SNR in bins of true SNR. use only galaxies with hsm_obs_sigma>x
    if 2 in list_plots:
        list_snr_centers = [10,15,20,25,30,40,50,60,70,80,90]
        list_snr_edges = plotstools.get_bins_edges(list_snr_centers)
        list_bias_result = []
        for isnr in range(1,len(list_snr_edges)):

            tag = 'snr=%.2f' % list_snr_centers[isnr-1]
            snr = list_snr_edges[isnr]
            selection_string = "select=  (cat_tru['snr'] > %2.4f) * (cat_tru['snr'] < %f) * (cat_tru['sf_hlr'] > %f) " % (
                                                list_snr_edges[isnr-1],
                                                list_snr_edges[isnr],
                                                size_cut)

            log.info(selection_string)
            bias_result, all_res, all_tru = get_mc(selection_string,tag)
            list_bias_result.append(bias_result)

        bias_results = np.concatenate(list_bias_result)
        m_mean = (bias_results['m1'] + bias_results['m2'])/2.
        m_mean_std = np.sqrt((bias_results['m1_std']**2 + bias_results['m2_std']**2))/np.sqrt(2.)


        pl.figure(figsize=figsize)
        pl.errorbar(list_snr_centers,bias_results['m1'],yerr=bias_results['m1_std'],fmt=':r',label='m1')
        pl.errorbar(list_snr_centers,bias_results['m2'],yerr=bias_results['m2_std'],fmt=':b',label='m2')
        pl.errorbar(list_snr_centers,m_mean,m_mean_std,fmt='k',label='<m>')

        pl.legend()
        pl.xlabel('SNR')
        pl.ylabel('multiplicative bias')
        plotstools.adjust_limits()
        plot_add_requirements()
        filename_fig = 'figs/bias_vs_snr.png' 
        pl.savefig(filename_fig)
        log.info('saved %s', filename_fig)

    # snr observed
    if 22 in list_plots:
        list_snr_centers = [10,15,20,25,30,40,50,60,70,80,90]
        list_snr_edges = plotstools.get_bins_edges(list_snr_centers)
        list_bias_result = []
        for isnr in range(1,len(list_snr_edges)):

            tag = 'snrobs=%.2f' % list_snr_centers[isnr-1]
            snr = list_snr_edges[isnr]
            selection_string = "select=  (cat_res['snr'] > %2.4f) * (cat_res['snr'] < %f) * (cat_tru['sf_hlr'] > %f) * (cat_tru['snr']) > %f" % (
                                                list_snr_edges[isnr-1],
                                                list_snr_edges[isnr],
                                                size_cut, snr_cut)

            log.info(selection_string)
            bias_result, all_res, all_tru = get_mc(selection_string,tag)
            list_bias_result.append(bias_result)

        bias_results = np.concatenate(list_bias_result)
        m_mean = (bias_results['m1'] + bias_results['m2'])/2.
        m_mean_std = np.sqrt((bias_results['m1_std']**2 + bias_results['m2_std']**2))/np.sqrt(2.)


        pl.figure(figsize=figsize)
        pl.errorbar(list_snr_centers,bias_results['m1'],yerr=bias_results['m1_std'],fmt=':r',label='m1')
        pl.errorbar(list_snr_centers,bias_results['m2'],yerr=bias_results['m2_std'],fmt=':b',label='m2')
        pl.errorbar(list_snr_centers,m_mean,m_mean_std,fmt='k',label='<m>')

        pl.legend()
        pl.xlabel('SNR')
        pl.ylabel('multiplicative bias')
        plotstools.adjust_limits()
        plot_add_requirements()
        filename_fig = 'figs/bias_vs_snrobs.png' 
        pl.savefig(filename_fig)
        log.info('saved %s', filename_fig)

    # --------------------------------------------------------------------------------------------------
    # plot m vs true size, using true FWHM measured from noise -free images, restrict SNR>x

    if 3 in list_plots:
        list_size_centers = np.linspace(0.5,1.5,10)
        list_size_edges = plotstools.get_bins_edges(list_size_centers)
        list_bias_result = []
        for isize in range(1,len(list_size_edges)):

            tag = 'hsmsize=%.2f' % list_size_centers[isize-1]
            size = list_size_edges[isize]
            selection_string = "select=     (cat_res['rgpp_rp_1']/cat_tru['psf_fwhm']  > %f) * (cat_res['rgpp_rp_1']/cat_tru['psf_fwhm'] < %f) * (cat_tru['sf_hlr']  < %f) * (cat_res['snr']>%f) " % (list_size_edges[isize-1],list_size_edges[isize],size_cut,snr_cut)

            log.info(selection_string)
            bias_result, all_res, all_tru = get_mc(selection_string,tag)
            list_bias_result.append(bias_result)

        bias_results = np.concatenate(list_bias_result)
        m_mean = (bias_results['m1'] + bias_results['m2'])/2.
        m_mean_std = np.sqrt((bias_results['m1_std']**2 + bias_results['m2_std']**2))/np.sqrt(2.)

        pl.figure(figsize=figsize)
        pl.errorbar(list_size_centers,bias_results['m1'],yerr=bias_results['m1_std'],fmt='r',label='m1')
        pl.errorbar(list_size_centers,bias_results['m2'],yerr=bias_results['m2_std'],fmt='b',label='m2')
        # pl.errorbar(list_size_centers,m_mean,m_mean_std,fmt='k',label='<m>')
        pl.legend()
        pl.xlabel('size')
        pl.ylabel('multiplicative bias')
        plotstools.adjust_limits()
        plot_add_requirements()
        filename_fig = 'figs/bias_vs_size.png' 
        pl.savefig(filename_fig)
        log.info('saved %s', filename_fig)

    # --------------------------------------------------------------------------------------------------
    # plot m vs FWHM size, using true FWHM measured from noise -free images, 
    
    if 4 in list_plots:
        list_size_centers = [1.,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6]
        list_size_edges = plotstools.get_bins_edges(list_size_centers)
        list_bias_result = []
        for isize in range(1,len(list_size_edges)):

            tag = 'fwhmratio=%.2f' % list_size_centers[isize-1]
            size = list_size_edges[isize]
            selection_string = "select= (cat_tru['hsm_obs_sigma']*2.355/cat_tru['psf_fwhm'] > %2.4f) * (cat_tru['hsm_obs_sigma']*2.355*0.27/cat_tru['psf_fwhm'] < %2.4f) * * (cat_tru['sf_hlr'] > %f) * (cat_tru['snr'] > %f) " % (list_size_edges[isize-1],list_size_edges[isize],size_cut,snr_cut)

            bias_result, all_res, all_tru = get_mc(selection_string,tag)
            list_bias_result.append(bias_result)

        bias_results = np.concatenate(list_bias_result)
        m_mean = (bias_results['m1'] + bias_results['m2'])/2.
        m_mean_std = np.sqrt((bias_results['m1_std']**2 + bias_results['m2_std']**2))/np.sqrt(2.)

        pl.figure(figsize=figsize)
        pl.errorbar(list_size_centers,bias_results['m1'],yerr=bias_results['m1_std'],fmt='r',label='m1')
        pl.errorbar(list_size_centers,bias_results['m2'],yerr=bias_results['m2_std'],fmt='b',label='m2')
        # pl.errorbar(list_size_centers,m_mean,m_mean_std,fmt='k',label='<m>')
        pl.legend()
        pl.xlabel('size')
        pl.ylabel('multiplicative bias')
        plotstools.adjust_limits()
        plot_add_requirements()
        filename_fig = 'figs/bias_vs_fwhmratio.png' 
        pl.savefig(filename_fig)
        log.info('saved %s', filename_fig)

    # --------------------------------------------------------------------------------------------------
    # plot m vs FWHM size, using true FWHM measured from noise -free images, 
    
    if 5 in list_plots:
        list_size_centers = [0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8]
        list_size_edges = plotstools.get_bins_edges(list_size_centers)
        list_bias_result = []
        for isize in range(1,len(list_size_edges)):

            tag = 'fwhmratio=%.2f' % list_size_centers[isize-1]
            size = list_size_edges[isize]
            selection_string = "select= (cat_res['rgpp_rp_1'] > %2.4f) * (cat_res['rgpp_rp_1'] < %2.4f) * (cat_tru['sf_hlr'] > %f) * (cat_tru['snr'] > %f) " % (list_size_edges[isize-1],list_size_edges[isize],size_cut,snr_cut)

            bias_result, all_res, all_tru = get_mc(selection_string,tag)
            list_bias_result.append(bias_result)

        bias_results = np.concatenate(list_bias_result)
        m_mean = (bias_results['m1'] + bias_results['m2'])/2.
        m_mean_std = np.sqrt((bias_results['m1_std']**2 + bias_results['m2_std']**2))/np.sqrt(2.)

        pl.figure(figsize=figsize)
        pl.errorbar(list_size_centers,bias_results['m1'],yerr=bias_results['m1_std'],fmt='r',label='m1')
        pl.errorbar(list_size_centers,bias_results['m2'],yerr=bias_results['m2_std'],fmt='b',label='m2')
        # pl.errorbar(list_size_centers,m_mean,m_mean_std,fmt='k',label='<m>')
        pl.legend()
        pl.xlabel('size')
        pl.ylabel('multiplicative bias')
        plotstools.adjust_limits()
        plot_add_requirements()
        filename_fig = 'figs/bias_vs_fwhmratio.png' 
        pl.savefig(filename_fig)
        log.info('saved %s', filename_fig)

    pl.show()


def show_example_galaxy_images():

    selection_string = "select= ( (cat_tru['fwhm']*0.27/cat_tru['psf_fwhm']) > %2.4f) * ((cat_tru['fwhm']*0.27/cat_tru['psf_fwhm']) < %2.4f) * (cat_res['snr']>40) " % (1.9 , 2)
    cols_res=['snr']
    cols_tru=['id','fwhm','psf_fwhm']  
    all_res,all_tru = get_selection_split(selection_string,cols_res,cols_tru)
    filename_meds1 = 'data/nbc2.meds.%03d.g%02d.noisefree.fits' % (0.,0.)
    filename_meds2 = 'data/nbc2.meds.%03d.g%02d.fits' % (0.,0.)
    id_gal = all_tru[0]['id'][2]
    import meds
    mm1=meds.MEDS(filename_meds1)
    mm2=meds.MEDS(filename_meds2)
    img1=mm1.get_cutout(id_gal,0)
    img2=mm2.get_cutout(id_gal,0)
    pl.subplot(1,2,1)
    pl.imshow(img1,interpolation='nearest'); 
    pl.subplot(1,2,2)
    pl.imshow(img2,interpolation='nearest'); 
    pl.suptitle(str(id_gal)+ ' in ' + filename_meds2 + '\n' + selection_string)
    pl.show()
    import pdb;pdb.set_trace()

    pass    

def kl_match():

    selection_string_sim = "select =  (cat_res['rgpp_rp_1']>-10) * (cat_tru['hsm_obs_sigma'] > 0) * (cat_tru['snr'] > 0)"
    selection_string_des = "select =  cat_res['rgpp_rp']>-10"

    print 'getting GREAT-DES selection'
    cat_res_all,cat_tru_all=get_selection(selection_string_sim, cols_tru=['fwhm','psf_fwhm', 'hsm_obs_sigma' ,'g1_true' , 'g2_true', 'snr' ], cols_res=['e1','e2','rgpp_rp_1','snr'] )
    print 'getting DES selection'
    cat_des_all = get_selection_DES(selection_string_des,cols=['e1','e2','rgpp_rp','snr'],n_files=200)

    n_grid=20
    grid_size = np.linspace(2,3,n_grid)
    grid_snr = np.linspace(5,10,n_grid)
    x,y = np.meshgrid(grid_size,grid_snr)
    kl_array = np.zeros([n_grid,n_grid])
    for isize,vsize in enumerate(grid_size):
        for isnr,vsnr in enumerate(grid_snr):

            select = (cat_tru_all['hsm_obs_sigma'] > vsize) * (cat_tru_all['snr'] > vsnr)
            cat_tru = cat_tru_all[select]
            cat_res = cat_res_all[select]
            cat_des = cat_des_all



            n_samples = 300000
            perm = np.random.permutation(len(cat_res))[:n_samples]
            cat_res = cat_res[perm,:]
            perm = np.random.permutation(len(cat_des))[:n_samples]
            cat_des = cat_des[perm,:]

            p = np.array(cat_des.tolist())
            q = np.array(cat_res.tolist())

            p  = p[~np.isnan(p).any(axis=1)]
            q  = q[~np.isnan(q).any(axis=1)]
            p  = p[~(p==0).any(axis=1)]
            q  = q[~(q==0).any(axis=1)]


            n_samples = len(p)
            perm = np.random.permutation(p.shape[0])[:n_samples]
            p = p[perm,:]
            perm = np.random.permutation(q.shape[0])[:n_samples]
            q = q[perm,:]

            import mathstools
            print 'getting kl' , isnr, isize
            kl = mathstools.get_kl_divergence_from_samples(p,q)
            kl_array[isize,isnr] = kl

            print vsize,vsnr,kl

    save_dict = {}
    save_dict['kl_array'] = kl_array
    save_dict['grid_size'] = grid_size
    save_dict['grid_snr'] = grid_snr
    import cPickle as pickle
    pickle.dump(save_dict,open('kl_array.pp2','w'),protocol=2)



def get_distributions():

    # selection_string_sim = "select =  (cat_res['rgpp_rp_1']>-10) * ((cat_tru['fwhm']*0.27/cat_tru['psf_fwhm']) > 1.37) * (cat_tru['snr'] > 4)"
    # selection_string_sim = "select =  (cat_res['rgpp_rp_1']>-10) * (cat_tru['hsm_obs_sigma']*2.355*0.27/cat_tru['psf_fwhm'] > 1.35) * (cat_tru['snr'] > 7.5)"
    selection_string_sim = "select =  (cat_res['rgpp_rp_1']>-10) * (cat_tru['sf_hlr'] > 0.65) * (cat_tru['snr'] > 10) "
    selection_string_des = "select =  cat_res['rgpp_rp']>-10"

    cat_res_all,cat_tru_all=get_selection(selection_string_sim, cols_tru=['fwhm','psf_fwhm', 'hsm_obs_sigma' ,'g1_true' , 'g2_true' ], cols_res=['e1','e2','gal_fwhm_1','snr', 'rgpp_rp_1' , 'radius' , 'flag'] )
    select_flag =  (cat_res_all['flag']&1==0)*(cat_res_all['flag']&4==0)*(cat_res_all['flag']&8==0)*(cat_res_all['flag']&64==0)*(cat_res_all['flag']&128==0)*(cat_res_all['flag']&256==0)
    # select_flag = True
    select_clean = (cat_res_all['snr'] > 0) * (cat_res_all['rgpp_rp_1'] > 0)
    cat_res=cat_res_all[select_flag*select_clean]
    cat_tru=cat_tru_all[select_flag*select_clean]
    print 'getting GREAT-DES selection n=', len(cat_tru_all), 'after_cuts' , len(cat_tru)

    if args.method!='im3shape': raise Exception('this function works only with method=im3shape')
    cat_des_all = get_selection_DES(selection_string_des,cols=['e1','e2','rgpp_rp','snr','radius','flag'],n_files=200)
    # select_flag = True
    select_flag =  (cat_des_all['flag']&1==0)*(cat_des_all['flag']&4==0)*(cat_des_all['flag']&8==0)*(cat_des_all['flag']&64==0)*(cat_des_all['flag']&128==0)*(cat_des_all['flag']&256==0)
    select_clean = (cat_des_all['snr'] > 0) * (cat_des_all['rgpp_rp'] > 0)
    cat_des=cat_des_all[select_flag*select_clean]
    print 'got DES galaxies n=' , len(cat_des_all) , ' , after cuts n=' , len(cat_des)

    great_des_e1 = cat_res['e1'] # - cat_tru['g1_true']
    great_des_e2 = cat_res['e2'] # - cat_tru['g2_true']

    cat_res_all,cat_tru_all=get_selection("select=(cat_res['rgpp_rp_1']>-10)", cols_tru=['hsm_obs_sigma' ], cols_res=['rgpp_rp_1'] )
    print 'selected galaxies:' , len(cat_res) , ' all_galaxies: ', len(cat_res_all)

    # cat_res = cat_res[~np.isinf(cat_res['rgpp_rp_1'])]
    # cat_res = cat_res[~np.isnan(cat_res['rgpp_rp_1'])]
    print selection_string_sim

    pl.figure(figsize=figsize)
    pl.subplot(2,2,1)
    pl.hist(great_des_e1, bins=np.linspace(-1,1,100),histtype='step',normed=True , label='GREAT-DES e1'      , color='r')
    pl.hist(great_des_e2, bins=np.linspace(-1,1,100),histtype='step',normed=True , label='GREAT-DES e2'      , color='m')
    pl.hist(cat_des['e1'], bins=np.linspace(-1,1,100),histtype='step',normed=True , label='im3shape-011-3 e1' , color='b')
    pl.hist(cat_des['e2'], bins=np.linspace(-1,1,100),histtype='step',normed=True , label='im3shape-011-3 e2' , color='c')
    pl.legend(mode='expand',ncol=2)
    ylim=list(pl.ylim()); ylim[1]*=1.5; pl.ylim(ylim)

    pl.subplot(2,2,2)
    pl.hist(cat_res['rgpp_rp_1'] ,bins=np.linspace(0,2,100),histtype='step',label='GREAT-DES rgpp_rp'      , normed=True, color='r') 
    pl.hist(cat_des['rgpp_rp'] ,bins=np.linspace(0,2,100),histtype='step',label='im3shape-011-4-r rgpp_rp' , normed=True, color='b') 
    pl.legend()
    ylim=list(pl.ylim()); ylim[1]*=1.5; pl.ylim(ylim)
    
    pl.subplot(2,2,3)
    pl.hist(cat_res['snr'] ,bins=np.linspace(1,100,100),histtype='step',label='GREAT-DES snr'      , normed=True, color='r') 
    pl.hist(cat_des['snr'] ,bins=np.linspace(1,100,100),histtype='step',label='im3shape-011-4-r snr' , normed=True, color='b') 
    pl.legend()

    pl.subplot(2,2,4)
    pl.hist(cat_res['radius'] ,bins=np.linspace(0,5,100),histtype='step',label='GREAT-DES radius'      , normed=True, color='r') 
    pl.hist(cat_des['radius'] ,bins=np.linspace(0,5,100),histtype='step',label='im3shape-011-4-r radius' , normed=True, color='b') 
    pl.legend()

    filename_fig = 'figs/histograms.match.png' 
    pl.savefig(filename_fig)
    log.info('saved %s', filename_fig)
    pl.show()

    # pl.figure()
    # select = (cat_des['snr'] < 100) * (cat_des['rgpp_rp'] < 2) * (cat_des['rgpp_rp'] > 0.5)
    # plotstools.plot_dist( np.concatenate( [cat_des['snr'][select][:,None], cat_des['rgpp_rp'][select][:,None]   ] , axis=1 ) , use_fraction=0.01   , contour=True, colormesh=False, color='b')
    # select = (cat_res['snr'] < 100) * (cat_res['rgpp_rp_1'] < 2) * (cat_res['rgpp_rp_1'] > 0.5)
    # plotstools.plot_dist( np.concatenate( [cat_res['snr'][select][:,None], cat_res['rgpp_rp_1'][select][:,None]   ] , axis=1 ) , use_fraction=0.01 , contour=True, colormesh=False, color='r')

    # plotstools.plot_dist( np.concatenate( [cat_res['snr'][:,None], cat_res['rgpp_rp_1'][:,None] ] , axis=1 ) )


    # pl.subplot(2,2,4)
    # pl.hist(cat_tru['hsm_obs_sigma'] ,bins=np.linspace(-2,10),histtype='step') 

    import pdb; pdb.set_trace()

def apply_calibration():

    size_cut = 0.65
    snr_cut = 10

    bins_snr = [10,15,20,25,30,40,50,60]
    bins_snr_edges = plotstools.get_bins_edges(bins_snr)
    list_bias_result = []
    for isnr in range(1,len(bins_snr_edges)):

        tag = 'snr=%.2f' % bins_snr[isnr-1]
        snr = bins_snr_edges[isnr]
        selection_string = "select=  (cat_res['snr'] > %2.4f) * (cat_res['snr'] < %f) * (cat_tru['sf_hlr'] > %f) * (cat_tru['snr'] > %f) " % (
                                            bins_snr_edges[isnr-1],
                                            bins_snr_edges[isnr],
                                            size_cut,snr_cut)

        log.info(selection_string)
        bias_result, all_res, all_tru = get_mc(selection_string,tag)
        list_bias_result.append(bias_result)

    bias_results = np.concatenate(list_bias_result)
    m_mean = (bias_results['m1'] + bias_results['m2'])/2.
    m_mean_std = np.sqrt((bias_results['m1_std']**2 + bias_results['m2_std']**2))/np.sqrt(2.)
    m_mean= np.append(m_mean,0)
    m_mean_std= np.append(m_mean_std,0)
    m_mean= np.insert(m_mean,0,0)
    m_mean_std= np.insert(m_mean_std,0,0)
    # print 'bins_snr' , bins_snr
    # print 'm_mean' , m_mean

    filelist = np.loadtxt('filelist_r.txt',dtype='a1024')
    n_files = len(filelist)
    list_results = []
    for filename_des in filelist[:n_files]:
        cat_res=tabletools.loadTable(filename_des,log=1)

        bin_membership = np.digitize(cat_res['snr'],bins=bins_snr_edges)
        column_m = m_mean[bin_membership]
        column_m_err = m_mean_std[bin_membership]

        select_flag =  (cat_res['flag']&1==0)*(cat_res['flag']&4==0)*(cat_res['flag']&8==0)*(cat_res['flag']&64==0)*(cat_res['flag']&128==0)*(cat_res['flag']&256==0)
        select_clean = (cat_res['snr'] > 0) * (cat_res['rgpp_rp'] > 0)

        if 'nbc_m' not in cat_res.dtype.names:
            cat_res = tabletools.appendColumn(rec=cat_res,arr=column_m,dtype='f8',name='nbc_m')
            cat_res = tabletools.appendColumn(rec=cat_res,arr=column_m_err,dtype='f8',name='nbc_m_err')
        else:
            cat_res['nbc_m'] = column_m
            cat_res['nbc_m_err'] = column_m_err

        # pl.plot(cat_res['snr'],cat_res['nbc_m'],'x')
        # pl.plot(bins_snr,m_mean[1:-1],'d')
        # pl.show()

        filename_des_calibrated = filename_des.replace('fits.gz','nbc.fits') 
        filename_des_calibrated = filename_des_calibrated.replace('/im3shape-011-4/','/im3shape-011-4-nbc/') 
        tabletools.saveTable(filename_des_calibrated,cat_res)



  

def main():

    global log , config , args

    description = 'Get statistics and plot results of noise bias calibration runs'
    parser = argparse.ArgumentParser(description=description, add_help=True)
    parser.add_argument('-v', '--verbosity', type=int, action='store', default=2, choices=(0, 1, 2, 3 ), help='integer verbosity level: min=0, max=3 [default=2]')
    parser.add_argument('-c', '--filename_config', default='nbc2_sva1.yaml',type=str, action='store', help='name of the yaml config file')
    parser.add_argument('-m', '--method', default='hsm',type=str, action='store', help='name of the yaml config file')
    parser.add_argument('-f', '--first', type=int, action='store', default=0, help='index of the first file to analyse')
    parser.add_argument('-n', '--num', type=int, action='store', default=-1, help='number of files to analyse, if -1 then =config[n_files]')

    # parser.add_argument('-a', '--actions', nargs='+' , type=str, action='store',  help='')
    
    args = parser.parse_args()
    # Parse the integer verbosity level from the command line args into a logging_level string
    logging_levels = { 0: logging.CRITICAL, 
                       1: logging.WARNING,
                       2: logging.INFO,
                       3: logging.DEBUG }
    logging_level = logging_levels[args.verbosity]
    

    log.info(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))

    config = yaml.load(open(args.filename_config))
    # get_distributions()
    # get_plots()
    # kl_match()
    # show_example_galaxy_images()
    apply_calibration()

 
    log.info(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))


if __name__ == '__main__':
    main()