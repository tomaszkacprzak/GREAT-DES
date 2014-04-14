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
log.addHandler(stream_handler)
log.propagate = False
       
req1_dg = 0.003
req2_dg = 0.02
req3_dg = 0.1
req1_nfail = 0.01
req2_nfail = 0.05

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

def get_shear_estimator(res):

    stats = np.zeros(1,dtype=dtype_stats)

    n_res = len(res)
    col_g1 = config['methods'][args.method]['col_g1']
    col_g2 = config['methods'][args.method]['col_g2']
    col_size = config['methods'][args.method]['col_size']
    flip_g1 =  config['methods'][args.method]['flip_g1']
    flip_g2 =  config['methods'][args.method]['flip_g2']
    select_successful = ((np.abs(res[col_g1] + 1j*res[col_g2]))<1) * (~np.isnan(res[col_size])) * (~np.isinf(res[col_size]))
    n_success = len(np.nonzero(select_successful)[0])
    n_fail = n_res - n_success
    res_success = res[select_successful]

    
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

    if stats['n_gals'] < 2: raise Exception('only one galaxy in sample')

    # check nans and infs
    for key in stats.dtype.names:
        if any(np.isinf(stats[key])) or any(np.isnan(stats[key])):
            import pdb; pdb.set_trace()


    return stats


def get_selection_DES(selection_string,cols,n_files=20):

    filelist_i = np.loadtxt('filelist_i.txt',dtype='a1024')
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



def get_selection_split(selection_string, cols_res, cols_tru):

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
                cat_tru = tabletools.loadTable(filename_tru,log=1)
            except:
                log.error('file %s not found' % filename_tru )
                continue

            for col in cols_tru:
                if col not in cat_tru.dtype.names:
                    raise Exception('column %s not found in truth catalog %s' % (col,filename_tru))

            try:
                cat_res = tabletools.loadTable(filename_res,log=1)
            except:
                log.error('file %s not found' % filename_res )
                continue
 
            if 'patch_g2' in config['methods'][args.method]:
                if config['methods'][args.method]['patch_g2']:
                    cat_res = tabletools.appendColumn(rec=cat_res, arr=cat_res['hsm_cor_g1'].copy(), name='hsm_cor_g2', dtype='f4')          
 

            for col in cols_res:
                if col not in cat_res.dtype.names:
                    raise Exception('column %s not found in results catalog %s' % (col,filename_res))


            exec selection_string
            if len(np.nonzero(select)[0]) < 1:
                raise Exception('select didnt give any results %s' % selection_string)

            selected_res = cat_res[select][cols_res]
            selected_tru = cat_tru[select][cols_tru]


            # selected_res = selected_res.astype( ['f4']*len(cols_res) )
            selected_tru = selected_tru.astype( dtype={ 'formats':['f4']*len(cols_tru) , 'names': cols_tru })

            list_all_res.append(selected_res)
            list_all_tru.append(selected_tru)
   
    n_total = 0
    for res in list_all_res:
        n_total += len(res)

        
    log.info('selected %d parts with average %d, total %d' % (len(list_all_res) , float(n_total)/float(len(list_all_res)), n_total) )
    return list_all_res, list_all_tru



def get_shear_results(selection_string):

    # these are the columns that we select for the results catalog (we are interested only in few, don't need all columns)
    col_g1 = config['methods'][args.method]['col_g1']
    col_g2 = config['methods'][args.method]['col_g2']
    col_size = config['methods'][args.method]['col_size']
    col_snr = config['methods'][args.method]['col_snr']
    cols_res = [col_g1, col_g2, col_size , col_snr]
    cols_tru = ['g1_true', 'g2_true' , 'snr','flux']

    list_all_res, list_all_tru = get_selection_split(selection_string, cols_res, cols_tru)
    n_split = len(list_all_tru)
    list_shears = []

    for il in range(n_split):

        try:
            stats=get_shear_estimator(list_all_res[il])
        except:
            log.error('not enough galaxies in sample %d' % il)
            continue
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
    
    stats = get_shear_results(selection_string)

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
    # d1 = np.std(g1_bias - g1_tru*m1,ddof=1)
    # d2 = np.std(g2_bias - g2_tru*m2,ddof=1)
    # [c1,m1,C1cm] = get_line_fit(g1_tru,g1_bias,np.ones_like(g1_tru)*d1)
    # [c2,m2,C2cm] = get_line_fit(g2_tru,g2_bias,np.ones_like(g2_tru)*d2)
    m1_std = np.sqrt(C1cm[1,1])
    m2_std = np.sqrt(C2cm[1,1])
    c1_std = np.sqrt(C1cm[0,0])
    c2_std = np.sqrt(C2cm[0,0])


    g_std = np.mean(stats['est_stdv_g1'])
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
    bias_result['g_std'] = g_std
    bias_result['n_fail'] = n_fail

    global fig_counter
    pl.figure()
    pl.errorbar(g1_tru,g1_bias,yerr=g1_err,fmt='b.',label='m1=%2.3f +/- %2.3f' % (m1,m1_std))
    pl.plot(g1_tru,g1_tru*m1 + c1,'b-')
    pl.errorbar(g2_tru,g2_bias,yerr=g2_err,fmt='r.',label='m1=%2.3f +/- %2.3f' % (m2,m2_std))
    pl.plot(g2_tru,g2_tru*m2 + c2,'r-')
    pl.title(tag)
    pl.legend()
    # pl.show()
    filename_fig = 'figs/fig.bias.%s.png' % (tag)
    pl.savefig(filename_fig)
    pl.close()
    log.info('saved %s' % filename_fig)
 
    # exec('func_get_shear='+config['methods'][args.method]['func_shear'])

    return bias_result


def get_plots():

    # --------------------------------------------------------------------------------------------------
    # list_bias_result = []
    # list_psf = [0.8,0.9,1.1,1.2,1.3]
    # for ipsf,vpsf in enumerate(list_psf):

    #     tag = 'psf_fwhm=%2.2f' % vpsf
    #     selection_string = "select=cat_tru['psf_fwhm']==%1.1f" % vpsf

    #     bias_result = get_mc(selection_string,tag)
    #     list_bias_result.append(bias_result)
    #     print vpsf , bias_result['m1'], bias_result['m1_std'] , bias_result['m2'], bias_result['m2_std']

    # bias_results = np.concatenate(list_bias_result)

    # pl.figure()
    # pl.errorbar(list_psf,bias_results['m1'],yerr=bias_results['m1_std'],fmt='r',label='m1')
    # pl.errorbar(list_psf,bias_results['m2'],yerr=bias_results['m2_std'],fmt='b',label='m2')
    # pl.legend()
    # pl.xlabel('PSF_FWHM')
    # pl.ylabel('multiplicative bias')
    # plotstools.adjust_limits()
    # filename_fig = 'figs/bias_vs_psf.png' 
    # pl.savefig(filename_fig)
    # log.info('saved %s', filename_fig)

    # --------------------------------------------------------------------------------------------------
    # list_snr_centers = [5,10,15,20,25,30,40,50,60,70,80,90,100]
    # list_snr_edges = [2.5,7.5,12.5,17.5,22.5,27.5,35,45,55,65,75,85,95,100]
    # list_bias_result = []
    # for isnr in range(1,len(list_snr_edges)):

    #     tag = 'snr=%.2f' % list_snr_centers[isnr-1]
    #     snr = list_snr_edges[isnr]
    #     selection_string = "select= ( cat_res['snr'] > %2.4f) * (cat_res['snr'] < %2.4f) " % (list_snr_edges[isnr-1],list_snr_edges[isnr])

    #     bias_result = get_mc(selection_string,tag)
    #     print list_snr_centers[isnr-1], bias_result['m1'], bias_result['m1_std']
    #     list_bias_result.append(bias_result)

    # bias_results = np.concatenate(list_bias_result)

    # pl.figure()
    # pl.errorbar(list_snr_centers,bias_results['m1'],yerr=bias_results['m1_std'],fmt='r',label='m1')
    # pl.errorbar(list_snr_centers,bias_results['m2'],yerr=bias_results['m2_std'],fmt='b',label='m2')
    # pl.legend()
    # pl.xlabel('SNR')
    # pl.ylabel('multiplicative bias')
    # plotstools.adjust_limits()
    # filename_fig = 'figs/bias_vs_snr.png' 
    # pl.savefig(filename_fig)
    # log.info('saved %s', filename_fig)
    # pl.show()

    # --------------------------------------------------------------------------------------------------

    # selection_string = "select = cat_tru['fwhm']*0.27/cat_tru['psf_fwhm']>1.2"


    list_size_centers = np.arange(1.5,4,0.25)
    list_size_edges = plotstools.get_bins_edges(list_size_centers)
    list_bias_result = []
    for isize in range(1,len(list_size_edges)):

        tag = 'size=%.2f' % list_size_centers[isize-1]
        size = list_size_edges[isize]
        selection_string = "select= ( cat_tru['hsm_obs_sigma'] > %2.4f) * (cat_tru['hsm_obs_sigma'] < %2.4f) * (cat_res['snr']>40) " % (list_size_edges[isize-1],list_size_edges[isize])

        bias_result = get_mc(selection_string,tag)
        print list_size_centers[isize-1], bias_result['m1'], bias_result['m1_std'] , bias_result['m2'], bias_result['m2_std']
        list_bias_result.append(bias_result)

    bias_results = np.concatenate(list_bias_result)

    pl.figure()
    pl.errorbar(list_size_centers,bias_results['m1'],yerr=bias_results['m1_std'],fmt='r',label='m1')
    pl.errorbar(list_size_centers,bias_results['m2'],yerr=bias_results['m2_std'],fmt='b',label='m2')
    pl.legend()
    pl.xlabel('size')
    pl.ylabel('multiplicative bias')
    plotstools.adjust_limits()
    filename_fig = 'figs/bias_vs_size.png' 
    pl.savefig(filename_fig)
    log.info('saved %s', filename_fig)
    pl.show()

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

    selection_string_sim = "select =  (cat_res['rgpp_rp_1']>-10) * (cat_tru['hsm_obs_sigma'] > 2.3) * (cat_tru['snr'] > 8)"
    selection_string_des = "select =  cat_res['rgpp_rp']>-10"

    print 'getting GREAT-DES selection'
    cat_res,cat_tru=get_selection(selection_string_sim, cols_tru=['fwhm','psf_fwhm', 'hsm_obs_sigma' ,'g1_true' , 'g2_true' ], cols_res=['e1','e2','gal_fwhm_1','snr', 'rgpp_rp_1'] )
    print 'getting DES selection'
    cat_des = get_selection_DES(selection_string_des,cols=['e1','e2','rgpp_rp','snr'],n_files=200)
    great_des_e1 = cat_res['e1'] #- cat_tru['g1_true']
    great_des_e2 = cat_res['e2'] #- cat_tru['g2_true']

    # cat_res = cat_res[~np.isinf(cat_res['rgpp_rp_1'])]
    # cat_res = cat_res[~np.isnan(cat_res['rgpp_rp_1'])]
    print cat_res.shape
    pl.figure()
    pl.subplot(2,2,1)
    pl.hist(great_des_e1, bins=np.linspace(-1,1,100),histtype='step',normed=True , label='GREAT-DES e1'      , color='r')
    pl.hist(great_des_e2, bins=np.linspace(-1,1,100),histtype='step',normed=True , label='GREAT-DES e2'      , color='m')
    pl.hist(cat_des['e1'], bins=np.linspace(-1,1,100),histtype='step',normed=True , label='im3shape-011-3 e1' , color='b')
    pl.hist(cat_des['e2'], bins=np.linspace(-1,1,100),histtype='step',normed=True , label='im3shape-011-3 e2' , color='c')
    pl.legend()
    pl.subplot(2,2,2)
    pl.hist(cat_res['rgpp_rp_1'] ,bins=np.linspace(0,2,100),histtype='step',label='GREAT-DES rgpp_rp'      , normed=True, color='r') 
    pl.hist(cat_des['rgpp_rp'] ,bins=np.linspace(0,2,100),histtype='step',label='im3shape-011-3 rgpp_rp' , normed=True, color='b') 
    pl.legend()
    pl.subplot(2,2,3)
    pl.hist(cat_res['snr'] ,bins=np.linspace(1,100,100),histtype='step',label='GREAT-DES snr'      , normed=True, color='r') 
    pl.hist(cat_des['snr'] ,bins=np.linspace(1,100,100),histtype='step',label='im3shape-011-3 snr' , normed=True, color='b') 
    pl.legend()
    pl.subplot(2,2,4)
    pl.hist(cat_tru['hsm_obs_sigma'] ,bins=np.linspace(-2,10),histtype='step') 
   
    pl.show()

    import pdb; pdb.set_trace()

  

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
    kl_match()


    log.info(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))


if __name__ == '__main__':
    main()