import numpy as np; import pylab as pl; import tabletools as tt
import  sys, logging, yaml, argparse, time, copy, itertools, fitting, warnings
warnings.simplefilter("once")
sys.path.append('/home/tomek/code/tktools')
import tabletools, plotstools
# from nbc2_dtypes import *
logging_level = logging.INFO; logger = logging.getLogger("nbc-v7"); logger.setLevel(logging_level)  
log_formatter = logging.Formatter("%(asctime)s %(name)s %(levelname)s   %(message)s", "%Y-%m-%d %H:%M:%S")
stream_handler = logging.StreamHandler(sys.stdout); stream_handler.setFormatter(log_formatter)
if logger.handlers == [] : logger.addHandler(stream_handler); logger.propagate = False
import nbc_v7_select

# info_vals = [1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536,131072,262144,524288,1048576]
info_vals = [1,2,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,65536,131072,262144,524288,1048576]
select_info = [ "((cat_res['info_flag']&%d)==0)&"%iv for iv in info_vals ]
select_info = (' '.join(select_info))[:-1]
# selection_string_sim = "select =  (cat_res['snr']>-1000000000) & (cat_tru['sf_hlr']>0.6) & (cat_tru['snr_mean_all'] > 9) "
selection_string_sim = "select =  (cat_res['error_flag']==0) & ( ( (cat_res['radius']<3.5) & (cat_res['stamp_size']==48) ) | (cat_res['stamp_size']!=48)  ) & (cat_tru['sf_hlr']> 0.2) & (cat_tru['cosmos_mag_auto'] < 23.2) & " + select_info 
selection_string_des = "select =  (cat_res['error_flag']==0) & ( ( (cat_res['radius']<3.5) & (cat_res['stamp_size']==48) ) | (cat_res['stamp_size']!=48)  ) & " + select_info 

def inv_snr_basis(x):

    n_points_x = x.shape[0]
    # X = np.concatenate( [ np.ones((n_points_x,1)), 1./x, 1./x**2, 1./x**4],axis=1 )
    X = np.concatenate( [ 1./x, 1./x**2, 1./x**4],axis=1 )
    return X
    


def plot_calibration():

    import fitsio
    filename_table_bias = 'bias_table.fits'
    bias_table = fitsio.read(filename_table_bias)
    dx=0.2

    colorscale=plotstools.get_colorscale(len(np.unique(bias_table['ipsf'])))


    for ipsf in np.unique(bias_table['ipsf']):

        select = bias_table['ipsf'] == ipsf
        bt1=bias_table[select]
        label = 'FWHM_RATIO=%2.1f-%2.1f'%(bt1['vpsf_min'][0],bt1['vpsf_max'][0])
        snr_mid = (bt1['vsnr_min']+bt1['vsnr_max'])/2.
        pl.errorbar(snr_mid+dx*ipsf,bt1['m']-1,yerr=bt1['std_m'],label=label,fmt='.',c=colorscale[ipsf-1])

        w, w_cov =fitting.fit(snr_mid, bt1['m']-1, s=bt1['std_m'],expand=inv_snr_basis)
        snr_pred =  np.linspace(snr_mid.min(),snr_mid.max(),1000)
        p, s = fitting.predict(snr_pred,w,w_cov,expand=inv_snr_basis)
        pl.plot(snr_pred,p,c=colorscale[ipsf-1])



    pl.axhline(1,c='k')
    pl.legend(mode='expand',ncol=2)
    ylim=list(pl.ylim()); ylim[1]*=1.; pl.ylim(ylim)


    pl.figure()
    for ipsf in np.unique(bias_table['ipsf']):

        select = bias_table['ipsf'] == ipsf
        bt1=bias_table[select]
        label = 'FWHM_RATIO=%2.1f-%2.1f'%(bt1['vpsf_min'][0],bt1['vpsf_max'][0])
        snr_mid = (bt1['vsnr_min']+bt1['vsnr_max'])/2.
        pl.errorbar(snr_mid+dx*ipsf,bt1['c1'],yerr=bt1['std_c1'],label=label,fmt='.')

    
    pl.legend()
    pl.show()


    import pdb; pdb.set_trace()

def plot_bias_vs_redshift():

    selection_string_sim = "select =  (cat_res['snr']>-1000000000) & (cat_tru['sf_hlr']>0.6) "
    selection_string_des = "select =  (cat_res['error_flag']==0) & ( ( (cat_res['radius']<3.5) & (cat_res['stamp_size']==48) ) | (cat_res['stamp_size']!=48)  ) & " + select_info 

    # res_sim,res_tru = nbc_v7_select.get_selection_sim(selection_string_sim,cols_res=['coadd_objects_id','ra_as','dec_as','e1','e2','snr','disc_A','bulge_A','mean_rgpp_rp','radius'],cols_tru=['snr','psf_e1','psf_e2','id_shear','cosmos_mag_auto','g1_true','g2_true'])





def get_mc(res_sim,res_tru,vpsf=0,vsnr=0):

    n_jack = 1
    mean_e1 = []; stdv_e1 = []; stdm_e1 = []; mean_e2 = []; stdv_e2 = []; stdm_e2 = []; true_e1 = []; true_e2 = [];

    for ig,vg in enumerate(config['shear']):

        select = res_tru['id_shear'] == ig
        sim_shear = res_sim[select]
        tru_shear = res_tru[select]

        for ij in range(n_jack):

            n_res = len(sim_shear)
            select = range(ij,n_res,n_jack)
            res_sim_select = sim_shear[select]
            tru_sim_select = tru_shear[select]
            n_select = len(res_sim_select)
            
            true_e1 += [ config['shear'][ig][0] ]
            mean_e1 += [ np.mean(res_sim_select['e1']) ]
            stdv_e1 += [ np.std(res_sim_select['e1'],ddof=1)  ]
            stdm_e1 += [ np.std(res_sim_select['e1'],ddof=1)/np.sqrt(n_select)     ]
            
            true_e2 += [ config['shear'][ig][1] ]
            mean_e2 += [ np.mean(res_sim_select['e2']) ]
            stdv_e2 += [ np.std(res_sim_select['e2'],ddof=1)  ]
            stdm_e2 += [ np.std(res_sim_select['e2'],ddof=1)/np.sqrt(n_select)     ]

    mean_e1 = np.array(mean_e1); stdv_e1 = np.array(stdv_e1); stdm_e1 = np.array(stdm_e1); mean_e2 = np.array(mean_e2); stdv_e2 = np.array(stdv_e2); stdm_e2 = np.array(stdm_e2); true_e1 = np.array(true_e1); true_e2 = np.array(true_e2);
    
    import fitting
    cc1,mm1,Ccm1=fitting.get_line_fit(true_e1,mean_e1,stdm_e1)
    cc2,mm2,Ccm2=fitting.get_line_fit(true_e2,mean_e2,stdm_e2)
    std_cc1 = np.sqrt(Ccm1[0,0])
    std_mm1 = np.sqrt(Ccm1[1,1])
    std_cc2 = np.sqrt(Ccm2[0,0])
    std_mm2 = np.sqrt(Ccm2[1,1])

    mm = (mm1+mm2)/2.
    std_mm = np.sqrt( (std_mm1**2 + std_mm2**2)/2. )

    cc = (cc1+cc2)/2.
    std_cc = np.sqrt( (std_cc1**2 + std_cc2**2)/2. )

    bins_psf_centers = np.unique(res_tru['psf_e1'])

    list_mean_e1 = []
    list_stdm_e1 = []
    list_mean_g1 = []
    
    list_mean_e2 = []
    list_stdm_e2 = []
    list_mean_g2 = []

    for ic in range(len(bins_psf_centers)):

        select = (res_tru['psf_e1'] > (bins_psf_centers[ic]-0.001)) & (res_tru['psf_e1'] < (bins_psf_centers[ic]+0.001))
        pmean_e1=np.mean(res_sim[select]['e1'])
        pstdm_e1=np.std(res_sim[select]['e1'],ddof=1)/np.sqrt(len(res_sim[select]))
        pmean_g1=np.mean(res_tru[select]['g1_true'])
        list_mean_e1.append(pmean_e1)
        list_stdm_e1.append(pstdm_e1)
        list_mean_g1.append(pmean_g1)

        select = (res_tru['psf_e2'] > (bins_psf_centers[ic]-0.001)) & (res_tru['psf_e2'] < (bins_psf_centers[ic]+0.001))
        pmean_e2=np.mean(res_sim[select]['e2'])
        pstdm_e2=np.std(res_sim[select]['e2'],ddof=1)/np.sqrt(len(res_sim[select]))
        pmean_g2=np.mean(res_tru[select]['g2_true'])
        list_mean_e2.append(pmean_e2)
        list_stdm_e2.append(pstdm_e2)
        list_mean_g2.append(pmean_g2)

    pmean_e1=np.array(list_mean_e1)
    pstdm_e1=np.array(list_stdm_e1)
    pmean_e2=np.array(list_mean_e2)
    pstdm_e2=np.array(list_stdm_e2)  
    pmean_g1=np.array(list_mean_g1)
    pmean_g2=np.array(list_mean_g2)
    
    pcc1,pmm1,pCcm1=fitting.get_line_fit(bins_psf_centers,pmean_e1,pstdm_e1)
    pcc2,pmm2,pCcm2=fitting.get_line_fit(bins_psf_centers,pmean_e2,pstdm_e2)
    std_pcc1 = np.sqrt(pCcm1[0,0])
    std_pmm1 = np.sqrt(pCcm1[1,1])
    std_pcc2 = np.sqrt(pCcm2[0,0])
    std_pmm2 = np.sqrt(pCcm2[1,1])

    print "--- m1 = % 2.5f +/- % 2.5f" % (mm1,std_mm1)
    print "--- m2 = % 2.5f +/- % 2.5f" % (mm2,std_mm2)
    print "--- c1 = % 2.5f +/- % 2.5f" % (cc1,std_cc1)
    print "--- c2 = % 2.5f +/- % 2.5f" % (cc2,std_cc2)
    print "--- m  = % 2.5f +/- % 2.5f" % (mm,std_mm)
    print "--- c  = % 2.5f +/- % 2.5f" % (cc,std_cc)
    print "--- pm1  = % 2.5f +/- % 2.5f" % (pmm1,std_pmm1)
    print "--- pc1  = % 2.5f +/- % 2.5f" % (pcc1,std_pcc1)
    print "--- pm2  = % 2.5f +/- % 2.5f" % (pmm2,std_pmm2)
    print "--- pc2  = % 2.5f +/- % 2.5f" % (pcc2,std_pcc2)

    pl.figure(figsize=(15,8))

    pl.subplot(2,2,1)
    pl.errorbar(true_e1,mean_e1-true_e1,yerr=stdm_e1,fmt='r.')
    pl.errorbar(true_e2,mean_e2-true_e2,yerr=stdm_e1,fmt='b.')
    pl.plot(true_e1,true_e1*0,'-k')
    pl.plot(true_e1,(mm1-1)*true_e1+cc1,'r-',label=r'sl=$% 2.4f \pm %2.4f$ int=$% 2.4f \pm % 2.4f$' % (mm1,std_mm1,cc1,std_cc1))
    pl.plot(true_e2,(mm2-1)*true_e2+cc2,'b-',label=r'sl=$% 2.4f \pm %2.4f$ int=$% 2.4f \pm % 2.4f$' % (mm2,std_mm2,cc2,std_cc2))
    pl.xlabel('g true')
    pl.ylabel('e true')
    ylim=list(pl.ylim()); ylim[1]*=2; pl.ylim(ylim)
    pl.legend(mode='expand')


    pl.subplot(2,2,2)
    pl.errorbar(bins_psf_centers,pmean_e1-pmean_g1,yerr=list_stdm_e1,fmt='rd') # ,label='e1 vs PSF e1'
    pl.errorbar(bins_psf_centers,pmean_e2-pmean_g2,yerr=list_stdm_e2,fmt='bd') # ,label='e2 vs PSF e2'
    pl.plot(bins_psf_centers,bins_psf_centers*pmm1+pcc1,'r-',label=r'sl=$% 2.4f \pm %2.4f$ int=$% 2.4f \pm %2.4f$' % (pmm1,std_pmm1,pcc1,std_pcc1))
    pl.plot(bins_psf_centers,bins_psf_centers*pmm2+pcc2,'b-',label=r'sl=$% 2.4f \pm %2.4f$ int=$% 2.4f \pm %2.4f$' % (pmm2,std_pmm2,pcc2,std_pcc2))
    # pl.plot(bins_psf_centers,pmean_g1,'ro',label='mean true g1')
    # pl.plot(bins_psf_centers,pmean_g2,'bo',label='mean true g2')
    pl.legend(mode='expand')
    pl.axhline(0,c='k')
    plotstools.adjust_limits()
    pl.xlabel('PSF e')
    pl.ylabel('<e>')
    pl.suptitle( 'SNR=%2.0f FWHM_RATIO=%2.1f' % (vsnr,vpsf) )
    ylim=list(pl.ylim()); ylim[1]*=2; pl.ylim(ylim)

    pl.subplot(2,2,3)
    pl.hist(res_sim['e1'],bins=np.linspace(-1,1,100),color='r',histtype='step')
    pl.hist(res_sim['e2'],bins=np.linspace(-1,1,100),color='m',histtype='step')

    pl.subplot(2,2,4)
    pl.hist(res_sim['mean_rgpp_rp'],bins=np.linspace(1,3),color='r',histtype='step')
    
    filename_fig = 'figs/bias.snr%2.2f.psf%2.2f.png' % (vsnr,vpsf)
    pl.savefig(filename_fig)
    logger.info('saved %s',filename_fig)
    pl.close()

    #TODO: make resampled error-bars

    return mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2

    




def get_calibration():

    selection_string_sim = "select =  (cat_res['snr']>-1000000000) & (cat_tru['sf_hlr']>0.6) "
    selection_string_des = "select =  (cat_res['error_flag']==0) & ( ( (cat_res['radius']<3.5) & (cat_res['stamp_size']==48) ) | (cat_res['stamp_size']!=48)  ) & " + select_info 

    res_sim,res_tru = nbc_v7_select.get_selection_sim(selection_string_sim,cols_res=['coadd_objects_id','ra_as','dec_as','e1','e2','snr','disc_A','bulge_A','mean_rgpp_rp','radius'],cols_tru=['snr','psf_e1','psf_e2','id_shear','cosmos_mag_auto','g1_true','g2_true'])
    
    # list_snr_centers = [5,10,15,20,25,30,40,50]
    list_snr_edges = [5,10,15,20,30,50,1000]
    list_snr_centers = plotstools.get_bins_centers(list_snr_edges)
    # list_psf_edges = [0.79,0.91,1.11,1.31]
    list_psf_edges = [1.0,1.2,1.4,1.6,2.0,100]
    list_psf_centers = plotstools.get_bins_centers(list_psf_edges)

    list_bias = []

    for ipsf in range(1,len(list_psf_edges)):
        for isnr in range(1,len(list_snr_edges)):

            select = (res_sim['snr'] > list_snr_edges[isnr-1]) & (res_sim['snr'] < list_snr_edges[isnr]) & (res_sim['mean_rgpp_rp'] > list_psf_edges[ipsf-1]) & (res_sim['mean_rgpp_rp'] < list_psf_edges[ipsf])
            print ipsf,isnr,len(np.nonzero(select)[0])
            res_sim_select = res_sim[select]
            res_tru_select = res_tru[select]

            vpsf_mid = list_psf_centers[ipsf-1]
            vsnr_mid = list_snr_centers[isnr-1]
            vpsf_min = list_psf_edges[ipsf-1]
            vsnr_min = list_snr_edges[isnr-1]
            vpsf_max = list_psf_edges[ipsf]
            vsnr_max = list_snr_edges[isnr]


            logger.info(selection_string_sim)
            mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2=get_mc(res_sim_select,res_tru_select,vpsf=vpsf_min,vsnr=vsnr_min)

            std_e = np.std(res_sim_select['e1'],ddof=1)
            list_bias.append( [ipsf,isnr,vpsf_min,vpsf_max,vsnr_min,vsnr_max,std_e,mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2] )


    arr_bias = tabletools.arr2rec(np.array(list_bias),dtype={'names': ["ipsf","isnr","vpsf_min","vpsf_max","vsnr_min","vsnr_max","std_e","m","std_m","c","std_c","m1","std_m1","m2","std_m2","c1","std_c1","c2","std_c2"], 'formats': ['i4']*2 + ['f8']*17 })

    filename_table_bias = 'bias_table.fits'
    import pyfits
    pyfits.writeto(filename_table_bias,arr_bias,clobber=True)
    logger.info('saved %s',filename_table_bias)

    import pdb; pdb.set_trace()



    #     print 'total number in all bins n=' , n_total
    #     bias_results = np.concatenate(list_bias_result)

    #     m_mean = (bias_results['m1'] + bias_results['m2'])/2.
    #     m_mean_std = np.sqrt((bias_results['m1_std']**2 + bias_results['m2_std']**2))/np.sqrt(2.)

    #     array_m[ipsf,:] = m_mean
    #     array_m_std[ipsf,:] = m_mean_std

    #     # pl.plot(list_snr_centers,bias_results['m1'],list_psf_colors[ipsf-1]+':',label='m1')
    #     # pl.plot(list_snr_centers,bias_results['m2'],list_psf_colors[ipsf-1]+'--',label='m2')
    #     pl.errorbar(list_snr_centers,m_mean,m_mean_std,label='PSF_FWHM in (%2.1f,%2.1f)'%(list_psf_edges[ipsf-1],list_psf_edges[ipsf]),fmt=list_psf_colors[ipsf-1]+'-' )

    # list_psf_colors = ['r' , 'g' , 'b']
    # pl.figure(figsize=figsize)
    # pl.legend()
    # pl.xlabel('SNR')
    # pl.ylabel('multiplicative shear bias')
    # # pl.title('PSF_FWHM in {%2.1f,%2.1f}'%(list_psf_edges[ipsf-1],list_psf_edges[ipsf]))
    # plotstools.adjust_limits()
    # plot_add_requirements()
    # pl.title(args.method)
    # filename_fig = 'figs/bias_vs_snr_vs_psf.png'
    # pl.savefig(filename_fig)
    # log.info('saved %s', filename_fig)

    # bias_result = { 'm_mean' : array_m, 'm_mean_std' : array_m_std}
    # filename_pickle = 'bias_m.pp2'
    # tabletools.savePickle(filename_pickle,bias_results)

def get_distributions():

    res_sim,res_tru = nbc_v7_select.get_selection_sim(selection_string_sim,cols_res=['coadd_objects_id','ra_as','dec_as','e1','e2','snr','disc_A','bulge_A','mean_rgpp_rp','radius'],cols_tru=['id','snr','psf_e1','psf_e2','cosmos_mag_auto'])
    # res_sim,res_tru = nbc_v7_select.get_selection_sim(selection_string_sim,cols_res=['ra_as','dec_as','e1','e2','snr','disc_A','bulge_A','rgpp_rp_1','radius'],cols_tru=['snr','psf_e1','psf_e2'])
    res_des = nbc_v7_select.get_selection_des(selection_string_des,cols=['ra_as','dec_as','e1','e2','snr','mean_psf_e1_sky','mean_psf_e2_sky','disc_A','bulge_A','mean_rgpp_rp','radius'],n_files=args.num)

    great_des_e1 = res_sim['e1'] #- cat_tru['g1_true']
    great_des_e2 = res_sim['e2'] #- cat_tru['g2_true']

    # res_sim = res_sim[~np.isinf(res_sim['rgpp_rp_1'])]
    # res_sim = res_sim[~np.isnan(res_sim['rgpp_rp_1'])]
    print selection_string_sim

    pl.figure(figsize=(15,5))
    pl.subplot(2,2,1)
    pl.hist(great_des_e1, bins=np.linspace(-1,1,100),histtype='step',normed=True , label='GREAT-DES e1'      , color='r')
    pl.hist(great_des_e2, bins=np.linspace(-1,1,100),histtype='step',normed=True , label='GREAT-DES e2'      , color='m')
    pl.hist(res_des['e1'], bins=np.linspace(-1,1,100),histtype='step',normed=True , label='im3shape-v7-r e1' , color='b')
    pl.hist(res_des['e2'], bins=np.linspace(-1,1,100),histtype='step',normed=True , label='im3shape-v7-r e2' , color='c')
    pl.legend(mode='expand',ncol=2)
    ylim=list(pl.ylim()); ylim[1]*=1.5; pl.ylim(ylim)

    pl.subplot(2,2,2)
    hsnr_res, _ , _= pl.hist(res_sim['snr'] ,bins=np.linspace(1,100,100),histtype='step',label='GREAT-DES snr'      , normed=True, color='r') 
    hsnr_des, _ , _= pl.hist(res_des['snr'] ,bins=np.linspace(1,100,100),histtype='step',label='im3shape-v7-r snr' , normed=True, color='b') 
    ylim=list(pl.ylim()); ylim[1]*=1.5; pl.ylim(ylim)
    pl.legend()

    pl.subplot(2,2,3)
    pl.hist(res_sim['mean_rgpp_rp'] ,bins=np.linspace(0,2,100),histtype='step',label='GREAT-DES rgpp_rp'      , normed=True, color='r') 
    pl.hist(res_sim['mean_rgpp_rp'] ,bins=np.linspace(0,2,100),histtype='step',label='GREAT-DES rgpp_rp'      , normed=True, color='r') 
    pl.hist(res_des['mean_rgpp_rp']   ,bins=np.linspace(0,2,100),histtype='step',label='im3shape-v7-r rgpp_rp' , normed=True, color='b') 
    pl.legend()
    ylim=list(pl.ylim()); ylim[1]*=1.5; pl.ylim(ylim)

    pl.subplot(2,2,4)
    pl.hist(res_sim['radius'] ,bins=np.linspace(0,4,100),histtype='step',label='GREAT-DES radius'     , normed=True, color='r') 
    pl.hist(res_des['radius'] ,bins=np.linspace(0,4,100),histtype='step',label='im3shape-v7-r radius' , normed=True, color='b') 
    pl.legend()
    ylim=list(pl.ylim()); ylim[1]*=1.5; pl.ylim(ylim)
    
    # pl.subplot(2,2,4)
    # pl.hist(res_sim['radius'] ,bins=np.linspace(0,5,100),histtype='step',label='GREAT-DES radius'      , normed=True, color='r') 
    # pl.hist(res_des['radius'] ,bins=np.linspace(0,5,100),histtype='step',label='im3shape-011-4-r radius' , normed=True, color='b') 
    # pl.legend()

    filename_fig = 'figs/histograms.match.png' 
    pl.savefig(filename_fig)
    logger.info('saved %s', filename_fig)

    # pl.figure()
    # pl.hist(res_des['FWHMPSF_IMAGE'],bins=100)
    # pl.show()


    list_snr_centers = [10,15,20,25,30,40,50,60,70,80,90]
    list_snr_edges = plotstools.get_bins_edges(list_snr_centers)

    for ib in range(1,len(list_snr_edges)):

        vb=list_snr_centers[ib-1]
        
        select1 = (list_snr_edges[ib-1] < res_sim['snr']) *  (res_sim['snr']< list_snr_edges[ib])
        select2 = (list_snr_edges[ib-1] < res_des['snr']) *  (res_des['snr']< list_snr_edges[ib])

        pl.figure()
        pl.subplot(1,2,1)
        pl.hist(res_sim['e1'][select1],bins=np.linspace(-1,1,100),histtype='step',normed=True , label='GREAT-DES' )
        pl.hist(res_des['e1'][select2],bins=np.linspace(-1,1,100),histtype='step',normed=True , label='im3shape-v7-r' )
        pl.xlabel('e1 SNR=%2.0f'%vb)   
        ylim=list(pl.ylim()); ylim[1]*=1.5; pl.ylim(ylim)
        pl.legend()

        pl.subplot(1,2,2)
        pl.hist(res_sim['e2'][select1],bins=np.linspace(-1,1,100),histtype='step',normed=True , label='GREAT-DES' )
        pl.hist(res_des['e2'][select2],bins=np.linspace(-1,1,100),histtype='step',normed=True , label='im3shape-v7-r' )
        pl.xlabel('e2 SNR=%2.0f'%vb)   
        ylim=list(pl.ylim()); ylim[1]*=1.5; pl.ylim(ylim)
        pl.legend()

        filename_fig = 'figs/histograms.e_vs_snr.%02d.png'  % ib
        pl.savefig(filename_fig)
        logger.info('saved %s', filename_fig)
        pl.close()

    # pl.figure()
    # pl.scatter(res_sim['snr'][:20000],cat_tru['snr_mean_all'][:20000],0.1)
    # pl.show()

    # pl.figure()
    # select = (res_des['snr'] < 100) * (res_des['rgpp_rp'] < 2) * (res_des['rgpp_rp'] > 0.5)
    # plotstools.plot_dist( np.concatenate( [res_des['snr'][select][:,None], res_des['rgpp_rp'][select][:,None]   ] , axis=1 ) , use_fraction=0.01   , contour=True, colormesh=False, color='b')
    # select = (res_sim['snr'] < 100) * (res_sim['rgpp_rp_1'] < 2) * (res_sim['rgpp_rp_1'] > 0.5)
    # plotstools.plot_dist( np.concatenate( [res_sim['snr'][select][:,None], res_sim['rgpp_rp_1'][select][:,None]   ] , axis=1 ) , use_fraction=0.01 , contour=True, colormesh=False, color='r')

    # plotstools.plot_dist( np.concatenate( [res_sim['snr'][:,None], res_sim['rgpp_rp_1'][:,None] ] , axis=1 ) )


    # pl.subplot(2,2,4)
    # pl.hist(cat_tru['hsm_obs_sigma'] ,bins=np.linspace(-2,10),histtype='step') 

    pl.figure(figsize=(15,5))
    pl.subplot(1,2,1)
    pl.hist(res_des['bulge_A'], bins=np.linspace(-200,400,100),histtype='step',normed=True , label='im3shape-v7-r bulge_A'  , color='r')
    pl.hist(res_sim['bulge_A'], bins=np.linspace(-200,400,100),histtype='step',normed=True , label='GREAT-DES bulge_A'      , color='b')
    pl.legend()
    ylim=list(pl.ylim()); ylim[1]*=1.2; pl.ylim(ylim)

    pl.subplot(1,2,2)
    pl.hist(res_des['disc_A'], bins=np.linspace(-2,2,100),histtype='step',normed=True , label='im3shape-v7-r disc_A'  , color='r')
    pl.hist(res_sim['disc_A'], bins=np.linspace(-2,2,100),histtype='step',normed=True , label='GREAT-DES disc_A'      , color='b')
    pl.legend()
    ylim=list(pl.ylim()); ylim[1]*=1.2; pl.ylim(ylim)

    pl.figure()
    pl.figure(figsize=(15,5))
    pl.subplot(1,2,1)
    pl.hist(res_des['dec_as'], bins=np.linspace(-200,400,100),histtype='step',normed=True , label='im3shape-v7-r bulge_A'  , color='r')
    pl.hist(res_sim['dec_as'], bins=np.linspace(-200,400,100),histtype='step',normed=True , label='GREAT-DES bulge_A'      , color='b')
    pl.legend()
    ylim=list(pl.ylim()); ylim[1]*=1.2; pl.ylim(ylim)

    pl.subplot(1,2,2)
    pl.hist(res_des['disc_A'], bins=np.linspace(-2,2,100),histtype='step',normed=True , label='im3shape-v7-r disc_A'  , color='r')
    pl.hist(res_sim['disc_A'], bins=np.linspace(-2,2,100),histtype='step',normed=True , label='GREAT-DES disc_A'      , color='b')
    pl.legend()
    ylim=list(pl.ylim()); ylim[1]*=1.2; pl.ylim(ylim)


    pl.figure()
    list_pos_centers = np.linspace(-0.03,0.03,10)
    list_pos_edges = plotstools.get_bins_edges(list_snr_centers)
    list_mean_e1 = []
    list_mean_e2 = []
    pos_type = 'ra_as'

    for ib in range(1,len(list_pos_edges)):

        print ib, len(list_pos_centers)
        
        select1 = (list_pos_edges[ib-1] < res_sim[pos_type]) *  (res_sim[pos_type]< list_pos_edges[ib])
        select2 = (list_pos_edges[ib-1] < res_des[pos_type]) *  (res_des[pos_type]< list_pos_edges[ib])

        list_mean_e1.append(np.mean(res_sim['e1'][select1]))
        list_mean_e2.append(np.mean(res_des['e2'][select2]))

    pl.show()
    import pdb; pdb.set_trace()
    list_mean_e1 = np.array(list_mean_e1)
    list_mean_e2 = np.array(list_mean_e2)

    pl.plot(list_pos_centers,list_mean_e1)
    pl.plot(list_pos_centers,list_mean_e2)

    pl.show()

    import pdb; pdb.set_trace()


def get_histograms():

    res_sim,res_tru = nbc_v7_select.get_selection_sim(selection_string_sim,cols_res=['e1','e2','snr','disc_A','bulge_A'],cols_tru=['snr','psf_e1','psf_e2'])
    res_des = nbc_v7_select.get_selection_des(selection_string_des,cols=['e1','e2','snr','mean_psf_e1_sky','mean_psf_e2_sky','disc_A','bulge_A'],n_files=args.num)

    snr_min = 10
    snr_max = 1000

    res = res_sim
    snr_min=bins_snr_edges[ibin-1]
    snr_max=bins_snr_edges[ibin]
    select = (res['snr']>snr_min) & (res['snr']<snr_max)
    res_sim = res_sim[select]

    res = res_des
    snr_min=bins_snr_edges[ibin-1]
    snr_max=bins_snr_edges[ibin]
    select = (res['snr']>snr_min) & (res['snr']<snr_max)
    res_des = res[select]
    
    pl.figure()
    pl.hist(res_sim['disc_A']);
    pl.show()    

def get_meane_vs_snr():

    res_sim,res_tru = nbc_v7_select.get_selection_sim(selection_string_sim,cols_res=['e1','e2','snr'],cols_tru=['snr','psf_e1','psf_e2'])
    res_des = nbc_v7_select.get_selection_des(selection_string_des,cols=['e1','e2','snr','mean_psf_e1_sky','mean_psf_e2_sky'],n_files=args.num)

    pl.figure()
    pl.subplot(2,2,1)
    pl.hist(res_tru['psf_e1'])
    pl.xlabel('SIM psf_e1')
    pl.subplot(2,2,2)
    pl.hist(res_tru['psf_e2'])
    pl.xlabel('SIM psf_e2')
    pl.subplot(2,2,3)
    pl.hist(res_des['mean_psf_e1_sky'])
    pl.xlabel('DES psf_e1')
    pl.subplot(2,2,4)
    pl.hist(res_des['mean_psf_e2_sky'])
    pl.xlabel('DES psf_e2')
    pl.show()


    logger.info('median PSF_e1 DES: %2.4f',np.median(res_des['mean_psf_e1_sky']))
    logger.info('median PSF_e2 DES: %2.4f',np.median(res_des['mean_psf_e2_sky']))
    logger.info('median PSF_e1 SIM: %2.4f',np.median(res_tru['psf_e1']))
    logger.info('median PSF_e2 SIM: %2.4f',np.median(res_tru['psf_e2']))

    bins_snr_edges = [0,5,10,15,20,30,40,60,100,1000]

    print 'total',len(res_des)

    list_bin_des = []
    list_bin_sim = []
    for ibin in range(1,len(bins_snr_edges)):

        # DES part
        res = res_des
        snr_min=bins_snr_edges[ibin-1]
        snr_max=bins_snr_edges[ibin]
        select = (res['snr']>snr_min) & (res['snr']<snr_max)
        nbin = len(np.nonzero(select)[0])
        mean_e1 = np.mean(res[select]['e1'])
        mean_e2 = np.mean(res[select]['e2'])
        stdm_e1 = np.std(res[select]['e1'],ddof=1)/np.sqrt(nbin)
        stdm_e2 = np.std(res[select]['e2'],ddof=1)/np.sqrt(nbin)
        mean_psfe1 = np.mean(res[select]['mean_psf_e1_sky'])
        mean_psfe2 = np.mean(res[select]['mean_psf_e2_sky'])
        stdm_psfe1 = np.std(res[select]['mean_psf_e1_sky'],ddof=1)
        stdm_psfe2 = np.std(res[select]['mean_psf_e2_sky'],ddof=1)

        list_bin_des.append( [ibin, nbin, snr_min, snr_max, mean_e1, stdm_e1,  mean_e2, stdm_e2 , mean_psfe1 , mean_psfe2 , stdm_psfe1 , stdm_psfe2] )

        # SIM part
        res = res_sim
        snr_min=bins_snr_edges[ibin-1]
        snr_max=bins_snr_edges[ibin]
        select = (res['snr']>snr_min) & (res['snr']<snr_max)
        nbin = len(np.nonzero(select)[0])
        mean_e1 = np.mean(res[select]['e1'])
        mean_e2 = np.mean(res[select]['e2'])
        stdm_e1 = np.std(res[select]['e1'],ddof=1)/np.sqrt(nbin)
        stdm_e2 = np.std(res[select]['e2'],ddof=1)/np.sqrt(nbin)
        mean_psfe1 = np.mean(res_tru[select]['psf_e1'])
        mean_psfe2 = np.mean(res_tru[select]['psf_e2'])
        stdm_psfe1 = np.std(res_tru[select]['psf_e1'],ddof=1)
        stdm_psfe2 = np.std(res_tru[select]['psf_e2'],ddof=1)

        list_bin_sim.append( [ibin, nbin, snr_min, snr_max, mean_e1, stdm_e1,  mean_e2, stdm_e2 , mean_psfe1 , mean_psfe2 , stdm_psfe1 , stdm_psfe2] )


    bin_des = tt.arr2rec(list_bin_des,dtype=zip(['ibin','nbin','snr_min','snr_max','mean_e1','stdm_e1','mean_e2','stdm_e2','mean_psfe1' , 'mean_psfe2','stdm_psfe1' , 'stdm_psfe2'],['float64']*len(list_bin_des[0])))  
    bin_sim = tt.arr2rec(list_bin_sim,dtype=zip(['ibin','nbin','snr_min','snr_max','mean_e1','stdm_e1','mean_e2','stdm_e2','mean_psfe1' , 'mean_psfe2','stdm_psfe1' , 'stdm_psfe2'],['float64']*len(list_bin_sim[0])))  

    import cPickle as pickle
    filename_pickle = 'plot_meane_vs_snr.pp2'
    pickle.dump((bin_des,bin_sim),open(filename_pickle,'w'))
    logger.info('saved %s',filename_pickle)

def plot_meane_vs_snr():

    import cPickle as pickle
    filename_pickle = 'plot_meane_vs_snr.pp2'
    bin_des,bin_sim = pickle.load(open(filename_pickle,'r'))
    logger.info('loaded %s',filename_pickle)

    cosmic_variance_stde = 0.0004
    dx=1

    logger.info('mean PSF_e1 DES: %2.4f',np.mean(bin_des['mean_psfe1']))
    logger.info('mean PSF_e2 DES: %2.4f',np.mean(bin_des['mean_psfe2']))
    logger.info('mean PSF_e1 SIM: %2.4f',np.mean(bin_sim['mean_psfe1']))
    logger.info('mean PSF_e2 SIM: %2.4f',np.mean(bin_sim['mean_psfe2']))


    pl.figure()
    pl.errorbar(bin_des['snr_max'],bin_des['mean_e1'],yerr=bin_des['stdm_e1'],fmt='.',label='e1')
    pl.errorbar(bin_des['snr_max']+dx,bin_des['mean_e2'],yerr=bin_des['stdm_e2'],fmt='.',label='e2')
    pl.axhline(0,c='k')
    # pl.yscale('symlog',linthreshy=0.0002)
    pl.xscale('log')
    pl.title('im3shape-v7-r, info_flags=0 (except SNR-related flags), radius<3.5 if stamp_size==48')
    pl.xlabel('SNR - im3shape (bin max)')
    pl.ylabel('mean(e)')
    pl.legend()

    pl.figure()
    pl.errorbar(bin_des['snr_max'],bin_des['mean_psfe1'],yerr=bin_des['stdm_psfe1'],fmt='.',label='PSF e1')
    pl.errorbar(bin_des['snr_max']+dx,bin_des['mean_psfe2'],yerr=bin_des['stdm_psfe2'],fmt='.',label='PSF e2')
    pl.xscale('log')
    pl.title('im3shape-v7-r, info_flags=0 (except SNR-related flags)')
    pl.xlabel('SNR - im3shape (bin max)')
    pl.ylabel('mean(PSF_e)')
    pl.legend()

    pl.figure()
    pl.errorbar(bin_sim['snr_max'],bin_sim['mean_e1'],yerr=bin_sim['stdm_e1'],fmt='.',label='e1')
    pl.errorbar(bin_sim['snr_max']+dx,bin_sim['mean_e2'],yerr=bin_sim['stdm_e2'],fmt='.',label='e2')
    # pl.yscale('symlog',linthreshy=0.0002)
    pl.axhline(0,c='k')
    pl.xscale('log')
    pl.title('GREAT-DES im3shape')
    pl.xlabel('SNR - im3shape (bin max)')
    pl.ylabel('mean(e)')
    pl.legend()

    pl.figure()
    pl.errorbar(bin_sim['snr_max'],bin_sim['mean_psfe1'],yerr=bin_sim['stdm_psfe1'],fmt='.',label='PSF e1')
    pl.errorbar(bin_sim['snr_max']+dx,bin_sim['mean_psfe2'],yerr=bin_sim['stdm_psfe2'],fmt='.',label='PSF e2')
    pl.xscale('log')
    pl.title('GREAT-DES')
    pl.xlabel('SNR - im3shape (bin max)')
    pl.ylabel('mean(PSF_e)')
    pl.legend()
    pl.show()


def get_PSF_leakage():

    res_sim,res_tru = nbc_v7_select.get_selection_sim(selection_string_sim,cols_res=['e1','e2','snr'],cols_tru=['snr','psf_e1','psf_e2','g1_true','g2_true'])
    # res_des = nbc_v7_select.get_selection_des(selection_string_des,cols=['e1','e2','snr','mean_psf_e1_sky','mean_psf_e2_sky'],n_files=args.num)

    bins_psf_centers = np.unique(res_tru['psf_e1'])
    print bins_psf_centers
    bins_psf_edges = plotstools.get_bins_edges(bins_psf_centers)

    list_mean_e1 = []
    list_stdm_e1 = []
    list_mean_g1 = []
    
    list_mean_e2 = []
    list_stdm_e2 = []
    list_mean_g2 = []

    for ic in range(len(bins_psf_centers)):

        select = (res_tru['psf_e1'] > (bins_psf_centers[ic]-0.001)) & (res_tru['psf_e1'] < (bins_psf_centers[ic]+0.001))
        print len(np.nonzero(select)[0])
        stdm_e1=np.std(res_sim[select]['e1'],ddof=1)/np.sqrt(len(res_sim[select]))
        mean_e1=np.mean(res_sim[select]['e1'])
        mean_g1=np.mean(res_tru[select]['g1_true'])
        list_mean_e1.append(mean_e1)
        list_stdm_e1.append(stdm_e1)
        list_mean_g1.append(mean_g1)

        select = (res_tru['psf_e2'] > (bins_psf_centers[ic]-0.001)) & (res_tru['psf_e2'] < (bins_psf_centers[ic]+0.001))
        stdm_e2=np.std(res_sim[select]['e2'],ddof=1)/np.sqrt(len(res_sim[select]))
        mean_e2=np.mean(res_sim[select]['e2'])
        mean_g2=np.mean(res_tru[select]['g2_true'])
        list_mean_e2.append(mean_e2)
        list_stdm_e2.append(stdm_e2)
        list_mean_g2.append(mean_g2)

    mean_e1=np.array(list_mean_e1)
    stdm_e1=np.array(list_stdm_e1)
    mean_e2=np.array(list_mean_e2)
    stdm_e2=np.array(list_stdm_e2)
    import fitting
    cc1,mm1,Ccm1=fitting.get_line_fit(bins_psf_centers,mean_e1,stdm_e1)
    cc2,mm2,Ccm2=fitting.get_line_fit(bins_psf_centers,mean_e2,stdm_e2)
    std_cc1 = np.sqrt(Ccm1[0,0])
    std_mm1 = np.sqrt(Ccm1[1,1])
    std_cc2 = np.sqrt(Ccm2[0,0])
    std_mm2 = np.sqrt(Ccm2[1,1])

    pl.figure()
    pl.errorbar(bins_psf_centers,list_mean_e1,yerr=list_stdm_e1,label='e1 vs PSF e1',fmt='rd')
    pl.errorbar(bins_psf_centers,list_mean_e2,yerr=list_stdm_e2,label='e2 vs PSF e2',fmt='bd')
    pl.plot(bins_psf_centers,bins_psf_centers*mm1+cc1,'r-',label=r'slope=$%2.4f \pm %2.4f$ interc=$%2.4f \pm %2.4f$' % (mm1,std_mm1,cc1,std_cc1))
    pl.plot(bins_psf_centers,bins_psf_centers*mm2+cc2,'b-',label=r'slope=$%2.4f \pm %2.4f$ interc=$%2.4f \pm %2.4f$' % (mm2,std_mm2,cc2,std_cc2))
    pl.plot(bins_psf_centers,list_mean_g1,'ro',label='mean true g1')
    pl.plot(bins_psf_centers,list_mean_g2,'bo',label='mean true g2')
    pl.legend(mode='expand')
    pl.axhline(0,c='k')
    plotstools.adjust_limits()
    pl.xlabel('PSF e')
    pl.ylabel('<e>')
    ylim=list(pl.ylim()); ylim[1]*=2; pl.ylim(ylim)
    pl.show()
    # show_plot()


    
    # import ipdb; ipdb.set_trace()
    from IPython import embed; embed()

# def show_plot(figure_id=None):    
#     if figure_id is not None:
#         fig = pl.figure(num=figure_id)
#     else:
#         fig = pl.gcf()


#     # wm = pl.get_current_fig_manager() 
#     # wm.window.attributes('-topmost', 1)
#     # wm.window.attributes('-topmost', 0)
#     # fig.canvas.manager.window.raise_()

#     # pl.show()
#     # pl.pause(1e-9)
#     # fig.canvas.manager.window.activateWindow()
#     # fig.canvas.manager.window.raise_()  


def main():

    valid_actions = ['plot_mc_vs_snr','plot_meane_vs_snr','get_meane_vs_snr','get_histograms','get_distributions','get_PSF_leakage','get_calibration','plot_calibration']

    global logger , config , args

    description = 'Get statistics and plot results of noise bias calibration runs'
    parser = argparse.ArgumentParser(description=description, add_help=True)
    parser.add_argument('-v', '--verbosity', type=int, action='store', default=2, choices=(0, 1, 2, 3 ), help='integer verbosity level: min=0, max=3 [default=2]')
    parser.add_argument('-c', '--filename_config', default='sva1.yaml',type=str, action='store', help='name of the yaml config file')
    parser.add_argument('-m', '--method', default='im3shape',type=str, action='store', help='name of the yaml config file')
    parser.add_argument('-f', '--first', type=int, action='store', default=0, help='index of the first file to analyse')
    parser.add_argument('-n', '--num', type=int, action='store', default=2, help='number of files to analyse, if -1 then =config[n_files]')
    parser.add_argument('-a', '--actions', nargs='+' ,default=None, type=str, action='store',  help='which actions to run, available: %s' % str(valid_actions) )
    
    args = parser.parse_args()
    logging_levels = { 0: logging.CRITICAL,1: logging.WARNING,2: logging.INFO,3: logging.DEBUG }
    logging_level = logging_levels[args.verbosity]; logger.setLevel(logging_level)  
    logger.info(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))

    config = yaml.load(open(args.filename_config))
    nbc_v7_select.config=config; nbc_v7_select.args = args;

    if args.actions==None:
        logger.error('no action specified, choose from %s' % valid_actions)
        return
    for action in valid_actions:
        if action in args.actions:
            logger.info('executing %s' % action)
            exec action+'()'
    for ac in args.actions:
        if ac not in valid_actions:
            print 'invalid action %s' % ac

 
    logger.info(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))


main()