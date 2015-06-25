import numpy as np; import pylab as pl; import tktools as tt;
import  sys, logging, yaml, argparse, time, copy, itertools, tktools, warnings, os, fitsio, pyfits;
warnings.simplefilter("once")
sys.path.append('/home/tomek/code/tktools');
sys.path.append('/Users/tomek/code/tktools');
logging_level = logging.INFO; logger = logging.getLogger("nbc-v7-stats"); logger.setLevel(logging_level)  
log_formatter = logging.Formatter("%(asctime)s %(name)s %(levelname)s   %(message)s", "%Y-%m-%d %H:%M:%S")
stream_handler = logging.StreamHandler(sys.stdout); stream_handler.setFormatter(log_formatter)
if logger.handlers == [] : logger.addHandler(stream_handler); logger.propagate = False
import nbc_v7_select
import plotstools, fitting

config = None
args = None

def get_shear_estimator(res,use_calibration=False,use_weights=False):

    e1 = res['e1']
    e2 = res['e2']
       
    if use_calibration:
        if ('nbc_m1' not in res.dtype.names) and ('nbc_m' not in res.dtype.names) :
            raise Exception("column nbc_m1 or nbc_m not found in res_des")
        if 'nbc_m' in res.dtype.names:
            nbc_m1 = res['nbc_m']
            nbc_m2 = res['nbc_m']
        elif 'nbc_m1' in res.dtype.names:
            nbc_m1 = res['nbc_m1']
            nbc_m2 = res['nbc_m2']
        nbc_c1 = res['nbc_c1']
        nbc_c2 = res['nbc_c2']     
    else:
        nbc_m1 = np.zeros(len(res))
        nbc_m2 = np.zeros(len(res))
        nbc_c1 = np.zeros(len(res))
        nbc_c2 = np.zeros(len(res))

    if use_weights:
        if 'w' not in res.dtype.names:
            raise Exception("column 'w' not found in res")
        # elif np.all(res['w']==1):
            # warnings.warn('using get_weights()')
            # res['w'] = add_weights.get_weights(res['snr'])
        else:
            warnings.warn("using weights from 'w' column")
            # warnings.warn("but actually replacing them with weight-snr!!!!!!")
            # res['w'] = add_weights.get_weights(res['snr'])
    else:
        warnings.warn("not using weights")
        res['w'] = np.ones(len(res))


    n_select = len(res)
    
    # mean_e1 += [ res_sim_select['e1'](np.mean(res_sim_select['e1'] - nbc_c1_sim_select)) / (1+np.mean(nbc_m1_sim_select))  ]
    mean_e1 =  np.sum(res['w']*(e1 - nbc_c1)) / np.sum(res['w']*(1+nbc_m1))  
    stdv_e1 =  np.std( (e1-nbc_c1)*(1+nbc_m1),ddof=1)  
    stdm_e1 =  stdv_e1/np.sqrt(n_select)     

    # mean_e2 += [ (np.mean(res['e2'] - nbc_c2_sim_select)) / (1+np.mean(nbc_m2_sim_select))  ]
    mean_e2 =  np.sum(res['w']*(e2 - nbc_c2)) / np.sum(res['w']*(1+nbc_m2))  
    stdv_e2 =  np.std( (e2-nbc_c2)*(1+nbc_m2),ddof=1)  
    stdm_e2 =  stdv_e2/np.sqrt(n_select)     

    mean_m1 = np.mean(nbc_m1)
    mean_m2 = np.mean(nbc_m2)

    logger.debug('get_shear_estimator: mean_e = %2.4f %2.4f' % (mean_e1,mean_e2))

    return mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2



def get_mc(res_sim,res_tru,res_des=None,filename_str=None,use_calibration=False,use_weights=False,resample_errors=False,correct_selection_bias=False):


    necessary_cols_res = ['mean_psf_e1_sky','mean_psf_e2_sky','w','nbc_m','nbc_c1','nbc_c2','nbc_alpha','e1','e2','mean_rgpp_rp','snr']
    necessary_cols_tru = ['sheared_e1','sheared_e2','intrinsic_e1','intrinsic_e2','id_shear','g1_true','g2_true','sf_hlr']

    logger.debug('getting necessary columns')
    res_sim = res_sim[necessary_cols_res].copy()
    res_tru = res_tru[necessary_cols_tru].copy()
    if res_des!=None: res_des = res_des[necessary_cols_res].copy()
    logger.debug('done!')


    n_jack = 1
    mean_e1 = []; stdv_e1 = []; stdm_e1 = []; mean_e2 = []; stdv_e2 = []; stdm_e2 = []; true_e1 = []; true_e2 = []; mean_mm = []; mean_ei1 = []; mean_ei2 = []; stdm_ei1 = []; stdm_ei2 = []; mean_es1 = []; mean_es2 = []; stdm_es1 = []; stdm_es2 = []; stdm_es1eo1 = []; stdm_es2eo2 = [];


    res_she = res_sim.copy()
    res_she['e1'] = res_tru['sheared_e1']
    res_she['e2'] = res_tru['sheared_e2']

    if correct_selection_bias:
        res_int = res_sim.copy()
        res_int['e1'] = res_tru['intrinsic_e1']
        res_int['e2'] = res_tru['intrinsic_e2']

    n_gals = len(res_sim)

    mean_e1_entire_sample, stdv_e1_entire_sample, stdm_e1_entire_sample, mean_e2_entire_sample, stdv_e2_entire_sample, stdm_e2_entire_sample, _, _ = get_shear_estimator(res_sim,use_calibration=use_calibration,use_weights=use_weights)
    # mean_e1i_entire_sample, stdv_e1i_entire_sample, stdm_e1i_entire_sample, mean_e2i_entire_sample, stdv_e2i_entire_sample, stdm_e2i_entire_sample, _, _ = get_shear_estimator(res_int,use_calibration=use_calibration,use_weights=use_weights)

    
    for ig,vg in enumerate(config['shear']):

        select = res_tru['id_shear'] == ig
    
        res_sim_select = res_sim[select]
        tru_sim_select = res_tru[select]
        if correct_selection_bias:
            res_int_select = res_int[select]
        res_she_select = res_she[select]

        n_select = len(res_sim_select)
        e1o_mean, e1o_stdv, e1o_stdm, e2o_mean, e2o_stdv, e2o_stdm, m1o_mean, m2o_mean = get_shear_estimator(res_sim_select,use_calibration=use_calibration,use_weights=use_weights)
        if correct_selection_bias:
            e1i_mean, e1i_stdv, e1i_stdm, e2i_mean, e2i_stdv, e2i_stdm, m1i_mean, m2i_mean = get_shear_estimator(res_int_select,use_calibration=False,use_weights=use_weights)
        e1s_mean, e1s_stdv, e1s_stdm, e2s_mean, e2s_stdv, e2s_stdm, m1s_mean, m2s_mean = get_shear_estimator(res_she_select,use_calibration=False,use_weights=use_weights)

        true_e1 += [ np.mean(tru_sim_select['g1_true']) ]
        true_e2 += [ np.mean(tru_sim_select['g2_true']) ]
        mean_e1  += [ e1o_mean ]
        stdv_e1  += [ e1o_stdv ]
        stdm_e1  += [ e1o_stdm ]
        mean_e2  += [ e2o_mean ]
        stdv_e2  += [ e2o_stdv ]
        stdm_e2  += [ e2o_stdm ]
        mean_es1 += [ e1s_mean ]
        mean_es2 += [ e2s_mean ]
        stdm_es1 += [ e1s_stdm ]
        stdm_es2 += [ e2s_stdm ]
        mean_mm  += [ m1o_mean ]

        if correct_selection_bias:
            mean_ei1 += [ e1i_mean ]
            mean_ei2 += [ e2i_mean ]
            stdm_ei1 += [ e1i_stdm ]
            stdm_ei2 += [ e2i_stdm ]
            covmat1 = np.cov(res_sim_select['e1'],tru_sim_select['intrinsic_e1'])
            covmat2 = np.cov(res_sim_select['e2'],tru_sim_select['intrinsic_e2'])
            import pdb; pdb.set_trace()
            stdm_es1eo1 += [ np.sqrt(covmat1[0,0]+covmat1[1,1]-2*covmat1[0,1])/np.sqrt(len(res_sim_select)) ]
            stdm_es2eo2 += [ np.sqrt(covmat2[0,0]+covmat2[1,1]-2*covmat2[0,1])/np.sqrt(len(res_sim_select)) ]

        # mean_e1 += [ res_sim_select['e1'](np.mean(res_sim_select['e1'] - nbc_c1_sim_select)) / (1+np.mean(nbc_m1_sim_select))  ]
        # mean_e1 += [ np.sum(res_sim_select['w']*(res_sim_select['e1'] - nbc_c1_sim_select)) / np.sum(res_sim_select['w']*(1+nbc_m1_sim_select))  ]
        # stdv_e1 += [ np.std(res_sim_select['e1'],ddof=1)  ]
        # stdm_e1 += [ np.std(res_sim_select['e1'],ddof=1)/np.sqrt(n_select)     ]
        
        # # mean_e2 += [ (np.mean(res_sim_select['e2'] - nbc_c2_sim_select)) / (1+np.mean(nbc_m2_sim_select))  ]
        # mean_e1 += [ np.sum(res_sim_select['w']*(res_sim_select['e2'] - nbc_c2_sim_select)) / np.sum(res_sim_select['w']*(1+nbc_m2_sim_select))  ]
        # stdv_e2 += [ np.std(res_sim_select['e2'],ddof=1)  ]
        # stdm_e2 += [ np.std(res_sim_select['e2'],ddof=1)/np.sqrt(n_select)     ]

        # selection bias 



    mean_e1 = np.array(mean_e1); stdv_e1 = np.array(stdv_e1); stdm_e1 = np.array(stdm_e1); mean_e2 = np.array(mean_e2); stdv_e2 = np.array(stdv_e2); stdm_e2 = np.array(stdm_e2); true_e1 = np.array(true_e1); true_e2 = np.array(true_e2);  mean_ei1 = np.array(mean_ei1);     mean_ei2 = np.array(mean_ei2);     stdm_ei1 = np.array(stdm_ei1);     stdm_ei2 = np.array(stdm_ei2);     mean_es1 = np.array(mean_es1);     mean_es2 = np.array(mean_es2);     stdm_es1 = np.array(stdm_es1);     stdm_es2 = np.array(stdm_es2);  stdm_es1eo1 = np.array(stdm_es1eo1); stdm_es2eo2 = np.array(stdm_es2eo2)
    

    if correct_selection_bias:

        e1c_mean = mean_e1 - mean_ei1 
        e2c_mean = mean_e2 - mean_ei2
        e1c_stdm = stdm_es1eo1
        e2c_stdm = stdm_es2eo2

        mean_e1 = e1c_mean
        mean_e2 = e2c_mean
        stdm_e1 = e1c_stdm
        stdm_e2 = e2c_stdm


    import fitting


    # fit line to measured shear
    cc1,mm1,Ccm1=fitting.get_line_fit(true_e1,mean_e1,stdm_e1)
    cc2,mm2,Ccm2=fitting.get_line_fit(true_e2,mean_e2,stdm_e2)
    if resample_errors:    
        t1=mm1*true_e1+cc1
        t2=mm2*true_e2+cc2
        res1=np.abs(t1-mean_e1)
        res2=np.abs(t2-mean_e2)
        err_e1 = np.mean(res1)
        err_e2 = np.mean(res2)
        cc1,mm1,Ccm1=fitting.get_line_fit(true_e1,mean_e1,stdm_e1*0+err_e1)
        cc2,mm2,Ccm2=fitting.get_line_fit(true_e2,mean_e2,stdm_e2*0+err_e2)
    std_cc1 = np.sqrt(Ccm1[0,0])
    std_mm1 = np.sqrt(Ccm1[1,1])
    std_cc2 = np.sqrt(Ccm2[0,0])
    std_mm2 = np.sqrt(Ccm2[1,1])

    # fit line to sheared intrinsic
    cc1s,mm1s,Ccm1s=fitting.get_line_fit(true_e1,mean_es1,stdm_es1)
    cc2s,mm2s,Ccm2s=fitting.get_line_fit(true_e2,mean_es2,stdm_es2)
    std_cc1s = np.sqrt(Ccm1s[0,0])
    std_mm1s = np.sqrt(Ccm1s[1,1])
    std_cc2s = np.sqrt(Ccm2s[0,0])
    std_mm2s = np.sqrt(Ccm2s[1,1])

    cc12,mm12,Ccm12=fitting.get_line_fit(true_e2,mean_e1-true_e1,stdm_e1)
    cc21,mm21,Ccm21=fitting.get_line_fit(true_e1,mean_e2-true_e2,stdm_e2)
    std_mm12 = np.sqrt(Ccm12[1,1])
    std_cc12 = np.sqrt(Ccm12[0,0])
    std_mm21 = np.sqrt(Ccm21[1,1])
    std_cc21 = np.sqrt(Ccm21[0,0])

    mm = (mm1+mm2)/2.
    std_mm = np.sqrt( (std_mm1**2 + std_mm2**2)/2. )
    cc = (cc1+cc2)/2.
    std_cc = np.sqrt( (std_cc1**2 + std_cc2**2)/2. )

    mms = (mm1s + mm2s)/2.
    std_mms = np.sqrt( (std_mm1s**2 + std_mm2s**2)/2. )
    ccs = (cc1s+cc2s)/2.
    std_ccs = np.sqrt( (std_cc1s**2 + std_cc2s**2)/2. )

    # leakage from SIM
    pcc1,pmm1,pcc2,pmm2,std_pcc1,std_pmm1,std_pcc2,std_pmm2,pmean_e1,pstdm_e1,pmean_e2,pstdm_e2,pmean_g1,pmean_g2,bins_psf_centers = get_PSF_leakage(res_sim,res_tru,use_calibration=use_calibration,use_weights=use_weights)
    pmm = (pmm1+pmm2)/2.
    pcc = (pcc1+pcc2)/2.
    std_pmm = np.sqrt( (std_pmm1**2 + std_pmm2**2)/2. )
    std_pcc = np.sqrt( (std_pcc1**2 + std_pcc2**2)/2. )

    # leakage from DES
    if res_des!=None: dcc1,dmm1,dcc2,dmm2,std_dcc1,std_dmm1,std_dcc2,std_dmm2,dmean_e1,dstdm_e1,dmean_e2,dstdm_e2,dmean_g1,dmean_g2,bins_psf_centers_d = get_PSF_leakage(res_des,None,use_calibration=use_calibration,use_weights=use_weights)


    # print out 

    logger.info('noise_bias_calibration=%d weights=%d correct_selection_bias=%d', use_calibration, use_weights, correct_selection_bias)
    logger.info( "--- SIM n_gals = %d", n_gals )        
    logger.info( "--- SIM mean_e1_all  = %2.5f +/- %2.5f" % (mean_e1_entire_sample, stdm_e1_entire_sample) )        
    logger.info( "--- SIM mean_e2_all  = %2.5f +/- %2.5f" % (mean_e2_entire_sample, stdm_e2_entire_sample) )        
    # logger.info( "--- SIM mean_e1i_all = %2.5f +/- %2.5f" % (mean_e1i_entire_sample, stdm_e1i_entire_sample) )        
    # logger.info( "--- SIM mean_e2i_all = %2.5f +/- %2.5f" % (mean_e2i_entire_sample, stdm_e2i_entire_sample) )        
    logger.info( "--- SIM m1  = % 2.5f +/- % 2.5f"   % (mm1,std_mm1)        )               
    logger.info( "--- SIM m2  = % 2.5f +/- % 2.5f"   % (mm2,std_mm2)        )               
    logger.info( "--- SIM c1  = % 2.5f +/- % 2.5f"   % (cc1,std_cc1)        )               
    logger.info( "--- SIM c2  = % 2.5f +/- % 2.5f"   % (cc2,std_cc2)        )               
    logger.info( "--- SIM m   = % 2.5f +/- % 2.5f"   % (mm,std_mm)          )             
    logger.info( "--- SIM c   = % 2.5f +/- % 2.5f"   % (cc,std_cc)          )             
    logger.info( "--- SIM m12 = % 2.5f +/- % 2.5f"   % (mm12,std_mm12)      )                  
    logger.info( "--- SIM m21 = % 2.5f +/- % 2.5f"   % (mm21,std_mm21)      )                  
    logger.info( "--- SIM m1s = % 2.5f +/- % 2.5f"   % (mm1s,std_mm1s)        )               
    logger.info( "--- SIM m2s = % 2.5f +/- % 2.5f"   % (mm2s,std_mm2s)        )               
    logger.info( "--- SIM c1s = % 2.5f +/- % 2.5f"   % (cc1s,std_cc1s)        )  
    if 'nbc_m' in res_sim.dtype.names:     logger.info( "--- SIM (this may be old) nbc_m = % 2.5f"            % (np.mean(res_sim['nbc_m'])  ))               
    if 'nbc_alpha' in res_sim.dtype.names: logger.info( "--- SIM (this may be old) nbc_a = % 2.5f"            % (np.mean(res_sim['nbc_alpha'])  ))               
    logger.info( "--- SIM c2s = % 2.5f +/- % 2.5f"   % (cc2s,std_cc2s)        )               
    logger.info( "--- SIM pm1 = % 2.5f +/- % 2.5f"  % (pmm1,std_pmm1)      )                   
    logger.info( "--- SIM pm2 = % 2.5f +/- % 2.5f"  % (pmm2,std_pmm2)      )                   
    logger.info( "--- SIM pc1 = % 2.5f +/- % 2.5f"  % (pcc1,std_pcc1)      )                   
    logger.info( "--- SIM pc2 = % 2.5f +/- % 2.5f"  % (pcc2,std_pcc2)      )                   
    if 'ra_as' in  res_sim.dtype.names: logger.info( "--- SIM x1 = % 2.4f +/- %2.4f" % ( np.mean(res_sim['ra_as']), np.std(res_sim['ra_as'],ddof=1)/np.sqrt(len(res_sim)) ) )
    if 'dec_as' in res_sim.dtype.names: logger.info( "--- SIM x2 = % 2.4f +/- %2.4f" % ( np.mean(res_sim['dec_as']), np.std(res_sim['dec_as'],ddof=1)/np.sqrt(len(res_sim)) ) )
    if res_des!=None: logger.info( "--- DES n_gals  = %d" % len(res_des) )
    if res_des!=None: logger.info( "--- DES pm1     = % 2.5f +/- % 2.5f" % (dmm1,std_dmm1) )
    if res_des!=None: logger.info( "--- DES pm2     = % 2.5f +/- % 2.5f" % (dmm2,std_dmm2) )
    if res_des!=None: logger.info( "--- DES pc1     = % 2.5f +/- % 2.5f" % (dcc1,std_dcc1) )
    if res_des!=None: logger.info( "--- DES pc2     = % 2.5f +/- % 2.5f" % (dcc2,std_dcc2) )
    if args.use_calibration:
        if res_des!=None: logger.info( "--- DES mean_e1 = % 2.5f +/- % 2.5f" % (np.mean(res_des['e1']-res_des['nbc_c1']),np.std(res_des['e1']-res_des['nbc_c1'],ddof=1)/np.sqrt(len(res_des))) )
        if res_des!=None: logger.info( "--- DES mean_e2 = % 2.5f +/- % 2.5f" % (np.mean(res_des['e2']-res_des['nbc_c2']),np.std(res_des['e2']-res_des['nbc_c2'],ddof=1)/np.sqrt(len(res_des))) )
    else:
        if res_des!=None: logger.info( "--- DES mean_e1 = % 2.5f +/- % 2.5f" % (np.mean(res_des['e1']),np.std(res_des['e1'],ddof=1)/np.sqrt(len(res_des))) )
        if res_des!=None: logger.info( "--- DES mean_e2 = % 2.5f +/- % 2.5f" % (np.mean(res_des['e2']),np.std(res_des['e2'],ddof=1)/np.sqrt(len(res_des))) )

    if filename_str!=None:
        pl.figure(figsize=(20,15))

        pl.subplot(3,3,1)
        pl.errorbar(true_e1,mean_e1-true_e1,yerr=stdm_e1,fmt='r.')
        pl.errorbar(true_e2,mean_e2-true_e2,yerr=stdm_e2,fmt='m.')
        if correct_selection_bias: pl.errorbar(true_e1,mean_ei1,yerr=stdm_ei1,fmt='b.')
        if correct_selection_bias: pl.errorbar(true_e2,mean_ei2,yerr=stdm_ei2,fmt='c.')
        pl.plot(true_e1,true_e1*0,'-k')
        pl.plot(true_e1,(mm1-1)*true_e1+cc1,    'r-',label=r'<e1> vs e1t sl=$% 2.4f \pm %2.4f$ int=$% 2.4f \pm % 2.4f$' % (mm1,std_mm1,cc1,std_cc1))
        pl.plot(true_e2,(mm2-1)*true_e2+cc2,    'm-',label=r'<e2> vs e2t sl=$% 2.4f \pm %2.4f$ int=$% 2.4f \pm % 2.4f$' % (mm2,std_mm2,cc2,std_cc2))
        pl.plot(true_e1,(mm1s-1)*true_e1+cc1s,  'b-',label=r'<e1> vs e1t sl=$% 2.4f \pm %2.4f$ int=$% 2.4f \pm % 2.4f$' % (mm1s,std_mm1s,cc1s,std_cc1s))
        pl.plot(true_e2,(mm2s-1)*true_e2+cc2s,  'c-',label=r'<e2> vs e2t sl=$% 2.4f \pm %2.4f$ int=$% 2.4f \pm % 2.4f$' % (mm2s,std_mm2s,cc2s,std_cc2s))
        pl.xlabel('g_true')
        pl.ylabel('<e> - g_true')
        ylim=list(pl.ylim()); ylim[1]*=2; pl.ylim(ylim)
        pl.legend(mode='expand',fontsize=8,framealpha=0.0,frameon=False)

        pl.subplot(3,3,2)
        pl.errorbar(true_e2,mean_e1-true_e1,yerr=stdm_e1,fmt='r.')
        pl.errorbar(true_e1,mean_e2-true_e2,yerr=stdm_e2,fmt='m.')
        pl.plot(true_e1,true_e1*0,'-k')
        pl.plot(true_e2,(mm12)*true_e2+cc12,'r-',label=r'<e1> vs e2t sl=$% 2.4f \pm %2.4f$ int=$% 2.4f \pm % 2.4f$' % (mm12,std_mm12,cc12,std_cc12))
        pl.plot(true_e1,(mm21)*true_e1+cc21,'m-',label=r'<e2> vs e1t sl=$% 2.4f \pm %2.4f$ int=$% 2.4f \pm % 2.4f$' % (mm21,std_mm21,cc21,std_cc21))
        pl.xlabel('g_true')
        pl.ylabel('<e> - g_true')
        ylim=list(pl.ylim()); ylim[1]*=2; pl.ylim(ylim)
        pl.legend(mode='expand',fontsize=8,framealpha=0.0,frameon=False)

        pl.subplot(3,3,3)
        pl.errorbar(bins_psf_centers,pmean_e1-pmean_g1,yerr=pstdm_e1,fmt='rd') # ,label='e1 vs PSF e1'
        pl.errorbar(bins_psf_centers,pmean_e2-pmean_g2,yerr=pstdm_e2,fmt='md') # ,label='e2 vs PSF e2'
        pl.plot(bins_psf_centers,bins_psf_centers*pmm1+pcc1,'r-',label=r'sl=$% 2.4f \pm %2.4f$ int=$% 2.4f \pm %2.4f$' % (pmm1,std_pmm1,pcc1,std_pcc1))
        pl.plot(bins_psf_centers,bins_psf_centers*pmm2+pcc2,'m-',label=r'sl=$% 2.4f \pm %2.4f$ int=$% 2.4f \pm %2.4f$' % (pmm2,std_pmm2,pcc2,std_pcc2))
        pl.legend(mode='expand',fontsize=8,framealpha=0.0,frameon=False)
        pl.axhline(0,c='k')
        pl.xlabel('PSF e')
        pl.ylabel('<e>')
        pl.suptitle( filename_str )
        ylim=list(pl.ylim()); ylim[1]*=1.5; pl.ylim(ylim)
        pl.xlim([-0.03,0.03])

        # pl.ylim([-0.005,0.01])



        try:
            pl.subplot(3,3,4)
            pl.hist(res_sim['e1'],bins=np.linspace(-1,1,100),color='r',histtype='step',normed=True,label='SIM')
            pl.hist(res_sim['e2'],bins=np.linspace(-1,1,100),color='m',histtype='step',normed=True)
            if res_des!= None: pl.hist(res_des['e1'],bins=np.linspace(-1,1,100),color='b',histtype='step',normed=True,label='DES')
            if res_des!= None: pl.hist(res_des['e2'],bins=np.linspace(-1,1,100),color='c',histtype='step',normed=True)       
            pl.legend(mode='expand',fontsize=8,framealpha=0.0,frameon=False)
            pl.xlabel('e')
            # ylim=list(pl.ylim()); ylim[1]*=1.5; pl.ylim(ylim)
        except Exception, errmsg:
            logger.error('plot 334 died %s', errmsg)

        try:
            pl.subplot(3,3,5)
            pl.hist(res_sim['snr'],bins=np.linspace(0,100,200),color='r',histtype='step',normed=True,label='SIM')
            if res_des!=None: pl.hist(res_des['snr'],bins=np.linspace(0,100,100),color='b',histtype='step',normed=True,label='DES')
            pl.legend(mode='expand',fontsize=8,framealpha=0.0,frameon=False)
            pl.xlabel('SNR')
            # ylim=list(pl.ylim()); ylim[1]*=2; pl.ylim(ylim)
        except Exception, errmsg:
            logger.error('plot 335 died %s', errmsg)

        try:
            pl.subplot(3,3,6)
            pl.hist(res_sim['mean_rgpp_rp'],bins=np.linspace(1,3,200),color='r',histtype='step',normed=True,label='SIM')
            if res_des!=None: pl.hist(res_des['mean_rgpp_rp'],bins=np.linspace(1,3),color='b',histtype='step',normed=True,label='DES')
            pl.legend(mode='expand',fontsize=8,framealpha=0.0,frameon=False)
            pl.xlabel('Rgp/Rp')
            # ylim=list(pl.ylim()); ylim[1]*=1.5; pl.ylim(ylim)
        except Exception, errmsg:
            logger.error('plot 336 died, %s', errmsg)

        try:
            pl.subplot(3,3,7)
            pl.hist(res_sim['radius'],bins=np.linspace(0,3,200),color='r',histtype='step',normed=True,label='SIM')
            pl.hist(res_tru['sf_hlr'],bins=np.linspace(0,3,200),color='b',histtype='step',normed=True,label='TRU')
            if res_des!=None: pl.hist(res_des['radius'],bins=np.linspace(0,3,200),color='b',histtype='step',normed=True,label='DES')
            pl.legend(mode='expand',fontsize=8,framealpha=0.0,frameon=False)
            pl.xlabel('radius')
            # ylim=list(pl.ylim()); ylim[1]*=1.5; pl.ylim(ylim)
        except Exception, errmsg:
            logger.error('plot 337 died %s', errmsg)

        try:
            pl.subplot(3,3,8)
            pl.hist(res_sim['ra_as'],1000,color='r',histtype='step',normed=True,label='SIM')
            pl.hist(res_sim['dec_as'],1000,color='b',histtype='step',normed=True,label='SIM')
            title = " x1 = % 2.4f +/- %2.4f" % ( np.mean(res_sim['ra_as']), np.std(res_sim['ra_as'],ddof=1)/np.sqrt(len(res_sim)) ) 
            title += "\n x2 = % 2.4f +/- %2.4f" % ( np.mean(res_sim['dec_as']), np.std(res_sim['dec_as'],ddof=1)/np.sqrt(len(res_sim)) ) 
            pl.title(title)
            pl.xlabel('ra_as/dec_as')
            # ylim=list(pl.ylim()); ylim[1]*=1.5; pl.ylim(ylim)
        except Exception, errmsg:
            logger.error('plot 338 died %s', errmsg)


        try:
            filename_fig = os.path.join(args.output_dir,'figs/bias.%s.png' % (filename_str))
        except:
            filename_fig = os.path.join('.','figs/bias.%s.png' % (filename_str))
            
        pl.savefig(filename_fig)
        logger.info('saved %s',filename_fig)
        pl.close()

    #TODO: make resampled error-bars

    result_list = (mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1_entire_sample,mean_e2_entire_sample,stdm_e1_entire_sample,stdm_e2_entire_sample,mms,std_mms,ccs,std_ccs,mm1s,std_mm1s,cc1s,std_cc1s,mm2s,std_mm2s,cc2s,std_cc2s)
    # result_str =  ['mm','std_mm','cc','std_cc','mm1','std_mm1','mm2','std_mm2','cc1','std_cc1','cc1','std_cc2','pmm','std_pmm','pcc','std_pcc','pmm1','std_pmm1','pmm2','std_pmm2','e1i_mean','e2i_mean','e1i_stdm','e2i_stdm','mms','std_mms','ccs','std_ccs','mm1s','std_mm1s','cc1s','std_cc1s','mm2s','std_mm2s','cc2s','std_cc2s']
    # result_dict = dict(zip(result_str,result_list))
    # return mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,e1i_mean,e2i_mean,e1i_stdm,e2i_stdm,mms,std_mms,ccs,std_ccs,

    return result_list



def get_PSF_leakage(res,res_tru=None,use_calibration=False,use_weights=False):

    # bins_psf_edges = np.linspace(-0.014,0.014,5)
    # bins_psf_centers = plotstools.get_bins_centers(bins_psf_edges)

    if args.method=='ngmix':
        # for NGMIX we use true PSF shape as the measured one is not provided in GREAT-DES results catalogs
        # so bin centers correspond to true PSF ellipticity
        bins_psf_centers = np.array([-0.02      , -0.00666667,  0.00666667,  0.02      ])
        bins_psf_edges = plotstools.get_bins_edges(bins_psf_centers)
    elif args.method=='im3shape':
        # there is a difference between measured shape in im3shape and in true PSF shape in sims -- 
        # so here we use the measured PSF shape in im3shape, as that is used in analysis in DES
        # there is also a sign flip on e1 for PSF
        if res_tru!=None:
            warnings.warn('using 5 PSF bins')
            bins_psf_edges = np.linspace(-0.014,0.014,5)
            # pl.hist(res['mean_psf_e1_sky'],np.linspace(-0.014,0.014,100)); pl.hist(res['mean_psf_e1_sky'],bins_psf_edges,histtype='step'); pl.show()
        else:
            bins_psf_edges = np.linspace(-0.01,0.01,50)
            warnings.warn('using more PSF bins')

        bins_psf_centers = plotstools.get_bins_centers(bins_psf_edges)
    else:
        raise Exception('invalid method name %s' % args.method)

    # pl.hist(res['mean_psf_e1_sky'],1000); pl.plot(bins_psf_edges,np.zeros_like(bins_psf_edges),'rd'); pl.show()

    logger.debug('bins_psf_edges: [ ' +  (' '.join(['% 5.4f']*len(bins_psf_edges)))%tuple(bins_psf_edges) + ' ]')


    list_mean_e1 = []
    list_stdm_e1 = []
    list_mean_g1 = []
    
    list_mean_e2 = []
    list_stdm_e2 = []
    list_mean_g2 = []

    list_psf_medians_e1 = []
    list_psf_medians_e2 = []
    for ic in range(1,len(bins_psf_edges)):

        select = (res['mean_psf_e1_sky'] > (bins_psf_edges[ic-1])) & (res['mean_psf_e1_sky'] < (bins_psf_edges[ic]))
        res_select=res[select]
        pmean_e1, pstdv_e1, pstdm_e1, pmean_e2, pstdv_e2, pstdm_e2, pmean_m1, pmean_m2 = get_shear_estimator(res_select,use_calibration=use_calibration,use_weights=use_weights)
        if res_tru!=None:
            pmean_g1=np.mean(res_tru[select]['g1_true'])
        else:
            pmean_g1=pmean_e1*0
        list_mean_e1.append(pmean_e1)
        list_stdm_e1.append(pstdm_e1)
        list_mean_g1.append(pmean_g1)
        median_psf_e1 = np.median(res_select['mean_psf_e1_sky'])
        list_psf_medians_e1.append(median_psf_e1)
        logger.debug('getting PSF leakage e1, bin %2d [% 5.4f % 5.4f], median=%2.4f, n_gals=%6d' % (ic, bins_psf_edges[ic-1], bins_psf_edges[ic], median_psf_e1, len(res_select)))

        select = (res['mean_psf_e2_sky'] > (bins_psf_edges[ic-1])) & (res['mean_psf_e2_sky'] < (bins_psf_edges[ic]))
        res_select=res[select]
        pmean_e1, pstdv_e1, pstdm_e1, pmean_e2, pstdv_e2, pstdm_e2, pmean_m1, pmean_m2 = get_shear_estimator(res_select,use_calibration=use_calibration,use_weights=use_weights)
        if res_tru!=None:
            pmean_g2=np.mean(res_tru[select]['g2_true'])
        else:
            pmean_g2=pmean_e2*0
        list_mean_e2.append(pmean_e2)
        list_stdm_e2.append(pstdm_e2)
        list_mean_g2.append(pmean_g2)     
        median_psf_e2 = np.median(res_select['mean_psf_e2_sky'])
        list_psf_medians_e2.append(median_psf_e2)
        logger.debug('getting PSF leakage e2, bin %2d [% 5.4f % 5.4f], median=%2.4f, n_gals=%6d' % (ic, bins_psf_edges[ic-1], bins_psf_edges[ic], median_psf_e2, len(res_select)))

        if np.abs(pmean_e1)>10:
            warnings.warn('np.abs(pmean_e1)>10')
            import pdb; pdb.set_trace()      

    pmean_e1=np.array(list_mean_e1)
    pstdm_e1=np.array(list_stdm_e1)
    pmean_e2=np.array(list_mean_e2)
    pstdm_e2=np.array(list_stdm_e2)  
    pmean_g1=np.array(list_mean_g1)
    pmean_g2=np.array(list_mean_g2)
    psf_medians_e1=np.array(list_psf_medians_e1)
    psf_medians_e2=np.array(list_psf_medians_e2)
  
    pcc1,pmm1,pCcm1=fitting.get_line_fit(psf_medians_e1,pmean_e1-pmean_g1,pstdm_e1)
    pcc2,pmm2,pCcm2=fitting.get_line_fit(psf_medians_e2,pmean_e2-pmean_g2,pstdm_e2)
    std_pcc1 = np.sqrt(pCcm1[0,0])
    std_pmm1 = np.sqrt(pCcm1[1,1])
    std_pcc2 = np.sqrt(pCcm2[0,0])
    std_pmm2 = np.sqrt(pCcm2[1,1])

    bins_psf_centers = np.array(list_psf_medians_e1)

    return pcc1,pmm1,pcc2,pmm2,std_pcc1,std_pmm1,std_pcc2,std_pmm2,pmean_e1,pstdm_e1,pmean_e2,pstdm_e2,pmean_g1,pmean_g2,bins_psf_centers

