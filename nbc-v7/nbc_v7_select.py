import numpy as np; import pylab as pl; import tktools as tt;
import  sys
sys.path.append('/Users/tomek/code/tktools')
import logging, yaml, argparse, time, copy, itertools, tabletools, warnings
warnings.simplefilter("once")
# from nbc2_dtypes import *
logging_level = logging.INFO; logger = logging.getLogger("nbc-v7"); logger.setLevel(logging_level)  
log_formatter = logging.Formatter("%(asctime)s %(name)s %(levelname)s   %(message)s", "%Y-%m-%d %H:%M:%S")
stream_handler = logging.StreamHandler(sys.stdout); stream_handler.setFormatter(log_formatter)
if logger.handlers == [] : logger.addHandler(stream_handler); logger.propagate = False

global means_snr
means_snr = None
global means_rgp
means_rgp = None

def rename_ngmix_cols(cat_res,cat_tru=None):

    if 'g' in cat_res.dtype.names:

        cat_res = tabletools.appendColumn(cat_res,'e1', cat_res['g'][:,0])
        cat_res = tabletools.appendColumn(cat_res,'e2', cat_res['g'][:,1])
        cat_res = tabletools.appendColumn(cat_res,'error_flag', cat_res['flags'])
        # cat_res = tabletools.appendColumn(cat_res,'info_flag', cat_res['arate'])
        cat_res = tabletools.appendColumn(cat_res,'snr', cat_res['s2n_r'])

    if 'arate' in cat_res.dtype.names:
        cat_res = tabletools.appendColumn(cat_res,'info_flag', cat_res['arate'])
    else:
        cat_res = tabletools.appendColumn(cat_res,'info_flag', np.ones(len(cat_res)),dtype=np.bool )

    if 'T' in cat_res.dtype.names:
        cat_res = tabletools.appendColumn(cat_res,'mean_rgpp_rp', cat_res['T'])
    if 'T' in cat_res.dtype.names:
        cat_res = tabletools.appendColumn(cat_res,'radius', cat_res['T'])

    if 'log_T' in cat_res.dtype.names:
        cat_res = tabletools.appendColumn(cat_res,'mean_rgpp_rp', cat_res['log_T'])
    if 'log_T' in cat_res.dtype.names:
        cat_res = tabletools.appendColumn(cat_res,'radius', cat_res['log_T'])

    elif 'ngmix_EXP_T' in cat_res.dtype.names:
        cat_res = tabletools.appendColumn(cat_res,'mean_rgpp_rp', cat_res['ngmix_EXP_T'])
    elif 'pars' in cat_res.dtype.names:
        cat_res = tabletools.appendColumn(cat_res,'mean_rgpp_rp', np.sqrt(np.abs(cat_res['pars'][:,4]/2.)))       
    else:
        cat_res = tabletools.appendColumn(cat_res,'mean_rgpp_rp', np.random.uniform(len(cat_res))*5 )

    if 'g_sens' in cat_res.dtype.names:

        # cat_res = tabletools.appendColumn(cat_res,'nbc_m1', (cat_res['g_sens'][:,0]-1) )
        # cat_res = tabletools.appendColumn(cat_res,'nbc_m2', (cat_res['g_sens'][:,1]-1) )
        cat_res = tabletools.appendColumn(cat_res,'nbc_c1', np.zeros(len(cat_res)))
        cat_res = tabletools.appendColumn(cat_res,'nbc_c2', np.zeros(len(cat_res)))
        # cat_res['nbc_c1'] = 0
        # cat_res['nbc_c2'] = 0

        warnings.warn('using mean g_sens')
        col_m = ( cat_res['g_sens'][:,0] + cat_res['g_sens'][:,1] ) / 2.
        cat_res = tabletools.appendColumn(cat_res,'nbc_m1', (col_m - 1) )
        cat_res = tabletools.appendColumn(cat_res,'nbc_m2', (col_m - 1) )

    # if 'ngmix_PSFREC_E_1' not in cat_res.dtype.names:

    #     try:
    #         cat_res = tabletools.appendColumn(cat_res,'mean_psf_e1_sky', cat_tru['psf_e1'])
    #         cat_res = tabletools.appendColumn(cat_res,'mean_psf_e2_sky', cat_tru['psf_e2'])
    #     except:

    #         pass


    if 'mean_psf_e1_sky'  not in cat_res.dtype.names:
        
        cat_res = tabletools.appendColumn(cat_res,'mean_psf_e1_sky', cat_res['psf_g'][:,0])
        cat_res = tabletools.appendColumn(cat_res,'mean_psf_e2_sky', cat_res['psf_g'][:,1])
        cat_res = tabletools.appendColumn(cat_res,'mean_psf_fwhm', cat_res['psf_T'])


    if 'ngmix_EXP_E_1' in cat_res.dtype.names:

        cat_res = tabletools.appendColumn(cat_res,'e1', cat_res['ngmix_EXP_E_1'])
        cat_res = tabletools.appendColumn(cat_res,'e2', cat_res['ngmix_EXP_E_2'])
        cat_res = tabletools.appendColumn(cat_res,'nbc_m1', (cat_res['ngmix_EXP_E_SENS_1']-1) )
        cat_res = tabletools.appendColumn(cat_res,'nbc_m2', (cat_res['ngmix_EXP_E_SENS_2']-1) )
        # cat_res = tabletools.appendColumn(cat_res,'nbc_c1', np.zeros(len(cat_res)))
        # cat_res = tabletools.appendColumn(cat_res,'nbc_c2', np.zeros(len(cat_res)))
        cat_res = tabletools.appendColumn(cat_res,'flags', cat_res['ngmix_FLAGS'])
        cat_res = tabletools.appendColumn(cat_res,'arate', cat_res['ngmix_EXP_ARATE'])
        cat_res = tabletools.appendColumn(cat_res,'T_s2n', cat_res['ngmix_EXP_T_S2N'])
        cat_res = tabletools.appendColumn(cat_res,'snr', cat_res['ngmix_EXP_S2N_W'])
        
        cat_res = tabletools.appendColumn(cat_res,'mean_psf_e1_sky', cat_res['ngmix_PSFREC_E_1'])
        cat_res = tabletools.appendColumn(cat_res,'mean_psf_e2_sky', cat_res['ngmix_PSFREC_E_2'])     

        cat_res['nbc_c1'] = 0
        cat_res['nbc_c2'] = 0



    if 'w' not in cat_res.dtype.names:

        if 'g_cov' in cat_res.dtype.names:

            e_cov_1_1 = cat_res['g_cov'][:,0]
            e_cov_2_2 = cat_res['g_cov'][:,1]
            e_cov_1_2 = cat_res['g_cov'][:,2]
            SN=0.32
            w = 1.0/( SN**2 + e_cov_1_1 + e_cov_2_2 + 2*e_cov_1_2 )
            cat_res = tabletools.appendColumn(cat_res,'w', w)
            warnings.warn('ngmix: using g_cov to create weight w')

        elif 'ngmix_EXP_E_COV_1_1' in cat_res.dtype.names:

            e_cov_1_1 = cat_res['ngmix_EXP_E_COV_1_1']
            e_cov_2_2 = cat_res['ngmix_EXP_E_COV_2_2']
            e_cov_1_2 = cat_res['ngmix_EXP_E_COV_1_2']
            SN=0.16
            w = 1.0/(2*SN**2 + e_cov_1_1 + e_cov_2_2 + 2*e_cov_1_2)
            cat_res = tabletools.appendColumn(cat_res,'w', w)

             

    cat_res = tabletools.appendColumn(cat_res,'coadd_objects_id', np.arange(len(cat_res)))
    cat_res = tabletools.appendColumn(cat_res,'ra_as' , cat_res['pars'][:,0])
    cat_res = tabletools.appendColumn(cat_res,'dec_as', cat_res['pars'][:,1])

    return cat_res

def rename_cols_truth(cat_tru):

    if 'shape_e1' in cat_tru.dtype.names:

        warnings.warn('adding true shape')

        e1,e2 = tt.distortion_to_shear(cat_tru['shape_e1'] , cat_tru['shape_e2'])
        e = e1+1j*e2
        e_rot = e*np.exp(2*cat_tru['rotation_angle']*1j)
        e1_rot = e_rot.real
        e2_rot = e_rot.imag


        warnings.warn('multiplying e1 intrinsic sign by %d' % config['methods'][args.method]['flip_g1'])
        e1_rot *= config['methods'][args.method]['flip_g1']

        e1_rot_sheared , e2_rot_sheared = tt.add_shear(e1_rot,e2_rot,cat_tru['g1_true'],cat_tru['g2_true'])

        # shape_e = (cat_tru['shape_e1'] + 1j*cat_tru['shape_e2'])*np.exp(1j*cat_tru['rotation_angle'])
        # xi = np.abs( shape_e )
        # ang  = np.angle( shape_e )
        # q = np.sqrt( (1-xi)/(1+xi) )
        # e = (1-q)/(1+q) * np.exp( ang * 1j)
        # e1 = e.real
        # e2 = e.imag

        cat_tru = tabletools.appendColumn(cat_tru,'intrinsic_e1',e1_rot)
        cat_tru = tabletools.appendColumn(cat_tru,'intrinsic_e2',e2_rot)
        cat_tru = tabletools.appendColumn(cat_tru,'sheared_e1',e1_rot_sheared)
        cat_tru = tabletools.appendColumn(cat_tru,'sheared_e2',e2_rot_sheared)

    elif 'intrinsic_g1' in  cat_tru.dtype.names:
        
        warnings.warn("e1 = -cat_tru['intrinsic_g1']")
        e1 = -cat_tru['intrinsic_g1'] 
        e2 =  cat_tru['intrinsic_g2'] 

        e = e1+1j*e2
        e_rot = e*np.exp(2*cat_tru['rotation_angle']*1j)
        e1_rot = e_rot.real
        e2_rot = e_rot.imag

        e1_rot_sheared , e2_rot_sheared = tt.add_shear(e1_rot,e2_rot,cat_tru['g1_true'],cat_tru['g2_true'])

        # shape_e = (cat_tru['shape_e1'] + 1j*cat_tru['shape_e2'])*np.exp(1j*cat_tru['rotation_angle'])
        # xi = np.abs( shape_e )
        # ang  = np.angle( shape_e )
        # q = np.sqrt( (1-xi)/(1+xi) )
        # e = (1-q)/(1+q) * np.exp( ang * 1j)
        # e1 = e.real
        # e2 = e.imag

        cat_tru = tabletools.appendColumn(cat_tru,'sheared_e1',e1_rot_sheared)
        cat_tru = tabletools.appendColumn(cat_tru,'sheared_e2',e2_rot_sheared)


    return cat_tru


def rename_im3shape_cols(cat_res,cat_tru=None,sim=None):

    global means_snr, means_rgp

    if means_snr==None:

        import cPickle as pickle
        filename_pickle = 'means_snr_rgp.cpickle'   
        pickle_dict = pickle.load(open(filename_pickle))
        means_num = pickle_dict['means_num']
        means_rgp = pickle_dict['means_rgp']/means_num
        means_snr = pickle_dict['means_snr']/means_num

    if 'im3shape_r_e1' in cat_res.dtype.names:

        cat_res = tabletools.appendColumn(cat_res,'e1', cat_res['im3shape_r_e1'])
        cat_res = tabletools.appendColumn(cat_res,'e2', cat_res['im3shape_r_e2'])
        cat_res = tabletools.appendColumn(cat_res,'ra_as', cat_res['im3shape_r_ra_as'])
        cat_res = tabletools.appendColumn(cat_res,'dec_as', cat_res['im3shape_r_dec_as'])
        cat_res = tabletools.appendColumn(cat_res,'radius', cat_res['im3shape_r_radius'])
        cat_res = tabletools.appendColumn(cat_res,'info_flag', cat_res['im3shape_r_info_flag'])
        cat_res = tabletools.appendColumn(cat_res,'error_flag', cat_res['im3shape_r_error_flag'])
        cat_res = tabletools.appendColumn(cat_res,'mean_psf_e1_sky', cat_res['im3shape_r_mean_psf_e1_sky'])
        cat_res = tabletools.appendColumn(cat_res,'mean_psf_e2_sky', cat_res['im3shape_r_mean_psf_e2_sky'])
        cat_res = tabletools.appendColumn(cat_res,'mean_rgpp_rp', cat_res['im3shape_r_mean_rgpp_rp'])
        cat_res = tabletools.appendColumn(cat_res,'mean_psf_fwhm', cat_res['im3shape_r_mean_psf_fwhm'])
        cat_res = tabletools.appendColumn(cat_res,'disc_flux', cat_res['im3shape_r_disc_flux'])
        cat_res = tabletools.appendColumn(cat_res,'bulge_flux', cat_res['im3shape_r_bulge_flux'])
        cat_res = tabletools.appendColumn(cat_res,'snr', cat_res['im3shape_r_snr'])


    if ('nbc_m' not in cat_res.dtype.names) & ( 'im3shape_r_w' in cat_res.dtype.names ) :
        
        cat_res = tabletools.appendColumn(cat_res,'w', cat_res['im3shape_r_w'])
        cat_res = tabletools.appendColumn(cat_res,'nbc_c1', cat_res['im3shape_r_nbc_c1'])
        cat_res = tabletools.appendColumn(cat_res,'nbc_c2', cat_res['im3shape_r_nbc_c2'])
        cat_res = tabletools.appendColumn(cat_res,'nbc_m', cat_res['im3shape_r_nbc_m'])


    if 'nbc_m' in cat_res.dtype.names:

        cat_res = tabletools.appendColumn(cat_res,'nbc_m1', cat_res['nbc_m'])
        cat_res = tabletools.appendColumn(cat_res,'nbc_m2', cat_res['nbc_m'])

    if 'w' not in cat_res.dtype.names:

        cat_res = tabletools.appendColumn(cat_res,'w', np.ones(len(cat_res)))
        warnings.warn('added weight column - all ones')

    # if cat_tru!=None:
    #     warnings.warn('using true psf e1 for mean_psf_e1_sky')
    #     if 'mean_psf_e1_sky' in cat_res.dtype.names:
    #         cat_res['mean_psf_e1_sky'] = cat_tru['psf_e1']
    #         cat_res['mean_psf_e2_sky'] = cat_tru['psf_e2']

    if 'nbc_m1' not in cat_res.dtype.names:
        cat_res = tabletools.appendColumn(cat_res,'nbc_c1', np.zeros(len(cat_res)))
        cat_res = tabletools.appendColumn(cat_res,'nbc_c2', np.zeros(len(cat_res)))
        cat_res = tabletools.appendColumn(cat_res,'nbc_m1', np.zeros(len(cat_res)))
        cat_res = tabletools.appendColumn(cat_res,'nbc_m2', np.zeros(len(cat_res)))


    if cat_tru!=None:
        if 'snr_avg' in config['cols_res']:

            col_snr = means_snr[cat_tru['id_cosmos'],cat_tru['id_psf']]
            col_rgp = means_rgp[cat_tru['id_cosmos'],cat_tru['id_psf']]
            cat_res = tabletools.appendColumn(cat_res,'snr_avg', col_snr)
            cat_res = tabletools.appendColumn(cat_res,'rgpp_rp_avg', col_rgp)



    return cat_res


def get_selection_des(selection_string,cols,n_files=30,get_calibrated=False):

    n_all = 0
    import glob
    if get_calibrated:
        filelist_i =glob.glob(config['filelist_des_calibrated'])
    else:
        filelist_i = glob.glob(config['filelist_des'])
    list_results = []
    logger.info('got %d DES files',len(filelist_i))
    n_files_loaded = 0
    n_gals_selected = 0
    for filename_des in filelist_i[:n_files]:
        try:
            cat_res=tt.load(filename_des,logger=0,remember=False)
        except Exception,errmsg:
            logger.error('file not loaded: %s, error: %s' % (filename_des,errmsg))
            continue

        if args.method == 'ngmix':
            cat_res = rename_ngmix_cols(cat_res)
        elif args.method == 'im3shape':
            cat_res = rename_im3shape_cols(cat_res)

        n_all+=len(cat_res)
        n_files_loaded+=1

        exec selection_string
        res_select = cat_res[select][cols]
        list_results.append(res_select)
        n_gals_selected += len(res_select)

        if n_files_loaded % 20 == 0:
            logger.info('loaded %d files, last %s, n_gals %d/%d %2.2f',n_files_loaded,filename_des,n_gals_selected,n_all,n_gals_selected/float(n_all))

    results = np.concatenate(list_results)
    logger.info('selected DES galaxies %d/%d %2.2f' , len(results),n_all,len(results)/float(n_all))

    if ('filepath_round_snr' in config) & ('round_snr' in cols):
        logger.info('adding round SNR')
        round_snr=tt.load(config['filepath_round_snr'])
        logger.info('running intersect')
        common_ids = np.intersect1d(round_snr['coadd_objects_id'],results['coadd_objects_id'])
        logger.info('running in1d')
        select_rsnr = np.in1d(common_ids,round_snr['coadd_objects_id'])
        logger.info('running in1d')
        select_results  = np.in1d(common_ids,results['coadd_objects_id'])
        logger.info('selectiong intersection')
        results = results[select_results]
        round_snr = round_snr[select_rsnr]
        logger.info('sorting')
        results.sort(order='coadd_objects_id')
        round_snr.sort(order='coadd_objects_id')
        logger.info('adding col')
        results = tt.add_col(results,'round_snr',round_snr['round_snr'])

    return results

def get_selection_sim(selection_string, cols_res, cols_tru, get_calibrated=False, short=False):

    if short==False:
        list_all_res, list_all_tru = get_selection_split(selection_string, cols_res, cols_tru, get_calibrated)
        all_tru = np.concatenate(list_all_tru)
        all_res = np.concatenate(list_all_res)
    else:
        if get_calibrated:
            results_filename_fmt = config['methods'][args.method]['filename_calibrated']  
        else:   
            results_filename_fmt = config['methods'][args.method]['filename_results']  



    return all_res, all_tru


def get_selection_split(selection_string, cols_res, cols_tru,get_calibrated=False):

    if get_calibrated:
        results_filename_fmt = config['methods'][args.method]['filename_calibrated']  
    else:   
        results_filename_fmt = config['methods'][args.method]['filename_results']  

    truth_filename_fmt = config['filename_truth']  
    list_shears = []
    list_all_res = []
    list_all_tru = []

    n_all_loaded=0
    n_all_selected=0

    ia=0
    n_missing=0

    for ig,vg in enumerate(config['shear']):

        list_results = []

        id_first = args.first
        id_last = id_first + args.num
        # if id_last > 200: 
        #     id_last=200;
        #     warnings.warn('hard coded max number of files 200')

        for ip in range(id_first,id_last):
         
            filename_tru = truth_filename_fmt % (ip,ig)
            filename_res = results_filename_fmt % (ip,ig)
            try:
                    cat_tru_all = tabletools.loadTable(filename_tru,log=1,remember=False)
                    cat_tru = cat_tru_all
                    logger.debug('loaded %05d galaxies from file: %s' % (len(cat_tru_all),filename_tru))

            except Exception,errmsg:
                logger.error('file %s : %s' % (filename_tru,errmsg) )
                continue


            try:
                    cat_res_all = tabletools.loadTable(filename_res,log=1,remember=False)
                    cat_res = cat_res_all
                    logger.debug('loaded %05d galaxies from file: %s' % (len(cat_res_all),filename_res))

            except Exception,errmsg:
                logger.debug('sth wrong with file %s errmsg %s' % (filename_res,errmsg) )
                n_missing+=1
                continue

 

            if ('e1' in cols_res) & (args.method=='im3shape'):
                # cat_res['e1'] = cat_res['e1']*config['methods'][args.method]['flip_g1']
                cat_tru['g1_true'] = -1*cat_tru['g1_true']
                warnings.warn('flipping g1 in truth cat for method %s'%args.method)


            if len(cat_tru) != len(cat_res):
                cat_tru=cat_tru[cat_res['coadd_objects_id']]

            if args.method == 'ngmix':
                cat_res = rename_ngmix_cols(cat_res,cat_tru)
            elif args.method == 'im3shape':
                cat_res = rename_im3shape_cols(cat_res,cat_tru)
            
            cat_tru = rename_cols_truth(cat_tru)

            for col in cols_tru:
                if col not in cat_tru.dtype.names:
                    raise Exception('column %s not found in truth catalog %s' % (col,filename_tru))
            for col in cols_res:
                if col not in cat_res.dtype.names:
                    raise Exception('column %s not found in results catalog %s' % (col,filename_res))

            n_all_loaded+=len(cat_res)
            try:
                exec selection_string
            except Exception,errmsg:
                print errmsg
                import pdb; pdb.set_trace()

            n_gal_after_cuts = len(np.nonzero(select)[0])
            n_all_selected += n_gal_after_cuts
            if  n_gal_after_cuts < 1:
                logger.debug('select didnt give any results %s, n_gal=%d' % (selection_string, len(cat_res)) )
                # import pdb; pdb.set_trace()
                # raise Exception('select didnt give any results %s' % selection_string)
            # else:
                # logger.debug('cuts %d/%d' , n_gal_after_cuts, len(cat_res))


            if len(cat_res) != len(cat_tru):
                logger.error('ip=%d ig=%d uneven number of rows %d %d' % (ip , ig , len(cat_res),len(cat_tru)))
                continue
            
            try:
                selected_res = cat_res[select][cols_res]
                selected_tru = cat_tru[select][cols_tru]
            except Exception,errmsg:
                logger.error(errmsg)
                import pdb;pdb.set_trace()

            
            # selected_res = selected_res.astype( ['f4']*len(cols_res) )
            selected_tru = selected_tru.astype( dtype={ 'formats':['f4']*len(cols_tru) , 'names': cols_tru })

            ia+=1
            list_all_res.append(selected_res)
            list_all_tru.append(selected_tru)

            if ia % 100 == 0:
                logger.info('loaded %5d files, last %s, n_gals_selected %d/%d %f' % (ia,filename_res,n_all_selected,n_all_loaded,float(n_all_selected)/float(n_all_loaded)))
  
    n_total = 0
    for res in list_all_res:
        n_total += len(res)

        
    logger.info('selected %d parts with average %d, total %d/%d loaded, n_missing_files %d' % (len(list_all_res) , float(n_total)/float(len(list_all_res)), n_total , n_all_loaded, n_missing) )
    return list_all_res, list_all_tru