import numpy as np; import pylab as pl
import  sys, logging, yaml, argparse, time, copy, itertools, fitting, warnings
warnings.simplefilter("once")
sys.path.append('/home/tomek/code/tktools')
import tabletools, plotstools
# from nbc2_dtypes import *
logging_level = logging.INFO; logger = logging.getLogger("nbc-v7"); logger.setLevel(logging_level)  
log_formatter = logging.Formatter("%(asctime)s %(name)s %(levelname)s   %(message)s", "%Y-%m-%d %H:%M:%S")
stream_handler = logging.StreamHandler(sys.stdout); stream_handler.setFormatter(log_formatter)
if logger.handlers == [] : logger.addHandler(stream_handler); logger.propagate = False

def rename_ngmix_cols(cat_res):

    if 'g' in cat_res.dtype.names:

        cat_res = tabletools.appendColumn(cat_res,'e1', cat_res['g'][:,0])
        cat_res = tabletools.appendColumn(cat_res,'e2', cat_res['g'][:,1])

    if 'g_sens' in cat_res.dtype.names:

        cat_res = tabletools.appendColumn(cat_res,'nbc_m1', (cat_res['g_sens'][:,0]-1) )
        cat_res = tabletools.appendColumn(cat_res,'nbc_m2', (cat_res['g_sens'][:,1]-1) )
        cat_res = tabletools.appendColumn(cat_res,'nbc_c1', np.zeros(len(cat_res)))
        cat_res = tabletools.appendColumn(cat_res,'nbc_c2', np.zeros(len(cat_res)))

    if 'ngmix009_EXP_E_1' in cat_res.dtype.names:

        cat_res = tabletools.appendColumn(cat_res,'e1', cat_res['ngmix009_EXP_E_1'])
        cat_res = tabletools.appendColumn(cat_res,'e2', cat_res['ngmix009_EXP_E_2'])
        cat_res = tabletools.appendColumn(cat_res,'nbc_m1', (cat_res['ngmix009_EXP_E_SENS_1']-1) )
        cat_res = tabletools.appendColumn(cat_res,'nbc_m2', (cat_res['ngmix009_EXP_E_SENS_2']-1) )
        cat_res = tabletools.appendColumn(cat_res,'nbc_c1', np.zeros(len(cat_res)))
        cat_res = tabletools.appendColumn(cat_res,'nbc_c2', np.zeros(len(cat_res)))
        cat_res = tabletools.appendColumn(cat_res,'flags', cat_res['ngmix009_FLAGS'])
        cat_res = tabletools.appendColumn(cat_res,'arate', cat_res['ngmix009_EXP_ARATE'])

    if 'w' not in cat_res.dtype.names:

        if 'g_cov' in cat_res.dtype.names:

            e_cov_1_1 = cat_res['g_cov'][:,0,0]
            e_cov_2_2 = cat_res['g_cov'][:,1,1]
            e_cov_1_2 = cat_res['g_cov'][:,0,1]
            SN=0.16
            w = 1.0/(2*SN**2 + e_cov_1_1 + e_cov_2_2 + 2*e_cov_1_2)
            cat_res = tabletools.appendColumn(cat_res,'w', w)

        elif 'ngmix009_EXP_E_COV_1_1' in cat_res.dtype.names:

            e_cov_1_1 = cat_res['ngmix009_EXP_E_COV_1_1']
            e_cov_2_2 = cat_res['ngmix009_EXP_E_COV_2_2']
            e_cov_1_2 = cat_res['ngmix009_EXP_E_COV_1_2']
            SN=0.16
            w = 1.0/(2*SN**2 + e_cov_1_1 + e_cov_2_2 + 2*e_cov_1_2)
            cat_res = tabletools.appendColumn(cat_res,'w', w)

              

    cat_res = tabletools.appendColumn(cat_res,'mean_rgpp_rp', np.random.uniform(1,3))
    cat_res = tabletools.appendColumn(cat_res,'mean_psf_e1_sky', np.random.uniform(-0.05,0.05))
    cat_res = tabletools.appendColumn(cat_res,'mean_psf_e2_sky', np.random.uniform(-0.05,0.05))
    cat_res = tabletools.appendColumn(cat_res,'snr', np.random.uniform(0.,100))
    cat_res = tabletools.appendColumn(cat_res,'coadd_objects_id', np.arange(len(cat_res)))

    return cat_res

def rename_im3shape_cols(cat_res):

    if 'nbc_m' in cat_res.dtype.names:

        cat_res = tabletools.appendColumn(cat_res,'nbc_m1', cat_res['nbc_m'])
        cat_res = tabletools.appendColumn(cat_res,'nbc_m2', cat_res['nbc_m'])

    if 'w' not in cat_res.dtype.names:

        cat_res = tabletools.appendColumn(cat_res,'w', np.ones(len(cat_res)))

    return cat_res


def get_selection_des(selection_string,cols,n_files=30,get_calibrated=False):

    n_all = 0
    if get_calibrated:
        filelist_i = np.loadtxt(config['filelist_des_calibrated'],dtype='a1024')
    else:
        filelist_i = np.loadtxt(config['filelist_des'],dtype='a1024')
    list_results = []
    for filename_des in filelist_i[:n_files]:
        cat_res=tabletools.loadTable(filename_des,log=logger,remember=False)
        if args.method == 'ngmix':
            cat_res = rename_ngmix_cols(cat_res)
        elif args.method == 'im3shape':
            cat_res = rename_im3shape_cols(cat_res)

        n_all+=len(cat_res)

        exec selection_string
        res_select = cat_res[select][cols]
        list_results.append(res_select)

    results = np.concatenate(list_results)
    logger.info('selected DES galaxies %d/%d %2.2f' , len(results),n_all,len(results)/float(n_all))

    return results

def get_selection_sim(selection_string, cols_res, cols_tru, get_calibrated=False):

    list_all_res, list_all_tru = get_selection_split(selection_string, cols_res, cols_tru, get_calibrated)
    all_tru = np.concatenate(list_all_tru)
    all_res = np.concatenate(list_all_res)

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

    ia=0
    n_missing=0

    for ig,vg in enumerate(config['shear']):
        
        list_results = []

        id_first = args.first
        id_last = id_first + args.num
        if id_last > 200: 
            id_last=200;
            warnings.warn('hard coded max number of files 200')

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

            for col in cols_tru:
                if col not in cat_tru.dtype.names:
                    raise Exception('column %s not found in truth catalog %s' % (col,filename_tru))

            try:
                    cat_res_all = tabletools.loadTable(filename_res,log=1,remember=False)
                    cat_res = cat_res_all
                    logger.debug('loaded %05d galaxies from file: %s' % (len(cat_res_all),filename_res))
                    if args.method == 'ngmix':
                        cat_res = rename_ngmix_cols(cat_res)
                    elif args.method == 'im3shape':
                        cat_res = rename_im3shape_cols(cat_res)


            except Exception,errmsg:
                logger.debug('sth wrong with file %s errmsg %s' % (filename_res,errmsg) )
                n_missing+=1
                continue
 
            for col in cols_res:
                if col not in cat_res.dtype.names:
                    raise Exception('column %s not found in results catalog %s' % (col,filename_res))

            if ('e1' in cols_res) & (args.method=='im3shape'):
                # cat_res['e1'] = cat_res['e1']*config['methods'][args.method]['flip_g1']
                cat_tru['g1_true'] = -1*cat_tru['g1_true']
                warnings.warn('flipping g1 in truth cat for method %s'%args.method)


            if len(cat_tru) != len(cat_res):
                cat_tru=cat_tru[cat_res['coadd_objects_id']]

            n_all_loaded+=len(cat_res)
            try:
                exec selection_string
            except Exception,errmsg:
                print errmsg
                import pdb; pdb.set_trace()

            if len(np.nonzero(select)[0]) < 1:
                logger.debug('select didnt give any results %s' % selection_string)
                # import pdb;pdb.set_trace()
                # raise Exception('select didnt give any results %s' % selection_string)


            if len(cat_res) != len(cat_tru):
                logger.error('ip=%d ig=%d uneven number of rows %d %d' % (ip , ig , len(cat_res),len(cat_tru)))
                continue
            
            try:
                selected_res = cat_res[select][cols_res]
                selected_tru = cat_tru[select][cols_tru]
            except:
                import pdb;pdb.set_trace()

            
            # selected_res = selected_res.astype( ['f4']*len(cols_res) )
            selected_tru = selected_tru.astype( dtype={ 'formats':['f4']*len(cols_tru) , 'names': cols_tru })

            ia+=1
            list_all_res.append(selected_res)
            list_all_tru.append(selected_tru)

            if ia % 100 == 0:
                logger.info('loaded %5d files, last %s' % (ia,filename_res))
  
    n_total = 0
    for res in list_all_res:
        n_total += len(res)

        
    logger.info('selected %d parts with average %d, total %d/%d loaded, n_missing_files %d' % (len(list_all_res) , float(n_total)/float(len(list_all_res)), n_total , n_all_loaded, n_missing) )
    return list_all_res, list_all_tru