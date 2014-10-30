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

def get_selection_des(selection_string,cols,n_files=30):

    filelist_i = np.loadtxt('filelist-v7-r.txt',dtype='a1024')
    list_results = []
    for filename_des in filelist_i[:n_files]:
        cat_res=tabletools.loadTable(filename_des,log=logger,remember=False)
        name_tile=filename_des.replace('/Users/tomek/data/DES/im3shape-v7/','').replace('-r.fits.gz','')       
        exec selection_string
        res_select = cat_res[select][cols]
        list_results.append(res_select)

    results = np.concatenate(list_results)
    logger.info('selected DES galaxies %d' , len(results))

    return results

def get_selection_sim(selection_string, cols_res, cols_tru):

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

    ia=0

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
                    logger.debug('loaded %05d galaxies from file: %s' % (len(res_tru_all),filename_tru))

            except:
                logger.error('file %s not found' % filename_tru )
                continue

            for col in cols_tru:
                if col not in cat_tru.dtype.names:
                    raise Exception('column %s not found in truth catalog %s' % (col,filename_tru))

            try:
                    cat_res_all = tabletools.loadTable(filename_res,log=1,remember=False)
                    cat_res = cat_res_all
                    logger.debug('loaded %05d galaxies from file: %s' % (len(cat_res_all),filename_res))

            except Exception,errmsg:
                logger.error('sth wrong with file %s errmsg %s' % (filename_res,errmsg) )
                continue
 
            for col in cols_res:
                if col not in cat_res.dtype.names:
                    raise Exception('column %s not found in results catalog %s' % (col,filename_res))

            if 'e1' in cols_res:
                cat_res['e1'] = cat_res['e1']*config['methods'][args.method]['flip_g1']
                warnings.warn('flipping g1 for method %s'%args.method)

            if len(cat_tru) != len(cat_res):
                cat_tru=cat_tru[cat_res['coadd_objects_id']]

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
  
    n_total = 0
    for res in list_all_res:
        n_total += len(res)

        
    logger.info('selected %d parts with average %d, total %d' % (len(list_all_res) , float(n_total)/float(len(list_all_res)), n_total) )
    return list_all_res, list_all_tru