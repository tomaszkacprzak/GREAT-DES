import numpy as np; import pylab as pl; import tktools as tt;
import  sys, logging, yaml, argparse, time, copy, itertools, tktools, warnings, os, fitsio, pyfits;
warnings.simplefilter("once")
sys.path.append('/home/tomek/code/tktools');
sys.path.append('/Users/tomek/code/tktools');
logging_level = logging.INFO; logger = logging.getLogger("nbc-v7"); logger.setLevel(logging_level)  
log_formatter = logging.Formatter("%(asctime)s %(name)s %(levelname)s   %(message)s", "%Y-%m-%d %H:%M:%S")
stream_handler = logging.StreamHandler(sys.stdout); stream_handler.setFormatter(log_formatter)
if logger.handlers == [] : logger.addHandler(stream_handler); logger.propagate = False
import nbc_v7_select, nbc_v7_stats
import plotstools, fitting

# info_vals = [1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536,131072,262144,524288,1048576]

# for im3shape-v8 catalogs: missing 4 which for v8 is the v7 connection
# snr > 10 32768
# rgpp_rp > 1.15 32

def get_params_covariance():

    res_sim, res_tru, res_des = load_selection()

    # abse= np.abs(res_sim['e1']+1j*res_sim['e2'])
    # X=np.concatenate([res_sim['snr'][:,None],res_sim['round_snr'][:,None],res_sim['e1'][:,None],res_sim['e2'][:,None],abse[:,None],res_sim['mean_rgpp_rp'][:,None],res_sim['radius'][:,None]],axis=1)
    # logger.info('calculating cov')
    # logger.info('shape of matrix %s',str(X.shape))
    # import pdb; pdb.set_trace()
    # C=np.cov(X.T)
    # pl.pcolormesh(np.log10(np.abs(C))); 
    # pl.colorbar()
    # pl.xticklabels()
    # pl.show()



def plot_add_requirements(level=0.01,target=1.0,mult=1.):

    corner = pl.xlim()[0]
    length = abs(pl.xlim()[1]) + abs(pl.xlim()[0])
    pl.gca().add_patch(pl.Rectangle(  (corner, target-mult*level), length , 2*mult*level , facecolor = '0.6' , edgecolor='none', alpha=0.5 ))


def polynomial_basis(x):

    xx = 1-x+0.05
    X = np.concatenate( [ 1/(xx)**0,1/(xx)**1,1./(xx)**2,1./(xx)**3] ,axis=1 ) 
    # X = np.concatenate( [ np.log(xx+0.1)*np.exp(-xx**2) ] ,axis=1 ) 

    return X



def inv_snr_basis(x):

    n_points_x = x.shape[0]
    # X = np.concatenate( [ np.ones((n_points_x,1)), 1./x, 1./x**2, 1./x**4],axis=1 )
    # X = np.concatenate( [  1/x, 1/x**1.5, 1./x**2, 1/x**2.5, 1./x**3, 1./x**4 ],axis=1 )
    # X = np.concatenate( [np.ones((n_points_x,1)),   1/x, 1./x**2, 1./x**3 ],axis=1 )
    # X = np.concatenate( [ np.ones((n_points_x,1)), 1/x, 1./x**1.5, 1./x**2, 1./x**2.5 ],axis=1 )
    
    # use this one for real calibration
    # X = np.concatenate( [  np.ones((n_points_x,1)), 1./x**1, 1/x**2, 1/x**3 ],axis=1 ) 

    # X = np.concatenate( [  1/x**1, 1/x**1.5, 1/x**2],axis=1 ) 
    basis_list = eval(config['inv_snr_basis_m'])
    X = np.concatenate( basis_list ,axis=1 ) 

    # use this one for figure
    # X = np.concatenate( [ 1./x**2,  1./x**3, 1/x**4],axis=1 ) 
    return X


def inv_snr_basis2(x):

    n_points_x = x.shape[0]
    # X = np.concatenate( [ np.ones((n_points_x,1)), 1./x, 1./x**2, 1./x**4],axis=1 )
    # X = np.concatenate( [ np.ones((n_points_x,1)), 1./x**2 ],axis=1 )
    # X = np.concatenate( [ 1./x, 1./x**1.25, 1./x**1.5, 1./x**2],axis=1 )
    # X = np.concatenate( [ 1/x**2, 1/x**4 ],axis=1 )
    # X = 1./x**2

    basis_list = eval(config['inv_snr_basis_a'])
    X = np.concatenate( basis_list ,axis=1 ) 

    # use this one for figure
    # X = np.concatenate( [ 1./x**2,  1./x**3, 1/x**4],axis=1 ) 
    return X


    return X

def apply_calibration_selection():

    res_sim, res_tru, res_des = load_selection()


    logger.info('applying calibration to DES')
    res_des_cal , calib_m, calib_a=get_calibration_columns(res_des)
    logger.info('done! mean(m) = % 2.4f mean(a) = % 2.4f n(m=0)=%d / %d %2.2f'  % (np.mean(calib_m),np.mean(calib_a),np.sum(calib_m==0),len(calib_m),np.sum(calib_m==0)/float(len(calib_m)) ))

    # entire sample
    logger.info('calculating bias for the entire sample - using cuts from final selection')
    logger.info(selection_string_final_des)
    cat_res = res_des
    exec selection_string_final_des
    select_des = select.copy()

    select_des = select_des & (~np.isnan(res_des['nbc_c1'])) & (~np.isnan(res_des['nbc_c2']))
    res_des_select = res_des_cal[select_des]
    logger.info('raw mean_e1 = %2.5f +/- %2.5f' % ( np.mean(res_des_select['e1']) , np.std( res_des_select['e1']) / np.sqrt(len(res_des_select)) ) )
    logger.info('raw mean_e2 = %2.5f +/- %2.5f' % ( np.mean(res_des_select['e2']) , np.std( res_des_select['e2']) / np.sqrt(len(res_des_select)) ) )

    logger.info('nbc mean_e1 = %2.5f +/- %2.5f' % ( np.mean(res_des_select['e1']-res_des_select['nbc_c1']) , np.std( res_des_select['e1']-res_des_select['nbc_c1']) / np.sqrt(len(res_des_select)) ) )
    logger.info('nbc mean_e2 = %2.5f +/- %2.5f' % ( np.mean(res_des_select['e2']-res_des_select['nbc_c2']) , np.std( res_des_select['e2']-res_des_select['nbc_c2']) / np.sqrt(len(res_des_select)) ) )
    logger.info('mean nbc_m  = %2.5f ' % ( np.mean( res_des_select['nbc_m'] ) ) )
    logger.info('mean nbc_a  = %2.5f ' % ( np.mean( res_des_select['nbc_alpha'] ) ) ) 
    logger.info('mean nbc_c1 = %2.5f ' % ( np.mean( res_des_select['nbc_c1'] ) ) )
    logger.info('mean nbc_c2 = %2.5f ' % ( np.mean( res_des_select['nbc_c2'] ) ) )

    import pdb; pdb.set_trace()
    logger.info('getting weights for DES')
    res_des_cal = get_weight_column(res_des_cal)
    logger.info('done! mean sigma = %2.3f' % (np.median(np.sqrt(1./res_des_cal['w']))))

    logger.info('applying calibration to SIM')
    res_sim_cal, calib_m, calib_a=get_calibration_columns(res_sim)
    logger.info('done! mean(m) = % 2.4f mean(a) = % 2.4f n(m=0)=%d / %d %2.2f'  % (np.mean(calib_m),np.mean(calib_a),np.sum(calib_m==0),len(calib_m),np.sum(calib_m==0)/float(len(calib_m))))

    z_bins = [ 0.3, 0.644, 0.901, 1.3 ]
    select = (res_tru['zphot'] > 0.3) & ((res_tru['zphot'] < 0.644)) & ((res_sim['snr'] > 15)) & ((res_sim['mean_rgpp_rp'] > 1.15))
    logger.info('zbin=1 mean nbc_a=%2.3f' % np.mean(res_sim_cal[select]['nbc_alpha']))
    select = (res_tru['zphot'] > 0.644) & ((res_tru['zphot'] < 0.901)) & ((res_sim['snr'] > 15)) & ((res_sim['mean_rgpp_rp'] > 1.15))
    logger.info('zbin=2 mean nbc_a=%2.3f' % np.mean(res_sim_cal[select]['nbc_alpha']))
    select = (res_tru['zphot'] > 0.901) & ((res_tru['zphot'] < 1.3)) & ((res_sim['snr'] > 15)) & ((res_sim['mean_rgpp_rp'] > 1.15))
    logger.info('zbin=3 mean nbc_a=%2.3f' % np.mean(res_sim_cal[select]['nbc_alpha']))

    logger.info('getting weights for SIM')
    res_sim_cal = get_weight_column(res_sim_cal)
    logger.info('done! mean sigma = %2.3f' % (np.median(np.sqrt(1./res_sim_cal['w']))))

    filename = '%s/res_sim.fits' % args.output_dir
    tktools.save(filename, res_sim_cal, clobber=True)
    filename = '%s/res_tru.fits' % args.output_dir
    tktools.save(filename, res_tru, clobber=True)
    filename = '%s/res_des.fits' % args.output_dir
    tktools.save(filename, res_des_cal, clobber=True)



    pl.figure()
    pl.hist(res_des_cal['nbc_m'],np.linspace(-0.5,1.1,100),histtype='step');
    pl.xlabel('nbc_m')
    pl.figure()
    pl.hist(res_des_cal['nbc_c1'],np.linspace(-0.01,0.01,100),histtype='step');
    pl.xlabel('nbc_c1')
    pl.figure()
    pl.hist(res_des_cal['nbc_c2'],np.linspace(-0.01,0.01,100),histtype='step');
    pl.xlabel('nbc_c2')
    pl.figure()
    pl.hist(res_des_cal['nbc_alpha'],np.linspace(-0.5,0.5,100),histtype='step');
    pl.xlabel('nbc_alpha')
    pl.figure()
    pl.hist2d(res_des_cal['snr'],res_des_cal['nbc_m'],[np.linspace(0,100,100),np.linspace(-0.4,0.1,100)]);
    pl.xlabel('snr')
    pl.xlabel('nbc_m')
    pl.figure()
    pl.hist2d(res_des_cal['snr'],res_des_cal['nbc_alpha'],[np.linspace(0,100,100),np.linspace(-0.4,0.1,100)]);
    pl.xlabel('snr')
    pl.xlabel('nbc_alpha')
    pl.show()





def apply_calibration_to_file():

    res_sim=tt.load(args.filename_to_calibrate)
    logger.info('applying calibration to %s',args.filename_to_calibrate)
    res_sim_cal, calib_m, calib_a=get_calibration_columns(res_sim)
    res_sim_cal = get_weight_column(res_sim_cal)
    filename_out = args.filename_to_calibrate.replace('.fits','.nbc.fits')
    tktools.save(filename_out, res_sim_cal, clobber='skip')

    # pl.figure()
    # pl.hist(res_sim_cal['nbc_m'],np.linspace(-0.5,1.1,100),histtype='step');
    # pl.xlabel('nbc_m')
    # pl.figure()
    # pl.hist(res_sim_cal['nbc_c1'],np.linspace(-0.01,0.01,100),histtype='step');
    # pl.xlabel('nbc_c1')
    # pl.figure()
    # pl.hist(res_sim_cal['nbc_c2'],np.linspace(-0.01,0.01,100),histtype='step');
    # pl.xlabel('nbc_c2')
    # pl.figure()
    # pl.hist(res_sim_cal['nbc_alpha'],np.linspace(-0.5,0.5,100),histtype='step');
    # pl.xlabel('nbc_alpha')
    # pl.show()

def apply_calibration_des():
    
    import glob
    filelist_r = glob.glob(config['filelist_des'])

    logger.info('found %d files in %s',len(filelist_r),config['filelist_des'])

    n_files = 10000
    n_all = 0
    n_calibrated = 0

    dirpath_calib = os.path.dirname(config['filelist_des_calibrated'])
    if not os.path.exists(dirpath_calib): 
        os.makedirs(dirpath_calib)
        logger.info('created directory: %s',dirpath_calib)
    else:
        logger.info('using directory: %s',dirpath_calib)

    logger.info('calibrating %s',os.path.dirname(filelist_r[0]))

    for filename_des in filelist_r[:n_files]:

        filename_tile = os.path.basename(filename_des)
        filename_calibrated = os.path.join(dirpath_calib,filename_tile)
        if os.path.isfile(filename_calibrated): 
            logger.info('file exists, skipping %s', filename_calibrated)
            continue
            
        # cat_res=tabletools.loadTable(filename_des,log=1)
        cat_res=pyfits.getdata(filename_des)
        n_all+=len(cat_res)
        
        cat_res, calib_m, calib_a=get_calibration_columns(cat_res)
        cat_res=get_weight_column(cat_res)

        # optionally remove calibration from galaxies which do not make the cut
        # exec selection_string_des
        # cat_res['nbc_m'][~select]=0
        # cat_res['nbc_c1'][~select]=0
        # cat_res['nbc_c2'][~select]=0

        fitsio.write(filename_calibrated,cat_res,clobber=False)
        logger.info('saved %s',filename_calibrated)

        # if np.any(select):
        #     mean_m  = np.mean(cat_res['nbc_m'][select])
        #     mean_c1 = np.mean(cat_res['nbc_c1'][select])
        #     mean_c2 = np.mean(cat_res['nbc_c2'][select])
        #     mean_a = np.mean(calib_a)
        #     logger.info('wrote %s with mean(m)=%2.4f mean(a)=%2.4f mean(c1)=%2.4f mean(c2)=%2.4f' % (filename_calibrated,mean_m,mean_a,mean_c1,mean_c2))
        # else:
        #     logger.info('wrote %s but there are no selected galaxies with select: %s' % (filename_calibrated,selection_string_des))
        #     mean_m = 0
        #     mean_c1 = 0
        #     mean_c2 = 0
        #     mean_a = 0
        # if any(~np.isfinite([mean_m,mean_c1,mean_c2,mean_a])):
        #     logger.error('not finite value in calibration')
        #     import pdb; pdb.set_trace()


        n_calibrated += len(cat_res)

    logger.info('calibrated %d galaxies out of %d' % (n_calibrated,n_all))


def apply_calibration_sim():

    n_files = 10000
    n_all = 0
    n_calibrated = 0
    n_calibrated_nonzero_tot = 0

    dirpath_calib = os.path.dirname(config['methods'][args.method]['filename_calibrated'])
    if not os.path.exists(dirpath_calib): 
        os.makedirs(dirpath_calib)
        logger.info('created directory: %s',dirpath_calib)
    else:
        logger.info('using directory: %s',dirpath_calib)

    logger.info('calibrating GREAT-DES')
    import glob
    filelist_des = glob.glob( os.path.dirname(config['methods'][args.method]['filename_results']) + '/nbc.meds*fits')
    logger.info('found %d results files' , len(filelist_des))
    for filename_des in filelist_des:


        filename_tile = os.path.basename(filename_des)
        filename_calibrated = os.path.join(dirpath_calib,filename_tile)
        if os.path.isfile(filename_calibrated): 
            logger.info('file exists, skipping %s', filename_calibrated)
            continue
        
        cat_res=tktools.load(filename_des,logger=1)
        filename_truth = filename_des.replace('.meds.','.truth.').replace('/main_cats/','/data/')
        cat_tru=tktools.load(filename_truth,logger=1)
        n_all+=len(cat_res)
        
        cat_res, m, a=get_calibration_columns(cat_res)
        cat_res=get_weight_column(cat_res)

        cat_tru = cat_tru[cat_res['coadd_objects_id']]

        exec "select="+selection_string_final_sim
        cat_res[~select]['nbc_m']=0
        cat_res[~select]['nbc_c1']=0
        cat_res[~select]['nbc_c2']=0
        n_calibrated_this = len(np.nonzero(select)[0])
        n_calibrated += n_calibrated_this

        tktools.save(filename_calibrated,cat_res,clobber=True)

        mean_m  = np.mean(cat_res['nbc_m'][select])
        mean_c1 = np.mean(cat_res['nbc_c1'][select])
        mean_c2 = np.mean(cat_res['nbc_c2'][select])
        mean_a = np.mean(a[select])

        select = (cat_res['nbc_m']>1e-5) | (cat_res['nbc_m']<-1e-5)
        n_calibrated_nonzero = len(np.nonzero(select)[0])
        n_calibrated_nonzero_tot += n_calibrated_nonzero
        logger.info('wrote %s with mean(m)=%2.4f mean(a)=%2.4f mean(c1)=%2.4f mean(c2)=%2.4f n_calibrated_nonzero=%d/%d/%d' % (filename_calibrated,mean_m,mean_a,mean_c1,mean_c2,n_calibrated_nonzero,n_calibrated_this,len(cat_res)))

    logger.info('calibrated %d galaxies out of %d' % (n_calibrated,n_all))


def get_bins_centers(bins_edges,constant_spacing=True):

    # ensure np array
    bins_edges=np.array(bins_edges)
    bins_centers=np.zeros(len(bins_edges)-1)

    for be in range(len(bins_edges)-1): bins_centers[be] = np.mean([bins_edges[be],bins_edges[be+1]])      

    return bins_centers

def add_col(rec, name, arr, dtype=None):
    import numpy
    arr = numpy.asarray(arr)
    if dtype is None:
        dtype = arr.dtype
    newdtype = numpy.dtype(rec.dtype.descr + [(name, dtype)])
    newrec = numpy.empty(rec.shape, dtype=newdtype)
    for field in rec.dtype.fields:
        newrec[field] = rec[field]
    newrec[name] = arr
    return newrec


def get_weight_column(cat):

    logger.info('applying weights ')
    filename_table = 'weights_interpolation_table.v2.cpickle'
    file_pickle = open(filename_table)
    import cPickle as pickle
    sigma_e, sigma_e_max, list_snr_edges, list_rgp_edges , vec_snr_hires, vec_rgp_hires, sigma_e_hires, n_gal, stdd0_e = pickle.load(file_pickle)
    logger.info('loaded %s'%filename_table)

    # set the noisy values to one high snr value
    # sigma_e[sigma_e<sigma_e_max] = sigma_e_max
    sigma_greater = sigma_e.copy()
    sigma_greater[stdd0_e>sigma_e]=stdd0_e[stdd0_e>sigma_e]
    sigma_e_max = min(sigma_greater.flatten())
    logger.info('sigma_e_max=%2.4f'%sigma_e_max)

    vec_snr = get_bins_centers(list_snr_edges)
    vec_rgp = get_bins_centers(list_rgp_edges)

    # logger.debug('hires')
    # vec_snr_hires = np.linspace(vec_snr.min(),vec_snr.max(),100)
    # vec_rgp_hires = np.linspace(vec_rgp.min(),vec_rgp.max(),100)

    import scipy.interpolate
    X1,X2 = np.meshgrid(vec_snr_hires,vec_rgp_hires)
    interp_method = 'linear'
    func_interp = scipy.interpolate.interp2d(vec_snr,vec_rgp,sigma_greater.T, kind=interp_method)
    sigma_e_hires = func_interp(vec_snr_hires,vec_rgp_hires)

    if 'snr'  in cat.dtype.names: 
        size_col = 'mean_rgpp_rp'
        snr_col = 'snr'
    else:
        snr_col = 'im3shape_r_snr'
        size_col = 'im3shape_r_mean_rgpp_rp'

    cat[snr_col][~np.isfinite(cat[snr_col])]=0
    cat[size_col][~np.isfinite(cat[size_col])]=0

    col_sig = scipy.interpolate.griddata((X1.flatten(),X2.flatten()),sigma_e_hires.flatten(),(cat[snr_col],cat[size_col]),method='linear',rescale=True,fill_value=0)

    col_sig[cat[snr_col]>list_snr_edges.max()] = sigma_e_max
    col_sig[cat[size_col]>list_rgp_edges.max()] = sigma_e_max
    col_sig[col_sig<sigma_e_max] = sigma_e_max

    col_w = 1./col_sig**2
    col_w[np.isnan(col_w)] = 0
    col_w[np.isinf(col_w)] = 0
    if np.any(np.isnan(col_w)): logger.warning('nans')
    if np.any(np.isinf(col_w)): logger.warning('infs')

    sigma_e_min = 0.24
    max_w = 1/sigma_e_min**2
    logger.info('setting maximum weight to %2.2f, sigma=%2.2f' % (max_w,sigma_e_min))
    col_w[col_w>max_w] = max_w

    if 'w' in cat.dtype.names:
        warnings.warn('overwriting existing column with weights')
        cat['w'] = col_w
        cat_added = cat
    else:
        warnings.warn('appending new column with weights')
        cat_added=add_col(cat,'w',col_w,'f4')

    return cat_added


def get_calibration_columns(res_des):

    if 'param_split1' in config:
        param_snr  = config['param_split1'] ; param_size = config['param_split2']
        if param_size == 'sf_hlr':
            param_size = 'radius'
            warnings.warn('using radius as nbc parameter')
        if param_size == 'rgpp_rp_avg':
            param_size = 'mean_rgpp_rp'
            warnings.warn('using mean_rgpp_rp as nbc parameter')
        if param_snr == 'snr_avg':
            param_snr = 'snr'
            warnings.warn('using snr as nbc parameter')

    else:
        param_size = 'mean_rgpp_rp' ; param_snr = 'snr'

    if (param_snr == 'round_snr') & (param_snr not in res_des.dtype.names):
        warnings.warn('workaround when there is no round_snr in DES results')
        param_snr = 'snr'

    import cPickle as pickle
    filename_bias_models = os.path.join(args.output_dir,'bias_models.cpickle')
    bm=pickle.load(open(filename_bias_models))
    warnings.warn('using %s' % filename_bias_models)


    X1,X2 = np.meshgrid(bm['model_snr'],bm['model_rgp'])
    x1=X1.flatten('C'); x2=X2.flatten('C'); y1=bm['model_m'].flatten('C'); y2=bm['model_a'].flatten('C')

    q1=np.array(res_des[param_snr])
    q2=np.array(res_des[param_size])

    # griddata_method = 'linear'
    # import scipy.interpolate
    # m = scipy.interpolate.griddata((x1,x2),y1,(q1,q2),method=griddata_method,fill_value=0)
    # a = scipy.interpolate.griddata((x1,x2),y2,(q1,q2),method=griddata_method,fill_value=0)

    logger.info('getting model predictions') 
    import model_multiplicative
    import model_alpha
    logger.info(model_alpha)
    logger.info(model_multiplicative)
    filename_table_bias = os.path.join(args.output_dir,'bias_table.fits')
    m = model_multiplicative.get_model_prediction(filename_table_bias,q1.flatten(),q2.flatten())
    a = model_alpha.get_model_prediction(filename_table_bias,q1.flatten(),q2.flatten())

    c1 = res_des['mean_psf_e1_sky']*a
    c2 = res_des['mean_psf_e2_sky']*a

    import tabletools

    logger.info('adding columns')
    res_des=tabletools.ensureColumn(res_des,'nbc_m' ,arr=m, dtype='f4')
    res_des=tabletools.ensureColumn(res_des,'nbc_c1',arr=c1,dtype='f4')
    res_des=tabletools.ensureColumn(res_des,'nbc_c2',arr=c2,dtype='f4')
    res_des=tabletools.ensureColumn(res_des,'nbc_alpha',arr=a,dtype='f4')

    res_des['nbc_m']=m
    res_des['nbc_c1']=c1
    res_des['nbc_c2']=c2
    res_des['nbc_alpha']=a

    logger.info('changing all calibration values outside the range to 0')
    select = (res_des['snr']<10) | (res_des['mean_rgpp_rp']<1.1)
    res_des[select]['nbc_m'] = 0
    res_des[select]['nbc_c1']= 0
    res_des[select]['nbc_c2']= 0
    res_des[select]['nbc_alpha']= 0


    plots = False
    if plots:

        import pylab as pl
        pl.figure()
        pl.scatter(q1[:10000],q2[:10000],s=20,c=m[:10000],lw=0)
        pl.xlim([0,bm['model_snr'].max()*2 ])
        pl.ylim([0,bm['model_rgp'].max()*2 ])
        pl.colorbar()
        pl.xlabel(param_snr)
        pl.ylabel(param_size)
        pl.title('nbc_m')

        pl.figure()
        pl.pcolormesh(bm['model_snr'],bm['model_rgp'],bm['model_m'])
        pl.xscale('log')
        pl.axis('tight')

        pl.figure()
        pl.scatter(q1[:10000],q2[:10000],s=20,c=a[:10000],lw=0)
        pl.xlim([0,bm['model_snr'].max()*2 ])
        pl.ylim([0,bm['model_rgp'].max()*2 ])
        pl.colorbar()
        pl.xlabel(param_snr)
        pl.ylabel(param_size)
        pl.title('nbc_alpha')

        pl.figure()
        pl.pcolormesh(bm['model_snr'],bm['model_rgp'],bm['model_a'])
        pl.xlabel(param_snr)
        pl.ylabel(param_size)
        pl.xscale('log')
        pl.axis('tight')

        pl.show(); import pdb; pdb.set_trace()

    return res_des, m, a

def plot_face_fig():

    import fitsio
    filename_table_bias = 'bias_table.004.fits'
    bias_table = fitsio.read(filename_table_bias)
    dx=0.5

    colorscale=plotstools.get_colorscale(len(np.unique(bias_table['ipsf'])))

    list_model_m = []
    list_arr_x = []
    list_arr_y = [] 

    n_upsample = 300
    max_snr = 300
 
    bt1=bias_table
    label = 'FWHM_RATIO=%2.2f-%2.2f'%(bt1['vpsf_min'][0],bt1['vpsf_max'][0])
    snr_mid = (bt1['vsnr_min']+bt1['vsnr_max'])/2.
    rgp_mid = (bt1['vpsf_min']+bt1['vpsf_max'])/2.
    warnings.warn('setting highest SNR bin to no bias')
    pl.errorbar(snr_mid,bt1['m']-1,yerr=bt1['std_m'],label=label,fmt='.:',c='r')
    pl.xlabel('SNR')
    pl.ylabel('multiplicative shear bias')

    snr_pred =  np.linspace(snr_mid.min(),max_snr,n_upsample)

    w, w_cov =fitting.fit(snr_mid, bt1['m']-1, s=bt1['std_m'],expand=inv_snr_basis)       
    p, s = fitting.predict(snr_pred,w,w_cov,expand=inv_snr_basis)
    pl.plot(snr_pred,p,c='r')

    pl.title('noise bias calibration - all sizes',fontsize=10)
    pl.axhline(0,c='k')
    pl.xticks(list_snr_edges)
    pl.grid()
    # ylim=list(pl.ylim()); ylim[1]*=1.; pl.ylim(ylim)
    pl.ylim([-0.5,0.1])
    pl.xlim([0,50])
    # pl.legend(mode='expand',ncol=2,loc='bottom')
    pl.show()



def get_bias_model():

    # pl.style.use('supermongo')

    import fitsio
    filename_table_bias = os.path.join(args.output_dir,'bias_table.fits')

    eps_m = 1e-5
    eps_a = 1e-2

    bias_table = tktools.load(filename_table_bias)
    dx=0.005

    interp_method = 'linear'

    n_psf_bins = len(np.unique(bias_table['ipsf']))
    colorscale=plotstools.get_colorscale(n_psf_bins,cmap_name='rainbow')
    # colorscale=plotstools.get_colorscale(n_psf_bins,cmap_name='gist_rainbow')

    list_model_m = []
    list_arr_x = []
    list_arr_y = [] 

    if 'param_split1' in config:
        param_snr  = config['param_split1'] ; param_size = config['param_split2']
    else:
        param_size = 'mean_rgpp_rp' ; param_snr = 'snr'


    n_upsample = 300
    max_snr = 300

    def log_separate(x,dx,n_sep):

        sep = []
        for i in range(len(x)):
            s=np.logspace(np.log10(x[i])-dx,np.log10(x[i])+dx,n_sep);
            sep.append(s)

        sep = np.array(sep)
        return sep

    def get_optimal_w(snr_mid, m, s, basis, n_gals=1):

        list_chi2 = []
        list_res_sig = []
        list_w = []
        list_w_cov = []
        eps_grid=np.linspace(-15,15,1000)
        for eps in eps_grid:
            w, w_cov =fitting.fit(snr_mid, m, s,expand=basis,eps=10**eps)
            p_mid, sigma = fitting.predict(snr_mid,w,w_cov,expand=basis)

            # chi2 =  np.sum( n_gals* (((m-p_mid)**2)/(s**2)) ) / np.sum(n_gals)
            # chi2 =  np.mean( (((m-p_mid)**2)/(s**2)) ) 
            # chi2 =  np.sum( (((m-p_mid)**2)/(s**2)) ) / float( len(snr_mid) - basis(snr_mid[:,None]).shape[1] - 1  )
            # chi2 =  np.sum( n_gals*(m-p_mid) ) / np.sum(n_gals)
            # chi2 =  np.mean( (((m-p_mid)**2)/(s**2)) )  - np.sum( n_gals*(m-p_mid) )/ float(np.sum(n_gals))
            chi2 =  np.mean( (((m-p_mid)**2)/(s**2)) ) 
            res_sig = np.sum( n_gals*(m-p_mid) )/ float(np.sum(n_gals))
            logger.debug('optimal w: eps=%2.2e chi2=%2.2e mean_res/sig=%2.4e' % (10**(eps),chi2,res_sig))
            list_chi2.append(chi2)
            list_w.append(w)
            list_w_cov.append(w_cov)
            list_res_sig.append(res_sig)

        arr_chi2 = np.array(list_chi2)

        select = np.argmin(np.abs(arr_chi2-1))

        w = list_w[select]
        w_cov = list_w_cov[select]
        eps = eps_grid[select]
        chi2 = list_chi2[select]
        res_sig = list_res_sig[select]

        logger.info('final optimal w: eps=%2.2e chi2=%2.2e res_sig=%2.5f w=%s' % (10**(eps),chi2,res_sig,np.array_str(w,precision=2)))

        return w, w_cov

    list_chi2 = []
    for ipsf in np.unique(bias_table['ipsf']):

        select = bias_table['ipsf'] == ipsf
        bt1=bias_table[select]
        label = '%s=%2.2f-%2.2f'%(param_size,bt1['vpsf_min'][0],bt1['vpsf_max'][0])
        # label = '%s=%2.2f' % (param_size,bt1['vpsf_mid'][0])
        # snr_mid = (bt1['vsnr_min']+bt1['vsnr_max'])/2.
        # rgp_mid = (bt1['vpsf_min']+bt1['vpsf_max'])/2.
        snr_mid = bt1['vsnr_mid']
        rgp_mid = bt1['vpsf_mid']
        logger.debug('rgp_mid = %2.4f' % rgp_mid[0])

        if (rgp_mid[0]<1.21) | (rgp_mid[0]>2):
            logger.warning('skipping bin 1<1.15 -- not including in pltos')
            continue

        import model_multiplicative
        # import model_multiplicative_alt as model_multiplicative
        pl.figure(1)
        pl.xlabel(param_snr)
        pl.ylabel('m')
        pl.errorbar(snr_mid,bt1['m']-1,yerr=bt1['std_m'],label=label,fmt='o',c=colorscale[ipsf-1])
        try:
            snr_pred =  np.linspace(snr_mid.min(),max_snr,n_upsample)
            # w, w_cov = get_optimal_w(snr_mid=snr_mid, m=bt1['m']-1, s=bt1['std_m'], basis=inv_snr_basis, n_gals=bt1['n_gals'])
            # w, w_cov =fitting.fit(snr_mid, bt1['m']-1, s=bt1['std_m'],expand=inv_snr_basis,eps=0.01*np.max(bt1['std_m'])**2)       
            # w, w_cov =fitting.fit(snr_mid, bt1['m']-1, s=bt1['std_m'],expand=inv_snr_basis,eps=eps_m)       
            p = model_multiplicative.get_model_prediction(filename_table_bias,snr_pred,np.ones_like(snr_pred)*np.mean(rgp_mid),plots=False)
            # p, s = fitting.predict(snr_pred,w,w_cov,expand=inv_snr_basis)
            pl.plot(snr_pred,p,c=colorscale[ipsf-1])
            # pl.plot(snr_pred,p,c=colorscale[ipsf-1],lw=np.sum(bt1['n_gals'])/1e6)
            list_model_m.append(p)
            list_arr_x.append(snr_pred)
            list_arr_y.append(np.ones_like(snr_pred)*rgp_mid[0])
        except Exception,errmsg:
            logger.error(errmsg)


        pl.figure(2)
        try:
            snr_mid_sep = log_separate(x=snr_mid,dx=dx,n_sep=n_psf_bins)
            p_mid = model_multiplicative.get_model_prediction(filename_table_bias,snr_mid,np.ones_like(snr_mid)*np.mean(rgp_mid),plots=False)
            # p_mid, s = fitting.predict(snr_mid,w,w_cov,expand=inv_snr_basis)
            pl.errorbar(snr_mid_sep[:,ipsf-1],(bt1['m']-1)-p_mid,yerr=bt1['std_m'],label=label,fmt='o',c=colorscale[ipsf-1])

            # chi2 = np.sum(   bt1['n_gals']*((((bt1['m']-1)-p_mid)**2)/(bt1['std_m']**2))   ) / np.sum( bt1['n_gals'] ) 
            chi2 = np.mean(  (((bt1['m']-1)-p_mid)**2)/(bt1['std_m']**2))  
            list_chi2.append(chi2)
            logger.info('%f %f' % ( ipsf,chi2))

        except Exception,errmsg:
            logger.error(errmsg)

    logger.info('chi2 m %f' , np.mean(np.array(list_chi2)))

    pl.figure(1,figsize=(8,5))
    # pl.title(selection_string_des,fontsize=10)
    pl.xscale('log')
    # pl.title('noise bias: multiplicative - im3shape - %s' % filename_table_bias)
    pl.axhline(0,c='k')
    pl.xticks(list_snr_edges)
    pl.grid()
    pl.legend(framealpha=0.0,frameon=False,loc='lower right',fontsize=14)
    pl.ylim([-0.35,0.15])
    pl.gca().xaxis.set_major_formatter(pl.matplotlib.ticker.ScalarFormatter())
    pl.xlim([12.5,80])
    # plot_add_requirements(level=0.03,target=0)
    # plot_add_requirements(level=0.01,target=0)
    # pl.tick_params(axis='both', which='major', labelsize=10)

    pl.figure(2)
    # pl.title(selection_string_des,fontsize=10)
    pl.xscale('log')
    # pl.title('noise bias: multiplicative - im3shape - %s' % filename_table_bias)
    pl.axhline(0,c='k')
    pl.xticks(list_snr_edges)
    pl.grid()
    pl.legend(framealpha=0.0,frameon=False,loc='lower right',fontsize=14)
    pl.ylim([-0.2,0.05])
    pl.gca().xaxis.set_major_formatter(pl.matplotlib.ticker.ScalarFormatter())
    pl.xlim([9,200])
    pl.ylabel('data - model')
    # plot_add_requirements(level=0.03,target=0)
    # plot_add_requirements(level=0.01,target=0)

    # -----------------------------------

    try:

        model_arr_m = np.array(list_model_m)
        model_arr_x = np.array(list_arr_x)
        model_arr_y = np.array(list_arr_y)
        vec_x = model_arr_x[0,:]
        vec_y = model_arr_y[:,0]
        vec_x_hires = np.linspace(vec_x.min(),vec_x.max(),n_upsample)
        vec_y_hires = np.linspace(vec_y.min(),vec_y.max(),n_upsample)
        import scipy.interpolate    
        func_interp=scipy.interpolate.interp2d(vec_x,vec_y,model_arr_m,kind=interp_method)
        model_arr_m_hires=func_interp(vec_x_hires,vec_y_hires)


        pl.figure()
        pl.pcolormesh(model_arr_x,model_arr_y,model_arr_m)
        pl.xticks(vec_x)
        pl.yticks(vec_y)
        pl.axis('tight')
        pl.xlabel(param_snr)
        pl.ylabel('Rgp/Rp')
        pl.title('multiplicative shear bias')
        
        # lin scale    
        pl.figure()
        pl.pcolormesh(vec_x_hires,vec_y_hires,model_arr_m_hires)
        pl.colorbar()
        # cs=pl.contour(vec_x_hires,vec_y_hires,model_arr_m_hires,[0,-0.05,0.10,-0.13])
        # pl.clabel(cs, fontsize=9, inline=1)
        pl.xscale('log')
        pl.xticks(list_snr_edges)
        pl.gca().xaxis.set_major_formatter(pl.matplotlib.ticker.ScalarFormatter())
        pl.yticks(list_psf_edges)
        pl.xlabel(param_snr)
        pl.ylabel(param_size)
        pl.xlim([12,100])
        pl.ylim([1.15,2])
        # pl.axis('tight')
        pl.title('multiplicative shear bias')

        # # now use the 2d model

        vec_x_2d = np.logspace(1,2,100)
        vec_y_2d = np.logspace(0.06,0.6,100)
        X,Y=np.meshgrid(vec_x_2d,vec_y_2d)
        import model_multiplicative
        model2d_vec_a_hires = model_multiplicative.get_model_prediction(filename_table_bias,X.flatten(),Y.flatten(),plots=True)
        model2d_arr_a_hires = np.reshape(model2d_vec_a_hires,[100,100])

        pl.figure()
        pl.pcolormesh(vec_x_2d,vec_y_2d,model2d_arr_a_hires)
        pl.axis('tight')
        pl.xlabel(param_snr)
        pl.ylabel('Rgp/Rp')
        pl.axis('tight')
        pl.title('multiplicaive - 2d model')
        pl.xscale('log')
        pl.ylim([1.1,3])
        pl.xlim([0,100])
        # pl.clim([-0.1,0.75])
        pl.colorbar()

    except Exception,errmsg:
        logger.error(errmsg)

    # now leakage

    list_model_a = []
    list_arr_x = []
    list_arr_y = []
    pl.figure()
    list_chi2 = []

   
    for ipsf in np.unique(bias_table['ipsf']):

        pl.figure(20)
        select = bias_table['ipsf'] == ipsf
        bt1=bias_table[select]
        label = '%s=%2.2f-%2.2f'%(param_size,bt1['vpsf_min'][0],bt1['vpsf_max'][0])
        snr_mid = bt1['vsnr_mid']
        rgp_mid = bt1['vpsf_mid']

        if (rgp_mid[0]<1.21) | (rgp_mid[0]>2):
            logger.warning('skipping bin 1<1.15 -- not including in pltos')
            continue

        # snr_mid = (bt1['vsnr_min']+bt1['vsnr_max'])/2.
        # rgp_mid = (bt1['vpsf_min']+bt1['vpsf_max'])/2.

        # warnings.warn('setting highest SNR bin to no bias')
        # bt1['m'][-1] = 0
        # bt1['std_m'][-1] = 0.001
        # snr_mid[-1] = 2000

        warnings.warn('setting highest SNR bin to no bias')
        pl.errorbar(snr_mid+dx*ipsf,bt1['pmm'],yerr=bt1['std_pmm'],label=label,fmt='o',c=colorscale[ipsf-1])


        # w, w_cov = get_optimal_w(snr_mid=snr_mid, m=bt1['pmm'], s=bt1['std_pmm'],basis=inv_snr_basis2,n_gals=bt1['n_gals'])
        # w, w_cov =fitting.fit(snr_mid, bt1['pmm'], s=bt1['std_pmm'],expand=inv_snr_basis2,eps=np.max(bt1['std_pmm'])**2)
        # w, w_cov =fitting.fit(snr_mid, bt1['pmm'], s=bt1['std_pmm'],expand=inv_snr_basis2,eps=eps_a)

        
        snr_pred =  np.linspace(snr_mid.min(),max_snr,n_upsample)

        import model_alpha
        # import model_alpha_alt as model_alpha
        # p, s = fitting.predict(snr_pred,w,w_cov,expand=inv_snr_basis2)
        p = model_alpha.get_model_prediction(filename_table_bias,snr_pred,np.ones_like(snr_pred)*np.mean(rgp_mid),plots=False)
        pl.plot(snr_pred,p,c=colorscale[ipsf-1])
        list_model_a.append(p)
        list_arr_x.append(snr_pred)
        list_arr_y.append(np.ones_like(snr_pred)*rgp_mid[0])

        pl.figure(21)
        try:
            snr_mid_sep = log_separate(x=snr_mid,dx=dx,n_sep=n_psf_bins)
            # p_mid, s = fitting.predict(snr_mid,w,w_cov,expand=inv_snr_basis2)
            p = model_alpha.get_model_prediction(filename_table_bias,snr_mid,np.ones_like(snr_mid)*np.mean(rgp_mid),plots=False)
            pl.errorbar(snr_mid_sep[:,ipsf-1],(bt1['pmm'])-p_mid,yerr=bt1['std_pmm'],label=label,fmt='o',c=colorscale[ipsf-1])

            # chi2 = np.mean((((bt1['pmm'])-p_mid)**2)/(bt1['std_pmm']**2))
            chi2 = np.sum( bt1['n_gals']*(((bt1['pmm'])-p_mid)**2)/(bt1['std_pmm']**2) ) / np.sum(bt1['n_gals'])
            list_chi2.append(chi2)
            logger.info('%d rgp_mid=%2.2f %f' % ( ipsf,rgp_mid[0],chi2))
        except Exception,errmsg:
            logger.error(errmsg)

    logger.info('chi2 alpha %f' , np.mean(np.array(list_chi2)))

    pl.figure(20,figsize=(8,5))
    # pl.title('noise bias: PSF leakage - im3shape %s' % filename_table_bias)
    # pl.legend(framealpha=0.0,frameon=False,loc='lower right',fontsize=14,mode='expand',ncol=2)
    pl.xscale('log')
    ylim=list(pl.ylim()); ylim[1]*=1.; pl.ylim(ylim)
    pl.axhline(0,c='k')
    pl.legend(framealpha=0.0,frameon=False,loc='upper right')
    pl.xlabel(param_snr)
    pl.ylabel(r'leakage $\alpha$')
    pl.xticks(list_snr_edges)
    pl.xlim([12.5,80])
    pl.ylim([-0.25,0.7])
    pl.grid()
    # plot_add_requirements(level=0.05,target=0)
    pl.gca().xaxis.set_major_formatter(pl.matplotlib.ticker.ScalarFormatter())

    
    # pl.imshow(model_arr_a,aspect='auto',interpolation='nearest'); pl.show()

    pl.figure(21)
    # pl.title(selection_string_des,fontsize=10)
    pl.xscale('log')
    # pl.title('noise bias: multiplicative - im3shape - %s' % filename_table_bias)
    pl.axhline(0,c='k')
    pl.xticks(list_snr_edges)
    pl.grid()
    pl.legend(framealpha=0.0,frameon=False,loc='lower right',fontsize=14)
    # pl.ylim([-0.2,0.05])
    pl.gca().xaxis.set_major_formatter(pl.matplotlib.ticker.ScalarFormatter())
    pl.xlim([9,200])
    pl.ylabel('data - model')
    # plot_add_requirements(level=0.03,target=0)
    # plot_add_requirements(level=0.01,target=0)


    try:

        model_arr_a = np.array(list_model_a)
        model_arr_x = np.array(list_arr_x)
        model_arr_y = np.array(list_arr_y)
        
        vec_x = model_arr_x[0,:]
        vec_y = model_arr_y[:,0]
        vec_x_hires = np.linspace(vec_x.min(),vec_x.max(),n_upsample)
        vec_y_hires = np.linspace(vec_y.min(),vec_y.max(),n_upsample)
        import scipy.interpolate    
        func_interp=scipy.interpolate.interp2d(vec_x,vec_y,model_arr_a,kind=interp_method)
        model_arr_a_hires=func_interp(vec_x_hires,vec_y_hires)

        pl.figure()
        pl.pcolormesh(model_arr_x,model_arr_y,model_arr_a)
        pl.xticks(vec_x)
        pl.yticks(vec_y)
        pl.xlabel(param_snr)
        pl.ylabel('Rgp/Rp')
        pl.axis('tight')
        pl.title('PSF leakage')
        pl.xscale('log')
        pl.colorbar()


        # lin scale
        pl.figure()
        pl.pcolormesh(vec_x_hires,vec_y_hires,model_arr_a_hires)
        pl.colorbar()
        cs=pl.contour(vec_x_hires,vec_y_hires,model_arr_a_hires,[0,0.05,0.10,0.15,0.20])
        pl.clabel(cs, fontsize=9, inline=1)
        pl.axis('tight')
        pl.xlabel(param_snr)
        pl.ylabel('Rgp/Rp')
        pl.axis('tight')
        pl.title('PSF leakage')
        pl.xscale('log')

        # # now use the 2d model


        vec_x_2d = np.logspace(0,2,100)
        vec_y_2d = np.logspace(0,0.3,100)
        X,Y=np.meshgrid(vec_x_2d,vec_y_2d)
        import model_alpha
        model2d_vec_a_hires = model_alpha.get_alpha_model_prediction(filename_table_bias,X.flatten(),Y.flatten(),plots=False)
        model2d_arr_a_hires = np.reshape(model2d_vec_a_hires,[100,100])

        pl.figure()
        pl.pcolormesh(vec_x_2d,vec_y_2d,model2d_arr_a_hires)
        pl.axis('tight')
        pl.xlabel(param_snr)
        pl.ylabel('Rgp/Rp')
        pl.axis('tight')
        pl.title('PSF leakage - 2d model')
        pl.xscale('log')
        pl.ylim([1.1,3])
        pl.xlim([0,100])
        # pl.clim([-0.1,0.75])
        pl.colorbar()




    except Exception,errmsg:
        logger.error(errmsg)

    # n_unique

    pl.figure()
    for ipsf in np.unique(bias_table['ipsf']):

        select = bias_table['ipsf'] == ipsf
        bt1=bias_table[select]
        # label = r'$R_{gp}/R_p=%2.2f-%2.2f$'%(bt1['vpsf_min'][0],bt1['vpsf_max'][0])
        snr_mid = (bt1['vsnr_min']+bt1['vsnr_max'])/2.
        rgp_mid = (bt1['vpsf_min']+bt1['vpsf_max'])/2.
        warnings.warn('setting highest SNR bin to no bias')
        # bt1['m'][-1] = 1
        # bt1['std_m'][-1] = 0.001
        # snr_mid[-1] = 2000
        # bt1['m'][-2] = 1
        # bt1['std_m'][-2] = 0.0001
        # pl.errorbar(snr_mid,bt1['m']-1,yerr=bt1['std_m'],label=label,fmt='.',c=colorscale[ipsf-1])
        pl.plot(snr_mid,bt1['n_unique'],c=colorscale[ipsf-1])
        pl.xlabel(param_snr)
        pl.ylabel('n_unique')


    # pl.title(selection_string_des,fontsize=10)
    pl.xscale('log')
    # pl.title('noise bias: multiplicative - im3shape - %s' % filename_table_bias)
    pl.axhline(0,c='k')
    pl.xticks(list_snr_edges)
    # pl.grid()
    pl.legend(framealpha=0.0,frameon=False,loc='lower right',fontsize=14)
    # pl.xlim([9,200])
    # pl.ylim([-0.5,0.15])
    pl.gca().xaxis.set_major_formatter(pl.matplotlib.ticker.ScalarFormatter())
    # pl.tick_params(axis='both', which='major', labelsize=10)

    # bulge_frac

    pl.figure()
    for ipsf in np.unique(bias_table['ipsf']):

        select = bias_table['ipsf'] == ipsf
        bt1=bias_table[select]
        # label = r'$R_{gp}/R_p=%2.2f-%2.2f$'%(bt1['vpsf_min'][0],bt1['vpsf_max'][0])
        snr_mid = (bt1['vsnr_min']+bt1['vsnr_max'])/2.
        rgp_mid = (bt1['vpsf_min']+bt1['vpsf_max'])/2.
        warnings.warn('setting highest SNR bin to no bias')
        # bt1['m'][-1] = 1
        # bt1['std_m'][-1] = 0.001
        # snr_mid[-1] = 2000
        # bt1['m'][-2] = 1
        # bt1['std_m'][-2] = 0.0001
        # pl.errorbar(snr_mid,bt1['m']-1,yerr=bt1['std_m'],label=label,fmt='.',c=colorscale[ipsf-1])
        pl.plot(snr_mid,bt1['bulge_fraction'],c=colorscale[ipsf-1])
        pl.xlabel(param_snr)
        pl.ylabel('bulge_fraction')


    # pl.title(selection_string_des,fontsize=10)
    pl.xscale('log')
    # pl.title('noise bias: multiplicative - im3shape - %s' % filename_table_bias)
    pl.axhline(0,c='k')
    pl.xticks(list_snr_edges)
    # pl.grid()
    pl.legend(framealpha=0.0,frameon=False,loc='lower right',fontsize=14)
    # pl.xlim([9,200])
    # pl.ylim([-0.5,0.15])
    pl.gca().xaxis.set_major_formatter(pl.matplotlib.ticker.ScalarFormatter())
    # pl.tick_params(axis='both', which='major', labelsize=10)

    # x1 x2

    pl.figure()
    true_center = -0.27
    for ipsf in np.unique(bias_table['ipsf']):

        select = bias_table['ipsf'] == ipsf
        bt1=bias_table[select]
        # label = r'$R_{gp}/R_p=%2.2f-%2.2f$'%(bt1['vpsf_min'][0],bt1['vpsf_max'][0])
        snr_mid = (bt1['vsnr_min']+bt1['vsnr_max'])/2.
        rgp_mid = (bt1['vpsf_min']+bt1['vpsf_max'])/2.
        warnings.warn('setting highest SNR bin to no bias')
        # bt1['m'][-1] = 1
        # bt1['std_m'][-1] = 0.001
        # snr_mid[-1] = 2000
        # bt1['m'][-2] = 1
        # bt1['std_m'][-2] = 0.0001
        pl.errorbar(snr_mid,bt1['x1'],yerr=bt1['std_x1'],label=label,fmt='x',c=colorscale[ipsf-1])
        pl.errorbar(snr_mid,bt1['x2'],yerr=bt1['std_x2'],label=label,fmt='+',c=colorscale[ipsf-1])
        # pl.errorbar(snr_mid,bt1['desx1'],yerr=bt1['std_desx1'],label=label,fmt='o',c=colorscale[ipsf-1])
        # pl.errorbar(snr_mid,bt1['desx2'],yerr=bt1['std_desx2'],label=label,fmt='d',c=colorscale[ipsf-1])
        pl.xlabel(param_snr)
        pl.ylabel('x')


    # pl.title(selection_string_des,fontsize=10)
    pl.xscale('log')
    # pl.title('noise bias: multiplicative - im3shape - %s' % filename_table_bias)
    pl.axhline(0,c='k')
    pl.xticks(list_snr_edges)
    # pl.grid()
    pl.legend(framealpha=0.0,frameon=False,loc='lower right',fontsize=14)
    # pl.xlim([9,200])
    # pl.ylim([-0.5,0.15])
    pl.gca().xaxis.set_major_formatter(pl.matplotlib.ticker.ScalarFormatter())
    # pl.tick_params(axis='both', which='major', labelsize=10)

    # ----------------------------------------------------------------------------------------------

    pl.figure()
    for ipsf in np.unique(bias_table['ipsf']):

        select = bias_table['ipsf'] == ipsf
        bt1=bias_table[select]
        label = '%s=%2.2f-%2.2f'%(param_size,bt1['vpsf_min'][0],bt1['vpsf_max'][0])
        snr_mid = (bt1['vsnr_min']+bt1['vsnr_max'])/2.
        rgp_mid = (bt1['vpsf_min']+bt1['vpsf_max'])/2.
        pl.errorbar(snr_mid,bt1['m1']-1,yerr=bt1['std_m1'],label=label,fmt='+:',c=colorscale[ipsf-1])
        pl.errorbar(snr_mid,bt1['m2']-1,yerr=bt1['std_m2'],label=label,fmt='x--',c=colorscale[ipsf-1])
        pl.xlabel(param_snr)
        pl.ylabel('m')

        # snr_pred =  np.linspace(snr_mid.min(),max_snr,n_upsample)
 
        # w, w_cov =fitting.fit(snr_mid, bt1['m']-1, s=bt1['std_m'],expand=inv_snr_basis)       
        # p, s = fitting.predict(snr_pred,w,w_cov,expand=inv_snr_basis)
        # pl.plot(snr_pred,p,c=colorscale[ipsf-1])
        # list_model_m.append(p)
        # list_arr_x.append(snr_pred)
        # list_arr_y.append(np.ones_like(snr_pred)*rgp_mid[0])


    # pl.title(selection_string_des,fontsize=10)
    pl.xscale('log')
    # pl.title('noise bias: multiplicative - im3shape - %s' % filename_table_bias)
    pl.axhline(0,c='k')
    pl.xticks(list_snr_edges)
    # pl.grid()
    pl.legend(framealpha=0.0,frameon=False,loc='lower right',fontsize=14)
    pl.xlim([9,200])
    pl.ylim([-0.5,0.15])
    pl.gca().xaxis.set_major_formatter(pl.matplotlib.ticker.ScalarFormatter())
    # pl.tick_params(axis='both', which='major', labelsize=10)


    # ----------------------------------------------------------------------------------------------

    pl.figure()

    for ipsf in np.unique(bias_table['ipsf']):

        select = bias_table['ipsf'] == ipsf
        bt1=bias_table[select]
        label = '%s=%2.2f-%2.2f'%(param_size,bt1['vpsf_min'][0],bt1['vpsf_max'][0])
        snr_mid = (bt1['vsnr_min']+bt1['vsnr_max'])/2.
        rgp_mid = (bt1['vpsf_min']+bt1['vpsf_max'])/2.
        snr_mid_sep = log_separate(x=snr_mid,dx=dx,n_sep=n_psf_bins)
        pl.errorbar(snr_mid_sep[:,ipsf-1],bt1['mms']-1,yerr=bt1['std_mms'],label=label,fmt='.',c=colorscale[ipsf-1])
        # pl.errorbar(snr_mid,bt1['mm1s']-1,yerr=bt1['std_mm1s'],fmt='+--',c=colorscale[ipsf-1])
        # pl.errorbar(snr_mid,bt1['mm2s']-1,yerr=bt1['std_mm2s'],fmt='x:',c=colorscale[ipsf-1])
        pl.xlabel(param_snr)
        pl.ylabel('m (sheared shape)')

        # snr_pred =  np.linspace(snr_mid.min(),max_snr,n_upsample)
 
        # w, w_cov =fitting.fit(snr_mid, bt1['m']-1, s=bt1['std_m'],expand=inv_snr_basis)       
        # p, s = fitting.predict(snr_pred,w,w_cov,expand=inv_snr_basis)
        # pl.plot(snr_pred,p,c=colorscale[ipsf-1])
        # list_model_m.append(p)
        # list_arr_x.append(snr_pred)
        # list_arr_y.append(np.ones_like(snr_pred)*rgp_mid[0])


    # pl.title(selection_string_des,fontsize=10)
    pl.xscale('log')
    # pl.title('noise bias: multiplicative - im3shape - %s' % filename_table_bias)
    pl.axhline(0,c='k')
    pl.xticks(list_snr_edges)
    # pl.grid()
    pl.legend(framealpha=0.0,frameon=False,loc='lower right',fontsize=14)
    pl.ylim([-0.05,0.05])
    pl.gca().xaxis.set_major_formatter(pl.matplotlib.ticker.ScalarFormatter())
    pl.xlim([9,200])
    plot_add_requirements(level=0.03,target=0)
    plot_add_requirements(level=0.01,target=0)

    # pl.tick_params(axis='both', which='major', labelsize=10)

    # -----------------------------------


    pl.show()

    filename_bias_models = os.path.join(args.output_dir,'bias_models.cpickle')
    import cPickle as pickle
    pickle_dict = {'model_m':model_arr_m_hires, 'model_a':model_arr_a_hires, 'model_snr':vec_x_hires, 'model_rgp':vec_y_hires}
    pickle.dump(pickle_dict,open(filename_bias_models,'w'),protocol=2)
    logger.info('wrote %s' % (filename_bias_models))




    import pdb; pdb.set_trace()

def load_selection():

        res_sim = tt.load( '%s/res_sim.fits' % args.output_dir)
        res_des = tt.load( '%s/res_des.fits' % args.output_dir)
        res_tru = tt.load( '%s/res_tru.fits' % args.output_dir)

        return res_sim, res_tru, res_des

def get_bias_vs_redshift():

    # cols_res=['coadd_objects_id','e1','e2', 'w', 'snr','mean_rgpp_rp','mean_psf_e1_sky','mean_psf_e2_sky']
    # cols_tru=['snr','psf_e1','psf_e2','id_shear','cosmos_mag_auto','g1_true','g2_true','zphot','id_cosmos']
    # cols_des=['coadd_objects_id','e1','e2', 'w', 'snr','mean_rgpp_rp','mean_psf_e1_sky','mean_psf_e2_sky']
    # # res_sim,res_tru = nbc_v7_select.get_selection_sim(selection_string_sim,cols_res=cols_res,cols_tru=cols_tru)
    # # res_des         = nbc_v7_select.get_selection_des(selection_string_des,cols=cols_des,n_files=10)
    # cols_res+= ['nbc_m1','nbc_m2','nbc_c1','nbc_c2']
    # cols_des+= ['nbc_m1','nbc_m2','nbc_c1','nbc_c2']
    # cal_sim,cal_tru = nbc_v7_select.get_selection_sim(selection_string_final_sim,cols_res=cols_res,cols_tru=cols_tru,get_calibrated=True)
    # cal_des         = nbc_v7_select.get_selection_des(selection_string_final_des,cols=cols_des,n_files=args.num,get_calibrated=True)
    # res_sim,res_tru = cal_sim,cal_tru
    # res_des         = cal_des        

    res_sim, res_tru, res_des = load_selection()

    min_sigma_ell = 0.24
    warnings.warn('setting max weight to %2.2f' % (1./min_sigma_ell**2) )
    res_sim['w'][ (1./np.sqrt(res_sim['w'])) < min_sigma_ell ] = 1./min_sigma_ell**2
    res_des['w'][ (1./np.sqrt(res_des['w'])) < min_sigma_ell ] = 1./min_sigma_ell**2

    logger.info('calculating bias for the entire sample - using cuts from final selection')
    logger.info(selection_string_final_sim)
    logger.info(selection_string_final_des)
    cat_res = res_sim
    cat_tru = res_tru
    exec selection_string_final_sim
    select_sim = select.copy()
    cat_res = res_des
    exec selection_string_final_des
    select_des = select.copy()

    res_sim_final = res_sim[select_sim]
    res_tru_final = res_tru[select_sim]
    res_des_final = res_des[select_des]


    logger.info("pre-selection %d final-selection %d" % (len(res_sim),len(res_sim_final)))
    correct_selection_bias = False
    logger.info('correct_selection_bias=%d',correct_selection_bias)


    logger.info('final selection SIM %d/%d %f',len(res_sim[select_sim]),len(res_sim), len(res_sim[select_sim])/float(len(res_sim)))
    logger.info('final selection DES %d/%d %f',len(res_des[select_des]),len(res_des), len(res_des[select_des])/float(len(res_des)))
     
    # entire sample
    filename_str = 'z.all.nonbc.%s' % args.method
    mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2,mms,std_mms,ccs,std_ccs,mm1s,std_mm1s,cc1s,std_cc1s,mm2s,std_mm2s,cc2s,std_cc2s=nbc_v7_stats.get_mc(res_sim_final,res_tru_final,res_des_final,use_calibration=False,use_weights=args.use_weights,filename_str=filename_str,correct_selection_bias=correct_selection_bias)

    filename_str = 'z.all.withnbc.%s' % args.method
    mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2,mms,std_mms,ccs,std_ccs,mm1s,std_mm1s,cc1s,std_cc1s,mm2s,std_mm2s,cc2s,std_cc2s=nbc_v7_stats.get_mc(res_sim_final,res_tru_final,res_des_final,use_calibration=True,use_weights=args.use_weights,filename_str=filename_str,correct_selection_bias=correct_selection_bias)


    z_bins = [ 0.3, 0.644, 0.901, 1.3 ]

    list_bias = []
    list_bias_calibr = []


    for ibin in range(1,len(z_bins)):

        logger.info('--------------------------------------------------')

        vbin = [z_bins[ibin-1],z_bins[ibin]]

        select = (res_tru_final['zphot'] > vbin[0]) & (res_tru_final['zphot'] < vbin[1])
        res_sim_select = res_sim_final[select]
        res_tru_select = res_tru_final[select]

        # dummy - ignore this, it's only here so that the code runs through
        res_des_select = res_des_final

        logger.info(selection_string_final_sim)
        filename_str = 'z%2.2f.%s' % ((vbin[0]+vbin[1])/2.,args.method)
        mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2,mms,std_mms,ccs,std_ccs,mm1s,std_mm1s,cc1s,std_cc1s,mm2s,std_mm2s,cc2s,std_cc2s=nbc_v7_stats.get_mc(res_sim_select,res_tru_select,None,use_calibration=False,use_weights=args.use_weights,filename_str=filename_str,correct_selection_bias=correct_selection_bias)

        std_e = np.std(res_sim_select['e1'],ddof=1)
        list_bias.append( [ibin,vbin[0],vbin[1],std_e,mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2,mms,std_mms,ccs,std_ccs,mm1s,std_mm1s,cc1s,std_cc1s,mm2s,std_mm2s,cc2s,std_cc2s] )

        filename_str = 'z%2.2f.%s.corr' % ((vbin[0]+vbin[1])/2.,args.method)
        mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2,mms,std_mms,ccs,std_ccs,mm1s,std_mm1s,cc1s,std_cc1s,mm2s,std_mm2s,cc2s,std_cc2s=nbc_v7_stats.get_mc(res_sim_select,res_tru_select,None,use_calibration=True,use_weights=args.use_weights,filename_str=filename_str,correct_selection_bias=correct_selection_bias)

        std_e = np.std(res_sim_select['e1'],ddof=1)
        list_bias_calibr.append( [ibin,vbin[0],vbin[1],std_e,mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2,mms,std_mms,ccs,std_ccs,mm1s,std_mm1s,cc1s,std_cc1s,mm2s,std_mm2s,cc2s,std_cc2s] )

    arr_bias = tktools.arr2rec(np.array(list_bias),dtype={'names': ["iz","z_min","z_max","std_e","m","std_m","c","std_c","m1","std_m1","m2","std_m2","c1","std_c1","c2","std_c2","pmm","std_pmm","pcc","std_pcc","pm1","std_pm1","pm2","std_pm2","mean_e1","mean_e2","stdm_e1","stdm_e2","mms","std_mms","ccs","std_ccs","mm1s","std_mm1s","cc1s","std_cc1s","mm2s","std_mm2s","cc2s","std_cc2s"], 'formats': ['i4']*1 + ['f8']*39 })
    arr_bias_calibr = tktools.arr2rec(np.array(list_bias_calibr),dtype={'names': ["iz","z_min","z_max","std_e","m","std_m","c","std_c","m1","std_m1","m2","std_m2","c1","std_c1","c2","std_c2","pmm","std_pmm","pcc","std_pcc","pm1","std_pm1","pm2","std_pm2","mean_e1","mean_e2","stdm_e1","stdm_e2","mms","std_mms","ccs","std_ccs","mm1s","std_mm1s","cc1s","std_cc1s","mm2s","std_mm2s","cc2s","std_cc2s"], 'formats': ['i4']*1 + ['f8']*39 })

    import cPickle as pickle
    
    if args.use_weights == True:
        filename_pickle = os.path.join(args.output_dir,'figdata.bias_vs_redshift.with-weights.pp2')
    else:
        filename_pickle = os.path.join(args.output_dir,'figdata.bias_vs_redshift.pp2')

    pickle.dump((arr_bias,arr_bias_calibr,z_bins),open(filename_pickle,'w'),protocol=2) 
    logger.info('wrote %s' % (filename_pickle))

    # now do jacknifes
    # n_jack = config['n_jack']
    # filename_cosmos = os.path.join(config['input']['real_catalog']['image_dir'],config['input']['real_catalog']['file_name'])
    # cat_cosmos = tktools.load(filename_cosmos)
    # col_jack = cat_cosmos[cal_tru['id_cosmos'].astype(np.int32)]['jack_id']
    # res_tru = tktools.add_col(rec=res_tru,arr=col_jack,name='jack_id',dtype='i4')

    # list_bias = []
    # list_bias_calibr = []

    # for ibin in range(1,len(z_bins)):

    #     list_jack = []
    #     list_jack_calibr = []
            
    #     for ijack in range(n_jack):

    #         vbin = [z_bins[ibin-1],z_bins[ibin]]

    #         select = (res_tru['zphot'] > vbin[0]) & (res_tru['zphot'] < vbin[1]) & (res_tru['jack_id'] != ijack) 
    #         res_sim_select = res_sim[select]
    #         res_tru_select = res_tru[select]

    #         select = (cal_tru['zphot'] > vbin[0]) & (cal_tru['zphot'] < vbin[1]) & (res_tru['jack_id'] != ijack) 
    #         cal_sim_select = cal_sim[select]
    #         cal_tru_select = cal_tru[select]

    #         # dummy - ignore this
    #         cal_des_select = cal_des
    #         res_des_select = res_des

    #         logger.info(selection_string_sim)
    #         filename_str = 'z%2.2f.%s.jack%02d' % ((vbin[0]+vbin[1])/2.,args.method,ijack)
    #         mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2=nbc_v7_stats.get_mc(res_sim_select,res_tru_select,None,use_calibration=False,use_weights=args.use_weights,filename_str=filename_str)

    #         std_e = np.std(res_sim_select['e1'],ddof=1)
    #         list_jack.append( [ibin,vbin[0],vbin[1],std_e,mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2] )
    
    #         filename_str = 'z%2.2f.%s.jack%02d.corr' % ((vbin[0]+vbin[1])/2.,args.method,ijack)
    #         mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2=nbc_v7_stats.get_mc(cal_sim_select,cal_tru_select,None,use_calibration=True,use_weights=args.use_weights,filename_str=filename_str)

    #         std_e = np.std(cal_sim_select['e1'],ddof=1)
    #         list_jack_calibr.append( [ibin,vbin[0],vbin[1],std_e,mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2] )


    #     arr_jack = tktools.arr2rec(np.array(list_jack),dtype={'names': ["iz","z_min","z_max","std_e","m","std_m","c","std_c","m1","std_m1","m2","std_m2","c1","std_c1","c2","std_c2","pmm","std_pmm","pcc","std_pcc","pm1","std_pm1","pm2","std_pm2"], 'formats': ['i4']*1 + ['f8']*23 })
    #     arr_jack_calibr = tktools.arr2rec(np.array(list_jack_calibr),dtype={'names': ["iz","z_min","z_max","std_e","m","std_m","c","std_c","m1","std_m1","m2","std_m2","c1","std_c1","c2","std_c2","pmm","std_pmm","pcc","std_pcc","pm1","std_pm1","pm2","std_pm2"], 'formats': ['i4']*1 + ['f8']*23 })

    #     import pdb; pdb.set_trace()



    import pdb; pdb.set_trace()



def plot_bias_vs_redshift():

    pl.style.use('supermongo')

    import cPickle as pickle
    if args.use_weights == True:
        title = 'with weights'
        filename_pickle = os.path.join(args.output_dir,'figdata.bias_vs_redshift.with-weights.pp2')
    else:
        title = 'without weights'
        filename_pickle = os.path.join(args.output_dir,'figdata.bias_vs_redshift.pp2')
    arr_bias,arr_bias_calibr,z_bins =pickle.load(open(filename_pickle))
    logger.info('opened %s' % (filename_pickle))
    dx=0.01

    pl.figure()
    pl.errorbar( (arr_bias['z_min']+arr_bias['z_max'])/2 - dx , arr_bias['m1']-1, yerr=arr_bias['std_m1'], fmt='rs', label='$m_1$ uncalibrated', markeredgecolor='none')
    pl.errorbar( (arr_bias['z_min']+arr_bias['z_max'])/2 - dx , arr_bias['m2']-1, yerr=arr_bias['std_m2'], fmt='md', label='$m_2$ uncalibrated', markeredgecolor='none')
    pl.errorbar( (arr_bias_calibr['z_min']+arr_bias_calibr['z_max'])/2 -dx, arr_bias_calibr['m1']-1, yerr=arr_bias_calibr['std_m1'], fmt='bs', label='$m_1$ calibrated', markeredgecolor='none')
    pl.errorbar( (arr_bias_calibr['z_min']+arr_bias_calibr['z_max'])/2 -dx, arr_bias_calibr['m2']-1, yerr=arr_bias_calibr['std_m2'], fmt='cd', label='$m_2$ calibrated', markeredgecolor='none')
    pl.errorbar( (arr_bias_calibr['z_min']+arr_bias_calibr['z_max'])/2 +dx, arr_bias_calibr['mm1s']-1, yerr=arr_bias_calibr['std_mm1s'], fmt='ks', label='$m_1$ selection bias', markeredgecolor='none')
    pl.errorbar( (arr_bias_calibr['z_min']+arr_bias_calibr['z_max'])/2 +dx, arr_bias_calibr['mm2s']-1, yerr=arr_bias_calibr['std_mm2s'], fmt='kd', label='$m_2$ selection bias', markeredgecolor='none')
    
    # pl.errorbar( (arr_bias_calibr['z_min']+arr_bias_calibr['z_max'])/2 , arr_bias_calibr['m']-1, yerr=arr_bias_calibr['std_m'], fmt='b.', label='$m$ with NBC')
    # pl.errorbar( (arr_bias['z_min']+arr_bias['z_max'])/2 , arr_bias['m']-1, yerr=arr_bias['std_m'], fmt='r.', label='$m$ no NBC')
    # pl.errorbar( (arr_bias_calibr['z_min']+arr_bias_calibr['z_max'])/2 , arr_bias_calibr['mms']-1, yerr=arr_bias_calibr['std_mms'], fmt='k.', label='$m$ selection bias')

    pl.xlabel('COSMOS photometric redshift')
    pl.ylabel('m')
    pl.xticks(z_bins)
    pl.yticks([-0.09,-0.06,-0.03,0,0.03])
    pl.grid()
    pl.ylim([-0.12,0.035])
    pl.text(0.35,0.02,'im3shape',fontsize=17)
    pl.axhline(0,c='k')
    # pl.legend(loc='lower left',framealpha=0.0,frameon=False,fontsize=15)
    pl.legend(ncol=3,mode='expand',loc='lower center',framealpha=0.0,frameon=False,fontsize=14)
    plot_add_requirements(level=0.03,target=0)
    # plot_add_requirements(level=0.01,target=1)
    filename_fig = os.path.join(args.output_dir,('m_vs_z.weights%d.%s' % (args.use_weights,args.fig_format)).replace('..','.'))
    pl.savefig(filename_fig)
    logger.info('wrote %s' % (filename_fig))
    # pl.title(title)

    pl.figure()
    pl.errorbar( (arr_bias['z_min']+arr_bias['z_max'])/2 , arr_bias['pm1'], yerr=arr_bias['std_pm1'], fmt='rs', label=r'$\alpha_1$ uncalibrated', markeredgecolor='none')
    pl.errorbar( (arr_bias['z_min']+arr_bias['z_max'])/2 , arr_bias['pm2'], yerr=arr_bias['std_pm2'], fmt='md', label=r'$\alpha_2$ uncalibrated', markeredgecolor='none')
    pl.errorbar( (arr_bias_calibr['z_min']+arr_bias_calibr['z_max'])/2 , arr_bias_calibr['pm1'], yerr=arr_bias_calibr['std_pm1'], fmt='bs', label=r'$\alpha_1$ calibrated', markeredgecolor='none')
    pl.errorbar( (arr_bias_calibr['z_min']+arr_bias_calibr['z_max'])/2 , arr_bias_calibr['pm2'], yerr=arr_bias_calibr['std_pm2'], fmt='cd', label=r'$\alpha_2$ calibrated', markeredgecolor='none')
    pl.xlabel('COSMOS photometric redshift')
    pl.ylabel(r'$\alpha$')
    pl.xticks(z_bins)
    pl.grid()
    pl.ylim([-0.125,0.175])
    pl.text(0.35,0.14,'im3shape',fontsize=17)
    pl.axhline(0,c='k')
    # pl.legend(ncol=2, mode='expand', loc='lower left',framealpha=0.0,frameon=False,fontsize=15)
    pl.legend(ncol=2, mode='expand', loc='lower center',framealpha=0.0,frameon=False,fontsize=14)
    plot_add_requirements(level=0.05,target=0)
    filename_fig = os.path.join(args.output_dir,('a_vs_z.weights%d.%s' % (args.use_weights,args.fig_format)).replace('..','.'))
    pl.savefig(filename_fig)
    logger.info('wrote %s' % (filename_fig))
    # pl.title(title)


    pl.figure() 
    pl.errorbar( (arr_bias['z_min']+arr_bias['z_max'])/2 , arr_bias['mean_e1'], yerr=arr_bias['stdm_e1'], fmt='r.', label=r'mean_e1 no NBC')
    pl.errorbar( (arr_bias['z_min']+arr_bias['z_max'])/2 , arr_bias['mean_e2'], yerr=arr_bias['stdm_e2'], fmt='m.', label=r'mean_e2 no NBC')
    pl.errorbar( (arr_bias_calibr['z_min']+arr_bias_calibr['z_max'])/2 , arr_bias_calibr['mean_e1'], yerr=arr_bias_calibr['stdm_e1'], fmt='b.', label=r'mean_e1 with NBC')
    pl.errorbar( (arr_bias_calibr['z_min']+arr_bias_calibr['z_max'])/2 , arr_bias_calibr['mean_e2'], yerr=arr_bias_calibr['stdm_e2'], fmt='c.', label=r'mean_e2 with NBC')
    pl.xlabel('COSMOS photometric redshift')
    pl.ylabel(r'mean_e')
    pl.xticks(z_bins)
    pl.grid()
    pl.axhline(0,c='k')
    pl.legend(ncol=2, mode='expand', loc='upper center',framealpha=0.0,frameon=False)
    # pl.ylim([-0.1,0.3])
    plot_add_requirements(level=0.0005,target=0)
    pl.title(title)
    filename_fig = os.path.join(args.output_dir,('meane_vs_z.weights%d.%s' % (args.use_weights,args.fig_format)).replace('..','.'))
    pl.savefig(filename_fig)
    logger.info('wrote %s' % (filename_fig))

    pl.figure() 
    pl.errorbar( (arr_bias['z_min']+arr_bias['z_max'])/2 , arr_bias['mms'], yerr=arr_bias['std_mms'], fmt='r.', label=r'true sheared shape')
    pl.errorbar( (arr_bias['z_min']+arr_bias['z_max'])/2 , arr_bias['mms'], yerr=arr_bias['std_mms'], fmt='m.', label=r'true sheared shape')
    pl.xlabel('z')
    pl.ylabel('m')
    pl.xticks(z_bins)
    pl.grid()
    pl.axhline(1,c='k')
    pl.legend(ncol=2, mode='expand', loc='upper center',framealpha=0.0,frameon=False)
    # pl.ylim([-0.1,0.3])
    plot_add_requirements(level=0.0005,target=0)
    pl.title(title)
    filename_fig = os.path.join(args.output_dir,('ms_vs_z.weights%d.%s'% (args.use_weights,args.fig_format)).replace('..','.'))
    pl.savefig(filename_fig)
    logger.info('wrote %s' % (filename_fig))

    pl.show()

    res_sim, res_tru, res_des = load_selection()

    logger.info('calculating bias for the entire sample - using cuts from final selection')
    logger.info(selection_string_final_sim)
    logger.info(selection_string_final_des)
    cat_res = res_sim
    cat_tru = res_tru
    exec selection_string_final_sim
    select_sim = select.copy()
    cat_res = res_des
    exec selection_string_final_des
    select_des = select.copy()

    res_sim_final = res_sim[select_sim]
    res_tru_final = res_tru[select_sim]
    res_des_final = res_des[select_des]

    cal_sim_final = res_sim[select_sim]
    cal_tru_final = res_tru[select_sim]
    cal_des_final = res_des[select_des]

    logger.info('final selection SIM %d/%d %f',len(res_sim[select_sim]),len(res_sim), len(res_sim[select_sim])/float(len(res_sim)))
    logger.info('final selection DES %d/%d %f',len(res_des[select_des]),len(res_des), len(res_des[select_des])/float(len(res_des)))
    pl.figure(101) 
    pl.figure(201)
    pl.figure(401)
    list_bulge_fraction = []
    for ibin in range(1,len(z_bins)):

        logger.info('--------------------------------------------------')

        vbin = [z_bins[ibin-1],z_bins[ibin]]

        select = (res_tru_final['zphot'] > vbin[0]) & (res_tru_final['zphot'] < vbin[1])
        res_sim_select = res_sim_final[select]
        res_tru_select = res_tru_final[select]

        select = (cal_tru_final['zphot'] > vbin[0]) & (cal_tru_final['zphot'] < vbin[1])
        cal_sim_select = cal_sim_final[select]
        cal_tru_select = cal_tru_final[select]

        # dummy - ignore this, it's only here so that the code runs through
        cal_des_select = cal_des_final
        res_des_select = res_des_final

        bulge_fraction = np.sum(res_sim_select['disc_flux']==0)/float(len(res_sim_select))
        list_bulge_fraction.append( ( vbin , bulge_fraction ) )

        pl.figure(101)
        pl.hist(cal_sim_select['snr'],bins=np.logspace(0,2,100),histtype='step',label='z=[%2.2f - %2.2f]'%(vbin[0],vbin[1]),color=tt.get_colorscale(len(z_bins))[ibin],normed=True)
        pl.xlabel('SNR')
        pl.figure(201)
        pl.hist(cal_sim_select['radius'],bins=np.linspace(0,3,100),histtype='step',label='z=[%2.2f - %2.2f]'%(vbin[0],vbin[1]),color=tt.get_colorscale(len(z_bins))[ibin],normed=True)
        pl.xlabel('radius')
        pl.figure(202)
        pl.hist(cal_sim_select['mean_rgpp_rp'],bins=np.linspace(1,3,100),histtype='step',label='z=[%2.2f - %2.2f]'%(vbin[0],vbin[1]),color=tt.get_colorscale(len(z_bins))[ibin],normed=True)


    pl.figure(101)
    pl.hist(res_sim_final['snr'],bins=np.logspace(0,2,100),histtype='step',label='z=all',color='k',normed=True)
    pl.legend(framealpha=0.0,frameon=False)
    pl.ylim([0,0.15])
    filename_fig = os.path.join(args.output_dir,('hist.snr_vs_z.weights%d.%s' % (args.use_weights,args.fig_format)).replace('..','.'))
    pl.savefig(filename_fig)
    logger.info('wrote %s' % (filename_fig))

    pl.figure(201)
    pl.hist(res_sim_final['radius'],bins=np.linspace(0,3,100),histtype='step',label='z=all',color='k',normed=True)
    pl.legend(framealpha=0.0,frameon=False)
    filename_fig = os.path.join(args.output_dir,('hist.radius_vs_z.weights%d.%s' % (args.use_weights,args.fig_format)).replace('..','.'))
    pl.savefig(filename_fig)
    logger.info('wrote %s' % (filename_fig))
    

    pl.figure(202)
    pl.hist(res_sim_final['mean_rgpp_rp'],bins=np.linspace(1,3,100),histtype='step',label='z=all',color='k',normed=True)
    pl.legend(framealpha=0.0,frameon=False)
    filename_fig = os.path.join(args.output_dir,('hist.meanrgpprp_vs_z.weights%d.%s' % (args.use_weights,args.fig_format)).replace('..','.'))
    pl.savefig(filename_fig)
    logger.info('wrote %s' % (filename_fig))

    pl.figure(301)
    pl.hist(res_tru_final['zphot'],bins=np.linspace(0,1.5,100),histtype='step',normed=True)
    filename_fig = os.path.join(args.output_dir,('hist.z_phot.weights%d.%s' % (args.use_weights,args.fig_format)).replace('..','.'))
    pl.savefig(filename_fig)
    logger.info('wrote %s' % (filename_fig))

    pl.figure(31)
    arr_bulge_fraction = np.array(list_bulge_fraction)
    pl.plot(arr_bulge_fraction[:,1])
    filename_fig = os.path.join(args.output_dir,('hist.bulge_fraction_vs_z.weights%d.%s' % (args.use_weights,args.fig_format)).replace('..','.'))
    pl.savefig(filename_fig)
    logger.info('wrote %s' % (filename_fig))

    pl.show()

    import pdb; pdb.set_trace()







def get_mc_vs_snr():

    # cols_res=['coadd_objects_id','e1','e2','w','snr','mean_rgpp_rp','mean_psf_e1_sky','mean_psf_e2_sky','nbc_m1', 'nbc_m2', 'nbc_c1', 'nbc_c2']
    # cols_tru=['snr','psf_e1','psf_e2','id_shear','cosmos_mag_auto','g1_true','g2_true']
    # cols_des=['coadd_objects_id','e1','e2','w','snr','mean_rgpp_rp','mean_psf_e1_sky','mean_psf_e2_sky','nbc_m1', 'nbc_m2', 'nbc_c1', 'nbc_c2']
    # if args.use_calibration:
    #     logger.info('testing calibration columns')
    #     cols_res += ['nbc_m1', 'nbc_m2', 'nbc_c1', 'nbc_c2']
    #     cols_des += ['nbc_m1', 'nbc_m2', 'nbc_c1', 'nbc_c2']
    # else:
    #     logger.info('measuring bias')
    # res_sim,res_tru = nbc_v7_select.get_selection_sim(selection_string_model_sim,cols_res=cols_res,cols_tru=cols_tru,get_calibrated=args.use_calibration)

    res_sim, res_tru, res_des, selection_string_model_sim, selection_string_model_des = load_selection()

    logger.info('calculating bias for the entire sample - using cuts from final selection')
    logger.info(selection_string_final_sim)
    logger.info(selection_string_final_des)
    cat_res = res_sim
    cat_tru = res_tru
    exec selection_string_final_sim
    select_sim = select.copy()
    cat_res = res_des
    exec selection_string_final_des
    select_des = select.copy()

    res_sim_final = res_sim[select_sim]
    res_tru_final = res_tru[select_sim]
    res_des_final = res_des[select_des]

    load_selection()
 
    list_snr_centers = plotstools.get_bins_centers(list_snr_edges)

    list_bias = []

    # entire sample
    logger.info('calculating bias for the entire sample')
    filename_str = 'all.%s' % args.method
    mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2,mms,std_mms,ccs,std_ccs,mm1s,std_mm1s,cc1s,std_cc1s,mm2s,std_mm2s,cc2s,std_cc2s=nbc_v7_stats.get_mc(res_sim_final,res_tru_final,None,use_calibration=args.use_calibration,use_weights=args.use_weights,filename_str=filename_str)

    for isnr in range(1,len(list_snr_edges)):

            # select = (res_sim['snr'] > list_snr_edges[isnr-1]) & (res_sim['snr'] < list_snr_edges[isnr]) 
            select = (res_sim_final['snr'] > list_snr_edges[isnr-1]) 
            n_gals_sim = len(np.nonzero(select)[0])
            logger.info('isnr=%2d snr=[%2.2f %2.2f] n_gals_sim=%d' % (isnr,list_snr_edges[isnr-1],list_snr_edges[isnr],n_gals_sim))

            res_sim_select = res_sim_final[select]
            res_tru_select = res_tru_final[select]

            vsnr_mid = list_snr_centers[isnr-1]

            filename_str = 'm_vs_snr.snr%2.2f' % (vsnr_mid)
            # mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2,mms,std_mms,ccs,std_ccs,mm1s,std_mm1s,cc1s,std_cc1s,mm2s,std_mm2s,cc2s,std_cc2s=nbc_v7_stats.get_mc(res_sim_select,res_tru_select,None,use_calibration=args.use_calibration,use_weights=args.use_weights,filename_str=filename_str)
            mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2,mms,std_mms,ccs,std_ccs,mm1s,std_mm1s,cc1s,std_cc1s,mm2s,std_mm2s,cc2s,std_cc2s=nbc_v7_stats.get_mc(res_des_select,None,None,use_calibration=args.use_calibration,use_weights=args.use_weights,filename_str=filename_str)

            m_mean  = np.mean(res_sim_select['nbc_m'])
            c1_mean = np.mean(res_sim_select['nbc_c1'])
            c2_mean = np.mean(res_sim_select['nbc_c2'])

            std_e = np.std(res_sim_select['e1'],ddof=1)
            list_bias.append( [isnr,n_gals_sim,vsnr_mid,std_e,mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,mean_e1,mean_e2,stdm_e1,stdm_e2,m_mean,c1_mean,c2_mean] )


    import tktools;
    arr_bias = tktools.arr2rec(np.array(list_bias),dtype={'names': ["isnr","n_gals","snr","std_e","m","std_m","c","std_c","m1","std_m1","m2","std_m2","c1","std_c1","c2","std_c2","pmm","std_pmm","pcc","std_pcc",'mean_e1','mean_e2','stdm_e1','stdm_e2','m_mean','c1_mean','c2_mean'], 'formats': ['i4']*2 + ['f8']*25 } )

    tktools.save('plot_mc_vs_snr.cpickle',[arr_bias,args.method,selection_string_model_sim,args.use_weights,args.use_calibration],clobber=True)

def plot_mc_vs_snr():

    [arr_bias,method,selection_string_model_sim,use_weights,use_calibration] = tktools.load('plot_mc_vs_snr.cpickle')

    pl.figure(figsize=(20,15))
    pl.errorbar(arr_bias['snr'],arr_bias['m'] ,yerr=arr_bias['std_m'],  fmt='k-', label='m')
    pl.errorbar(arr_bias['snr'],arr_bias['m1'],yerr=arr_bias['std_m'],  fmt='r:', label='m1')
    pl.errorbar(arr_bias['snr'],arr_bias['m2'],yerr=arr_bias['std_m'],  fmt='m:', label='m2')
    pl.axhline(1,c='k')
    pl.xscale('log')
    pl.xlabel('snr')
    pl.ylabel('multiplicative bias m')
    pl.xticks(arr_bias['snr'])
    pl.gca().xaxis.set_major_formatter(pl.matplotlib.ticker.ScalarFormatter())
    plot_add_requirements()
    pl.xlim([5,250])
    pl.title('%s\n%s\n weights=%d calib/sens=%d' % (method,selection_string_model_sim,use_weights,use_calibration))
    pl.xlabel('snr')
    pl.ylabel('mean m bias')
    pl.legend(framealpha=0.0,frameon=False,loc='lower right')
    filename_fig = os.path.join(args.output_dir,'figs/plot_mc_vs_snr.m.weights%d.nbc%d.png'% (use_weights,use_calibration))
    pl.savefig(filename_fig)
    logger.info('wrote %s' % (filename_fig))
    
    pl.figure(figsize=(20,15))
    pl.errorbar(arr_bias['snr'],arr_bias['pmm'],yerr=arr_bias['std_pmm'],  fmt='k-' , label='alpha')
    pl.axhline(0,c='k')
    pl.xscale('log')
    pl.xticks(arr_bias['snr'])
    pl.gca().xaxis.set_major_formatter(pl.matplotlib.ticker.ScalarFormatter())
    plot_add_requirements()
    pl.xlim([5,250])
    pl.xlabel('snr')
    pl.ylabel('mean alpha bias')
    pl.title('%s\n%s\n weights=%d calib/sens=%d' % (method,selection_string_model_sim,use_weights,use_calibration))
    pl.legend(framealpha=0.0,frameon=False,loc='lower right')
    filename_fig = os.path.join(args.output_dir,'figs/plot_mc_vs_snr.alpha.weights%d.nbc%d.png'% (use_weights,use_calibration))
    pl.savefig(filename_fig)
    logger.info('wrote %s' % (filename_fig))

    pl.figure(figsize=(20,15))
    pl.axhline(0,c='k')
    pl.errorbar(arr_bias['snr'],arr_bias['mean_e1'],yerr=arr_bias['stdm_e1'],  fmt='-' , label='mean_e1')
    pl.errorbar(arr_bias['snr'],arr_bias['mean_e2'],yerr=arr_bias['stdm_e2'],  fmt='-' , label='mean_e2')
    pl.xscale('log')
    pl.xticks(arr_bias['snr'])
    pl.gca().xaxis.set_major_formatter(pl.matplotlib.ticker.ScalarFormatter())
    plot_add_requirements()
    pl.xlim([5,250])
    pl.xlabel('snr')
    pl.ylabel('mean_e (applied shear subtracted)')
    pl.title('%s\n%s\n weights=%d calib/sens=%d' % (method,selection_string_model_sim,use_weights,use_calibration))
    pl.legend(framealpha=0.0,frameon=False,loc='lower right')
    filename_fig = os.path.join(args.output_dir,'figs/plot_mc_vs_snr.meane.weights%d.nbc%d.png'% (use_weights,use_calibration))
    pl.savefig(filename_fig)
    logger.info('wrote %s' % (filename_fig))


    pl.figure(figsize=(20,15))
    pl.axhline(0,c='k')
    pl.plot(arr_bias['snr'],arr_bias['m_mean'], '-x' , label='mean_m')
    pl.xscale('log')
    pl.xticks(arr_bias['snr'])
    pl.gca().xaxis.set_major_formatter(pl.matplotlib.ticker.ScalarFormatter())
    pl.xlim([5,250])
    pl.title('%s\n%s\n weights=%d calib/sens=%d' % (method,selection_string_model_sim,use_weights,use_calibration))
    pl.legend(framealpha=0.0,frameon=False,loc='lower right')
    pl.xlabel('snr')
    pl.ylabel('mean m correction')
    filename_fig = os.path.join(args.output_dir,'figs/plot_mc_vs_snr.m_applied.weights%d.nbc%d.png'% (use_weights,use_calibration))
    pl.savefig(filename_fig)
    logger.info('wrote %s' % (filename_fig))

    pl.figure(figsize=(20,15))
    pl.axhline(0,c='k')
    pl.plot(arr_bias['snr'],arr_bias['c1_mean'], '-x' , label='mean_c1')
    pl.plot(arr_bias['snr'],arr_bias['c2_mean'], '-x' , label='mean_c2')
    pl.xscale('log')
    pl.xticks(arr_bias['snr'])
    pl.gca().xaxis.set_major_formatter(pl.matplotlib.ticker.ScalarFormatter())
    pl.xlim([5,250])
    pl.title('%s\n%s\n weights=%d calib/sens=%d' % (method,selection_string_model_sim,use_weights,use_calibration))
    pl.legend(framealpha=0.0,frameon=False,loc='lower right')
    pl.xlabel('snr')
    pl.ylabel('mean m correction')
    filename_fig = os.path.join(args.output_dir,'figs/plot_mc_vs_snr.c_applied.weights%d.nbc%d.png' % (use_weights,use_calibration))
    pl.savefig(filename_fig)
    logger.info('wrote %s' % (filename_fig))

    pl.show()

    import pdb; pdb.set_trace()

def save_selection():
    
    cols_res=config['cols_res']
    cols_tru=config['cols_tru']
    cols_des=config['cols_des']

    logger.info('testing calibration columns')
    cols_res += ['nbc_m', 'nbc_alpha', 'nbc_c1', 'nbc_c2']
    cols_des += ['nbc_m', 'nbc_alpha', 'nbc_c1', 'nbc_c2']

    logger.info('loading SIM results with following selection:')
    logger.info(selection_string_store_sim)
    res_sim,res_tru = nbc_v7_select.get_selection_sim(selection_string_store_sim,cols_res=cols_res,cols_tru=cols_tru)
    logger.info('loading DES results with following selection:')
    logger.info(selection_string_store_des)

    res_des = nbc_v7_select.get_selection_des(selection_string_store_des,cols=cols_des,n_files=args.n_des_files)

    if 'bulge_flux' in res_sim.dtype.names: 
        bulge_fraction_sim = np.sum(res_sim['bulge_flux']!=0)/float(len(res_sim))
        bulge_fraction_des = np.sum(res_des['bulge_flux']!=0)/float(len(res_des))
        logger.info('bulge_fraction SIM=%2.4f',bulge_fraction_sim)
        logger.info('bulge_fraction DES=%2.4f',bulge_fraction_des)

    tt.save(os.path.join(args.output_dir,'res_sim.fits'),res_sim,clobber=True)
    tt.save(os.path.join(args.output_dir,'res_tru.fits'),res_tru,clobber=True)
    tt.save(os.path.join(args.output_dir,'res_des.fits'),res_des,clobber=True)
    

def get_gold_cut(res_des):

    filename_flatcats='/Users/tomek/data/DES/sva1_gold_wl_flatcats_v4_info.fits'
    gold=tktools.load(filename_flatcats)
    select_gold = (gold['sva1_gold_flags']==0)&(gold['sva1_spte_flags']==0)&(gold['sva1_gold_mag_flags']==0)
    intersect = np.intersect1d(res_des['coadd_objects_id'],gold['coadd_objects_id'][select_gold])
    select_gold = np.in1d(res_des['coadd_objects_id'], intersect)

    return select_gold

def get_radec_split():

    res_sim, res_tru, res_des, _, _ = load_selection()

    cosmos=tktools.load('/Users/tomek/data/COSMOS/COSMOS_23.5_training_sample/real_galaxy_catalog_23.5.fits')

    import kmeans_radec
    from kmeans_radec import KMeans, kmeans_sample
    X = np.concatenate([cosmos['RA'].astype(np.float)[:,np.newaxis],cosmos['DEC'].astype(np.float)[:,np.newaxis]],axis=1)

    n_patches = 144
    logger.info('running kmeans')
    km = kmeans_sample(X, n_patches, maxiter=10, tol=1.0e-5); print("converged?",km.converged)
    logger.info('done')

    labels = km.labels[res_tru['id_cosmos'].astype(np.int)]

    pl.figure()
    pl.scatter(cosmos['RA'],cosmos['DEC'],c=km.labels)

    pl.figure()
    list_m_mean = []
    list_m_stdm = []
    list_mm_mean = []
    list_mm_stdm = []
    list_mc_mean = []
    list_mc_stdm = []

    for ip in range(n_patches):
        logger.info('patch %d',ip)
        select = labels==ip
        res_sim_select = res_sim[select] 
        res_tru_select = res_tru[select] 
        mean_m, stdm_m = np.mean(res_sim_select['nbc_m']), np.std(res_sim_select['nbc_m'],ddof=1)/np.sqrt(float(len(res_sim_select)))
        list_m_mean.append(mean_m)
        list_m_stdm.append(stdm_m)

        filename_str = 'patch%03d.nocal.png' % ip
        mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2,mms,std_mms,ccs,std_ccs,mm1s,std_mm1s,cc1s,std_cc1s,mm2s,std_mm2s,cc2s,std_cc2s=nbc_v7_stats.get_mc(res_sim_select,res_tru_select,None,use_calibration=False,use_weights=True,filename_str=filename_str)
        list_mm_mean.append(mm)
        list_mm_stdm.append(std_mm)

        filename_str = 'patch%03d.cal.png' % ip
        mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2,mms,std_mms,ccs,std_ccs,mm1s,std_mm1s,cc1s,std_cc1s,mm2s,std_mm2s,cc2s,std_cc2s=nbc_v7_stats.get_mc(res_sim_select,res_tru_select,None,use_calibration=True,use_weights=True,filename_str=filename_str)
        list_mc_mean.append(mm)
        list_mc_stdm.append(std_mm)

        pl.subplot(1,2,1)
        pl.hist(res_sim_select['snr'],bins=np.logspace(1,2,100),histtype='step',normed=True)
        pl.xscale('log')
        pl.subplot(1,2,2)
        pl.hist(res_sim_select['mean_rgpp_rp'],bins=np.linspace(1,2,100),histtype='step',normed=True)


    pl.figure()
    pl.errorbar(range(n_patches),np.array(list_mm_mean)-1, fmt='x', yerr=list_mm_stdm,c='b',label='calibrated')
    pl.errorbar(range(n_patches),np.array(list_mc_mean)-1, fmt='+', yerr=list_mc_stdm,c='r',label='uncalibrated')
    pl.errorbar(range(n_patches),np.array(list_m_mean),    fmt='.', yerr=list_m_stdm, c='g',label='correction (no weights used)')
    pl.axhline(0,color='k')
    pl.xlim([-1,n_patches+1])
    pl.ylim([-0.2,0.1])
    plot_add_requirements(level=0.01,target=0)
    pl.xlabel('patch number')
    pl.ylabel('m')
    pl.legend(framealpha=0.0,frameon=False, mode='expand',ncol=2,loc='lower right')
    pl.show()


    import pdb; pdb.set_trace()











    # n_side = 10
    # mean_ra = np.linspace([res_tru['ra'].max(),res_tru['ra'].min()],n_side)
    # mean_de = np.linspace([res_tru['de'].max(),res_tru['de'].min()],n_side)

    # select_n = (res_tru['ra'] > mean_ra) & (res_tru['de'] > mean_de) 
    # select_s = (res_tru['ra'] < mean_ra) & (res_tru['de'] < mean_de) 
    # select_e = 
    # select_w = 

    # m_n, std_m_n = np.mean(res_sim[select_n]['nbc_m']), np.std(res_sim[select_n]['nbc_m'],ddof=1)/float(len(res_sim[select_n]))
    # m_s, std_m_s = np.mean(res_sim[select_s]['nbc_m']), np.std(res_sim[select_s]['nbc_m'],ddof=1)/float(len(res_sim[select_w]))
    # m_e, std_m_e = np.mean(res_sim[select_e]['nbc_m']), np.std(res_sim[select_e]['nbc_m'],ddof=1)/float(len(res_sim[select_e]))
    # m_w, std_m_w = np.mean(res_sim[select_w]['nbc_m']), np.std(res_sim[select_w]['nbc_m'],ddof=1)/float(len(res_sim[select_w]))

    # print m_n, m_s, m_e, m_w


    # pl.show()

def test_shape_noise():


    res_sim, res_tru, res_des = load_selection()

    xi = res_tru['shape_e1'] + res_tru['shape_e2']*1j
    q = np.sqrt( (1 - np.abs(xi))/(1 + np.abs(xi)) )
    e = (1-q)/(1+q) * np.exp(1j*np.angle(xi))
    e1 = e.real
    e2 = e.imag
    e1,e2 = tt.distortion_to_shear(res_tru['shape_e1'] , res_tru['shape_e2'])
    e_rot = e*np.exp(2*res_tru['rotation_angle']*1j)
    e1_rot = e_rot.real
    e2_rot = e_rot.imag

    ng_g1 = res_sim['e1']
    ng_g2 = res_sim['e2']

    pl.figure()
    pl.scatter(xi[:10000].real, e1[:10000], c=np.abs(e)[:10000] ); 
    pl.xlabel('xi1')
    pl.ylabel('e1, color=abs(e)')
    pl.colorbar()

    pl.figure()
    pl.scatter(xi[:10000].imag, e2[:10000], c=np.abs(e)[:10000] ); 
    pl.xlabel('xi2')
    pl.ylabel('e2, color=abs(e)')
    pl.colorbar()

    pl.figure()
    pl.scatter(  ng_g1[:10000], e_rot.real[:10000] ); 
    pl.xlabel('g ngmix')
    pl.ylabel('e_rot')

    pl.figure()
    pl.scatter(  ng_g2[:10000], e_rot.imag[:10000] ); 
    pl.xlabel('g ngmix')
    pl.ylabel('e_rot')

    # cut here
    select = res_sim['snr'] < 4000
    n_s4000 = np.sum(select)

    np.mean( ng_g1 )
    np.mean( ng_g2 )

    np.mean( e1_rot )
    np.mean( e2_rot )

    filename_str = 'ng_high_snr_s4000_ngmix'
    filename_config = 'nbc-sersics.yaml'
    mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2,mms,std_mms,ccs,std_ccs,mm1s,std_mm1s,cc1s,std_cc1s,mm2s,std_mm2s,cc2s,std_cc2s = get_mc(res_sim[select],res_tru[select],res_des,use_calibration=False,use_weights=False,filename_str=filename_str)
  
    res_sim['e1'] = res_tru['sheared_e1']
    res_sim['e2'] = res_tru['sheared_e2']

    filename_str = 'ng_high_snr_s4000_intrinsic'
    mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2,mms,std_mms,ccs,std_ccs,mm1s,std_mm1s,cc1s,std_cc1s,mm2s,std_mm2s,cc2s,std_cc2s = get_mc(res_sim[select],res_tru[select],res_des,use_calibration=False,use_weights=False,filename_str=filename_str)

    mm_s4000 = mm

    select = res_sim['snr'] > 4000
    n_g4000 = np.sum(select)

    np.mean( ng_g1 )
    np.mean( ng_g2 )

    np.mean( e1_rot )
    np.mean( e2_rot )

    filename_str = 'ng_high_snr_g4000_ngmix'
    filename_config = 'nbc-sersics.yaml'
    mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2,mms,std_mms,ccs,std_ccs,mm1s,std_mm1s,cc1s,std_cc1s,mm2s,std_mm2s,cc2s,std_cc2s = get_mc(res_sim[select],res_tru[select],res_des,use_calibration=False,use_weights=False,filename_str=filename_str)
  
    res_sim['e1'] = res_tru['sheared_e1']
    res_sim['e2'] = res_tru['sheared_e2']

    filename_str = 'ng_high_snr_g4000_intrinsic'
    mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2,mms,std_mms,ccs,std_ccs,mm1s,std_mm1s,cc1s,std_cc1s,mm2s,std_mm2s,cc2s,std_cc2s = get_mc(res_sim[select],res_tru[select],res_des,use_calibration=False,use_weights=False,filename_str=filename_str)

    mm_g4000 = mm

    print mm_g4000, n_g4000
    print mm_s4000, n_s4000

    n_total = len(res_sim)

    print (mm_g4000/float(n_total))*n_g4000 + (mm_s4000/float(n_total))*n_s4000

    pl.show()



def get_selection_bias():

    res_sim, res_tru, res_des = load_selection()

    select = np.isfinite(res_tru['sheared_e1']) & np.isfinite(res_tru['sheared_e2'])

    res_sim = res_sim[select]
    res_tru = res_tru[select]

    # entire sample
    logger.info('calculating bias for the entire sample - using cuts from model selection')
    logger.info(selection_string_model_sim)
    logger.info(selection_string_model_des)
    cat_res = res_sim
    cat_tru = res_tru
    exec selection_string_model_sim
    select_sim = select.copy()
    cat_res = res_des
    exec selection_string_model_des
    select_des = select.copy()
    filename_str = 'all.%s' % args.method
    res_sim = res_sim[select_sim]
    res_tru = res_tru[select_sim]
    res_des = res_des[select_des]

    # res_sim = res_sim[::10]
    # res_tru = res_tru[::10]
    # res_sim['e1'] = res_tru['sheared_e1']
    # res_sim['e2'] = res_tru['sheared_e2']
    # res_sim['mean_psf_e1_sky'] = res_tru['psf_e1']
    # res_sim['mean_psf_e2_sky'] = res_tru['psf_e2']

    snr_cuts = np.linspace(10,30,6)
    rgp_cuts = np.linspace(1.15,1.3,6)
    hlr_cuts = np.linspace(0.1,1.1,2)

    list_snr_mms = []
    list_rgp_mms = []
    list_hlr_mms = []

    default_min_rgpprp = 1.15
    default_min_snr = 1.15

    for i1,v1 in enumerate(snr_cuts):

            select = (res_sim['snr']>v1) & (res_sim['mean_rgpp_rp']>default_min_rgpprp) 
            res_sim_final = res_sim[select]
            res_tru_final = res_tru[select]
            
            filename_str = 'selection_bias.snr>%2.2f' % v1
            mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2,mms,std_mms,ccs,std_ccs,mm1s,std_mm1s,cc1s,std_cc1s,mm2s,std_mm2s,cc2s,std_cc2s=nbc_v7_stats.get_mc(res_sim_final,res_tru_final,res_des,use_calibration=args.use_calibration,use_weights=args.use_weights,filename_str=filename_str,correct_selection_bias=False)

            list_snr_mms.append([v1,mms,std_mms,mm1s,std_mm1s,mm2s,std_mm2s])
            # list_snr_mms.append([v1,mm,std_mm,mm1,std_mm1,mm2,std_mm2])

    for i1,v1 in enumerate(rgp_cuts):

            select = (res_sim['mean_rgpp_rp']>v1) & (res_sim['mean_rgpp_rp']>default_min_snr) 
            res_sim_final = res_sim[select]
            res_tru_final = res_tru[select]

            filename_str = 'selection_bias.mean_rgpp_rp>%2.2f' % v1
            mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2,mms,std_mms,ccs,std_ccs,mm1s,std_mm1s,cc1s,std_cc1s,mm2s,std_mm2s,cc2s,std_cc2s=nbc_v7_stats.get_mc(res_sim_final,res_tru_final,res_des,use_calibration=args.use_calibration,use_weights=args.use_weights,filename_str=filename_str,correct_selection_bias=False)

            list_rgp_mms.append([v1,mms,std_mms,mm1s,std_mm1s,mm2s,std_mm2s])
            # list_rgp_mms.append([v1,mm,std_mm,mm1,std_mm1,mm2,std_mm2])

    for i1,v1 in enumerate(hlr_cuts):

            select = (res_tru['sf_hlr']>v1) 
            res_sim_final = res_sim[select]
            res_tru_final = res_tru[select]

            filename_str = 'selection_bias.truehlr>%2.2f' % v1
            mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2,mms,std_mms,ccs,std_ccs,mm1s,std_mm1s,cc1s,std_cc1s,mm2s,std_mm2s,cc2s,std_cc2s=nbc_v7_stats.get_mc(res_sim_final,res_tru_final,res_des,use_calibration=args.use_calibration,use_weights=args.use_weights,filename_str=filename_str,correct_selection_bias=False)

            list_hlr_mms.append([v1,mms,std_mms,mm1s,std_mm1s,mm2s,std_mm2s])
            # list_hlr_mms.append([v1,mm,std_mm,mm1,std_mm1,mm2,std_mm2])

    snr_mms = tktools.arr2rec(np.array(list_snr_mms),dtype={'names': ['snr','mms','mms_std','mm1s','mm1s_std','mm2s','mm2s_std'], 'formats': ['f8']*7 })
    rgp_mms = tktools.arr2rec(np.array(list_rgp_mms),dtype={'names': ['rgp','mms','mms_std','mm1s','mm1s_std','mm2s','mm2s_std'], 'formats': ['f8']*7 })
    hlr_mms = tktools.arr2rec(np.array(list_hlr_mms),dtype={'names': ['hlr','mms','mms_std','mm1s','mm1s_std','mm2s','mm2s_std'], 'formats': ['f8']*7 })

    pickle_dict = {'mms_snr':snr_mms, 'mms_rgp':rgp_mms, 'mms_hlr':hlr_mms}

    import cPickle as pickle
    filename_pickle = os.path.join(args.output_dir,'selection_bias_vs_cuts_table.nbc%d.weights%d.cpickle' % (args.use_calibration, args.use_weights))
    pickle.dump(pickle_dict,open(filename_pickle,'w'),protocol=2)
    logger.info('wrote %s' % (filename_pickle))
    



def plot_selection_bias():

    filename_pickle = os.path.join(args.output_dir,'selection_bias_vs_cuts_table.nbc%d.weights%d.cpickle' % (args.use_calibration, args.use_weights))
    import cPickle as pickle
    pickle_dict = pickle.load(open(filename_pickle))
    logger.info('loaded %s', filename_pickle)
    snr_mms = pickle_dict['mms_snr']
    rgp_mms = pickle_dict['mms_rgp']
    hlr_mms = pickle_dict['mms_hlr']

    import pdb; pdb.set_trace()


    pl.figure()
    pl.errorbar(snr_mms['snr'],snr_mms['mms'],yerr=snr_mms['mms_std'],fmt='.')
    pl.errorbar(snr_mms['snr'],snr_mms['mm1s'],yerr=snr_mms['mm1s_std'],fmt='+:')
    pl.errorbar(snr_mms['snr'],snr_mms['mm2s'],yerr=snr_mms['mm2s_std'],fmt='x--')
    pl.xlabel('round snr lower cut')
    pl.ylabel('multiplicative bias (from true sheared shape)')
    # pl.ylabel('multiplicative bias m from im3shape shapes')
    # pl.xlim([-1,21])
    # pl.ylim([0.9,1.1])
    pl.axhline(1,color='k')
    pl.title('selection biases when cutting on SNR')
    plot_add_requirements(level=0.03,target=1)

    pl.figure()
    pl.errorbar(rgp_mms['rgp'],rgp_mms['mms'],yerr=rgp_mms['mms_std'],fmt='.')
    # pl.errorbar(rgp_mms['rgp'],rgp_mms['mm1s'],yerr=rgp_mms['mm1s_std'],fmt='+:')
    # pl.errorbar(rgp_mms['rgp'],rgp_mms['mm2s'],yerr=rgp_mms['mm2s_std'],fmt='x--')
    pl.xlabel('radius lower cut')
    pl.ylabel('multiplicative bias (from true sheared shape)')
    # pl.ylabel('multiplicative bias m from im3shape shapes')
    pl.axhline(1,color='k')
    # pl.xlim([0.9,1.6])
    # pl.ylim([0.9,1.2])
    pl.title('selection biases when cutting on Rgp/Rp')
    plot_add_requirements(level=0.03,target=1)

    pl.figure()
    pl.errorbar(hlr_mms['hlr'],hlr_mms['mms'],yerr=hlr_mms['mms_std'],fmt='.')
    # pl.errorbar(hlr_mms['rgp'],hlr_mms['mm1s'],yerr=hlr_mms['mm1s_std'],fmt='+:')
    # pl.errorbar(hlr_mms['rgp'],hlr_mms['mm2s'],yerr=hlr_mms['mm2s_std'],fmt='x--')
    pl.xlabel('true hlr lower cut')
    pl.ylabel('multiplicative bias (from true sheared shape)')
    # pl.ylabel('multiplicative bias m from im3shape shapes')
    pl.axhline(1,color='k')
    pl.xlim([0.,1.2])
    pl.ylim([0.9,1.1])
    pl.title('selection biases when cutting on hlr')
    plot_add_requirements(level=0.03,target=1)
    
    pl.show()


    import pdb; pdb.set_trace()


    pass

def plot_mc_vs_true_params():

       # pl.style.use('supermongo')

    import fitsio
    if args.use_calibration:
        filename_table_bias = os.path.join(args.output_dir,'get_mc_vs_true_params.calibrated.fits')
    else:
        filename_table_bias = os.path.join(args.output_dir,'get_mc_vs_true_params.fits')


    bias_table = tktools.load(filename_table_bias)

    n_psf_bins = len(np.unique(bias_table['ipsf']))
    colorscale=plotstools.get_colorscale(n_psf_bins,cmap_name='rainbow')

    param_snr = 'cosmos_mag_auto'
    param_size = 'sf_hlr'

    n_upsample = 200

    for ipsf in np.unique(bias_table['ipsf']):

        select = bias_table['ipsf'] == ipsf
        bt1=bias_table[select]
        label = '%s=%2.2f-%2.2f'%(param_size,bt1['vpsf_min'][0],bt1['vpsf_max'][0])
        snr_mid = bt1['vsnr_mid']
        rgp_mid = bt1['vpsf_mid']
        logger.debug('rgp_mid = %2.4f' % rgp_mid[0])

        pl.figure(1)
        pl.xlabel(param_snr)
        pl.ylabel('m')
        pl.errorbar(snr_mid,bt1['m']-1,yerr=bt1['std_m'],label=label,fmt='o',c=colorscale[ipsf-1])
        
        # snr_pred =  np.linspace(bt1['vsnr_min'].min(),bt1['vsnr_max'].max(),n_upsample)
        # fitdata =fitting.fit2(snr_mid, bt1['m']-1, s=bt1['std_m'],expand=polynomial_basis)       
        # p, s = fitting.predict2(snr_pred,fitdata)
        # # pl.plot(snr_pred,p,c=colorscale[ipsf-1],lw=np.sum(bt1['n_gals'])/1e6)
        # pl.plot(snr_pred,p,c=colorscale[ipsf-1])

    pl.figure(1)
    pl.axhline(0,c='k')
    pl.grid()
    pl.legend(framealpha=0.0,frameon=False,loc='lower right',fontsize=14)
    pl.gca().xaxis.set_major_formatter(pl.matplotlib.ticker.ScalarFormatter())
    plot_add_requirements(level=0.03,target=0)
    plot_add_requirements(level=0.01,target=0)

    pl.show()

    import pdb; pdb.set_trace()

def plot_mc_vs_psf_size():


    res_sim, res_tru, res_des = load_selection()  

    pl.style.use('supermongo')
    filename_data = os.path.join(args.output_dir,'figdata.get_mc_vs_psf_size.weights%d.cpickle' % (args.use_weights))
    import cPickle as pickle
    pickle_dict = pickle.load(open(filename_data))
    arr_bias = pickle_dict['arr_bias']
    arr_bias_calibr = pickle_dict['arr_bias_calibr'] 

    pl.figure()
    pl.errorbar( arr_bias['psf_min'], arr_bias['m1'], yerr=arr_bias['std_m1'], fmt='r.', label='$m_1$ no NBC')
    pl.errorbar( arr_bias['psf_min'], arr_bias['m2'], yerr=arr_bias['std_m2'], fmt='m.', label='$m_2$ no NBC')
    pl.errorbar( arr_bias_calibr['psf_min'] , arr_bias_calibr['m1'], yerr=arr_bias_calibr['std_m1'], fmt='b.', label='$m_1$ with NBC')
    pl.errorbar( arr_bias_calibr['psf_min'] , arr_bias_calibr['m2'], yerr=arr_bias_calibr['std_m2'], fmt='c.', label='$m_2$ with NBC')
    pl.xlabel('z')
    pl.ylabel('m')
    pl.grid()
    pl.ylim([0.9,1.05])
    pl.axhline(1,c='k')
    pl.legend(ncol=2, mode='expand', loc='lower center',framealpha=0.0,frameon=False)
    plot_add_requirements(level=0.03,target=1)
    # plot_add_requirements(level=0.01,target=1)
    filename_fig = os.path.join(args.output_dir,('m_vs_psffwhm.weights%d.%s' % (args.use_weights,args.fig_format)).replace('..','.'))
    pl.savefig(filename_fig)
    logger.info('wrote %s' % (filename_fig))
    # pl.title(title)

    pl.figure()
    pl.errorbar( arr_bias['psf_min'] , arr_bias['pm1'], yerr=arr_bias['std_pm1'], fmt='r.', label=r'$\alpha_1$ no NBC')
    pl.errorbar( arr_bias['psf_min'] , arr_bias['pm2'], yerr=arr_bias['std_pm2'], fmt='m.', label=r'$\alpha_2$ no NBC')
    pl.errorbar( arr_bias_calibr['psf_min'] , arr_bias_calibr['pm1'], yerr=arr_bias_calibr['std_pm1'], fmt='b.', label=r'$\alpha_1$ with NBC')
    pl.errorbar( arr_bias_calibr['psf_min'] , arr_bias_calibr['pm2'], yerr=arr_bias_calibr['std_pm2'], fmt='c.', label=r'$\alpha_2$ with NBC')
    pl.xlabel('z')
    pl.ylabel(r'$\alpha$')
    # pl.xticks(z_bins)
    pl.grid()
    pl.axhline(0,c='k')
    pl.legend(ncol=2, mode='expand', loc='lower center',framealpha=0.0,frameon=False)
    pl.ylim([-0.1,0.2])
    plot_add_requirements(level=0.05,target=0)
    filename_fig = os.path.join(args.output_dir,('a_vs_psffwhm.weights%d.%s' % (args.use_weights,args.fig_format)).replace('..','.'))
    pl.savefig(filename_fig)
    logger.info('wrote %s' % (filename_fig))
    # pl.title(title)

    pl.show()

def get_mc_vs_psf_size():

    res_sim, res_tru, res_des = load_selection()

    # list_psf_edges = np.array([ 0.8,  1.0,  1.2,  1.4]) - 0.05
    list_psf_edges = np.array([ 0.75,  1.05,  1.25,  1.45])

    param1 = 'snr'
    param2 = 'mean_psf_fwhm'

    # model
    logger.info('calculating bias for the model - using cuts from model selection')
    logger.info(selection_string_final_sim)
    cat_res = res_sim; cat_tru = res_tru
    exec selection_string_final_sim
    select_final_sim = select.copy()

    filename_str = 'all.%s' % args.method
    logger.info('final selection SIM %d/%d %f',len(res_sim[select_final_sim]),len(res_sim), len(res_sim[select_final_sim])/float(len(res_sim)))
    res_sim_model = res_sim[select_final_sim]
    res_tru_model = res_tru[select_final_sim]

    list_bias = []
    list_bias_calibr = []
  
    warnings.warn('using %s and %s' % (param1,param2) )

    isnr=0
    for ipsf in range(1,len(list_psf_edges)):

            select = ( res_tru_model['psf_fwhm'] < list_psf_edges[ipsf] ) & ( res_tru_model['psf_fwhm'] > list_psf_edges[ipsf-1] )
            n_gals_sim = len(np.nonzero(select)[0])
            if n_gals_sim < 100:
                import pdb; pdb.set_trace()
            res_sim_select = res_sim_model[select]
            res_tru_select = res_tru_model[select]

            vpsf_mid = np.median( res_tru_select['psf_fwhm'] )
            vsnr_mid = np.median(res_sim_select[param1])
            vpsf_min = list_psf_edges[ipsf-1]
            vsnr_min = 0
            vpsf_max = list_psf_edges[ipsf]
            vsnr_max = 10000

            n_unique = len(np.unique(res_tru_select['id_cosmos'].astype(np.float)))
            if 'disc_flux' in res_sim_select.dtype.names: 
                bulge_fraction = np.sum(res_sim_select['disc_flux']==0)/float(len(res_sim_select))
            else:
                bulge_fraction = 0


            if ('ra_as'  in res_sim_select) & ('dec_as' in res_sim_select): 
                x1 = np.mean(res_sim_select['ra_as']);
                x2 = np.mean(res_sim_select['dec_as']);
                std_x1 = np.std(res_sim_select['ra_as'],ddof=1)/np.sqrt(len(res_sim_select));
                std_x2 = np.std(res_sim_select['dec_as'],ddof=1)/np.sqrt(len(res_sim_select));
                desx1 = np.mean(res_des_select['ra_as']);
                desx2 = np.mean(res_des_select['dec_as']);
                std_desx1 = np.std(res_des_select['ra_as'],ddof=1)/np.sqrt(len(res_des_select));
                std_desx2 = np.std(res_des_select['dec_as'],ddof=1)/np.sqrt(len(res_des_select));
            else:
                x1 = np.random.randn()
                x2 = np.random.randn()
                std_x1 = np.random.randn()
                std_x2 = np.random.randn()
                desx1 = np.random.randn()
                desx2 = np.random.randn()
                std_desx1 = np.random.randn()
                std_desx2 = np.random.randn()

            ibin = ipsf
            vbin = [vpsf_mid,vpsf_mid]

            current_psf_fwhm = np.median( res_tru_select['psf_fwhm'])
            logger.info('------------------- psf_fwhm %2.2f' % current_psf_fwhm )
            logger.info('ipsf=%2d isnr=%2d snr_mid=%2.2f size_mid=%2.2f mag=[%2.2f %2.2f] hlr=[%2.2f %2.2f] n_gals_sim=%d n_unique=%d bulge_fraction=%2.2f' % (ipsf,isnr,vsnr_mid,vpsf_mid,vsnr_min,vsnr_max,vpsf_min,vpsf_max,n_gals_sim,n_unique,bulge_fraction))
            std_e = np.std(res_sim_select['e1'],ddof=1)
            logger.info(selection_string_final_sim)

            filename_str = 'psffwhm=%2.2f.nonbc' % (current_psf_fwhm)

            mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2,mms,std_mms,ccs,std_ccs,mm1s,std_mm1s,cc1s,std_cc1s,mm2s,std_mm2s,cc2s,std_cc2s=nbc_v7_stats.get_mc(res_sim_select,res_tru_select,None,use_calibration=False,use_weights=args.use_weights,filename_str=filename_str,correct_selection_bias=config['correct_selection_bias'])
            list_bias.append( [ibin,vbin[0],vbin[1],std_e,mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2,mms,std_mms,ccs,std_ccs,mm1s,std_mm1s,cc1s,std_cc1s,mm2s,std_mm2s,cc2s,std_cc2s] )

            filename_str = 'psffwhm=%2.2f.nbc' % (current_psf_fwhm)
            mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2,mms,std_mms,ccs,std_ccs,mm1s,std_mm1s,cc1s,std_cc1s,mm2s,std_mm2s,cc2s,std_cc2s=nbc_v7_stats.get_mc(res_sim_select,res_tru_select,None,use_calibration=True ,use_weights=args.use_weights,filename_str=filename_str,correct_selection_bias=config['correct_selection_bias'])

            list_bias_calibr.append( [ibin,vbin[0],vbin[1],std_e,mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2,mms,std_mms,ccs,std_ccs,mm1s,std_mm1s,cc1s,std_cc1s,mm2s,std_mm2s,cc2s,std_cc2s] )
            

    arr_bias = tktools.arr2rec(np.array(list_bias),dtype={'names': ["ipsf","psf_min","psf_max","std_e","m","std_m","c","std_c","m1","std_m1","m2","std_m2","c1","std_c1","c2","std_c2","pmm","std_pmm","pcc","std_pcc","pm1","std_pm1","pm2","std_pm2","mean_e1","mean_e2","stdm_e1","stdm_e2","mms","std_mms","ccs","std_ccs","mm1s","std_mm1s","cc1s","std_cc1s","mm2s","std_mm2s","cc2s","std_cc2s"], 'formats': ['i4']*1 + ['f8']*39 })
    arr_bias_calibr = tktools.arr2rec(np.array(list_bias_calibr),dtype={'names': ["ipsf","psf_min","psf_max","std_e","m","std_m","c","std_c","m1","std_m1","m2","std_m2","c1","std_c1","c2","std_c2","pmm","std_pmm","pcc","std_pcc","pm1","std_pm1","pm2","std_pm2","mean_e1","mean_e2","stdm_e1","stdm_e2","mms","std_mms","ccs","std_ccs","mm1s","std_mm1s","cc1s","std_cc1s","mm2s","std_mm2s","cc2s","std_cc2s"], 'formats': ['i4']*1 + ['f8']*39 })

    filename_table_bias = os.path.join(args.output_dir,'figdata.get_mc_vs_psf_size.weights%d.cpickle' % (args.use_weights))
    import cPickle as pickle
    pickle_dict = { 'arr_bias' : arr_bias, 'arr_bias_calibr' : arr_bias_calibr }
    pickle.dump(pickle_dict,open(filename_table_bias,'w'))
    logger.info('saved %s',filename_table_bias)

    import pdb; pdb.set_trace()



def get_mc_vs_true_params():

    res_sim, res_tru, res_des = load_selection()  
    
    list_hlr_edges = [0,3]
    list_mag_edges = [17,19,21,21.5,21.75,22,22.1,22.2,22.3,22.4,22.5,22.6,22.7,22.8,22.9,23,24]

    param_snr = 'cosmos_mag_auto'
    param_size = 'sf_hlr'

    # model
    logger.info('calculating bias for the model - using cuts from model selection')
    logger.info(selection_string_final_sim)
    cat_res = res_sim
    cat_tru = res_tru
    exec selection_string_final_sim
    select_final_sim = select.copy()

    filename_str = 'all.%s' % args.method
    logger.info('final selection SIM %d/%d %f',len(res_sim[select_final_sim]),len(res_sim), len(res_sim[select_final_sim])/float(len(res_sim)))
    # res_sim_model = res_sim[select_final_sim]
    # res_tru_model = res_tru[select_final_sim]
    res_sim_model = res_sim
    res_tru_model = res_tru

    list_bias = []

    correct_selection_bias = False
     
    warnings.warn('using %s and %s' % (param_snr,param_size) )

    for ipsf in range(1,len(list_hlr_edges)):
        for isnr in range(1,len(list_mag_edges)):

            select = (res_tru_model[param_snr] > list_mag_edges[isnr-1]) & (res_tru_model[param_snr] < list_mag_edges[isnr]) & (res_tru_model[param_size] > list_hlr_edges[ipsf-1]) & (res_tru_model[param_size] < list_hlr_edges[ipsf])
            n_gals_sim = len(np.nonzero(select)[0])
            if n_gals_sim < 100:
                import pdb; pdb.set_trace()
            res_sim_select = res_sim_model[select]
            res_tru_select = res_tru_model[select]

            vpsf_mid = np.median(res_tru_select[param_size])
            vsnr_mid = np.median(res_tru_select[param_snr])
            vpsf_min = list_hlr_edges[ipsf-1]
            vsnr_min = list_mag_edges[isnr-1]
            vpsf_max = list_hlr_edges[ipsf]
            vsnr_max = list_mag_edges[isnr]

            n_unique = len(np.unique(res_tru_select['id_cosmos'].astype(np.float)))
            if 'disc_flux' in res_sim_select.dtype.names: 
                bulge_fraction = np.sum(res_sim_select['disc_flux']==0)/float(len(res_sim_select))
            else:
                bulge_fraction = 0


            if ('ra_as'  in res_sim_select) & ('dec_as' in res_sim_select): 
                x1 = np.mean(res_sim_select['ra_as']);
                x2 = np.mean(res_sim_select['dec_as']);
                std_x1 = np.std(res_sim_select['ra_as'],ddof=1)/np.sqrt(len(res_sim_select));
                std_x2 = np.std(res_sim_select['dec_as'],ddof=1)/np.sqrt(len(res_sim_select));
                desx1 = np.mean(res_des_select['ra_as']);
                desx2 = np.mean(res_des_select['dec_as']);
                std_desx1 = np.std(res_des_select['ra_as'],ddof=1)/np.sqrt(len(res_des_select));
                std_desx2 = np.std(res_des_select['dec_as'],ddof=1)/np.sqrt(len(res_des_select));
            else:
                x1 = np.random.randn()
                x2 = np.random.randn()
                std_x1 = np.random.randn()
                std_x2 = np.random.randn()
                desx1 = np.random.randn()
                desx2 = np.random.randn()
                std_desx1 = np.random.randn()
                std_desx2 = np.random.randn()

            logger.info('ipsf=%2d isnr=%2d snr_mid=%2.2f size_mid=%2.2f mag=[%2.2f %2.2f] hlr=[%2.2f %2.2f] n_gals_sim=%d n_unique=%d bulge_fraction=%2.2f' % (ipsf,isnr,vsnr_mid,vpsf_mid,vsnr_min,vsnr_max,vpsf_min,vpsf_max,n_gals_sim,n_unique,bulge_fraction))

            logger.info(selection_string_final_sim)
            filename_str = 'mag=[%2.2f,%2.2f].hlr=[%2.2f,%2.2f]' % (vsnr_min,vsnr_max,vpsf_min,vpsf_max)
            mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2,mms,std_mms,ccs,std_ccs,mm1s,std_mm1s,cc1s,std_cc1s,mm2s,std_mm2s,cc2s,std_cc2s=nbc_v7_stats.get_mc(res_sim_select,res_tru_select,None,use_calibration=args.use_calibration,use_weights=args.use_weights,filename_str=filename_str,correct_selection_bias=correct_selection_bias)

            std_e = np.std(res_sim_select['e1'],ddof=1)
            list_bias.append( [ipsf,isnr,n_gals_sim,n_unique,vpsf_mid,vsnr_mid,vpsf_min,vpsf_max,vsnr_min,vsnr_max,std_e,mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,bulge_fraction,x1,std_x1,x2,std_x2,desx1,std_desx1,desx2,std_desx2,mms,std_mms,ccs,std_ccs,mm1s,std_mm1s,cc1s,std_cc1s,mm2s,std_mm2s,cc2s,std_cc2s] )


    arr_bias = tktools.arr2rec(np.array(list_bias),dtype={'names': ["ipsf","isnr","n_gals",'n_unique',"vpsf_mid","vsnr_mid","vpsf_min","vpsf_max","vsnr_min","vsnr_max","std_e","m","std_m","c","std_c","m1","std_m1","m2","std_m2","c1","std_c1","c2","std_c2","pmm","std_pmm","pcc","std_pcc","bulge_fraction","x1","std_x1","x2","std_x2","desx1","std_desx1","desx2","std_desx2",'mms','std_mms','ccs','std_ccs','mm1s','std_mm1s','cc1s','std_cc1s','mm2s','std_mm2s','cc2s','std_cc2s'], 'formats': ['i4']*4 + ['f8']*44 })

    if args.use_calibration:
        filename_table_bias = os.path.join(args.output_dir,'get_mc_vs_true_params.calibrated.fits')
    else:
        filename_table_bias = os.path.join(args.output_dir,'get_mc_vs_true_params.fits')
    import pyfits
    pyfits.writeto(filename_table_bias,arr_bias,clobber=True)
    logger.info('saved %s',filename_table_bias)

    import pdb; pdb.set_trace()


def get_calibration():


    res_sim, res_tru, res_des = load_selection()


    if 'param_split1' in config:
        param_snr  = config['param_split1'] 
        param_size = config['param_split2']

        if (param_snr not in res_sim.dtype.names) and (param_snr in res_tru.dtype.names):
            logger.info('found %s in truth catalog, adding to results catalog' % param_snr)
            res_sim = add_col(res_sim,param_size,res_tru[param_size])

        if (param_size not in res_sim.dtype.names) and (param_size in res_tru.dtype.names):
            logger.info('found %s in truth catalog, adding to results catalog' % param_size)
            res_sim = res_sim[config['cols_res']]
            res_des = res_des[config['cols_des']]
            res_tru = res_tru[config['cols_tru']]
            res_sim = add_col(res_sim,param_size,res_tru[param_size])

        if (param_size not in res_des.dtype.names):
            res_des = add_col(res_des,param_size,np.random.uniform(low=list_psf_edges[0],high=list_psf_edges[-1],size=len(res_des)))
            res_des['e1'] = np.random.uniform(low=-1,high=1,size=len(res_des))
            res_des['e2'] = np.random.uniform(low=-1,high=1,size=len(res_des))

        if (param_snr not in res_des.dtype.names):
            res_des = add_col(res_des,param_snr,np.random.uniform(low=list_snr_edges[0],high=list_snr_edges[-1],size=len(res_des)))

            # pl.figure()
            # for ipsf in range(1,len(list_psf_edges)):
            #     select = (res_sim[param_size] > list_psf_edges[ipsf-1]) & (res_sim[param_size] < list_psf_edges[ipsf])
            #     print np.sum(select)
            #     pl.hist(res_sim[select]['radius'],bins=np.logspace(-1,0.3,100),label=str(ipsf),histtype='step')
            #     pl.legend(framealpha=0.0,frameon=False)
            #     pl.xscale('log')
            # pl.show()

            warnings.warn('added %s to res_sim' % param_size)
    else:
        param_size = 'mean_rgpp_rp' ; param_snr = 'snr'

    list_snr_centers = plotstools.get_bins_centers(list_snr_edges)
    list_psf_centers = plotstools.get_bins_centers(list_psf_edges)

    list_bias = []
    
    # entire sample
    logger.info('calculating bias for the entire sample - using cuts from final selection')
    logger.info(selection_string_final_sim)
    logger.info(selection_string_final_des)
    cat_res = res_sim
    cat_tru = res_tru
    exec selection_string_final_sim
    select_sim = select.copy()
    cat_res = res_des
    exec selection_string_final_des
    select_des = select.copy()
    filename_str = 'all.%s' % args.method


    # filename_flatcats='/Users/tomek/data/DES/sva1_gold_wl_flatcats_v4_info.fits'
    # gold=tktools.load(filename_flatcats)
    # select_gold = (gold['sva1_gold_flags']==0)&(gold['sva1_spte_flags']==0)&(gold['sva1_gold_mag_flags']==0)
    # intersect = np.intersect1d(res_des['coadd_objects_id'],gold['coadd_objects_id'][select_gold])
    # select_gold = np.in1d(res_des['coadd_objects_id'], intersect)
    # select_des = select_des*select_gold

    res_sim_final = res_sim[select_sim]
    res_tru_final = res_tru[select_sim]
    res_des_final = res_des[select_des]

    logger.info('final selection SIM %d/%d %f',len(res_sim[select_sim]),len(res_sim), len(res_sim[select_sim])/float(len(res_sim)))
    logger.info('final selection DES %d/%d %f',len(res_des[select_des]),len(res_des), len(res_des[select_des])/float(len(res_des)))
  
    # select_sim = get_subselection(res_sim_final,res_tru_final,res_des)
    # res_sim_final = res_sim_final[select_sim]
    # res_tru_final = res_tru_final[select_sim]

    mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2,mms,std_mms,ccs,std_ccs,mm1s,std_mm1s,cc1s,std_cc1s,mm2s,std_mm2s,cc2s,std_cc2s=nbc_v7_stats.get_mc(res_sim_final,res_tru_final,res_des_final,use_calibration=args.use_calibration,use_weights=args.use_weights,filename_str=filename_str,correct_selection_bias=config['correct_selection_bias'])
    final_mm = mms

    del(res_sim_final)
    del(res_tru_final)
    del(res_des_final)

    # model
    logger.info('calculating bias for the model - using cuts from model selection')
    logger.info(selection_string_model_sim)
    logger.info(selection_string_model_des)
    cat_res = res_sim
    cat_tru = res_tru
    exec selection_string_model_sim
    select_model_sim = select.copy()
    cat_res = res_des
    exec selection_string_model_des
    select_model_des = select.copy()
    filename_str = 'all.%s' % args.method
    logger.info('model selection SIM %d/%d %f',len(res_sim[select_model_sim]),len(res_sim), len(res_sim[select_model_sim])/float(len(res_sim)))
    logger.info('model selection DES %d/%d %f',len(res_des[select_model_des]),len(res_des), len(res_des[select_model_des])/float(len(res_des)))
    res_sim_model = res_sim[select_model_sim]
    res_tru_model = res_tru[select_model_sim]
    res_des_model = res_des[select_model_des]


    warnings.warn('using %s and %s' % (param_snr,param_size) )

    import pdb; pdb.set_trace()
    for ipsf in range(1,len(list_psf_edges)):
        for isnr in range(1,len(list_snr_edges)):

            select = (res_sim_model[param_snr] > list_snr_edges[isnr-1]) & (res_sim_model[param_snr] < list_snr_edges[isnr]) & (res_sim_model[param_size] > list_psf_edges[ipsf-1]) & (res_sim_model[param_size] < list_psf_edges[ipsf])
            n_gals_sim = len(np.nonzero(select)[0])
            if n_gals_sim < 100:
                import pdb; pdb.set_trace()
            res_sim_select = res_sim_model[select]
            res_tru_select = res_tru_model[select]

            select = (res_des_model[param_snr] > list_snr_edges[isnr-1]) & (res_des_model[param_snr] < list_snr_edges[isnr]) & (res_des_model[param_size] > list_psf_edges[ipsf-1]) & (res_des_model[param_size] < list_psf_edges[ipsf])
            n_gals_des = len(np.nonzero(select)[0])
            res_des_select = res_des_model[select]

            vpsf_mid = np.median(res_sim_select[param_size])
            vsnr_mid = np.median(res_sim_select[param_snr])
            vpsf_min = list_psf_edges[ipsf-1]
            vsnr_min = list_snr_edges[isnr-1]
            vpsf_max = list_psf_edges[ipsf]
            vsnr_max = list_snr_edges[isnr]

            # if param_size == 'sf_hlr':
            #     warnings.warn('adjusting bin edges and mid ')
            #     vpsf_mid = np.median(res_sim_select['radius'])
            #     vpsf_max = vpsf_mid + np.std(res_sim_select['radius'],ddof=1)
            #     vpsf_min = vpsf_mid - np.std(res_sim_select['radius'],ddof=1)

            n_unique = len(np.unique(res_tru_select['id_cosmos'].astype(np.float)))
            if 'disc_flux' in res_sim_select.dtype.names: 
                bulge_fraction = np.sum(res_sim_select['disc_flux']==0)/float(len(res_sim_select))
            else:
                bulge_fraction = 0


            if ('ra_as'  in res_sim_select) & ('dec_as' in res_sim_select): 
                x1 = np.mean(res_sim_select['ra_as']);
                x2 = np.mean(res_sim_select['dec_as']);
                std_x1 = np.std(res_sim_select['ra_as'],ddof=1)/np.sqrt(len(res_sim_select));
                std_x2 = np.std(res_sim_select['dec_as'],ddof=1)/np.sqrt(len(res_sim_select));
                desx1 = np.mean(res_des_select['ra_as']);
                desx2 = np.mean(res_des_select['dec_as']);
                std_desx1 = np.std(res_des_select['ra_as'],ddof=1)/np.sqrt(len(res_des_select));
                std_desx2 = np.std(res_des_select['dec_as'],ddof=1)/np.sqrt(len(res_des_select));
            else:
                x1 = np.random.randn()
                x2 = np.random.randn()
                std_x1 = np.random.randn()
                std_x2 = np.random.randn()
                desx1 = np.random.randn()
                desx2 = np.random.randn()
                std_desx1 = np.random.randn()
                std_desx2 = np.random.randn()

            # print res_tru_select['psf_e1'].max(),res_tru_select['psf_e2'].max(),res_tru_select['psf_e1'].min(),res_tru_select['psf_e2'].min()
            logger.info('-------------------------------------')
            logger.info('ipsf=%2d isnr=%2d %s=%2.2f %s=%2.2f snr=[%2.2f %2.2f] rgprp=[%2.2f %2.2f] n_gals_sim=%d n_gals_des=%d n_unique=%d bulge_fraction=%2.2f' % (ipsf,isnr,param_snr,vsnr_mid,param_size,vpsf_mid,vsnr_min,vsnr_max,vpsf_min,vpsf_max,n_gals_sim,n_gals_des,n_unique,bulge_fraction))

            logger.info(selection_string_model_sim)
            filename_str = 'snr=[%2.2f,%2.2f].psf=[%2.2f,%2.2f]' % (vsnr_min,vsnr_max,vpsf_min,vpsf_max)
            mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2,mms,std_mms,ccs,std_ccs,mm1s,std_mm1s,cc1s,std_cc1s,mm2s,std_mm2s,cc2s,std_cc2s=nbc_v7_stats.get_mc(res_sim_select,res_tru_select,None,use_calibration=args.use_calibration,use_weights=args.use_weights,filename_str=filename_str,correct_selection_bias=config['correct_selection_bias'])

            std_e = np.std(res_sim_select['e1'],ddof=1)
            list_bias.append( [ipsf,isnr,n_gals_sim,n_unique,vpsf_mid,vsnr_mid,vpsf_min,vpsf_max,vsnr_min,vsnr_max,std_e,mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,bulge_fraction,x1,std_x1,x2,std_x2,desx1,std_desx1,desx2,std_desx2,mms,std_mms,ccs,std_ccs,mm1s,std_mm1s,cc1s,std_cc1s,mm2s,std_mm2s,cc2s,std_cc2s] )



    arr_bias = tktools.arr2rec(np.array(list_bias),dtype={'names': ["ipsf","isnr","n_gals",'n_unique',"vpsf_mid","vsnr_mid","vpsf_min","vpsf_max","vsnr_min","vsnr_max","std_e","m","std_m","c","std_c","m1","std_m1","m2","std_m2","c1","std_c1","c2","std_c2","pmm","std_pmm","pcc","std_pcc","bulge_fraction","x1","std_x1","x2","std_x2","desx1","std_desx1","desx2","std_desx2",'mms','std_mms','ccs','std_ccs','mm1s','std_mm1s','cc1s','std_cc1s','mm2s','std_mm2s','cc2s','std_cc2s'], 'formats': ['i4']*4 + ['f8']*44 })

    expected_mm = np.sum(arr_bias['mms']*arr_bias['n_gals']) / np.sum(arr_bias['n_gals'])
    logger.info('final_mms = %2.5f sum_mms = %2.5f' % (final_mm,expected_mm))

    if args.use_calibration:
        filename_table_bias = os.path.join(args.output_dir,'bias_table.calibrated.fits')
    else:
        filename_table_bias = os.path.join(args.output_dir,'bias_table.fits')
    import pyfits
    pyfits.writeto(filename_table_bias,arr_bias,clobber=True)
    logger.info('saved %s',filename_table_bias)

    import pdb; pdb.set_trace()


def get_match_weights(res_sim,res_des):

    res_sim = add_col(res_sim, 'eabs', np.abs(res_sim['e1']+1j*res_sim['e2']) )
    res_des = add_col(res_des, 'eabs', np.abs(res_des['e1']+1j*res_des['e2']) )
    res_sim = add_col(res_sim, 'log10_snr', np.log10(res_sim['snr']) )
    res_des = add_col(res_des, 'log10_snr', np.log10(res_des['snr']) )

    X = np.concatenate([res_sim['e1'][:,np.newaxis],res_sim['e2'][:,np.newaxis],res_sim['eabs'][:,np.newaxis],res_sim['snr'][:,np.newaxis],res_sim['mean_rgpp_rp'][:,np.newaxis],res_sim['radius'][:,np.newaxis],res_sim['disc_flux'][:,np.newaxis],res_sim['bulge_flux'][:,np.newaxis]],axis=1)
    n_use = 500000
    select = np.random.permutation(X.shape[0])[:n_use]
    Xs = X[select,:]
    res_sim_select = res_sim[select]


    params = ['e1','e2','eabs','mean_rgpp_rp','log10_snr','radius','bulge_flux','disc_flux']
    ranges = {'e1': [-1,1], 'e2': [-1,1], 'eabs': [0,0.9], 'mean_rgpp_rp' : [1.1,2], 'log10_snr' : [0.5,2.3] , 'radius' : [0,5], 'bulge_flux' : [-10,20], 'disc_flux' : [-10,10] }

    n_bins=40
    list_bins = []
    list_H = []
    list_hd = []
    
    for ip,vp in enumerate(params):

        bins = np.linspace(ranges[vp][0],ranges[vp][1],n_bins)
        digitized = np.digitize(res_sim_select[vp], bins)
        H = np.zeros([len(Xs),n_bins+1])
        for i in range(H.shape[0]): H[i,digitized[i]]=1
        hd, bd = np.histogram(res_des[vp],bins,normed=True)
        list_bins.append(bins)
        list_H.append(H[:,1:-1].T)
        hd /= np.sum(hd)
        list_hd.append(hd)

    H = np.concatenate(list_H,axis=0)
    H /= np.sum(H.sum(axis=1))

    h = np.concatenate(list_hd,axis=0)
    h /= np.sum(h)
    t = h - H.sum(axis=1)# wtf?

    pl.figure()
    pl.plot(t)
    pl.title('t')
    pl.xlabel('histogram diff vector t')

    pl.figure()
    pl.plot(h)
    pl.title('h')
    pl.xlabel('histogram vector h')

    pl.figure()
    pl.plot(H.sum(axis=1))
    pl.title('H.sum(axis=1)')
    pl.xlabel('histogram vector H.sum(ax=1)')

    n_bins_min = n_bins - 1

    logger.debug('calculating weights')
    from sklearn.linear_model import Ridge
    sigma_regularisation = 2e-10
    regression = Ridge(alpha = sigma_regularisation, fit_intercept=False, normalize=False, copy_X=False, max_iter=None, tol=0.00001, solver='auto')

    
    # regression.fit(H,t)
    regression.fit(H[2*n_bins_min:,:],t[2*n_bins_min:])
    weight_match = regression.coef_

    pl.figure()
    pl.hist(weight_match,np.linspace(-2,2,100))
    pl.xlabel('weights value')

    list_hp = []
    for ip,vp in enumerate(params):
        hp,bp = np.histogram(res_sim_select[vp],bins=list_bins[ip],weights=(weight_match+1),normed=True)
        list_hp.append(hp)

    w = np.ones(H.shape[1])
    
    list_hw = []
    for ip,vp in enumerate(params):

        pl.figure()
        hw = np.dot(H[ip*n_bins_min:(ip+1)*n_bins_min,:],w)
        hw /= np.sum(hw)
        list_hw.append(hw)
        list_hd[ip] /= np.sum(list_hd[ip])
        list_hp[ip] /= np.sum(list_hp[ip])
        pl.plot(list_bins[ip][:-1],list_hw[ip],'+m',label='COSMOS'); 
        pl.plot(list_bins[ip][:-1],list_hd[ip],'or',label='data',markerfacecolor='None' )
        pl.plot(list_bins[ip][:-1],list_hp[ip],'xb',label='match-weighted')
        pl.xlabel(vp)
        pl.legend()


    pl.show()

    import pdb; pdb.set_trace()

    return np.ones(len(res_sim))

def get_subselection(res_sim,res_tru,res_des,param='snr',param_min=0,param_max=300,id_max=200):

    pl.figure()
    # hsnr_res, _ , _= pl.hist(res_sim['snr'] ,bins=np.linspace(5,300,300),histtype='step',label='SIM snr'      , normed=True, color='r') 
    # hsnr_des, _ , _= pl.hist(res_des['snr'] ,bins=np.linspace(5,300,300),histtype='step',label='%s snr' % config['methods'][args.method]['label'] , normed=True, color='b') 
    hsnr_des, bins_snr , _= pl.hist(res_des[param] ,bins=np.logspace(param_min,param_max,300),histtype='step',label='%s' % config['methods'][args.method]['label'] , normed=True, color='b') 
    hsnr_res, bins_snr , _= pl.hist(res_sim[param] ,bins=np.logspace(param_min,param_max,300),histtype='step',label='SIM' , normed=True, color='r') 
    ylim=list(pl.ylim()); ylim[1]*=2; pl.ylim(ylim)
    pl.xscale('log')
    pl.legend(framealpha=0.0,frameon=False)
    pl.xlabel('SNR')
    pl.ylim(1e-3,2e-1)
    pl.yscale('log')

    prob = hsnr_des / hsnr_res
    p_max = prob[id_max]
    print bins_snr[id_max]
    print p_max
    prob /= p_max
    prob[np.isfinite(prob)==False] = 0
    prob = np.concatenate([ np.array([0]), prob, np.array([1]) ])
    digitized = np.digitize(res_sim[param],bins_snr)

    prob_col = prob[digitized]
    select = prob_col > np.random.uniform(size=len(prob_col)) 

    pl.figure()
    pl.plot(bins_snr[1:],prob[1:-1])
    pl.axhline(1,color='k')

    pl.figure()
    # hsnr_res, _ , _= pl.hist(res_sim['snr'] ,bins=np.linspace(5,300,300),histtype='step',label='SIM snr'      , normed=True, color='r') 
    # hsnr_des, _ , _= pl.hist(res_des['snr'] ,bins=np.linspace(5,300,300),histtype='step',label='%s snr' % config['methods'][args.method]['label'] , normed=True, color='b') 
    hsnr_des, bins_snr , _= pl.hist(res_des[param] ,bins=np.logspace(param_min,param_max,300),histtype='step',label='%s' % config['methods'][args.method]['label'] , normed=True, color='b') 
    hsnr_res, bins_snr , _= pl.hist(res_sim[param][select] ,bins=np.logspace(param_min,param_max,300),histtype='step',label='SIM' , normed=True, color='r') 
    ylim=list(pl.ylim()); ylim[1]*=2; pl.ylim(ylim)
    pl.xscale('log')
    pl.legend(framealpha=0.0,frameon=False)
    pl.xlabel('SNR')
    pl.ylim(1e-3,2e-1)
    # pl.yscale('log')


    logger.info('subselect to match SNR distr, selected %d/%d %f', np.sum(select),len(select), np.sum(select)/float(len(select)))

    return select


def plot_distributions():

    if args.method=='im3shape':

        col_snr = 'snr'
        col_size = 'mean_rgpp_rp'

    elif args.method=='ngmix':
        col_snr = 's2n_r'
        col_size = 'T_r'

    else:
        raise Exception('unknown method %s',args.method)


    pl.style.use('supermongo')

    lagel_great_des = 'GREAT-DES simulation'

    res_sim, res_tru, res_des = load_selection()


    great_des_e1 = res_sim['e1'] - res_tru['g1_true']
    great_des_e2 = res_sim['e2'] - res_tru['g2_true']

    fig_format = args.fig_format


    logger.info('calculating bias for the entire sample - using cuts from final selection')
    logger.info(selection_string_final_sim)
    logger.info(selection_string_final_des)
    
    cat_res = res_sim
    cat_tru = res_tru
    exec selection_string_final_sim
    select_sim = select.copy()
    cat_res = res_des
    exec selection_string_final_des
    select_des = select.copy()
    filename_str = 'all.%s' % args.method

    try:
        bulge_fraction_des = np.sum(res_des['bulge_flux']!=0)/float(len(res_des))
        bulge_fraction_sim = np.sum(res_sim['bulge_flux']!=0)/float(len(res_sim))
        logger.info('bulge_fraction SIM=%2.4f',bulge_fraction_sim)
        logger.info('bulge_fraction DES=%2.4f',bulge_fraction_des)
    except Exception,errmsg:
        logger.error(errmsg)


    # select_gold = get_gold_cut(res_des)
    # select_des = select_gold & select_des

    res_sim = res_sim[select_sim]
    res_tru = res_tru[select_sim]
    res_des = res_des[select_des]
    logger.info('final selection SIM %d/%d %f',len(res_sim),len(select_sim), len(res_sim)/float(len(select_sim)))
    logger.info('final selection DES %d/%d %f',len(res_des),len(select_des), len(res_des)/float(len(select_des)))

    try:
        bulge_fraction_sim = np.sum(res_sim['bulge_flux']!=0)/float(len(res_sim))
        bulge_fraction_des = np.sum(res_des['bulge_flux']!=0)/float(len(res_des))
        logger.info('bulge_fraction SIM=%2.4f',bulge_fraction_sim)
        logger.info('bulge_fraction DES=%2.4f',bulge_fraction_des)
    except Exception,errmsg:
        logger.error(errmsg)

    try:

        list_snr_centers = [15,20,25,30,40,50,100,200,500,1000]
        list_snr_edges = plotstools.get_bins_edges(list_snr_centers)

        pl.figure()
        list_bulge_fraction_des = []
        list_bulge_fraction_sim = []
        list_n_gals = []

        for ib in range(1,len(list_snr_edges)):

            vb=list_snr_centers[ib-1]
            
            select_sim = (list_snr_edges[ib-1] < res_sim['snr']) *  (res_sim['snr'] < list_snr_edges[ib])
            select_des = (list_snr_edges[ib-1] < res_des['snr']) *  (res_des['snr'] < list_snr_edges[ib])

            res_sim_select = res_sim[select_sim]
            res_des_select = res_des[select_des]

            bulge_fraction_des = np.sum(res_des_select['bulge_flux']!=0)/float(len(res_des_select))
            bulge_fraction_sim = np.sum(res_sim_select['bulge_flux']!=0)/float(len(res_sim_select))
            list_n_gals.append( float(len(res_des_select)) )
            logger.info('bulge_fraction SIM=%2.4f',bulge_fraction_sim)
            logger.info('bulge_fraction DES=%2.4f',bulge_fraction_des)

            list_bulge_fraction_sim.append(bulge_fraction_des)
            list_bulge_fraction_des.append(bulge_fraction_sim)
            
        calib_relative = 0.3
        arr_bulge_fraction_sim = np.array(list_bulge_fraction_sim)
        arr_bulge_fraction_des = np.array(list_bulge_fraction_des)
        arr_n_gals = np.array(list_n_gals)
        pl.figure()
        pl.plot(list_snr_centers,arr_bulge_fraction_sim)
        pl.plot(list_snr_centers,arr_bulge_fraction_des)
        pl.plot(list_snr_centers, calib_relative*(arr_bulge_fraction_des-arr_bulge_fraction_sim)*arr_n_gals/float(len(res_des)) )
        pl.title('0.5*(bulge_fraction_des(snr)-bulge_fraction_sim(snr))*n_gals(snr)/n_total_gals')
        pl.xscale('log')
        pl.show()

        np.sum( calib_relative*(arr_bulge_fraction_des-arr_bulge_fraction_sim)*arr_n_gals )/np.sum(arr_n_gals)


            

    except Exception,errmsg:
        logger.error(errmsg)

    import pdb; pdb.set_trace()



    # res_sim,res_tru,res_des = get_subselection(res_sim,res_tru,res_des,param='radius',param_min=0,param_max=1)
    # mpl.rcParams['xtick.labelsize'] = label_size 
    pl.figure()
    pl.hist(great_des_e1, bins=np.linspace(-1,1,200), histtype='step',normed=True , label=r'%s $e_1$' % lagel_great_des, color='r')
    pl.hist(res_des['e1'], bins=np.linspace(-1,1,200),histtype='step',normed=True , label=r'%s $e_1$' % config['methods'][args.method]['label'] , color='b')
    pl.legend(framealpha=0.0,frameon=False,loc='upper left')
    ylim=list(pl.ylim()); ylim[1]*=1.1; pl.ylim(ylim)
    pl.xlabel(r'$e_1$')
    filename_fig = os.path.join(args.output_dir,'figs/dist-DES-SIM.ell1%s'%fig_format)
    pl.savefig(filename_fig)
    logger.info('wrote %s' % (filename_fig))

    pl.figure()
    pl.hist(great_des_e2, bins=np.linspace(-1,1,200), histtype='step',normed=True , label=r'%s $e_2$' % lagel_great_des, color='r')
    pl.hist(res_des['e2'], bins=np.linspace(-1,1,200),histtype='step',normed=True , label=r'%s $e_2$' % config['methods'][args.method]['label'] , color='b')
    pl.legend(framealpha=0.0,frameon=False,loc='upper left')    
    ylim=list(pl.ylim()); ylim[1]*=1.1; pl.ylim(ylim)
    pl.xlabel(r'$e_2$')
    filename_fig = os.path.join(args.output_dir,'figs/dist-DES-SIM.ell2%s'%fig_format)
    pl.savefig(filename_fig)
    logger.info('wrote %s' % (filename_fig))

    pl.figure()
    pl.hist2d(great_des_e1,great_des_e2,bins=np.linspace(-1,1,100))
    pl.legend(framealpha=0.0,frameon=False, mode='expand',ncol=2)
    pl.xlabel(r'$e_1$')
    pl.ylabel(r'$e_2$')
    filename_fig = os.path.join(args.output_dir,'figs/dist-DES-SIM.ell2d%s'%fig_format)
    pl.savefig(filename_fig)
    logger.info('wrote %s' % (filename_fig))

    pl.figure(figsize=(4,4))
    # hsnr_res, _ , _= pl.hist(res_sim['snr'] ,bins=np.linspace(5,300,300),histtype='step',label='SIM snr'      , normed=True, color='r') 
    # hsnr_des, _ , _= pl.hist(res_des['snr'] ,bins=np.linspace(5,300,300),histtype='step',label='%s snr' % config['methods'][args.method]['label'] , normed=True, color='b') 
    hsnr_des, _ , _= pl.hist(res_des[col_snr] ,bins=np.logspace(0,3,100),histtype='step',label='%s' % config['methods'][args.method]['label'] , normed=True, color='b') 
    hsnr_res, _ , _= pl.hist(res_sim[col_snr] ,bins=np.logspace(0,3,100),histtype='step',label='%s' % lagel_great_des , normed=True, color='r') 
    # ylim=list(pl.ylim()); ylim[1]*=2; pl.ylim(ylim)
    pl.xlim([5,1200])
    pl.xscale('log')
    pl.legend(framealpha=0.0,frameon=False)
    pl.xlabel('SNR')
    filename_fig = os.path.join(args.output_dir,'figs/dist-DES-SIM.snr-logscale%s'%fig_format)
    pl.savefig(filename_fig)
    logger.info('wrote %s' % (filename_fig))
    
    pl.figure()
    hsnr_res, _ , _= pl.hist(res_sim[col_snr] ,bins=np.linspace(0,100,300),histtype='step',label='%s S/N' % lagel_great_des , normed=True, color='r') 
    hsnr_des, _ , _= pl.hist(res_des[col_snr] ,bins=np.linspace(0,100,300),histtype='step',label='%s S/N' % config['methods'][args.method]['label'] , normed=True, color='b') 
    # ylim=list(pl.ylim()); ylim[1]*=2; pl.ylim(ylim)
    pl.xlim([5,100])
    pl.legend(framealpha=0.0,frameon=False)
    pl.xlabel('S/N')
    filename_fig = os.path.join(args.output_dir,'figs/dist-DES-SIM.snr%s'%fig_format)
    pl.savefig(filename_fig)
    logger.info('wrote %s' % (filename_fig))

    pl.figure()
    # pl.hist( ( res_sim[col_size]) ,bins=np.linspace(1,2.5,200),histtype='step',label='%s' % lagel_great_des      , normed=True, color='r') 
    # pl.hist( ( res_des[col_size]) ,bins=np.linspace(1,2.5,200),histtype='step',label='%s' % config['methods'][args.method]['label'] , normed=True, color='b') 
    pl.hist( ( res_sim[col_size]) ,bins=np.linspace(1.2,2,200),histtype='step',label='%s $R_{gpp}/R_{p}$' % lagel_great_des      , normed=True, color='r') 
    pl.hist( ( res_des[col_size]) ,bins=np.linspace(1.2,2,200),histtype='step',label='%s $R_{gpp}/R_{p}$' % config['methods'][args.method]['label'] , normed=True, color='b') 
    pl.legend(framealpha=0.0,frameon=False)
    # ylim=list(pl.ylim()); ylim[1]*=1.2; pl.ylim(ylim)
    pl.xlabel(r'$R_{gpp}/R_{p}$')
    filename_fig = os.path.join(args.output_dir,'figs/dist-DES-SIM.size%s'%fig_format)
    pl.savefig(filename_fig)
    logger.info('wrote %s' % (filename_fig))

    pl.figure()
    pl.hist2d(res_sim['snr'],res_sim['nbc_alpha'],[np.linspace(0,100,100),np.linspace(-0.1,0.5,100)]);
    pl.xlabel('snr')
    pl.ylabel('nbc_alpha')

    pl.figure()
    pl.hist2d(res_sim['snr'],res_sim['nbc_m'],[np.linspace(0,100,100),np.linspace(-0.4,0.1,100)]);
    pl.xlabel('snr')
    pl.ylabel('nbc_m')


    pl.figure()
    xs=np.concatenate([res_sim['snr'][:,None][:50000],res_sim['mean_rgpp_rp'][:,None][:50000]],axis=1)
    logger.info('making kde %s',xs.shape[0])
    points1 = np.logspace(1,2,100)
    points2 = np.linspace(1.1,2,100)
    P1,P2 = np.meshgrid(points1,points2)
    points = np.concatenate([P1.flatten()[:,None],P2.flatten()[:,None]],axis=1)
    P1,P2 = np.meshgrid(points1[1:],points2[1:])
    points_centers = np.concatenate([P1.flatten()[:,None],P2.flatten()[:,None]],axis=1)
    ha,_,_ = np.histogram2d(res_sim['snr'], res_sim['mean_rgpp_rp'], bins = [points1,points2])
    import scipy.interpolate
    import pdb; pdb.set_trace()
    # z = scipy.interpolate.griddata(points_centers,ha.T.flatten(),xs,fill_value=0,method='linear')
    # pl.scatter(xs[:,0], xs[:,1], c=z, s=20, edgecolor='')
    # pl.xscale('log')
    pl.xlim([15,100])
    pl.pcolormesh(points1[1:],points2[1:],ha)
    pl.xticks([15,20,30,50,100])
    pl.gca().xaxis.set_major_formatter(pl.matplotlib.ticker.ScalarFormatter())
    pl.ylim([1.15,2])
    pl.xlabel(r'SNR')
    pl.ylabel(r'$R_{gpp}/R_{p}$')
    pl.title('SNR - Rgpp/Rp distribution of galaxies in sim \n (color=galaxy number density)')
    pl.colorbar()



    try:
        pl.figure()
        pl.hist(res_sim['radius'] ,bins=np.linspace(0,2,200),histtype='step',label='%s' % lagel_great_des     , normed=True, color='r') 
        pl.hist(res_des['radius'] ,bins=np.linspace(0,2,200),histtype='step',label='%s' % config['methods'][args.method]['label'] , normed=True, color='b') 
        # pl.hist(res_tru['sf_hlr']*0.7, bins=np.linspace(0,6,200),histtype='step',normed=True , label='TRU'      , color='c')
        pl.legend(framealpha=0.0,frameon=False)
        # ylim=list(pl.ylim()); ylim[1]*=1.2; pl.ylim(ylim)
        pl.xlabel('radius [arcmin]')
        filename_fig = os.path.join(args.output_dir,'figs/dist-DES-SIM.radius%s'%fig_format)
        pl.savefig(filename_fig)
        logger.info('wrote %s' % (filename_fig))

        pl.figure()
        pl.hist(res_sim['bulge_flux'], bins=np.linspace(-10,20,200),histtype='step',normed=True , label='%s' % lagel_great_des      , color='r')
        pl.hist(res_des['bulge_flux'], bins=np.linspace(-10,20,200),histtype='step',normed=True , label='%s' % config['methods'][args.method]['label']  , color='b')
        pl.legend(framealpha=0.0,frameon=False)
        # ylim=list(pl.ylim()); ylim[1]*=1.2; pl.ylim(ylim)
        pl.xlabel('bulge flux')
        filename_fig = os.path.join(args.output_dir,'figs/dist-DES-SIM.bulgeflux%s'%fig_format)
        pl.savefig(filename_fig)
        logger.info('wrote %s' % (filename_fig))

        pl.figure()
        pl.hist(res_sim['disc_flux'], bins=np.linspace(-10,10,200),histtype='step',normed=True , label='%s' % lagel_great_des      , color='r')
        pl.hist(res_des['disc_flux'], bins=np.linspace(-10,10,200),histtype='step',normed=True , label='%s' % config['methods'][args.method]['label']  , color='b')
        pl.legend(framealpha=0.0,frameon=False)
        # ylim=list(pl.ylim()); ylim[1]*=1.2; pl.ylim(ylim)
        pl.xlabel('disc flux')
        filename_fig = os.path.join(args.output_dir,'figs/dist-DES-SIM.discflux%s'%fig_format)
        pl.savefig(filename_fig)
        logger.info('wrote %s' % (filename_fig))

        pl.figure()
        pl.hist(res_sim['dec_as'], bins=np.linspace(-0.5,0.5,100),histtype='step',normed=True , label='%s' % lagel_great_des      , color='r')
        pl.hist(res_des['dec_as'], bins=np.linspace(-0.5,0.5,100),histtype='step',normed=True , label='%s' % config['methods'][args.method]['label']  , color='b')
        pl.legend()
        # ylim=list(pl.ylim()); ylim[1]*=1.2; pl.ylim(ylim)
        filename_fig = os.path.join(args.output_dir,'figs/dist-DES-SIM.dec%s'%fig_format)
        pl.savefig(filename_fig)
        logger.info('wrote %s' % (filename_fig))

        pl.figure()    
        pl.hist(res_sim['ra_as'], bins=np.linspace(-0.5,0.5,100),histtype='step',normed=True , label='%s' % lagel_great_des     , color='r')
        pl.hist(res_des['ra_as'], bins=np.linspace(-0.5,0.5,100),histtype='step',normed=True , label='%s' % config['methods'][args.method]['label']  , color='b')
        pl.legend()
        # ylim=list(pl.ylim()); ylim[1]*=1.2; pl.ylim(ylim)
        filename_fig = os.path.join(args.output_dir,'figs/dist-DES-SIM.ra%s'%fig_format)
        pl.savefig(filename_fig)
        logger.info('wrote %s' % (filename_fig))

        pl.figure()    
        pl.hist(res_tru['sf_hlr'], bins=np.linspace(0,5,100),histtype='step',normed=True , label='%s' % lagel_great_des      , color='r')
        pl.xlabel(r'$sf_hlr$')
        pl.legend()
        # ylim=list(pl.ylim()); ylim[1]*=1.2; pl.ylim(ylim)
        filename_fig = os.path.join(args.output_dir,'figs/dist-SIM.sf_hlr%s'%fig_format)
        pl.savefig(filename_fig)
        logger.info('wrote %s' % (filename_fig))

        pl.figure()    
        pl.hist(res_tru['cosmos_mag_auto'], bins=np.linspace(15,25,100),histtype='step',normed=True , label='%s' % lagel_great_des      , color='r')
        pl.xlabel(r'\verb|mag_auto|')
        pl.legend()
        # ylim=list(pl.ylim()); ylim[1]*=1.2; pl.ylim(ylim)
        filename_fig = os.path.join(args.output_dir,'figs/dist-SIM.mag_auto%s'%fig_format)
        pl.savefig(filename_fig)
        logger.info('wrote %s' % (filename_fig))

        pl.figure()    



        pl.suptitle(selection_string_final_sim + '\n' + selection_string_final_des)

        pl.show()
        import pdb; pdb.set_trace()

        list_snr_centers = [10,20,25,30,40,50,60,70,80,90]
        list_snr_edges = plotstools.get_bins_edges(list_snr_centers)

        for ib in range(1,len(list_snr_edges)):

            vb=list_snr_centers[ib-1]
            
            select1 = (list_snr_edges[ib-1] < res_sim['snr']) *  (res_sim['snr']< list_snr_edges[ib])
            select2 = (list_snr_edges[ib-1] < res_des['snr']) *  (res_des['snr']< list_snr_edges[ib])

            pl.figure()
            pl.subplot(1,2,1)
            pl.hist(res_sim['e1'][select1],bins=np.linspace(-1,1,100),histtype='step',normed=True , label='SIM' )
            pl.hist(res_des['e1'][select2],bins=np.linspace(-1,1,100),histtype='step',normed=True , label=config['methods'][args.method]['label'] )
            pl.xlabel('e1 SNR=%2.0f'%vb)   
            ylim=list(pl.ylim()); ylim[1]*=1.5; pl.ylim(ylim)
            pl.legend()

            pl.subplot(1,2,2)
            pl.hist(res_sim['e2'][select1],bins=np.linspace(-1,1,100),histtype='step',normed=True , label='SIM' )
            pl.hist(res_des['e2'][select2],bins=np.linspace(-1,1,100),histtype='step',normed=True , label=config['methods'][args.method]['label'] )
            pl.xlabel('e2 SNR=%2.0f'%vb)   
            ylim=list(pl.ylim()); ylim[1]*=1.5; pl.ylim(ylim)
            pl.legend()

            filename_fig = os.path.join(args.output_dir,'figs/histograms.e_vs_snr.%02d%s'% (ib,fig_format))
            pl.savefig(filename_fig)
            logger.info('saved %s', filename_fig)

        # bulge fraction

        print 'aaa'
        import pdb; pdb.set_trace()

    except Exception,errmsg:
        logger.error(errmsg)




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

    # list_pos_centers = np.linspace(-0.03,0.03,10)
    # list_pos_edges = plotstools.get_bins_edges(list_snr_centers)
    # list_mean_e1 = []
    # list_mean_e2 = []
    # pos_type = 'ra_as'

    # for ib in range(1,len(list_pos_edges)):

    #     print ib, len(list_pos_centers)
        
    #     select1 = (list_pos_edges[ib-1] < res_sim[pos_type]) *  (res_sim[pos_type]< list_pos_edges[ib])
    #     select2 = (list_pos_edges[ib-1] < res_des[pos_type]) *  (res_des[pos_type]< list_pos_edges[ib])

    #     list_mean_e1.append(np.mean(res_sim['e1'][select1]))
    #     list_mean_e2.append(np.mean(res_des['e2'][select2]))

    # pl.show()

    # list_mean_e1 = np.array(list_mean_e1)
    # list_mean_e2 = np.array(list_mean_e2)

    # pl.plot(list_pos_centers,list_mean_e1)
    # pl.plot(list_pos_centers,list_mean_e2)

    pl.show()

    import pdb; pdb.set_trace()


def get_histograms():

    res_sim,res_tru = nbc_v7_select.get_selection_sim(selection_string_sim,cols_res=['e1','e2','snr','disc_A','bulge_A'],cols_tru=['snr','psf_e1','psf_e2'])
    res_des = nbc_v7_select.get_selection_des(selection_string_des,cols=['e1','e2','snr','mean_psf_e1_sky','mean_psf_e2_sky','disc_A','bulge_A'],n_files=args.num)

    for ipsf in range(1,len(list_psf_edges)):
        for isnr in range(1,len(list_snr_edges)):

            select_sim = (res_sim['snr'] > list_snr_edges[isnr-1]) & (res_sim['snr'] < list_snr_edges[isnr]) & (res_sim['mean_rgpp_rp'] > list_psf_edges[ipsf-1]) & (res_sim['mean_rgpp_rp'] < list_psf_edges[ipsf])
            select_des = (res_des['snr'] > list_snr_edges[isnr-1]) & (res_des['snr'] < list_snr_edges[isnr]) & (res_des['mean_rgpp_rp'] > list_psf_edges[ipsf-1]) & (res_des['mean_rgpp_rp'] < list_psf_edges[ipsf])
            print ipsf,isnr,len(np.nonzero(select)[0])
            res_sim_select = res_sim[select_sim]
            res_tru_select = res_tru[select_tru]
            res_des_select = res_des[select_des]

            pl.figure()
            pl.hist(res_sim_des['e1'],c='r',bins=np.linspace(0,100),histtype='step',normed=True)
            pl.hist(res_sim_des['e2'],c='m',bins=np.linspace(0,100),histtype='step',normed=True)
            pl.hist(res_sim_sim['e1'],c='b',bins=np.linspace(0,100),histtype='step',normed=True)
            pl.hist(res_sim_sim['e2'],c='c',bins=np.linspace(0,100),histtype='step',normed=True)
            pl.title('SNR \in [%2.0f %2.0f], Rgpp/Rp \in [%2.2f %2.2f]')

            filename_fig = 'hist-e.snr%2.0f'



def get_meane_vs_size():

    res_sim, res_tru, res_des = load_selection()    
    

    # bins_snr_edges = [0,5,10,15,20,30,40,60,100,1000]    
    bins_snr_edges = [10,1000,2000]    
    if args.method=='im3shape':
        bins_rgp_edges = np.linspace(1,2,5)
        label_size = 'mean_rgpp_rp'
        use_xlim = [0.9,2.2]
    elif args.method=='ngmix':
        bins_rgp_edges = np.linspace(0,3,10)
        # raise Exception('implement ngmix ')
        label_size = 'pars[:,4]'
        use_xlim = [-0.1,2.2]
    else:
        raise Exception('unknown method %s',args.method)

    # import pdb; pdb.set_trace()

    bins_hlr_edges = np.linspace(0.2,1.5,5)

    list_results1 = []
    list_results2 = []
    list_results3 = []
    list_results4 = []
    list_results5 = []
    list_results6 = []
    list_results7 = []
    list_results8 = []

    list_results_hlr1 = []
    list_results_hlr2 = []
    list_results_hlr3 = []
    list_results_hlr4 = []
    list_results_hlr5 = []
    list_results_hlr6 = []
    list_results_hlr7 = []
    list_results_hlr8 = []


    for isnr in range(1,len(bins_snr_edges)):
        
        list_results1.append( [] );
        list_results2.append( [] );
        list_results3.append( [] );
        list_results4.append( [] );
        list_results5.append( [] );
        list_results6.append( [] );
        list_results7.append( [] );
        list_results8.append( [] );
       
        for irgp in range(1,len(bins_rgp_edges)):

            select1 = (res_sim['snr'] > bins_snr_edges[isnr-1]) & (res_sim['snr'] < bins_snr_edges[isnr])
            select2 = (res_sim['mean_rgpp_rp'] > bins_rgp_edges[irgp-1]) & (res_sim['mean_rgpp_rp'] < bins_rgp_edges[irgp])

            select3 = np.isclose(res_tru['psf_e1'],-0.02)
            select = select1 & select2 & select3          
            # mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True )
            mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = nbc_v7_stats.get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True, res_tru=res_tru[select] )
            list_results1[isnr-1].append( [mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2] )

            select3 = np.isclose(res_tru['psf_e1'], 0.02)
            select = select1 & select2 & select3          
            # mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = nbc_v7_stats.get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True )
            mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = nbc_v7_stats.get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True, res_tru=res_tru[select] )
            list_results2[isnr-1].append( [mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2] )

            select3 = np.isclose(res_tru['psf_e2'],-0.02)
            select = select1 & select2 & select3          
            # mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = nbc_v7_stats.get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True )
            mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = nbc_v7_stats.get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True, res_tru=res_tru[select] )
            list_results3[isnr-1].append( [mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2] )

            select3 = np.isclose(res_tru['psf_e2'], 0.02)
            select = select1 & select2 & select3          
            # mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = nbc_v7_stats.get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True )
            mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = nbc_v7_stats.get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True, res_tru=res_tru[select] )
            list_results4[isnr-1].append( [mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2] )


            select3 = np.isclose(res_tru['psf_e1'],-0.006666666666)
            select = select1 & select2 & select3          
            # mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = nbc_v7_stats.get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True )
            mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = nbc_v7_stats.get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True, res_tru=res_tru[select] )
            list_results5[isnr-1].append( [mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2] )

            select3 = np.isclose(res_tru['psf_e1'], 0.006666666666)
            select = select1 & select2 & select3          
            # mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = nbc_v7_stats.get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True )
            mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = nbc_v7_stats.get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True, res_tru=res_tru[select] )
            list_results6[isnr-1].append( [mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2] )

            select3 = np.isclose(res_tru['psf_e2'],-0.006666666666)
            select = select1 & select2 & select3          
            # mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = nbc_v7_stats.get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True )
            mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = nbc_v7_stats.get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True, res_tru=res_tru[select] )
            list_results7[isnr-1].append( [mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2] )

            select3 = np.isclose(res_tru['psf_e2'],  0.006666666666)
            select = select1 & select2 & select3          
            # mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = nbc_v7_stats.get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True )
            mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = nbc_v7_stats.get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True, res_tru=res_tru[select] )
            list_results8[isnr-1].append( [mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2] )

        # sf_hlr
        list_results_hlr1.append( [] );
        list_results_hlr2.append( [] );
        list_results_hlr3.append( [] );
        list_results_hlr4.append( [] );
        list_results_hlr5.append( [] );
        list_results_hlr6.append( [] );
        list_results_hlr7.append( [] );
        list_results_hlr8.append( [] );

        for ihlr in range(1,len(bins_hlr_edges)):

            select1 = (res_sim['snr'] > bins_snr_edges[isnr-1]) & (res_sim['snr'] < bins_snr_edges[isnr])
            select2 = (res_tru['sf_hlr'] > bins_hlr_edges[ihlr-1]) & (res_tru['sf_hlr'] < bins_hlr_edges[ihlr])

            select3 = np.isclose(res_tru['psf_e1'],-0.02)
            select = select1 & select2 & select3          
            # mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True )
            mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = nbc_v7_stats.get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True, res_tru=res_tru[select] )
            list_results_hlr1[isnr-1].append( [mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2] )

            select3 = np.isclose(res_tru['psf_e1'], 0.02)
            select = select1 & select2 & select3          
            # mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = nbc_v7_stats.get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True )
            mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = nbc_v7_stats.get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True, res_tru=res_tru[select] )
            list_results_hlr2[isnr-1].append( [mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2] )

            select3 = np.isclose(res_tru['psf_e2'],-0.02)
            select = select1 & select2 & select3          
            # mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = nbc_v7_stats.get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True )
            mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = nbc_v7_stats.get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True, res_tru=res_tru[select] )
            list_results_hlr3[isnr-1].append( [mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2] )

            select3 = np.isclose(res_tru['psf_e2'], 0.02)
            select = select1 & select2 & select3          
            # mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = nbc_v7_stats.get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True )
            mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = nbc_v7_stats.get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True, res_tru=res_tru[select] )
            list_results_hlr4[isnr-1].append( [mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2] )


            select3 = np.isclose(res_tru['psf_e1'],-0.006666666666)
            select = select1 & select2 & select3          
            # mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = nbc_v7_stats.get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True )
            mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = nbc_v7_stats.get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True, res_tru=res_tru[select] )
            list_results_hlr5[isnr-1].append( [mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2] )

            select3 = np.isclose(res_tru['psf_e1'], 0.006666666666)
            select = select1 & select2 & select3          
            # mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = nbc_v7_stats.get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True )
            mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = nbc_v7_stats.get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True, res_tru=res_tru[select] )
            list_results_hlr6[isnr-1].append( [mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2] )

            select3 = np.isclose(res_tru['psf_e2'],-0.006666666666)
            select = select1 & select2 & select3          
            # mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = nbc_v7_stats.get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True )
            mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = nbc_v7_stats.get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True, res_tru=res_tru[select] )
            list_results_hlr7[isnr-1].append( [mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2] )

            select3 = np.isclose(res_tru['psf_e2'],  0.006666666666)
            select = select1 & select2 & select3          
            # mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = nbc_v7_stats.get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True )
            mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = nbc_v7_stats.get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True, res_tru=res_tru[select] )
            list_results_hlr8[isnr-1].append( [mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2] )

            print bins_snr_edges[isnr], bins_rgp_edges[irgp]

    means1 = np.array(list_results1)
    means2 = np.array(list_results2)
    means3 = np.array(list_results3)
    means4 = np.array(list_results4)
    means5 = np.array(list_results5)
    means6 = np.array(list_results6)
    means7 = np.array(list_results7)
    means8 = np.array(list_results8)

    means_hlr1 = np.array(list_results_hlr1)
    means_hlr2 = np.array(list_results_hlr2)
    means_hlr3 = np.array(list_results_hlr3)
    means_hlr4 = np.array(list_results_hlr4)
    means_hlr5 = np.array(list_results_hlr5)
    means_hlr6 = np.array(list_results_hlr6)
    means_hlr7 = np.array(list_results_hlr7)
    means_hlr8 = np.array(list_results_hlr8)


    pl.figure()
    pl.errorbar(bins_rgp_edges[1:],means1[0,:,0],yerr=means1[0,:,2],label = 'mean(e1) psf_e1=% 2.2f' % -0.02)
    pl.errorbar(bins_rgp_edges[1:],means2[0,:,0],yerr=means2[0,:,2],label = 'mean(e1) psf_e1=% 2.2f' %  0.02)
    pl.errorbar(bins_rgp_edges[1:],means3[0,:,3],yerr=means3[0,:,5],label = 'mean(e2) psf_e2=% 2.2f' % -0.02)
    pl.errorbar(bins_rgp_edges[1:],means4[0,:,3],yerr=means4[0,:,5],label = 'mean(e2) psf_e2=% 2.2f' %  0.02)
    pl.axhline(0,color='k')
    pl.xlabel(label_size)
    pl.ylabel('e')
    pl.legend(framealpha=0.0,frameon=False, loc='lower right')
    pl.xlim(use_xlim)
    pl.title('%s' % args.method)

    pl.figure()
    pl.errorbar(bins_hlr_edges[1:],means_hlr1[0,:,0],yerr=means_hlr1[0,:,2],label = 'mean(e1) psf_e1=% 2.2f' % -0.02)
    pl.errorbar(bins_hlr_edges[1:],means_hlr2[0,:,0],yerr=means_hlr2[0,:,2],label = 'mean(e1) psf_e1=% 2.2f' %  0.02)
    pl.errorbar(bins_hlr_edges[1:],means_hlr3[0,:,3],yerr=means_hlr3[0,:,5],label = 'mean(e2) psf_e2=% 2.2f' % -0.02)
    pl.errorbar(bins_hlr_edges[1:],means_hlr4[0,:,3],yerr=means_hlr4[0,:,5],label = 'mean(e2) psf_e2=% 2.2f' %  0.02)
    pl.axhline(0,color='k')
    pl.xlabel('COSMOS HLR [arcmin]')
    pl.ylabel('e')
    pl.legend(framealpha=0.0,frameon=False, loc='lower right')
    pl.xlim([0.4,1.6])
    pl.title('%s' % args.method)

    pl.figure()
    pl.hist(res_sim['mean_rgpp_rp'],bins=np.linspace(0,2,100))
    pl.xlabel(label_size)

    pl.figure()
    pl.hist(res_tru['sf_hlr'],bins=np.linspace(0,2,100))
    pl.xlabel('COSMOS hlr')

    pl.show()

            
    import pdb; pdb.set_trace()






    

def get_meane_vs_snr():

    res_sim, res_tru, res_des = load_selection()

    import pdb; pdb.set_trace()

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


    bin_des = arraytools.arr2rec(list_bin_des,dtype=zip(['ibin','nbin','snr_min','snr_max','mean_e1','stdm_e1','mean_e2','stdm_e2','mean_psfe1' , 'mean_psfe2','stdm_psfe1' , 'stdm_psfe2'],['float64']*len(list_bin_des[0])))  
    bin_sim = arraytools.arr2rec(list_bin_sim,dtype=zip(['ibin','nbin','snr_min','snr_max','mean_e1','stdm_e1','mean_e2','stdm_e2','mean_psfe1' , 'mean_psfe2','stdm_psfe1' , 'stdm_psfe2'],['float64']*len(list_bin_sim[0])))  

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


def main():

    valid_actions = ['plot_mc_vs_psf_size','get_mc_vs_psf_size','get_mc_vs_psf_size','plot_mc_vs_true_params','get_mc_vs_true_params','apply_calibration_to_file','get_params_covariance','get_selection_bias','plot_selection_bias','get_radec_split','apply_calibration_selection','save_selection','get_mc_vs_snr','plot_mc_vs_snr','plot_meane_vs_snr','get_meane_vs_snr','get_histograms','get_distributions','plot_distributions','get_PSF_leakage','get_calibration','get_bias_model','apply_calibration_sim','apply_calibration_des','plot_bias_vs_redshift','plot_face_fig','get_bias_vs_redshift','get_jacknife_regions','get_meane_vs_size','test_shape_noise']

    global logger , config , args

    description = 'Get statistics and plot results of noise bias calibration runs'
    parser = argparse.ArgumentParser(description=description, add_help=True)
    parser.add_argument('-v', '--verbosity', type=int, action='store', default=2, choices=(0, 1, 2, 3 ), help='integer verbosity level: min=0, max=3 [default=2]')
    parser.add_argument('-c', '--filename_config', default='sva1.yaml',type=str, action='store', help='name of the yaml config file')
    parser.add_argument('-m', '--method', default='im3shape',type=str, action='store', help='im3shape or ngmix')
    parser.add_argument('-f', '--first', type=int, action='store', default=0, help='index of the first file to analyse')
    parser.add_argument('-n', '--num', type=int, action='store', default=2, help='number of files to analyse, if -1 then =config[n_files]')
    parser.add_argument('-a', '--actions', nargs='+' ,default=None, type=str, action='store',  help='which actions to run, available: %s' % str(valid_actions) )
    parser.add_argument('--use_calibration', default=False, action='store_true', help='if to apply calibration columns')
    parser.add_argument('--use_weights', default=False, action='store_true', help='if to apply weights')
    parser.add_argument('--n_des_files', default=460, type=int, action='store', help='number of DES files to read')
    parser.add_argument('--fig_format', default='.png', type=str, action='store', help='format of the figure files')
    parser.add_argument('-o','--output_dir', default='.', type=str, action='store', help='dir to store results')
    parser.add_argument('--filename_to_calibrate', default='res_sim.fits', type=str, action='store', help='used with apply_calibration_to_file, filename to be calibrated')



    args = parser.parse_args()
    logging_levels = { 0: logging.CRITICAL,1: logging.WARNING,2: logging.INFO,3: logging.DEBUG }
    logging_level = logging_levels[args.verbosity]; logger.setLevel(logging_level)  
    logger.info(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))

    config = yaml.load(open(args.filename_config))
    nbc_v7_select.config=config; nbc_v7_select.args = args; 
    nbc_v7_stats.config=config; nbc_v7_stats.args = args; 

    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)
        os.makedirs(os.path.join(args.output_dir,'figs'))
        logger.info('created dir %s' % (args.output_dir))



    import shutil
    if os.path.abspath(args.output_dir) != os.path.abspath(os.path.dirname(args.filename_config)):
        filename_config_copy = os.path.join(args.output_dir,os.path.basename(args.filename_config))
        shutil.copyfile(args.filename_config,filename_config_copy)
        logger.info('copied config to %s' % (filename_config_copy))

    global selection_string_store_sim;
    global selection_string_store_des;
    global selection_string_model_sim;
    global selection_string_model_des;
    global selection_string_final_sim;
    global selection_string_final_des;
    selection_string_store_sim = config['selection_string_store_sim']
    selection_string_store_des = config['selection_string_store_des']
    selection_string_model_sim = config['selection_string_model_sim']
    selection_string_model_des = config['selection_string_model_des']
    selection_string_final_sim = config['selection_string_final_sim']
    selection_string_final_des = config['selection_string_final_des']

    global list_snr_edges
    global list_psf_edges
    list_snr_edges = config['list_snr_edges']
    list_psf_edges = config['list_psf_edges']

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

if __name__=='__main__':
    main()
