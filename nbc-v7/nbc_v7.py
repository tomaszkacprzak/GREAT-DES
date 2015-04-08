import numpy as np; import pylab as pl; import tktools as tt;
import  sys, logging, yaml, argparse, time, copy, itertools, tktools, warnings, os, fitsio, pyfits;
warnings.simplefilter("once")
sys.path.append('/home/tomek/code/tktools')
logging_level = logging.INFO; logger = logging.getLogger("nbc-v7"); logger.setLevel(logging_level)  
log_formatter = logging.Formatter("%(asctime)s %(name)s %(levelname)s   %(message)s", "%Y-%m-%d %H:%M:%S")
stream_handler = logging.StreamHandler(sys.stdout); stream_handler.setFormatter(log_formatter)
if logger.handlers == [] : logger.addHandler(stream_handler); logger.propagate = False
import nbc_v7_select
import plotstools, fitting

# info_vals = [1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536,131072,262144,524288,1048576]

# for im3shape-v8 catalogs: missing 4 which for v8 is the v7 connection
# snr > 10 32768
# rgpp_rp > 1.15 32

def plot_add_requirements(level=0.01,target=1.0,mult=1.):

    corner = pl.xlim()[0]
    length = abs(pl.xlim()[1]) + abs(pl.xlim()[0])
    pl.gca().add_patch(pl.Rectangle(  (corner, target-mult*level), length , 2*mult*level , facecolor = '0.9' , edgecolor='w' ))


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


    logger.info('applying calibration to SIM')
    res_sim_cal, calib_m, calib_a=get_calibration_columns(res_sim)
    res_sim_cal = get_weight_column(res_sim_cal)

    logger.info('applying calibration to DES')
    res_des_cal , calib_m, calib_a=get_calibration_columns(res_des)
    res_des_cal = get_weight_column(res_des_cal)

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
    pl.show()


    filename = '%s/res_sim.fits' % args.output_dir
    tktools.save(filename, res_sim_cal, clobber=True)
    filename = '%s/res_tru.fits' % args.output_dir
    tktools.save(filename, res_tru, clobber=True)
    filename = '%s/res_des.fits' % args.output_dir
    tktools.save(filename, res_des_cal, clobber=True)










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
    filename_table = 'weights_interpolation_table.cpickle'
    file_pickle = open(filename_table)
    import cPickle as pickle
    sigma_e, sigma_e_max, list_snr_edges, list_rgp_edges , vec_snr_hires, vec_rgp_hires, sigma_e_hires, n_gal, stdd0_e = pickle.load(file_pickle)

    # set the noisy values to one high snr value
    # sigma_e[sigma_e<sigma_e_max] = sigma_e_max
    sigma_greater = sigma_e.copy()
    sigma_greater[stdd0_e>sigma_e]=stdd0_e[stdd0_e>sigma_e]
    sigma_e_max = min(sigma_greater.flatten())

    vec_snr = get_bins_centers(list_snr_edges)
    vec_rgp = get_bins_centers(list_rgp_edges)
    vec_snr_hires = np.linspace(vec_snr.min(),vec_snr.max(),500)
    vec_rgp_hires = np.linspace(vec_rgp.min(),vec_rgp.max(),500)

    import scipy.interpolate
    X1,X2 = np.meshgrid(vec_snr_hires,vec_rgp_hires)

    logger.debug('hires')
    import scipy.interpolate
    func_interp = scipy.interpolate.interp2d(vec_snr,vec_rgp,sigma_greater.T, kind='cubic')
    sigma_e_hires = func_interp(vec_snr_hires,vec_rgp_hires)

    if 'snr'  in cat.dtype.names: 
        size_col = 'mean_rgpp_rp'
        snr_col = 'snr'
    else:
        snr_col = 'im3shape_r_snr'
        size_col = 'im3shape_r_mean_rgpp_rp'

    cat[snr_col][~np.isfinite(cat[snr_col])]=0
    cat[size_col][~np.isfinite(cat[size_col])]=0

    col_sig = scipy.interpolate.griddata((X1.flatten(),X2.flatten()),sigma_e_hires.flatten(),(cat[snr_col],cat[size_col]),method='nearest',rescale=True)

    col_sig[cat[snr_col]>list_snr_edges.max()] = sigma_e_max
    col_sig[cat[size_col]>list_rgp_edges.max()] = sigma_e_max
    col_sig[col_sig<sigma_e_max] = sigma_e_max

    col_w = 1./col_sig**2
    col_w[np.isnan(col_w)] = 0
    col_w[np.isinf(col_w)] = 0
    if np.any(np.isnan(col_w)): logger.warning('nans')
    if np.any(np.isinf(col_w)): logger.warning('infs')

    if 'w' in cat.dtype.names:
        warnings.warn('overwriting existing column with weights')
        cat['w'] = col_w
        cat_added = cat
    else:
        warnings.warn('appending new column with weights')
        cat_added=add_col(cat,'w',col_w,'f4')

    return cat_added



def get_calibration_columns(res_des):

    import cPickle as pickle
    filename_bias_models = os.path.join(args.output_dir,'bias_models.cpickle')
    bm=pickle.load(open(filename_bias_models))
    warnings.warn('using %s' % filename_bias_models)

    X1,X2 = np.meshgrid(bm['model_snr'],bm['model_rgp'])
    x1=X1.flatten('C'); x2=X2.flatten('C'); y1=bm['model_m'].flatten('C'); y2=bm['model_a'].flatten('C')

    if 'snr'  in res_des.dtype.names: 
        q1=np.array(res_des['snr'])
        q2=np.array(res_des['mean_rgpp_rp'])
    else:
        q1=np.array(res_des['im3shape_r_snr'])
        q2=np.array(res_des['im3shape_r_mean_rgpp_rp'])

    griddata_method = 'linear'

    import scipy.interpolate
    m = scipy.interpolate.griddata((x1,x2),y1,(q1,q2),method=griddata_method,fill_value=0)
    a = scipy.interpolate.griddata((x1,x2),y2,(q1,q2),method=griddata_method,fill_value=0)
 
    if 'mean_psf_e1_sky' in res_des.dtype.names:
        c1 = res_des['mean_psf_e1_sky']*a
        c2 = res_des['mean_psf_e2_sky']*a
    else:
        c1 = res_des['im3shape_r_mean_psf_e1_sky']*a
        c2 = res_des['im3shape_r_mean_psf_e2_sky']*a
    # warnings.warn('adding snr-based weights')
    # w = add_weights.get_weights(res_des['snr'])

    import tabletools
    res_des=tabletools.ensureColumn(res_des,'nbc_m' ,arr=m, dtype='f4')
    res_des=tabletools.ensureColumn(res_des,'nbc_c1',arr=c1,dtype='f4')
    res_des=tabletools.ensureColumn(res_des,'nbc_c2',arr=c2,dtype='f4')
    res_des=tabletools.ensureColumn(res_des,'nbc_alpha',arr=a,dtype='f4')
    # res_des=tabletools.ensureColumn(res_des,'w',arr=,dtype='f4')  

    plots = False
    if plots:

        import pylab as pl
        pl.figure()
        pl.scatter(q1,q2,s=20,c=m,lw=0)
        pl.xlim([0,bm['model_snr'].max()*2 ])
        pl.ylim([0,bm['model_rgp'].max()*2 ])
        pl.colorbar()
        pl.xlabel('SNR')
        pl.ylabel('Rgp/Rp')
        pl.title('nbc_m')

        pl.figure()
        pl.pcolormesh(bm['model_snr'],bm['model_rgp'],bm['model_m'])
        pl.axis('tight')

        pl.figure()
        pl.scatter(q1,q2,s=20,c=a,lw=0)
        pl.xlim([0,bm['model_snr'].max()*2 ])
        pl.ylim([0,bm['model_rgp'].max()*2 ])
        pl.colorbar()
        pl.xlabel('SNR')
        pl.ylabel('Rgp/Rp')
        pl.title('nbc_alpha')

        pl.figure()
        pl.pcolormesh(bm['model_snr'],bm['model_rgp'],bm['model_a'])
        pl.xlabel('SNR')
        pl.ylabel('Rgp/Rp')
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
    bias_table = tktools.load(filename_table_bias)
    dx=0.05

    colorscale=plotstools.get_colorscale(len(np.unique(bias_table['ipsf'])),cmap_name='gist_rainbow')

    list_model_m = []
    list_arr_x = []
    list_arr_y = [] 

    n_upsample = 300
    max_snr = 300

    for ipsf in np.unique(bias_table['ipsf']):

        select = bias_table['ipsf'] == ipsf
        bt1=bias_table[select]
        label = r'$R_{gp}/R_p=%2.2f-%2.2f$'%(bt1['vpsf_min'][0],bt1['vpsf_max'][0])
        snr_mid = (bt1['vsnr_min']+bt1['vsnr_max'])/2.
        rgp_mid = (bt1['vpsf_min']+bt1['vpsf_max'])/2.
        # warnings.warn('setting highest SNR bin to no bias')
        # bt1['m'][-1] = 1
        # bt1['std_m'][-1] = 0.001
        # snr_mid[-1] = 2000
        # bt1['m'][-2] = 1
        # bt1['std_m'][-2] = 0.0001
        pl.errorbar(snr_mid,bt1['m']-1,yerr=bt1['std_m'],label=label,fmt='.',c=colorscale[ipsf-1])
        pl.xlabel('SNR')
        pl.ylabel('m')

        snr_pred =  np.linspace(snr_mid.min(),max_snr,n_upsample)
 
        w, w_cov =fitting.fit(snr_mid, bt1['m']-1, s=bt1['std_m'],expand=inv_snr_basis)       
        p, s = fitting.predict(snr_pred,w,w_cov,expand=inv_snr_basis)
        pl.plot(snr_pred,p,c=colorscale[ipsf-1])
        list_model_m.append(p)
        list_arr_x.append(snr_pred)
        list_arr_y.append(np.ones_like(snr_pred)*rgp_mid[0])


    # pl.title(selection_string_des,fontsize=10)
    pl.xscale('log')
    # pl.title('noise bias: multiplicative - im3shape - %s' % filename_table_bias)
    pl.axhline(0,c='k')
    pl.xticks(list_snr_edges[::2])
    # pl.grid()
    pl.legend(framealpha=0.0,frameon=False,loc='lower right',fontsize=14)
    pl.xlim([9,200])
    pl.ylim([-0.5,0.15])
    pl.gca().xaxis.set_major_formatter(pl.matplotlib.ticker.ScalarFormatter())
    # pl.tick_params(axis='both', which='major', labelsize=10)

    model_arr_m = np.array(list_model_m)
    model_arr_x = np.array(list_arr_x)
    model_arr_y = np.array(list_arr_y)
    vec_x = model_arr_x[0,:]
    vec_y = model_arr_y[:,0]
    vec_x_hires = np.linspace(vec_x.min(),vec_x.max(),n_upsample)
    vec_y_hires = np.linspace(vec_y.min(),vec_y.max(),n_upsample)
    import scipy.interpolate    
    func_interp=scipy.interpolate.interp2d(vec_x,vec_y,model_arr_m,kind='linear')
    model_arr_m_hires=func_interp(vec_x_hires,vec_y_hires)

    pl.figure()
    pl.pcolormesh(model_arr_x,model_arr_y,model_arr_m)
    pl.xticks(vec_x)
    pl.yticks(vec_y)
    pl.axis('tight')
    pl.xlabel('SNR')
    pl.ylabel('Rgp/Rp')
    pl.title('multiplicative shear bias')
    
    # lin scale    
    pl.figure()
    pl.pcolormesh(vec_x_hires,vec_y_hires,model_arr_m_hires)
    pl.colorbar()
    cs=pl.contour(vec_x_hires,vec_y_hires,model_arr_m_hires,[0,-0.05,0.10,-0.13])
    pl.clabel(cs, fontsize=9, inline=1)
    pl.xlabel('SNR')
    pl.ylabel('Rgp/Rp')
    pl.axis('tight')
    pl.title('multiplicative shear bias')

    # now leakage

    list_model_a = []
    list_arr_x = []
    list_arr_y = []
    pl.figure()
    for ipsf in np.unique(bias_table['ipsf']):

        select = bias_table['ipsf'] == ipsf
        bt1=bias_table[select]
        label = r'$R_{gp}/R_p=%2.2f-%2.2f$'%(bt1['vpsf_min'][0],bt1['vpsf_max'][0])
        snr_mid = (bt1['vsnr_min']+bt1['vsnr_max'])/2.
        rgp_mid = (bt1['vpsf_min']+bt1['vpsf_max'])/2.

        # warnings.warn('setting highest SNR bin to no bias')
        # bt1['m'][-1] = 0
        # bt1['std_m'][-1] = 0.001
        # snr_mid[-1] = 2000

        warnings.warn('setting highest SNR bin to no bias')
        pl.errorbar(snr_mid+dx*ipsf,bt1['pmm'],yerr=bt1['std_pmm'],label=label,fmt='.',c=colorscale[ipsf-1])

        w, w_cov =fitting.fit(snr_mid, bt1['pmm'], s=bt1['std_pmm'],expand=inv_snr_basis2)
        
        snr_pred =  np.linspace(snr_mid.min(),max_snr,n_upsample)

        p, s = fitting.predict(snr_pred,w,w_cov,expand=inv_snr_basis2)
        pl.plot(snr_pred,p,c=colorscale[ipsf-1])
        list_model_a.append(p)
        list_arr_x.append(snr_pred)
        list_arr_y.append(np.ones_like(snr_pred)*rgp_mid[0])


    # pl.title('noise bias: PSF leakage - im3shape %s' % filename_table_bias)
    # pl.legend(framealpha=0.0,frameon=False,loc='lower right',fontsize=14,mode='expand',ncol=2)
    pl.xscale('log')
    ylim=list(pl.ylim()); ylim[1]*=1.; pl.ylim(ylim)
    pl.axhline(0,c='k')
    # pl.legend(framealpha=0.0,frameon=False,loc='lower right',mode='expand',ncol=2)
    pl.xlabel('SNR')
    pl.ylabel(r'leakage $\alpha$')
    pl.xticks(list_snr_edges[::2])
    pl.xlim([9,200])
    pl.ylim([-1,1])
    # pl.grid()
    pl.gca().xaxis.set_major_formatter(pl.matplotlib.ticker.ScalarFormatter())

    
    # pl.imshow(model_arr_a,aspect='auto',interpolation='nearest'); pl.show()

    model_arr_a = np.array(list_model_a)
    model_arr_x = np.array(list_arr_x)
    model_arr_y = np.array(list_arr_y)
    vec_x = model_arr_x[0,:]
    vec_y = model_arr_y[:,0]
    vec_x_hires = np.linspace(vec_x.min(),vec_x.max(),n_upsample)
    vec_y_hires = np.linspace(vec_y.min(),vec_y.max(),n_upsample)
    import scipy.interpolate    
    func_interp=scipy.interpolate.interp2d(vec_x,vec_y,model_arr_a,kind='linear')
    model_arr_a_hires=func_interp(vec_x_hires,vec_y_hires)

    pl.figure()
    pl.pcolormesh(model_arr_x,model_arr_y,model_arr_a)
    pl.xticks(vec_x)
    pl.yticks(vec_y)
    pl.xlabel('SNR')
    pl.ylabel('Rgp/Rp')
    pl.axis('tight')
    pl.title('PSF leakage')

    # lin scale
    pl.figure()
    pl.pcolormesh(vec_x_hires,vec_y_hires,model_arr_a_hires)
    pl.colorbar()
    cs=pl.contour(vec_x_hires,vec_y_hires,model_arr_a_hires,[0,0.05,0.10,0.15,0.20])
    pl.clabel(cs, fontsize=9, inline=1)
    pl.axis('tight')
    pl.xlabel('SNR')
    pl.ylabel('Rgp/Rp')
    pl.axis('tight')
    pl.title('PSF leakage')


    filename_bias_models = os.path.join(args.output_dir,'bias_models.cpickle')
    import cPickle as pickle
    pickle_dict = {'model_m':model_arr_m_hires, 'model_a':model_arr_a_hires, 'model_snr':vec_x_hires, 'model_rgp':vec_y_hires}
    pickle.dump(pickle_dict,open(filename_bias_models,'w'),protocol=2)
    logger.info('wrote %s' % (filename_bias_models))

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
        pl.xlabel('SNR')
        pl.ylabel('n_unique')


    # pl.title(selection_string_des,fontsize=10)
    pl.xscale('log')
    # pl.title('noise bias: multiplicative - im3shape - %s' % filename_table_bias)
    pl.axhline(0,c='k')
    pl.xticks(list_snr_edges[::2])
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
        pl.xlabel('SNR')
        pl.ylabel('bulge_fraction')


    # pl.title(selection_string_des,fontsize=10)
    pl.xscale('log')
    # pl.title('noise bias: multiplicative - im3shape - %s' % filename_table_bias)
    pl.axhline(0,c='k')
    pl.xticks(list_snr_edges[::2])
    # pl.grid()
    pl.legend(framealpha=0.0,frameon=False,loc='lower right',fontsize=14)
    # pl.xlim([9,200])
    # pl.ylim([-0.5,0.15])
    pl.gca().xaxis.set_major_formatter(pl.matplotlib.ticker.ScalarFormatter())
    # pl.tick_params(axis='both', which='major', labelsize=10)

    pl.show()
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
     
    # entire sample
    filename_str = 'all.nonbc.%s' % args.method
    mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2=get_mc(res_sim_final,res_tru_final,res_des_final,use_calibration=False,use_weights=args.use_weights,filename_str=filename_str)

    filename_str = 'all.withnbc.%s' % args.method
    mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2=get_mc(res_sim_final,res_tru_final,res_des_final,use_calibration=True,use_weights=args.use_weights,filename_str=filename_str)


    z_bins = [ 0.3, 0.644, 0.901, 1.3 ]

    list_bias = []
    list_bias_calibr = []


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

        logger.info(selection_string_final_sim)
        filename_str = 'z%2.2f.%s' % ((vbin[0]+vbin[1])/2.,args.method)
        mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2=get_mc(res_sim_select,res_tru_select,None,use_calibration=False,use_weights=args.use_weights,filename_str=filename_str)

        std_e = np.std(res_sim_select['e1'],ddof=1)
        list_bias.append( [ibin,vbin[0],vbin[1],std_e,mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2] )

        filename_str = 'z%2.2f.%s.corr' % ((vbin[0]+vbin[1])/2.,args.method)
        mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2=get_mc(cal_sim_select,cal_tru_select,None,use_calibration=True,use_weights=args.use_weights,filename_str=filename_str)

        std_e = np.std(cal_sim_select['e1'],ddof=1)
        list_bias_calibr.append( [ibin,vbin[0],vbin[1],std_e,mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2] )

    arr_bias = tktools.arr2rec(np.array(list_bias),dtype={'names': ["iz","z_min","z_max","std_e","m","std_m","c","std_c","m1","std_m1","m2","std_m2","c1","std_c1","c2","std_c2","pmm","std_pmm","pcc","std_pcc","pm1","std_pm1","pm2","std_pm2","mean_e1","mean_e2","stdm_e1","stdm_e2"], 'formats': ['i4']*1 + ['f8']*27 })
    arr_bias_calibr = tktools.arr2rec(np.array(list_bias_calibr),dtype={'names': ["iz","z_min","z_max","std_e","m","std_m","c","std_c","m1","std_m1","m2","std_m2","c1","std_c1","c2","std_c2","pmm","std_pmm","pcc","std_pcc","pm1","std_pm1","pm2","std_pm2","mean_e1","mean_e2","stdm_e1","stdm_e2"], 'formats': ['i4']*1 + ['f8']*27 })

    import cPickle as pickle
    
    if args.use_weights == True:
        filename_pickle = 'figdata.bias_vs_redshift.with-weights.pp2'
    else:
        filename_pickle = 'figdata.bias_vs_redshift.pp2'

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
    #         mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2=get_mc(res_sim_select,res_tru_select,None,use_calibration=False,use_weights=args.use_weights,filename_str=filename_str)

    #         std_e = np.std(res_sim_select['e1'],ddof=1)
    #         list_jack.append( [ibin,vbin[0],vbin[1],std_e,mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2] )
    
    #         filename_str = 'z%2.2f.%s.jack%02d.corr' % ((vbin[0]+vbin[1])/2.,args.method,ijack)
    #         mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2=get_mc(cal_sim_select,cal_tru_select,None,use_calibration=True,use_weights=args.use_weights,filename_str=filename_str)

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
        filename_pickle = 'figdata.bias_vs_redshift.with-weights.pp2'
    else:
        title = 'without weights'
        filename_pickle = 'figdata.bias_vs_redshift.pp2'
    arr_bias,arr_bias_calibr,z_bins =pickle.load(open(filename_pickle))
    logger.info('opened %s' % (filename_pickle))

    pl.figure()
    pl.errorbar( (arr_bias['z_min']+arr_bias['z_max'])/2 , arr_bias['m1'], yerr=arr_bias['std_m1'], fmt='r.', label='$m_1$ no NBC')
    pl.errorbar( (arr_bias['z_min']+arr_bias['z_max'])/2 , arr_bias['m2'], yerr=arr_bias['std_m2'], fmt='m.', label='$m_2$ no NBC')
    pl.errorbar( (arr_bias_calibr['z_min']+arr_bias_calibr['z_max'])/2 , arr_bias_calibr['m1'], yerr=arr_bias_calibr['std_m1'], fmt='b.', label='$m_1$ with NBC')
    pl.errorbar( (arr_bias_calibr['z_min']+arr_bias_calibr['z_max'])/2 , arr_bias_calibr['m2'], yerr=arr_bias_calibr['std_m2'], fmt='c.', label='$m_2$ with NBC')
    pl.xlabel('z')
    pl.ylabel('m')
    pl.xticks(z_bins)
    pl.grid()
    pl.ylim([0.8,1.05])
    pl.axhline(1,c='k')
    pl.legend(ncol=2, mode='expand', loc='lower center',framealpha=0.0,frameon=False)
    plot_add_requirements(level=0.03,target=1)
    # pl.title(title)

    pl.figure()
    pl.errorbar( (arr_bias['z_min']+arr_bias['z_max'])/2 , arr_bias['pm1'], yerr=arr_bias['std_pm1'], fmt='r.', label=r'$\alpha_1$ no NBC')
    pl.errorbar( (arr_bias['z_min']+arr_bias['z_max'])/2 , arr_bias['pm2'], yerr=arr_bias['std_pm2'], fmt='m.', label=r'$\alpha_2$ no NBC')
    pl.errorbar( (arr_bias_calibr['z_min']+arr_bias_calibr['z_max'])/2 , arr_bias_calibr['pm1'], yerr=arr_bias_calibr['std_pm1'], fmt='b.', label=r'$\alpha_1$ with NBC')
    pl.errorbar( (arr_bias_calibr['z_min']+arr_bias_calibr['z_max'])/2 , arr_bias_calibr['pm2'], yerr=arr_bias_calibr['std_pm2'], fmt='c.', label=r'$\alpha_2$ with NBC')
    pl.xlabel('z')
    pl.ylabel(r'$\alpha$')
    pl.xticks(z_bins)
    pl.grid()
    pl.axhline(0,c='k')
    pl.legend(ncol=2, mode='expand', loc='lower center',framealpha=0.0,frameon=False)
    pl.ylim([-0.2,0.3])
    plot_add_requirements(level=0.025,target=0)
    # pl.title(title)


    pl.figure() 
    pl.errorbar( (arr_bias['z_min']+arr_bias['z_max'])/2 , arr_bias['mean_e1'], yerr=arr_bias['stdm_e1'], fmt='r.', label=r'mean_e1 no NBC')
    pl.errorbar( (arr_bias['z_min']+arr_bias['z_max'])/2 , arr_bias['mean_e2'], yerr=arr_bias['stdm_e2'], fmt='m.', label=r'mean_e2 no NBC')
    pl.errorbar( (arr_bias_calibr['z_min']+arr_bias_calibr['z_max'])/2 , arr_bias_calibr['mean_e1'], yerr=arr_bias_calibr['stdm_e1'], fmt='b.', label=r'mean_e1 with NBC')
    pl.errorbar( (arr_bias_calibr['z_min']+arr_bias_calibr['z_max'])/2 , arr_bias_calibr['mean_e2'], yerr=arr_bias_calibr['stdm_e2'], fmt='c.', label=r'mean_e2 with NBC')
    pl.xlabel('z')
    pl.ylabel(r'mean_e')
    pl.xticks(z_bins)
    pl.grid()
    pl.axhline(0,c='k')
    pl.legend(ncol=2, mode='expand', loc='upper center',framealpha=0.0,frameon=False)
    # pl.ylim([-0.1,0.3])
    plot_add_requirements(level=0.0005,target=0)
    pl.title(title)
    pl.show()


    import pdb; pdb.set_trace()


def get_shear_estimator(res,use_calibration=False,use_weights=False,res_tru=None):

    if res_tru==None:
        e1 = res['e1']
        e2 = res['e2']
    else:
        gt = res_tru['g1_true'] + 1j*res_tru['g2_true']
        eo = res['e1'] + 1j*res['e2']
        ei = (gt-eo)/(eo*gt.conjugate()-1)
        e1,e2 = ei.real, ei.imag
        
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

    return mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2



def get_mc(res_sim,res_tru,res_des=None,filename_str='default',use_calibration=False,use_weights=False,resample_errors=False):

    n_jack = 1
    mean_e1 = []; stdv_e1 = []; stdm_e1 = []; mean_e2 = []; stdv_e2 = []; stdm_e2 = []; true_e1 = []; true_e2 = []; mean_mm = [];



    for ig,vg in enumerate(config['shear']):

        select = res_tru['id_shear'] == ig
        sim_shear = res_sim[select]
        tru_shear = res_tru[select]

        for ij in range(n_jack):

            n_res = len(sim_shear)
            select = range(ij,n_res,n_jack)
            res_sim_select = sim_shear[select]
            tru_sim_select = tru_shear[select]
            # nbc_c1_sim_select = nbc_c1_sim[select]
            # nbc_c2_sim_select = nbc_c2_sim[select]
            # nbc_m1_sim_select = nbc_m1_sim[select]
            # nbc_m2_sim_select = nbc_m2_sim[select]
            
            n_select = len(res_sim_select)
            e1_mean, e1_stdv, e1_stdm, e2_mean, e2_stdv, e2_stdm, m1_mean, m2_mean = get_shear_estimator(res_sim_select,use_calibration=use_calibration,use_weights=use_weights)
            true_e1 += [ np.mean(tru_sim_select['g1_true']) ]
            true_e2 += [ np.mean(tru_sim_select['g2_true']) ]

            mean_e1 += [ e1_mean ]
            stdv_e1 += [ e1_stdv ]
            stdm_e1 += [ e1_stdm ]
            mean_e2 += [ e2_mean ]
            stdv_e2 += [ e2_stdv ]
            stdm_e2 += [ e2_stdm ]

            mean_mm += [ m1_mean ]

            # mean_e1 += [ res_sim_select['e1'](np.mean(res_sim_select['e1'] - nbc_c1_sim_select)) / (1+np.mean(nbc_m1_sim_select))  ]
            # mean_e1 += [ np.sum(res_sim_select['w']*(res_sim_select['e1'] - nbc_c1_sim_select)) / np.sum(res_sim_select['w']*(1+nbc_m1_sim_select))  ]
            # stdv_e1 += [ np.std(res_sim_select['e1'],ddof=1)  ]
            # stdm_e1 += [ np.std(res_sim_select['e1'],ddof=1)/np.sqrt(n_select)     ]
            
            # # mean_e2 += [ (np.mean(res_sim_select['e2'] - nbc_c2_sim_select)) / (1+np.mean(nbc_m2_sim_select))  ]
            # mean_e1 += [ np.sum(res_sim_select['w']*(res_sim_select['e2'] - nbc_c2_sim_select)) / np.sum(res_sim_select['w']*(1+nbc_m2_sim_select))  ]
            # stdv_e2 += [ np.std(res_sim_select['e2'],ddof=1)  ]
            # stdm_e2 += [ np.std(res_sim_select['e2'],ddof=1)/np.sqrt(n_select)     ]

            

    mean_e1 = np.array(mean_e1); stdv_e1 = np.array(stdv_e1); stdm_e1 = np.array(stdm_e1); mean_e2 = np.array(mean_e2); stdv_e2 = np.array(stdv_e2); stdm_e2 = np.array(stdm_e2); true_e1 = np.array(true_e1); true_e2 = np.array(true_e2);
    
    import fitting
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

    cc12,mm12,Ccm12=fitting.get_line_fit(true_e2,mean_e1-true_e1,stdm_e1)
    cc21,mm21,Ccm21=fitting.get_line_fit(true_e1,mean_e2-true_e2,stdm_e2)
    std_cc1 = np.sqrt(Ccm1[0,0])
    std_mm1 = np.sqrt(Ccm1[1,1])
    std_cc2 = np.sqrt(Ccm2[0,0])
    std_mm2 = np.sqrt(Ccm2[1,1])
    std_mm12 = np.sqrt(Ccm12[1,1])
    std_cc12 = np.sqrt(Ccm12[0,0])
    std_mm21 = np.sqrt(Ccm21[1,1])
    std_cc21 = np.sqrt(Ccm21[0,0])

    mm = (mm1+mm2)/2.
    std_mm = np.sqrt( (std_mm1**2 + std_mm2**2)/2. )

    cc = (cc1+cc2)/2.
    std_cc = np.sqrt( (std_cc1**2 + std_cc2**2)/2. )

    e1i_mean, e1i_stdv, e1i_stdm, e2i_mean, e2i_stdv, e2i_stdm, m1i_mean, m2i_mean = get_shear_estimator(res_sim,use_calibration=use_calibration,use_weights=use_weights,res_tru=res_tru)

    # now c
    pcc1,pmm1,pcc2,pmm2,std_pcc1,std_pmm1,std_pcc2,std_pmm2,pmean_e1,pstdm_e1,pmean_e2,pstdm_e2,pmean_g1,pmean_g2,bins_psf_centers = get_PSF_leakage(res_sim,res_tru,use_calibration=use_calibration,use_weights=use_weights)
    logger.info( "--- SIM m1 = % 2.5f +/- % 2.5f"    % (mm1,std_mm1)        )               
    logger.info( "--- SIM m2 = % 2.5f +/- % 2.5f"    % (mm2,std_mm2)        )               
    logger.info( "--- SIM c1 = % 2.5f +/- % 2.5f"    % (cc1,std_cc1)        )               
    logger.info( "--- SIM c2 = % 2.5f +/- % 2.5f"    % (cc2,std_cc2)        )               
    logger.info( "--- SIM m  = % 2.5f +/- % 2.5f"    % (mm,std_mm)          )             
    logger.info( "--- SIM c  = % 2.5f +/- % 2.5f"    % (cc,std_cc)          )             
    logger.info( "--- SIM m12 = % 2.5f +/- % 2.5f"   % (mm12,std_mm12)      )                  
    logger.info( "--- SIM m21 = % 2.5f +/- % 2.5f"   % (mm21,std_mm21)      )                  
    logger.info( "--- SIM pm1  = % 2.5f +/- % 2.5f"  % (pmm1,std_pmm1)      )                   
    logger.info( "--- SIM pc1  = % 2.5f +/- % 2.5f"  % (pcc1,std_pcc1)      )                   
    logger.info( "--- SIM pm2  = % 2.5f +/- % 2.5f"  % (pmm2,std_pmm2)      )                   
    logger.info( "--- SIM pc2  = % 2.5f +/- % 2.5f"  % (pcc2,std_pcc2)      )                   
    logger.info( "--- SIM mean_e1 = %2.5f +/- %2.5f" % (e1i_mean, e1i_stdm) )        
    logger.info( "--- SIM mean_e2 = %2.5f +/- %2.5f" % (e2i_mean, e2i_stdm) )        

    if res_des!=None: dcc1,dmm1,dcc2,dmm2,std_dcc1,std_dmm1,std_dcc2,std_dmm2,dmean_e1,dstdm_e1,dmean_e2,dstdm_e2,dmean_g1,dmean_g2,bins_psf_centers_d = get_PSF_leakage(res_des,None,use_calibration=use_calibration,use_weights=use_weights)
    if res_des!=None: logger.info( "--- DES pm1  = % 2.5f +/- % 2.5f" % (dmm1,std_dmm1) )
    if res_des!=None: logger.info( "--- DES pm2  = % 2.5f +/- % 2.5f" % (dmm2,std_dmm2) )
    if res_des!=None: logger.info( "--- DES pc1  = % 2.5f +/- % 2.5f" % (dcc1,std_dcc1) )
    if res_des!=None: logger.info( "--- DES pc2  = % 2.5f +/- % 2.5f" % (dcc2,std_dcc2) )

    pmm = (pmm1+pmm2)/2.
    pcc = (pcc1+pcc2)/2.
    std_pmm = np.sqrt( (std_pmm1**2 + std_pmm2**2)/2. )
    std_pcc = np.sqrt( (std_pcc1**2 + std_pcc2**2)/2. )

    # now calculate mean ellipticity, no shear




    pl.figure(figsize=(15,8))

    pl.subplot(3,3,1)
    pl.errorbar(true_e1,mean_e1-true_e1,yerr=stdm_e1,fmt='r.')
    pl.errorbar(true_e2,mean_e2-true_e2,yerr=stdm_e2,fmt='m.')
    pl.plot(true_e1,true_e1*0,'-k')
    pl.plot(true_e1,(mm1-1)*true_e1+cc1,  'r-',label=r'<e1> vs e1t sl=$% 2.4f \pm %2.4f$ int=$% 2.4f \pm % 2.4f$' % (mm1,std_mm1,cc1,std_cc1))
    pl.plot(true_e2,(mm2-1)*true_e2+cc2,  'm-',label=r'<e2> vs e2t sl=$% 2.4f \pm %2.4f$ int=$% 2.4f \pm % 2.4f$' % (mm2,std_mm2,cc2,std_cc2))
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


    filename_fig = os.path.join(args.output_dir,'figs/bias.%s.png' % (filename_str))
    pl.savefig(filename_fig)
    logger.info('saved %s',filename_fig)
    pl.close()

    #TODO: make resampled error-bars

    return mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,e1i_mean,e2i_mean,e1i_stdm,e2i_stdm

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
    mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2=get_mc(res_sim_final,res_tru_final,None,use_calibration=args.use_calibration,use_weights=args.use_weights,filename_str=filename_str)

    for isnr in range(1,len(list_snr_edges)):

            # select = (res_sim['snr'] > list_snr_edges[isnr-1]) & (res_sim['snr'] < list_snr_edges[isnr]) 
            select = (res_sim_final['snr'] > list_snr_edges[isnr-1]) 
            n_gals_sim = len(np.nonzero(select)[0])
            logger.info('isnr=%2d snr=[%2.2f %2.2f] n_gals_sim=%d' % (isnr,list_snr_edges[isnr-1],list_snr_edges[isnr],n_gals_sim))

            res_sim_select = res_sim_final[select]
            res_tru_select = res_tru_final[select]

            vsnr_mid = list_snr_centers[isnr-1]

            filename_str = 'm_vs_snr.snr%2.2f' % (vsnr_mid)
            # mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2=get_mc(res_sim_select,res_tru_select,None,use_calibration=args.use_calibration,use_weights=args.use_weights,filename_str=filename_str)
            mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2=get_mc(res_des_select,None,None,use_calibration=args.use_calibration,use_weights=args.use_weights,filename_str=filename_str)

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

    if args.use_calibration:
        logger.info('testing calibration columns')
        cols_res += ['nbc_m1', 'nbc_m2', 'nbc_c1', 'nbc_c2']
        cols_des += ['nbc_m1', 'nbc_m2', 'nbc_c1', 'nbc_c2']
    else:
        logger.info('will measure bias, no calibration applied')

    logger.info('loading SIM results with following selection:')
    logger.info(selection_string_model_sim)
    res_sim,res_tru = nbc_v7_select.get_selection_sim(selection_string_model_sim,cols_res=cols_res,cols_tru=cols_tru,get_calibrated=args.use_calibration)
    logger.info('loading DES results with following selection:')
    logger.info(selection_string_model_des)

    res_des = nbc_v7_select.get_selection_des(selection_string_model_des,cols=cols_des,n_files=args.n_des_files,get_calibrated=args.use_calibration)

    bulge_fraction_sim = np.sum(res_sim['bulge_flux']!=0)/float(len(res_sim))
    bulge_fraction_des = np.sum(res_des['bulge_flux']!=0)/float(len(res_des))
    logger.info('bulge_fraction SIM=%2.4f',bulge_fraction_sim)
    logger.info('bulge_fraction DES=%2.4f',bulge_fraction_des)

    tt.save(os.path.join(args.output_dir,'res_sim.fits'),res_sim,clobber=True)
    tt.save(os.path.join(args.output_dir,'res_des.fits'),res_des,clobber=True)
    tt.save(os.path.join(args.output_dir,'res_tru.fits'),res_tru,clobber=True)

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
        mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2=get_mc(res_sim_select,res_tru_select,None,use_calibration=False,use_weights=True,filename_str=filename_str)
        list_mm_mean.append(mm)
        list_mm_stdm.append(std_mm)

        filename_str = 'patch%03d.cal.png' % ip
        mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2=get_mc(res_sim_select,res_tru_select,None,use_calibration=True,use_weights=True,filename_str=filename_str)
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



def get_calibration():


    res_sim, res_tru, res_des = load_selection()
 
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

    filename_flatcats='/Users/tomek/data/DES/sva1_gold_wl_flatcats_v4_info.fits'
    gold=tktools.load(filename_flatcats)
    select_gold = (gold['sva1_gold_flags']==0)&(gold['sva1_spte_flags']==0)&(gold['sva1_gold_mag_flags']==0)
    intersect = np.intersect1d(res_des['coadd_objects_id'],gold['coadd_objects_id'][select_gold])
    select_gold = np.in1d(res_des['coadd_objects_id'], intersect)
    select_des = select_des*select_gold

    res_sim_final = res_sim[select_sim]
    res_tru_final = res_tru[select_sim]
    res_des_final = res_des[select_des]

    logger.info('final selection SIM %d/%d %f',len(res_sim[select_sim]),len(res_sim), len(res_sim[select_sim])/float(len(res_sim)))
    logger.info('final selection DES %d/%d %f',len(res_des[select_des]),len(res_des), len(res_des[select_des])/float(len(res_des)))

    # model
    # logger.info('calculating bias for the model - using cuts from model selection')
    # logger.info(selection_string_model_sim)
    # logger.info(selection_string_model_des)
    # cat_res = res_sim
    # cat_tru = res_tru
    # exec selection_string_model_sim
    # select_sim = select.copy()
    # cat_res = res_des
    # exec selection_string_model_des
    # select_des = select.copy()
    # filename_str = 'all.%s' % args.method
    # logger.info('model selection SIM %d/%d %f',len(res_sim[select_sim]),len(res_sim), len(res_sim[select_sim])/float(len(res_sim)))
    # logger.info('model selection DES %d/%d %f',len(res_des[select_des]),len(res_des), len(res_des[select_des])/float(len(res_des)))
    res_sim_model = res_sim
    res_tru_model = res_tru
    res_des_model = res_des
   
    # select_sim = get_subselection(res_sim_final,res_tru_final,res_des)
    # res_sim_final = res_sim_final[select_sim]
    # res_tru_final = res_tru_final[select_sim]

    mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2=get_mc(res_sim_final,res_tru_final,res_des_final,use_calibration=args.use_calibration,use_weights=args.use_weights,filename_str=filename_str)

    import pdb; pdb.set_trace()

    for ipsf in range(1,len(list_psf_edges)):
        for isnr in range(1,len(list_snr_edges)):

            select = (res_sim_model['snr'] > list_snr_edges[isnr-1]) & (res_sim_model['snr'] < list_snr_edges[isnr]) & (res_sim_model['mean_rgpp_rp'] > list_psf_edges[ipsf-1]) & (res_sim_model['mean_rgpp_rp'] < list_psf_edges[ipsf])
            n_gals_sim = len(np.nonzero(select)[0])
            if n_gals_sim < 100:
                import pdb; pdb.set_trace()
            res_sim_select = res_sim_model[select]
            res_tru_select = res_tru_model[select]

            select = (res_des_model['snr'] > list_snr_edges[isnr-1]) & (res_des_model['snr'] < list_snr_edges[isnr]) & (res_des_model['mean_rgpp_rp'] > list_psf_edges[ipsf-1]) & (res_des_model['mean_rgpp_rp'] < list_psf_edges[ipsf])
            n_gals_des = len(np.nonzero(select)[0])
            res_des_select = res_des_model[select]

            vpsf_mid = list_psf_centers[ipsf-1]
            vsnr_mid = list_snr_centers[isnr-1]
            vpsf_min = list_psf_edges[ipsf-1]
            vsnr_min = list_snr_edges[isnr-1]
            vpsf_max = list_psf_edges[ipsf]
            vsnr_max = list_snr_edges[isnr]

            n_unique = len(np.unique(res_tru_select['id_cosmos'].astype(np.float)))
            bulge_fraction = np.sum(res_sim_select['disc_flux']==0)/float(len(res_sim_select))

            # print res_tru_select['psf_e1'].max(),res_tru_select['psf_e2'].max(),res_tru_select['psf_e1'].min(),res_tru_select['psf_e2'].min()

            logger.info('ipsf=%2d isnr=%2d snr=[%2.2f %2.2f] rgprp=[%2.2f %2.2f] n_gals_sim=%d n_gals_des=%d n_unique=%d bulge_fraction=%2.2f' % (ipsf,isnr,vsnr_min,vsnr_max,vpsf_min,vpsf_max,n_gals_sim,n_gals_des,n_unique,bulge_fraction))

            logger.info(selection_string_model_sim)
            filename_str = 'snr%2.2f.psf%2.2f' % (vsnr_mid,vpsf_mid)
            mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2=get_mc(res_sim_select,res_tru_select,None,use_calibration=args.use_calibration,use_weights=args.use_weights,filename_str=filename_str)

            std_e = np.std(res_sim_select['e1'],ddof=1)
            list_bias.append( [ipsf,isnr,n_gals_sim,n_unique,vpsf_min,vpsf_max,vsnr_min,vsnr_max,std_e,mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,bulge_fraction] )



    arr_bias = tktools.arr2rec(np.array(list_bias),dtype={'names': ["ipsf","isnr","n_gals",'n_unique',"vpsf_min","vpsf_max","vsnr_min","vsnr_max","std_e","m","std_m","c","std_c","m1","std_m1","m2","std_m2","c1","std_c1","c2","std_c2","pmm","std_pmm","pcc","std_pcc","bulge_fraction"], 'formats': ['i4']*4 + ['f8']*22 })

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

    res_sim, res_tru, res_des = load_selection()

    bulge_fraction_des = np.sum(res_des['bulge_flux']!=0)/float(len(res_des))
    bulge_fraction_sim = np.sum(res_sim['bulge_flux']!=0)/float(len(res_sim))
    logger.info('bulge_fraction SIM=%2.4f',bulge_fraction_sim)
    logger.info('bulge_fraction DES=%2.4f',bulge_fraction_des)

    great_des_e1 = res_sim['e1'] #- cat_tru['g1_true']
    great_des_e2 = res_sim['e2'] #- cat_tru['g2_true']

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

    # select_gold = get_gold_cut(res_des)
    # select_des = select_gold & select_des

    res_sim = res_sim[select_sim]
    res_tru = res_tru[select_sim]
    res_des = res_des[select_des]
    logger.info('final selection SIM %d/%d %f',len(res_sim),len(res_sim), len(res_sim)/float(len(res_sim)))
    logger.info('final selection DES %d/%d %f',len(res_des),len(res_des), len(res_des)/float(len(res_des)))

    bulge_fraction_sim = np.sum(res_sim['bulge_flux']!=0)/float(len(res_sim))
    bulge_fraction_des = np.sum(res_des['bulge_flux']!=0)/float(len(res_des))
    logger.info('bulge_fraction SIM=%2.4f',bulge_fraction_sim)
    logger.info('bulge_fraction DES=%2.4f',bulge_fraction_des)


    # res_sim,res_tru,res_des = get_subselection(res_sim,res_tru,res_des,param='radius',param_min=0,param_max=1)
    pl.figure()
    pl.hist(great_des_e1, bins=np.linspace(-1,1,100), histtype='step',normed=True , label=r'SIM $e_1$'      , color='r')
    pl.hist(great_des_e2, bins=np.linspace(-1,1,100), histtype='step',normed=True , label=r'SIM $e_2$'      , color='m')
    pl.hist(res_des['e1'], bins=np.linspace(-1,1,100),histtype='step',normed=True , label=r'%s $e_1$' % config['methods'][args.method]['label'] , color='b')
    pl.hist(res_des['e2'], bins=np.linspace(-1,1,100),histtype='step',normed=True , label=r'%s $e_2$' % config['methods'][args.method]['label'] , color='c')
    pl.legend(framealpha=0.0,frameon=False, mode='expand',ncol=2)
    ylim=list(pl.ylim()); ylim[1]*=1.2; pl.ylim(ylim)
    pl.xlabel(r'$e_i$')
    filename_fig = os.path.join(args.output_dir,'figs/dist-DES-SIM.ell%s'%fig_format)
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

    pl.figure()
    # hsnr_res, _ , _= pl.hist(res_sim['snr'] ,bins=np.linspace(5,300,300),histtype='step',label='SIM snr'      , normed=True, color='r') 
    # hsnr_des, _ , _= pl.hist(res_des['snr'] ,bins=np.linspace(5,300,300),histtype='step',label='%s snr' % config['methods'][args.method]['label'] , normed=True, color='b') 
    hsnr_des, _ , _= pl.hist(res_des['snr'] ,bins=np.logspace(0,3,300),histtype='step',label='%s' % config['methods'][args.method]['label'] , normed=True, color='b') 
    hsnr_res, _ , _= pl.hist(res_sim['snr'] ,bins=np.logspace(0,3,300),histtype='step',label='SIM' , normed=True, color='r') 
    ylim=list(pl.ylim()); ylim[1]*=2; pl.ylim(ylim)
    pl.xlim([5,1200])
    pl.xscale('log')
    pl.legend(framealpha=0.0,frameon=False)
    pl.xlabel('SNR')
    filename_fig = os.path.join(args.output_dir,'figs/dist-DES-SIM.snr-logscale%s'%fig_format)
    pl.savefig(filename_fig)
    logger.info('wrote %s' % (filename_fig))
    
    pl.figure()
    hsnr_des, _ , _= pl.hist(res_des['snr'] ,bins=np.linspace(0,200,300),histtype='step',label='%s' % config['methods'][args.method]['label'] , normed=True, color='b') 
    hsnr_res, _ , _= pl.hist(res_sim['snr'] ,bins=np.linspace(0,200,300),histtype='step',label='SIM' , normed=True, color='r') 
    ylim=list(pl.ylim()); ylim[1]*=2; pl.ylim(ylim)
    pl.xlim([5,200])
    pl.legend(framealpha=0.0,frameon=False)
    pl.xlabel('SNR')
    filename_fig = os.path.join(args.output_dir,'figs/dist-DES-SIM.snr%s'%fig_format)
    pl.savefig(filename_fig)
    logger.info('wrote %s' % (filename_fig))

    pl.figure()
    pl.hist( ( res_sim['mean_rgpp_rp']) ,bins=np.linspace(1,3,200),histtype='step',label='SIM'      , normed=True, color='r') 
    pl.hist( ( res_des['mean_rgpp_rp']) ,bins=np.linspace(1,3,200),histtype='step',label='%s' % config['methods'][args.method]['label'] , normed=True, color='b') 
    pl.legend(framealpha=0.0,frameon=False)
    ylim=list(pl.ylim()); ylim[1]*=1.2; pl.ylim(ylim)
    pl.xlabel(r'$R_{gpp}/R_{p}$')
    filename_fig = os.path.join(args.output_dir,'figs/dist-DES-SIM.size%s'%fig_format)
    pl.savefig(filename_fig)
    logger.info('wrote %s' % (filename_fig))

    pl.figure()
    pl.hist(res_sim['radius'] ,bins=np.linspace(0,6,200),histtype='step',label='SIM'     , normed=True, color='r') 
    pl.hist(res_des['radius'] ,bins=np.linspace(0,6,200),histtype='step',label='%s' % config['methods'][args.method]['label'] , normed=True, color='b') 
    pl.hist(res_tru['sf_hlr'], bins=np.linspace(0,6,200),histtype='step',normed=True , label='TRU'      , color='c')
    pl.legend(framealpha=0.0,frameon=False)
    ylim=list(pl.ylim()); ylim[1]*=1.2; pl.ylim(ylim)
    pl.xlabel('radius [arcmin]')
    filename_fig = os.path.join(args.output_dir,'figs/dist-DES-SIM.radius%s'%fig_format)
    pl.savefig(filename_fig)
    logger.info('wrote %s' % (filename_fig))

    pl.figure()
    pl.hist(res_sim['bulge_flux'], bins=np.linspace(-10,20,200),histtype='step',normed=True , label='SIM'      , color='r')
    pl.hist(res_des['bulge_flux'], bins=np.linspace(-10,20,200),histtype='step',normed=True , label='%s' % config['methods'][args.method]['label']  , color='b')
    pl.legend(framealpha=0.0,frameon=False)
    ylim=list(pl.ylim()); ylim[1]*=1.2; pl.ylim(ylim)
    pl.xlabel('bulge flux')
    filename_fig = os.path.join(args.output_dir,'figs/dist-DES-SIM.bulgeflux%s'%fig_format)
    pl.savefig(filename_fig)
    logger.info('wrote %s' % (filename_fig))


    pl.figure()
    pl.hist(res_sim['disc_flux'], bins=np.linspace(-10,10,200),histtype='step',normed=True , label='SIM'      , color='r')
    pl.hist(res_des['disc_flux'], bins=np.linspace(-10,10,200),histtype='step',normed=True , label='%s' % config['methods'][args.method]['label']  , color='b')
    pl.legend(framealpha=0.0,frameon=False)
    ylim=list(pl.ylim()); ylim[1]*=1.2; pl.ylim(ylim)
    pl.xlabel('disc flux')
    filename_fig = os.path.join(args.output_dir,'figs/dist-DES-SIM.discflux%s'%fig_format)
    pl.savefig(filename_fig)
    logger.info('wrote %s' % (filename_fig))

    pl.figure()
    pl.hist(res_sim['dec_as'], bins=np.linspace(-0.5,0.5,100),histtype='step',normed=True , label='SIM'      , color='r')
    pl.hist(res_des['dec_as'], bins=np.linspace(-0.5,0.5,100),histtype='step',normed=True , label='%s' % config['methods'][args.method]['label']  , color='b')
    pl.legend()
    ylim=list(pl.ylim()); ylim[1]*=1.2; pl.ylim(ylim)
    filename_fig = os.path.join(args.output_dir,'figs/dist-DES-SIM.dec%s'%fig_format)
    pl.savefig(filename_fig)
    logger.info('wrote %s' % (filename_fig))

    pl.figure()    
    pl.hist(res_sim['ra_as'], bins=np.linspace(-0.5,0.5,100),histtype='step',normed=True , label='SIM'      , color='r')
    pl.hist(res_des['ra_as'], bins=np.linspace(-0.5,0.5,100),histtype='step',normed=True , label='%s' % config['methods'][args.method]['label']  , color='b')
    pl.legend()
    ylim=list(pl.ylim()); ylim[1]*=1.2; pl.ylim(ylim)
    filename_fig = os.path.join(args.output_dir,'figs/dist-DES-SIM.ra%s'%fig_format)
    pl.savefig(filename_fig)
    logger.info('wrote %s' % (filename_fig))

    pl.figure()    
    pl.hist(res_tru['sf_hlr'], bins=np.linspace(0,5,100),histtype='step',normed=True , label='SIM'      , color='r')
    pl.xlabel('sf_hlr')
    pl.legend()
    ylim=list(pl.ylim()); ylim[1]*=1.2; pl.ylim(ylim)
    filename_fig = os.path.join(args.output_dir,'figs/dist-SIM.sf_hlr%s'%fig_format)
    pl.savefig(filename_fig)
    logger.info('wrote %s' % (filename_fig))

    pl.figure()    
    pl.hist(res_tru['cosmos_mag_auto'], bins=np.linspace(15,25,100),histtype='step',normed=True , label='SIM'      , color='r')
    pl.xlabel('mag_auto')
    pl.legend()
    ylim=list(pl.ylim()); ylim[1]*=1.2; pl.ylim(ylim)
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
    
    res_sim,res_tru = nbc_v7_select.get_selection_sim(selection_string_sim,cols_res=['e1','e2','snr','mean_rgpp_rp','nbc_m1','nbc_m2','nbc_c1','nbc_c2','w'],cols_tru=['snr','psf_e1','psf_e2','g1_true','g2_true','sf_hlr'])

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
            mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True, res_tru=res_tru[select] )
            list_results1[isnr-1].append( [mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2] )

            select3 = np.isclose(res_tru['psf_e1'], 0.02)
            select = select1 & select2 & select3          
            # mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True )
            mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True, res_tru=res_tru[select] )
            list_results2[isnr-1].append( [mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2] )

            select3 = np.isclose(res_tru['psf_e2'],-0.02)
            select = select1 & select2 & select3          
            # mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True )
            mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True, res_tru=res_tru[select] )
            list_results3[isnr-1].append( [mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2] )

            select3 = np.isclose(res_tru['psf_e2'], 0.02)
            select = select1 & select2 & select3          
            # mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True )
            mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True, res_tru=res_tru[select] )
            list_results4[isnr-1].append( [mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2] )


            select3 = np.isclose(res_tru['psf_e1'],-0.006666666666)
            select = select1 & select2 & select3          
            # mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True )
            mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True, res_tru=res_tru[select] )
            list_results5[isnr-1].append( [mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2] )

            select3 = np.isclose(res_tru['psf_e1'], 0.006666666666)
            select = select1 & select2 & select3          
            # mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True )
            mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True, res_tru=res_tru[select] )
            list_results6[isnr-1].append( [mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2] )

            select3 = np.isclose(res_tru['psf_e2'],-0.006666666666)
            select = select1 & select2 & select3          
            # mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True )
            mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True, res_tru=res_tru[select] )
            list_results7[isnr-1].append( [mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2] )

            select3 = np.isclose(res_tru['psf_e2'],  0.006666666666)
            select = select1 & select2 & select3          
            # mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True )
            mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True, res_tru=res_tru[select] )
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
            mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True, res_tru=res_tru[select] )
            list_results_hlr1[isnr-1].append( [mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2] )

            select3 = np.isclose(res_tru['psf_e1'], 0.02)
            select = select1 & select2 & select3          
            # mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True )
            mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True, res_tru=res_tru[select] )
            list_results_hlr2[isnr-1].append( [mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2] )

            select3 = np.isclose(res_tru['psf_e2'],-0.02)
            select = select1 & select2 & select3          
            # mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True )
            mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True, res_tru=res_tru[select] )
            list_results_hlr3[isnr-1].append( [mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2] )

            select3 = np.isclose(res_tru['psf_e2'], 0.02)
            select = select1 & select2 & select3          
            # mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True )
            mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True, res_tru=res_tru[select] )
            list_results_hlr4[isnr-1].append( [mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2] )


            select3 = np.isclose(res_tru['psf_e1'],-0.006666666666)
            select = select1 & select2 & select3          
            # mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True )
            mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True, res_tru=res_tru[select] )
            list_results_hlr5[isnr-1].append( [mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2] )

            select3 = np.isclose(res_tru['psf_e1'], 0.006666666666)
            select = select1 & select2 & select3          
            # mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True )
            mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True, res_tru=res_tru[select] )
            list_results_hlr6[isnr-1].append( [mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2] )

            select3 = np.isclose(res_tru['psf_e2'],-0.006666666666)
            select = select1 & select2 & select3          
            # mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True )
            mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True, res_tru=res_tru[select] )
            list_results_hlr7[isnr-1].append( [mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2] )

            select3 = np.isclose(res_tru['psf_e2'],  0.006666666666)
            select = select1 & select2 & select3          
            # mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True )
            mean_e1, stdv_e1, stdm_e1, mean_e2, stdv_e2, stdm_e2, mean_m1, mean_m2 = get_shear_estimator( res_sim[select], use_weights=True, use_calibration=True, res_tru=res_tru[select] )
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

    res_sim,res_tru = nbc_v7_select.get_selection_sim(selection_string_final_sim,cols_res=['e1','e2','snr'],cols_tru=['snr','psf_e1','psf_e2'])
    res_des = nbc_v7_select.get_selection_des(selection_string_final_des,cols=['e1','e2','snr','mean_psf_e1_sky','mean_psf_e2_sky'],n_files=100)

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
        else:
            bins_psf_edges = np.linspace(-0.01,0.01,50)
            warnings.warn('using more PSF bins')

        bins_psf_centers = plotstools.get_bins_centers(bins_psf_edges)
    else:
        raise Exception('invalid method name %s' % args.method)

    # pl.hist(res['mean_psf_e1_sky'],1000); pl.plot(bins_psf_edges,np.zeros_like(bins_psf_edges),'rd'); pl.show()

    logger.debug('bins_psf_edges: [ ' +  (' '.join(['% 5.2f']*len(bins_psf_edges)))%tuple(bins_psf_edges) + ' ]')


    list_mean_e1 = []
    list_stdm_e1 = []
    list_mean_g1 = []
    
    list_mean_e2 = []
    list_stdm_e2 = []
    list_mean_g2 = []

    for ic in range(1,len(bins_psf_edges)):

        select = (res['mean_psf_e1_sky'] > (bins_psf_edges[ic-1])) & (res['mean_psf_e1_sky'] < (bins_psf_edges[ic]))
        res_select=res[select]

        logger.debug('getting PSF leakage e1, bin %2d [% 5.2f % 5.2f], n_gals=%6d' % (ic, bins_psf_edges[ic-1], bins_psf_edges[ic], len(res_select)))

        pmean_e1, pstdv_e1, pstdm_e1, pmean_e2, pstdv_e2, pstdm_e2, pmean_m1, pmean_m2 = get_shear_estimator(res_select,use_calibration=use_calibration,use_weights=use_weights)

        if res_tru!=None:
            pmean_g1=np.mean(res_tru[select]['g1_true'])
        else:
            pmean_g1=pmean_e1*0
        list_mean_e1.append(pmean_e1)
        list_stdm_e1.append(pstdm_e1)
        list_mean_g1.append(pmean_g1)

        select = (res['mean_psf_e2_sky'] > (bins_psf_edges[ic-1])) & (res['mean_psf_e2_sky'] < (bins_psf_edges[ic]))
        res_select=res[select]

        logger.debug('getting PSF leakage e2, bin %2d [% 5.2f % 5.2f], n_gals=%6d' % (ic, bins_psf_edges[ic-1], bins_psf_edges[ic], len(res_select)))

        pmean_e1, pstdv_e1, pstdm_e1, pmean_e2, pstdv_e2, pstdm_e2, pmean_m1, pmean_m2 = get_shear_estimator(res_select,use_calibration=use_calibration,use_weights=use_weights)

        if res_tru!=None:
            pmean_g2=np.mean(res_tru[select]['g2_true'])
        else:
            pmean_g2=pmean_e2*0
        list_mean_e2.append(pmean_e2)
        list_stdm_e2.append(pstdm_e2)
        list_mean_g2.append(pmean_g2)     

        if np.abs(pmean_e1)>10:
            warnings.warn('np.abs(pmean_e1)>10')
            import pdb; pdb.set_trace()      

    pmean_e1=np.array(list_mean_e1)
    pstdm_e1=np.array(list_stdm_e1)
    pmean_e2=np.array(list_mean_e2)
    pstdm_e2=np.array(list_stdm_e2)  
    pmean_g1=np.array(list_mean_g1)
    pmean_g2=np.array(list_mean_g2)
    
    pcc1,pmm1,pCcm1=fitting.get_line_fit(bins_psf_centers,pmean_e1-pmean_g1,pstdm_e1)
    pcc2,pmm2,pCcm2=fitting.get_line_fit(bins_psf_centers,pmean_e2-pmean_g2,pstdm_e2)
    std_pcc1 = np.sqrt(pCcm1[0,0])
    std_pmm1 = np.sqrt(pCcm1[1,1])
    std_pcc2 = np.sqrt(pCcm2[0,0])
    std_pmm2 = np.sqrt(pCcm2[1,1])

    return pcc1,pmm1,pcc2,pmm2,std_pcc1,std_pmm1,std_pcc2,std_pmm2,pmean_e1,pstdm_e1,pmean_e2,pstdm_e2,pmean_g1,pmean_g2,bins_psf_centers

def main():

    valid_actions = ['get_radec_split','apply_calibration_selection','save_selection','get_mc_vs_snr','plot_mc_vs_snr','plot_meane_vs_snr','get_meane_vs_snr','get_histograms','get_distributions','plot_distributions','get_PSF_leakage','get_calibration','get_bias_model','apply_calibration_sim','apply_calibration_des','plot_bias_vs_redshift','plot_face_fig','get_bias_vs_redshift','get_jacknife_regions','get_meane_vs_size']

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



    args = parser.parse_args()
    logging_levels = { 0: logging.CRITICAL,1: logging.WARNING,2: logging.INFO,3: logging.DEBUG }
    logging_level = logging_levels[args.verbosity]; logger.setLevel(logging_level)  
    logger.info(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))

    config = yaml.load(open(args.filename_config))
    nbc_v7_select.config=config; nbc_v7_select.args = args; 

    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)
        os.makedirs(os.path.join(args.output_dir,'figs'))
        logger.info('created dir %s' % (args.output_dir))



    import shutil
    if os.path.abspath(args.output_dir) != os.path.abspath(os.path.dirname(args.filename_config)):
        filename_config_copy = os.path.join(args.output_dir,os.path.basename(args.filename_config))
        shutil.copyfile(args.filename_config,filename_config_copy)
        logger.info('copied config to %s' % (filename_config_copy))

    global selection_string_model_sim;
    global selection_string_model_des;
    global selection_string_final_sim;
    global selection_string_final_des;
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


main()