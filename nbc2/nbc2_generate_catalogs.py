import os 
import matplotlib as mpl
if 'DISPLAY' not in os.environ:
    mpl.use('agg')
    print 'using backend ' , mpl.get_backend()
import sys, logging, yaml, argparse, time, meds, pyfits, plotstools, tabletools, nbc2_dtypes, galsim, copy, galsim.des, warnings, subprocess
import pylab as pl
import numpy as np
from astropy.stats import sigma_clip
warnings.simplefilter('once')

logging_level = logging.INFO
log = logging.getLogger("nbc2_generate_catalogs") 
log.setLevel(logging_level)  
log_formatter = logging.Formatter("%(asctime)s %(name)s %(levelname)s   %(message)s", "%Y-%m-%d %H:%M:%S")
stream_handler = logging.StreamHandler(sys.stdout)
stream_handler.setFormatter(log_formatter)
log.addHandler(stream_handler)

bins_fwhm_centers = np.array([0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3])
bins_ell_centers  = np.array([-0.02, 0.0, 0.02])
bins_snr_centers  = np.linspace(1,100,100)


bins_fwhm = plotstools.get_bins_edges( bins_fwhm_centers )
bins_ell  = plotstools.get_bins_edges( bins_ell_centers )
bins_snr  = plotstools.get_bins_edges( bins_snr_centers )

DES_PIXEL_SIZE = 0.27

def plot_meds():
    

    filename_cat= 'DES0555-5957_i_cat.fits'
    filename_meds = 'DES0555-5957-i-meds-011.fits.fz'

    medsobj = meds.MEDS(filename_meds)
    sexcat = pyfits.getdata(filename_cat)

    medcat = medsobj._cat

    print 'n_objects meds   ' , len(medcat)
    print 'n_objects sexcat ' , len(sexcat)

    select = sexcat['FLUXERR_AUTO'] > 1e-3
    
    snr = sexcat['FLUX_AUTO'][select]/sexcat['FLUXERR_AUTO'][select]

    pl.figure()
    pl.hist(snr,bins=np.linspace(0,100))

    select = (snr<50) * (snr>49)

    high_snrs = np.nonzero(select)[0]
    print 'len high_snrs' , len(high_snrs) 
 
    obj_id = high_snrs[2]
    # obj_id = 6

    print 'ncutout',medsobj._cat['ncutout'][obj_id]
    print 'obj_id' , obj_id

    mosaic = medsobj.get_mosaic(obj_id)
    pl.figure()
    pl.imshow(mosaic,interpolation='nearest')    

    pl.figure()
    pl.plot(mosaic.flatten())

    print 'np.std(mosaic.flatten())' , np.std(mosaic.flatten(),ddof=1)
    print 'median_absolute_deviation' , median_absolute_deviation(mosaic.flatten())


 
    clipped_mosaic,_ = sigma_clip(mosaic.flatten(),sig=2)
    print 'np.std(clipped_mosaic)' , np.std(clipped_mosaic,ddof=1)
    import pdb; pdb.set_trace()

    pl.figure()
    pl.plot(clipped_mosaic.flatten())

    pl.figure()
    pl.hist(mosaic.flatten(),bins=100)

    pl.show()



    import pdb; pdb.set_trace()

def get_noise_level():

    
    filename_cat= 'DES0555-5957_i_cat.fits'
    filename_meds = 'DES0555-5957-i-meds-011.fits.fz'

    medsobj = meds.MEDS(filename_meds)
    sexcat = pyfits.getdata(filename_cat)

    list_res = []
    
    for iobj,vobj in enumerate(sexcat):

        n_cutouts = medsobj._cat['ncutout'][iobj]

        if n_cutouts < 1:
            print 'no cutouts'
        else:

            try:
                mosaic = medsobj.get_mosaic(iobj)
                snr = sexcat['FLUX_AUTO'][iobj]/sexcat['FLUXERR_AUTO'][iobj]
                sd = get_std(mosaic)
                list_res.append( [sd, snr, n_cutouts] )
            except:
                print 'error'
                pass

        if iobj % 100 == 0: print iobj

    arr_res = np.array(list_res)

    pl.hist(arr_res[:,0] ,100)
    pl.show()

    print 'median noise_sigma' , np.median(arr_res[:,0])

    import pdb; pdb.set_trace()



def get_std(mosaic):

    mos=mosaic.flatten()

    min_pixel = min(mos)

    select = mos < -min_pixel
    mos_select=mos[select]

    n_sigma = 3
    clipped_mosaic,_ = sigma_clip(mos_select,sig=n_sigma)

    return np.std(clipped_mosaic,ddof=1)

def get_psf_snr_dist():

    files_im3 = np.loadtxt('filelist_im3.txt',dtype='a1024')
   
    # init histograms
    cat=pyfits.getdata(files_im3[0])

    log.info('using %d im3shape tiles' % len(cat))


    hist_fwhm_all = np.ones_like(bins_fwhm_centers)
    hist_e1_all   = np.ones_like(bins_ell_centers)
    hist_e2_all   = np.ones_like(bins_ell_centers)
    hist_snr_all  = np.ones_like(bins_snr_centers)

    # for fi , fv in enumerate(files_im3[:10]):
    for fi , fv in enumerate(files_im3):
        try:
            cat=pyfits.getdata(fv)
        except:
            print 'opening failed' , fv
            continue

        n_gals_total = len(cat)

        list_fwhm = []
        list_e1 = []
        list_e2 = []
        for i in range(1,10):
            field_name = 'psf_fwhm_%d' % i
            list_fwhm.extend(cat[field_name])
            field_name = 'psf_e1_%d' % i
            list_e1.extend(cat[field_name])
            field_name = 'psf_e2_%d' % i
            list_e2.extend(cat[field_name])

        list_snr = cat['snr']

        list_fwhm = np.array(list_fwhm)*DES_PIXEL_SIZE
        list_e1 = np.array(list_e1)
        list_e2 = np.array(list_e2)
        
        select= (np.array(list_fwhm) > 0) * (np.abs(np.array(list_e1) + 1j*np.array(list_e2)) < 1) 
        list_fwhm = list_fwhm[select]
        list_e1 = list_e1[select]
        list_e2 = list_e2[select]
        
        select= (np.array(list_snr)>0)
        list_snr = list_snr[select]

        hist_fwhm,_=pl.histogram(list_fwhm,bins=bins_fwhm)
        hist_e1,_=pl.histogram(list_e1,bins=bins_ell)
        hist_e2,_=pl.histogram(list_e2,bins=bins_ell)
        hist_snr,_=pl.histogram(list_snr,bins=bins_snr)

        hist_fwhm_all += hist_fwhm
        hist_e1_all   += hist_e1
        hist_e2_all   += hist_e2
        hist_snr_all  += hist_snr

        log.info('%3d %50s %8d %8d' , fi, fv, n_gals_total, len(list_fwhm))

        del(cat)


   
    pl.figure()
    pl.clf()

    pl.figure()
    bins_fwhm_centered = plotstools.get_bins_centers(bins_fwhm)
    pl.plot(bins_fwhm_centered,hist_fwhm_all,'+-')
    pl.xlabel('exposure fwhm [arsec]')
    filename_fig = 'psf_fwhm_dist.png'
    pl.savefig(filename_fig)
    print 'saved' , filename_fig

    pl.figure()
    bins_e1_centered = plotstools.get_bins_centers(bins_ell)
    pl.plot(bins_e1_centered,hist_e1_all,'+-')
    pl.xlabel('psf e1')
    filename_fig = 'psf_e1_dist.png'
    pl.savefig(filename_fig)
    print 'saved' , filename_fig

    pl.figure()
    bins_e2_centered = plotstools.get_bins_centers(bins_ell)
    pl.plot(bins_e2_centered,hist_e2_all,'+-')
    pl.xlabel('psf e2')
    filename_fig = 'psf_e2_dist.png'
    pl.savefig(filename_fig)
    print 'saved' , filename_fig

    pl.figure()
    bins_snr_centered = plotstools.get_bins_centers(bins_snr)
    pl.plot(bins_snr_centered,hist_snr_all,'+-')
    pl.xlabel('snr')
    pl.xscale('log')
    filename_fig = 'snr_dist.png'
    pl.savefig(filename_fig)
    print 'saved' , filename_fig

   
    filename_psf = 'psf_fwhm_dist.txt'
    hist_fwhm_all /= sum(hist_fwhm_all)
    data_save = np.array([bins_fwhm_centered,hist_fwhm_all]).T
    np.savetxt(filename_psf,data_save,header='# bin_center_fwhm prob_fwhm')
    print 'saved' , filename_psf

    filename_psf = 'psf_ell_dist.txt'
    hist_e1_all /= sum(hist_e1_all)
    hist_e2_all /= sum(hist_e2_all)
    data_save = np.array([bins_e1_centered,hist_e1_all,hist_e2_all]).T
    np.savetxt(filename_psf,data_save,header='# bin_center_e e1 e2')
    print 'saved' , filename_psf

    filename_snr = 'snr_dist.txt'
    hist_snr_all /= sum(hist_snr_all)
    data_save = np.array([bins_snr_centered,hist_snr_all]).T
    np.savetxt(filename_snr,data_save,header='# bin_center_snr snr')
    print 'saved' , filename_snr

    

def median_absolute_deviation(data):

    import numpy as np
    return np.median(np.absolute(data - np.median(data)))

def get_psf_images():


    filename_lores = 'nbc2.psf.lores.fits' 
    filename_hires = 'nbc2.psf.hires.fits' 
    filename_field = 'nbc2.psf.field.fits' 

    if os.path.isfile(filename_lores): os.remove(filename_lores)
    if os.path.isfile(filename_hires): os.remove(filename_hires)
    if os.path.isfile(filename_field): os.remove(filename_field)

    orig_pixel_scale = config['pixel_scale']
    orig_image_size = config['cutout_size']

    # make a master PSF config copy
    config_psf = copy.deepcopy(config)

    filename_cat = 'psf_key.fits'
    config_psf['input']['catalog']['file_name'] = filename_cat

    # make a dummy galaxy to simplify PSF generation
    config_psf['gal'] = {}
    config_psf['gal']['type'] = 'Exponential'
    config_psf['gal']['half_light_radius'] = 3

    n_psfs = len(bins_fwhm_centers)*len(bins_ell_centers)**2
    
    iall = 0
    for ifwhm,fwhm in enumerate(bins_fwhm_centers):
        for ie1,e1 in enumerate(bins_ell_centers):
            for ie2,e2 in enumerate(bins_ell_centers):
                                                          
                log.debug('getting single PSF at the pixel scale of a galaxy')
                
                config_copy1=copy.deepcopy(config_psf)
                config_copy1['psf']['fwhm'] = fwhm
                config_copy1['psf']['ellip']['g1'] = e1
                config_copy1['psf']['ellip']['g2'] = e2
                img_gal,img_psf,_,_  = galsim.config.BuildImages(config=config_copy1,image_num=0,obj_num=0,make_psf_image=True,nimages=1)   
                img_psf = img_psf[0]
                img_psf = img_psf[galsim.BoundsI(1, orig_image_size, 1, orig_image_size)]
                pyfits.append(filename_lores,img_psf.array)
                # filename_lores = 'nbc2.psf.lores.%03d.fits' % iall
                # img_psf.write(filename_lores)
                                                              
                # now the hires PSF, centered in the middle
                log.debug('getting single PSF at high resolution')                 
                config_copy2=copy.deepcopy(config_psf)
                n_sub = config['upsampling']
                n_pad = config['padding']
                n_pix_hires = (orig_image_size + n_pad) * n_sub
                pixel_scale_hires = float(config_copy2['image']['pixel_scale']) / float(n_sub)
                config_copy2['image']['pixel_scale'] = pixel_scale_hires
                config_copy2['image']['size'] = n_pix_hires
                config_copy2['psf']['fwhm'] = fwhm
                config_copy2['psf']['ellip']['g1'] = e1
                config_copy2['psf']['ellip']['g2'] = e2            
                img_gal,img_psf,_,_ = galsim.config.BuildImages(config=config_copy2,image_num=0,obj_num=0,make_psf_image=True,nimages=1)      
                img_psf = img_psf[0]
                img_psf = img_psf[galsim.BoundsI(1, int(n_pix_hires), 1, int(n_pix_hires))]             
                pyfits.append(filename_hires,img_psf.array)
                # filename_hires = 'nbc2.psf.hires.%03d.fits' % iall
                # img_psf.write(filename_hires)

                # now field
                log.debug('getting low res PSF in a field')
                config_copy3=copy.deepcopy(config_psf) 
                config_copy3['image']['pixel_scale'] = orig_pixel_scale
                config_copy3['image']['stamp_size'] = orig_image_size
                config_copy3['image']['type'] = 'Tiled'
                config_copy3['image']['nx_tiles'] = 10
                config_copy3['image']['ny_tiles'] = 10
                config_copy3['image']['stamp_xsize'] = orig_image_size
                config_copy3['image']['stamp_ysize'] = orig_image_size
                config_copy3['psf']['fwhm'] = fwhm
                config_copy3['psf']['ellip']['g1'] = e1
                config_copy3['psf']['ellip']['g2'] = e2
                if 'size' in config_copy3['image']:    del(config_copy3['image']['size'])
                img_gal,img_psf,_,_ = galsim.config.BuildImage(config=config_copy3,image_num=0,obj_num=0,make_psf_image=True)    
                pyfits.append(filename_field,img_psf.array)
                # filename_field = 'nbc2.psf.field.%03d.fits' % iall
                # img_psf.write(filename_field)

                iall += 1
                log.info('generated %3d / %3d psfs fwhm=%2.2f e1=% 2.2f e2=% 2.2f  ' , iall, n_psfs , fwhm, e1, e2)



    hdus_lores=pyfits.open(filename_lores)
    hdus_hires=pyfits.open(filename_hires)
    hdus_field=pyfits.open(filename_field)

    log.info('finished writing %s with %d hdus' , filename_lores , len(hdus_lores))
    log.info('finished writing %s with %d hdus' , filename_hires , len(hdus_hires))
    log.info('finished writing %s with %d hdus' , filename_field , len(hdus_field))
    



def get_meds(noise=True):

    # get the start and end index of files

    id_start = args.first
    if args.num == -1:
        id_last = config['n_files']
    else:
        id_last = id_start + args.num

    for ip in range(config['n_files']):
    
        for ig,vg in enumerate(config['shear']):

            filename_cat = 'nbc2.truth.%03d.g%02d.fits' % (ip,ig)
            if noise: filename_meds = 'nbc2.meds.%03d.g%02d.fits' % (ip,ig)
            else: filename_meds = 'nbc2.meds.%03d.g%02d.noisefree.fits' % (ip,ig)

            config_copy = copy.deepcopy(config)
            if noise==False: del(config_copy['image']['noise'])
            config_copy['input']['catalog']['file_name'] = filename_cat
            config_copy['output']['file_name'] = filename_meds

            log.info('getting %s, noise=%s' % (filename_meds , noise ))
            galsim.config.Process(config_copy)
            fpack(filename_meds)
    
    log.info('done all meds')

def fpack(filename):

    filename_fz = filename + '.fz'
    if os.path.isfile(filename_fz):
        os.remove(filename_fz)

    cmd=['fpack' , '-t' , '10240,1' , filename]
    subprocess.call(cmd)
    os.remove(filename)
    os.rename(filename_fz,filename)
    log.debug('compressed file %s ...' % filename)


def get_psf_index(ifwhm,ie1,ie2):

    filename_key = 'psf_key.fits'
    psf_table = tabletools.loadTable(filename_key)
    key_all , key_fwhm , key_e1 , key_e2  = psf_table['id_psf'] , psf_table['id_psf_fwhm'] , psf_table['id_psf_e1'] , psf_table['id_psf_e2']

    indices = np.ones_like(ifwhm)
    for ii,vv in enumerate(ifwhm):
        select = (ifwhm[ii]==key_fwhm) * (ie1[ii]==key_e1) * (ie2[ii]==key_e2)
        nz = np.nonzero(select)[0][0]
        indices[ii] = key_all[nz] 

    return indices

def get_psf_key():

    n_psfs = len(bins_fwhm_centers)*len(bins_ell_centers)**2
    psf_key = np.zeros([n_psfs,7])
    iall = 0
    for ifwhm,fwhm in enumerate(bins_fwhm_centers):
        for ie1,e1 in enumerate(bins_ell_centers):
            for ie2,e2 in enumerate(bins_ell_centers):
                                                          
                psf_key[iall,:] = iall,ifwhm,ie1,ie2,fwhm,e1,e2
                iall+=1

    filename_key = 'psf_key.fits'
    tabletools.saveTable(filename_key,psf_key,dtype=nbc2_dtypes.dtype_psfkey)
    log.info('saved %s' , filename_key )

def get_truth_catalogs():

    filename_psf_fwhm = 'psf_fwhm_dist.txt'
    filename_psf_ell = 'psf_ell_dist.txt'
    filename_snr = 'snr_dist.txt'

    bins_psf_fwhm,prob_psf_fwhm=np.loadtxt(filename_psf_fwhm).T
    bins_psf_e,prob_psf_e1,prob_psf_e2=np.loadtxt(filename_psf_ell).T
    # prob_psf_e1 = np.ones_like(bins_psf_e)/float(len(bins_psf_e)) , np.ones_like(bins_psf_e)/float(len(bins_psf_e))
    bins_snr,prob_snr=np.loadtxt(filename_snr).T

    n_shears = len(config['shear'])
    n_gals = int(float(config['n_gals_per_file'])) 
    n_pairs = n_gals/2 

    filename_cosmos_catalog = os.path.join(config['input']['real_catalog']['dir'],config['input']['real_catalog']['file_name'])
    cosmos_catalog = pyfits.getdata(filename_cosmos_catalog)
    n_cosmos_gals = len(cosmos_catalog)
    log.info('opened %s with %d images' , filename_cosmos_catalog, n_cosmos_gals)

    for ip in range(config['n_files']):
    
        for ig,vg in enumerate(config['shear']):

            filename_cat = 'nbc2.truth.%03d.g%02d.fits' % (ip,ig)
            
            ids = np.arange(n_gals)[:,None]

            cosmos_ids = np.random.choice(n_cosmos_gals,size=n_pairs)
            rotation_angle = np.random.uniform(low=0,high=2*np.pi,size=n_pairs)          
            # obj_snr = np.random.choice(a=bins_snr,size=n_pairs,p=prob_snr)
            # shear_ids = np.random.choice(n_shears,size=n_pairs)
            
            # get pairs
            # obj_snr = np.concatenate([obj_snr[:,None] , obj_snr[:,None]] , axis=1) ; obj_snr = obj_snr.flatten()[:,None]
            cosmos_ids = np.concatenate([cosmos_ids[:,None] , cosmos_ids[:,None]] , axis=1) ; cosmos_ids = cosmos_ids.flatten()[:,None]
            rotation_angle = np.concatenate([rotation_angle[:,None] , rotation_angle[:,None] + np.pi/2.] , axis=1) ; rotation_angle = rotation_angle.flatten()[:,None] # rotate one of the pair nby 90 deg
            # shear_ids = np.concatenate([shear_ids[:,None] , shear_ids[:,None]] , axis=1) ; shear_ids = shear_ids.flatten()[:,None]
            # shear_g1 = np.concatenate([shear_g1[:,None] , shear_g1[:,None]] , axis=1) ; shear_g1 = shear_g1.flatten()[:,None]
            # shear_g2 = np.concatenate([shear_g2[:,None] , shear_g2[:,None]] , axis=1) ; shear_g2 = shear_g2.flatten()[:,None]

            shear_ids = np.ones_like(ids) * ig
            shear_g1  = np.ones_like(ids) * vg[0]
            shear_g2  = np.ones_like(ids) * vg[1]

            # create pairs
            
            psf_fwhm_ids = np.random.choice(a=len(bins_psf_fwhm),size=n_gals,p=prob_psf_fwhm)
            psf_e1_ids = np.random.choice(a=len(bins_psf_e),size=n_gals,p=prob_psf_e1)
            psf_e2_ids = np.random.choice(a=len(bins_psf_e),size=n_gals,p=prob_psf_e2)
            psf_fwhm = bins_psf_fwhm[psf_fwhm_ids]
            psf_e1 = bins_psf_e[psf_e1_ids]
            psf_e2 = bins_psf_e[psf_e2_ids]
            psf_ids = get_psf_index(psf_fwhm_ids ,psf_e1_ids,psf_e2_ids)

            psf_fwhm=psf_fwhm[:,None]
            psf_e1=psf_e1[:,None]
            psf_e2=psf_e2[:,None]
            psf_ids=psf_ids[:,None]
            obj_snr = np.ones_like(psf_ids)
            obj_flux = np.ones_like(psf_ids)

            # 'names' : [ 'id' ,  'id_cosmos' , 'id_shear' , 'id_psf' ,  'g1_true' , 'g2_true' ,  'snr' ,  'psf_fwhm' , 'psf_e1' , 'psf_e2' , 'rotation_angle']  ,
            catalog = np.concatenate([ids,cosmos_ids,shear_ids,psf_ids,shear_g1,shear_g2,obj_snr,obj_flux,psf_fwhm,psf_e1,psf_e2,rotation_angle],axis=1)
            catalog = tabletools.array2recarray(catalog,dtype=nbc2_dtypes.dtype_truth)

            tabletools.saveTable(filename_cat,catalog)

def get_snr_in_truth_table():


    noise_std = config['des_pixel_noise_sigma']

    for ip in range(config['n_files']):

        all_snr=[]
   
        for ig,vg in enumerate(config['shear']):

            list_normsq = []

            filename_cat = 'nbc2.truth.%03d.g%02d.fits' % (ip,ig)
            filename_meds = 'nbc2.meds.%03d.g%02d.noisefree.fits' % (ip,ig)

            log.info('using %s and %s' , filename_meds, filename_cat)

            noisless_gals = meds.MEDS(filename_meds)
            n_gals = len(noisless_gals._cat)
            cat = tabletools.loadTable(filename_cat)

            for ig in range(n_gals):

                img = noisless_gals.get_cutout(ig,0)
                normsq= np.sum( img.flatten()**2 )
                snr = np.sqrt(normsq)/noise_std
                flux = np.sum(img.flatten())
                cat[ig]['snr'] = snr
                cat[ig]['flux'] = flux

                                
                if ig % 1000 == 0: print 'getting norm of galaxy ' , ig

            tabletools.saveTable(filename_cat, cat)
            all_snr+=cat['snr'].tolist()

        pl.hist(all_snr,bins=np.linspace(0,100,200),histtype='step')
        filename_fig = 'truth_table_snr.png'
        pl.savefig(filename_fig)
        log.info('saved %s' , filename_fig)
        # pl.show()




def get_snr_from_cosmos():

    # filename_noiseless = 'nbc2.meds.%03d.g%02d.fits' % (0,0)
    filename_noiseless = 'noisless_galaxies.meds.fits'
    noisless_gals = meds.MEDS(filename_noiseless)

    n_gals = len(noisless_gals._cat)

    filename_snr = 'snr_dist.txt'   
    b_des, h_des = np.loadtxt(filename_snr,unpack=True)
    h_des = h_des/sum(h_des)

    list_normsq = []

    for ni in range(n_gals):

        img = noisless_gals.get_cutout(ni,0)
        normsq= np.sum( img.flatten()**2 )
        list_normsq.append(normsq)
        if ni % 1000 == 0: print 'getting norm of galaxy ' , ni

    from scipy.optimize import fmin
    arr_normsq = np.array(list_normsq)

    noise_std = 16.7552909851
    min_snr_to_use=15

    print 'using' , sum(h_des[b_des>min_snr_to_use])/sum(h_des)

    def get_kl(x):
        
        snr = np.sqrt(arr_normsq*(x**2)) / noise_std     
        h_sim, b_sim = pl.histogram(snr,bins=plotstools.get_bins_edges(b_des[b_des>min_snr_to_use],constant_spacing=True),normed=True)
        h_sim = h_sim/sum(h_sim)

        h_des_use = h_des[b_des>min_snr_to_use]
        KL_divergence = -sum(  h_des_use * np.log(h_sim)  ) + sum( h_des_use * np.log(h_des_use) )

        print 'KL_divergence=' , KL_divergence , 'scale=' , x
        return KL_divergence

    def get_kl_sigma(sigma):
        
        snr = np.sqrt(arr_normsq) / sigma     
        h_sim, b_sim = pl.histogram(snr,bins=plotstools.get_bins_edges(b_des[b_des>min_snr_to_use],constant_spacing=True),normed=True)
        h_sim = h_sim/sum(h_sim)
        # pl.plot(b_des,h_sim)
        # pl.show()

        h_des_use = h_des[b_des>min_snr_to_use]
        KL_divergence = -sum(  h_des_use * np.log(h_sim)  ) + sum( h_des_use * np.log(h_des_use) )

        print 'KL_divergence=' , KL_divergence , 'sigma=' , sigma
        return KL_divergence


    # best_scale = fmin(get_kl,100.)
    best_sigma = fmin(get_kl_sigma,0.2)
    # best_sigma=0.18


    # snr = np.sqrt(arr_normsq*(best_scale**2)) / noise_std     
    snr = np.sqrt(arr_normsq) / best_sigma     
    # snr = np.sqrt(arr_normsq) / 0.1
    h_sim, b_sim = pl.histogram(snr,bins=plotstools.get_bins_edges(b_des,constant_spacing=True),normed=True)
    h_sim = h_sim/sum(h_sim)

    best_test = 0.3
    snr = np.sqrt(arr_normsq) / best_test     
    # snr = np.sqrt(arr_normsq) / 0.1
    h_sim2, b_sim = pl.histogram(snr,bins=plotstools.get_bins_edges(b_des,constant_spacing=True),normed=True)
    h_sim2 = h_sim2/sum(h_sim2)


    des_scale_pixels = noise_std / best_sigma
    print 'des images noise sigma=' , noise_std
    print 'use noise_sigma=' , best_sigma
    print 'post multiply images by scale=' , des_scale_pixels

    
    pl.plot(b_des,h_sim,'rx-',label='GREAT-DES')
    pl.plot(b_des,h_sim2,'cx-',label='GREAT-DES -test')
    pl.plot(b_des,h_des,'bo-',label='im3shape 011 v3')
    pl.xlabel('SNR')
    # pl.xscale('log')
    pl.legend()
    filename_fig = 'match_GREATDES_and_IM3011v3.png'
    pl.savefig(filename_fig)
    print 'saved' , filename_fig

    return best_sigma, des_scale_pixels , noise_std

def preview_result():

    import meds
    m=meds.MEDS('nbc2.meds.000.g05.fits.fz')
    for i in range(100): print tru['snr'][i];pl.imshow(m.get_cutout(i,0),interpolation='nearest');pl.show()


def main():


    global log , config , args , cat

    description = 'Get input catalogs for GREAT-DES NBC2 simulation'
    parser = argparse.ArgumentParser(description=description, add_help=True)
    parser.add_argument('-v', '--verbosity', type=int, action='store', default=2, choices=(0, 1, 2, 3 ), help='integer verbosity level: min=0, max=3 [default=2]')
    parser.add_argument('-c', '--filename_config', type=str, action='store', default='nbc2.yaml', help='name of the config file')
    parser.add_argument('-f', '--first', type=int, action='store', default=0, help='index of the first file to create')
    parser.add_argument('-n', '--num', type=int, action='store', default=-1, help='number of files to create, if -1 then =config[n_files]')

    args = parser.parse_args()
    logging_levels = { 0: logging.CRITICAL, 
                       1: logging.WARNING,
                       2: logging.INFO,
                       3: logging.DEBUG }
    logging_level = logging_levels[args.verbosity]
    log.setLevel(logging_level)  
    
    log.info(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))

    config = yaml.load(open(args.filename_config))
    # plot_meds()
    # get_noise_level() # 16.7552914481
    # get_snr_from_cosmos()
    
    get_psf_snr_dist();     
    get_psf_key()
    get_psf_images()
    get_truth_catalogs()
    get_meds(noise=False)
    get_snr_in_truth_table()
    get_meds(noise=True)
    
    log.info(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))

main()