                # now the hires PSF, centers in the middle
                # log.debug('getting single PSF at high resolution')                 
                # config_copy2=copy.deepcopy(config_psf)
                # n_sub = config['upsampling']
                # n_pad = config['padding']
                # n_pix_hires = (orig_image_size + n_pad) * n_sub
                # pixel_scale_hires = float(config_copy2['image']['pixel_scale']) / float(n_sub)
                # config_copy2['image']['pixel_scale'] = pixel_scale_hires
                # config_copy2['image']['size'] = n_pix_hires
                # config_copy2['psf']['fwhm'] = fwhm
                # config_copy2['psf']['ellip']['g1'] = e1
                # config_copy2['psf']['ellip']['g2'] = e2            

                # # no pixel convolut
                # config_copy2['pix'] = {}
                # config_copy2['pix']['type'] = 'Pixel'
                # config_copy2['pix']['xw'] = orig_pixel_scale

                # img_gal,img_psf,_,_ = galsim.config.BuildImages(config=config_copy2,image_num=0,obj_num=0,make_psf_image=True,nimages=1)      
                # img_psf = img_psf[0]
                # img_psf = img_psf[galsim.BoundsI(1, int(n_pix_hires), 1, int(n_pix_hires))]             
                # pyfits.append(filename_hires,img_psf.array)

import os 
import matplotlib as mpl
if 'DISPLAY' not in os.environ:
    mpl.use('agg')
    print 'using backend ' , mpl.get_backend()
import sys, logging, yaml, argparse, time, meds, pyfits, plotstools, tabletools, galsim, copy, galsim.des, warnings, subprocess
import pylab as pl
import numpy as np
from nbc_dtypes import *
warnings.simplefilter('once')

logging_level = logging.INFO
log = logging.getLogger("nbc2_gener") 
log.setLevel(logging_level)  
log_formatter = logging.Formatter("%(asctime)s %(name)s %(levelname)s   %(message)s", "%Y-%m-%d %H:%M:%S")
stream_handler = logging.StreamHandler(sys.stdout)
stream_handler.setFormatter(log_formatter)
log.addHandler(stream_handler)
log.propagate = False

DES_PIXEL_SIZE = 0.27
ACS_PIXEL_SCALE = 0.03
BOX_SIZES=[32,48,64,96,128]
GRID_SNR=np.linspace(1,100,100)

def image_array_to_galsim(array):

    nx,ny = array.shape
    img_gs = galsim.ImageD(ny,nx)
    for ix in range(nx):
        for iy in range(ny):
            img_gs.setValue(iy+1,ix+1,float(array[ix][iy]))

    return img_gs


def get_fwhm(image, fwxm=0.5, upsampling=1):
    """
    @ brief Computes the FWHM of an i3 image. Per default, it computes the
    radial averaged profile for determining the FWHM. Alternatively,
    it computes the FWHMs of the profiles along the x and y axes in +
    and - direction (total of 4) and returns their average as an
    estimator of the FWHM.    
    @param image np.array with pixels
    """
    # compute weighted moments to get centroid
    x0,y0=np.unravel_index(np.argmax(image), image.shape)
    # x0 = numpy.floor(moments.x0)
    # y0 = numpy.floor(moments.y0)

    profile_x = image[int(x0), :]
    profile_y = image[:, int(y0)]

    max_val = image[int(x0), int(y0)]
    cut_val = max_val * fwxm

    fwhms = []
    for i in range(4):
        if i == 0:
            profile = profile_x[int(y0)::]
            dc0 = int(x0) - x0
        if i == 1:
            profile = profile_x[0:int(y0)+1][::-1]
            dc0 = -int(x0) + x0
        if i == 2:
            profile = profile_y[int(x0)::]
            dc0 = int(y0) - y0
        if i == 3:
            profile = profile_y[0:int(x0)+1][::-1]
            dc0 = -int(y0) + y0
     
        diff = abs(profile - cut_val)
     
        # fhwm code from Tomek
        f1 = 0.
        f2 = 0.
        x1 = 0
        x2 = 0
        
        x1 = np.argmin(diff)
        f1 = profile[x1]
     
        if( f1 < cut_val ):  x2 = x1+1
        else:       x2 = x1-1
        f2 = profile[x2];
     
        a = (f1-f2)/(x1 - x2)
        b = f1 - a*x1;
        x3 = (cut_val - b)/a;
     
        fwhms.append(2.* (dc0 + x3))

        fwhm =  np.mean(np.array(fwhms))/upsampling

    return fwhm

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
    from astropy.stats import sigma_clip
    clipped_mosaic,_ = sigma_clip(mos_select,sig=n_sigma)

    return np.std(clipped_mosaic,ddof=1)

def get_params_dist_from_DES():

    files_im3 = np.loadtxt('filelist_im3.txt',dtype='a1024')
   
    # init histograms
    cat=pyfits.getdata(files_im3[0])

    log.info('using %d im3shape tiles' % len(cat))


    hist_fwhm_all = np.ones_like(config['bins_fwhm_centers'])
    hist_e1_all   = np.ones_like(config['bins_ell_centers'])
    hist_e2_all   = np.ones_like(config['bins_ell_centers'])
    hist_snr_all = np.ones_like(GRID_SNR)
    hist_box_all = np.ones_like(BOX_SIZES)


    if 'n_im3shape_results_files' in config:

        if config['n_im3shape_results_files'] == -1:
            n_files_use = len(files_im3)
        else:
            n_files_use = config['n_im3shape_results_files']
    else: 
        n_files_use = config['n_im3shape_results_files']
    
    for fi , fv in enumerate(files_im3[:n_files_use]):
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

        list_fwhm = np.array(list_fwhm)*DES_PIXEL_SIZE
        list_e1 = np.array(list_e1)
        list_e2 = np.array(list_e2)
        
        select= (np.array(list_fwhm) > 0) * (np.abs(np.array(list_e1) + 1j*np.array(list_e2)) < 1) 
        list_fwhm = list_fwhm[select]
        list_e1 = list_e1[select]
        list_e2 = list_e2[select]
        
        hist_fwhm,_=pl.histogram(list_fwhm,bins=config['bins_fwhm'])
        hist_e1,_=pl.histogram(list_e1,bins=config['bins_ell'])
        hist_e2,_=pl.histogram(list_e2,bins=config['bins_ell'])
        hist_snr,_=pl.histogram(cat['snr'],bins=config['bins_snr'])
        hist_box,_=pl.histogram(cat['stamp_size'],bins=config['bins_box'])

        hist_snr_all += hist_snr
        hist_box_all += hist_box
        hist_fwhm_all += hist_fwhm
        hist_e1_all   += hist_e1
        hist_e2_all   += hist_e2

        log.info('%3d %50s %8d %8d' , fi, fv, n_gals_total, len(list_fwhm))

        del(cat)


   
    pl.figure()
    pl.clf()

    pl.figure()
    pl.plot(config['bins_fwhm_centers'],hist_fwhm_all,'+-')
    pl.xlabel('exposure fwhm [arsec]')
    filename_fig = 'psf_fwhm_dist.png'
    pl.savefig(filename_fig)
    print 'saved' , filename_fig

    pl.figure()
    pl.plot(config['bins_ell_centers'],hist_e1_all,'+-')
    pl.xlabel('psf e1')
    filename_fig = 'psf_e1_dist.png'
    pl.savefig(filename_fig)
    print 'saved' , filename_fig

    pl.figure()
    pl.plot(config['bins_ell_centers'],hist_e2_all,'+-')
    pl.xlabel('psf e2')
    filename_fig = 'psf_e2_dist.png'
    pl.savefig(filename_fig)
    print 'saved' , filename_fig

    pl.figure()
    pl.plot(BOX_SIZES,hist_box_all/float(sum(hist_box_all)),'+-')
    pl.xlabel('box_size')
    filename_fig = 'boxsize_dist.png'
    pl.savefig(filename_fig)
    print 'saved' , filename_fig

    pl.figure()
    pl.plot(GRID_SNR,hist_snr_all,'+-')
    pl.xlabel('snr')
    filename_fig = 'snr_dist.png'
    pl.savefig(filename_fig)
    print 'saved' , filename_fig
   
    filename_snr = 'snr_dist.txt'
    hist_snr_all /= sum(hist_snr_all)
    data_save = np.array([GRID_SNR,hist_snr_all]).T
    np.savetxt(filename_snr,data_save,header='# bin_center_snr prob_snr')
    print 'saved' , filename_snr
   
    filename_psf = 'psf_fwhm_dist.txt'
    hist_fwhm_all /= sum(hist_fwhm_all)
    data_save = np.array([config['bins_fwhm_centers'],hist_fwhm_all]).T
    np.savetxt(filename_psf,data_save,header='# bin_center_fwhm prob_fwhm')
    print 'saved' , filename_psf

    filename_psf = 'psf_ell_dist.txt'
    hist_e1_all /= sum(hist_e1_all)
    hist_e2_all /= sum(hist_e2_all)
    data_save = np.array([config['bins_ell_centers'],hist_e1_all,hist_e2_all]).T
    np.savetxt(filename_psf,data_save,header='# bin_center_e e1 e2')
    print 'saved' , filename_psf




def get_psf_dist():
    
    p_fwhm = np.ones_like(config['bins_fwhm_centers'])  / float(len(config['bins_fwhm_centers']))
    p_ell  = np.ones_like(config['bins_ell_centers'])  / float(len(config['bins_ell_centers']))

    filename_psf = 'psf_fwhm_dist.txt'

    data_save = np.array([config['bins_fwhm_centers'],p_fwhm]).T
    np.savetxt(filename_psf,data_save,header='# bin_center_fwhm prob_fwhm')
    print 'saved' , filename_psf

    filename_psf = 'psf_ell_dist.txt'
    data_save = np.array([config['bins_ell_centers'],p_ell,p_ell]).T
    np.savetxt(filename_psf,data_save,header='# bin_center_e e1 e2')
    print 'saved' , filename_psf


def median_absolute_deviation(data):

    import numpy as np
    return np.median(np.absolute(data - np.median(data)))

def get_psf_images():



    filename_lores = os.path.join(args.out_dir,'nbc.psf.lores.fits')
    filename_hires = os.path.join(args.out_dir,'nbc.psf.hires.fits')
    filename_field = os.path.join(args.out_dir,'nbc.psf.field.fits')
    
    if os.path.isfile(filename_lores): os.remove(filename_lores); log.info('removed existing %s',filename_lores)
    if os.path.isfile(filename_hires): os.remove(filename_hires); log.info('removed existing %s',filename_hires)
    if os.path.isfile(filename_field): os.remove(filename_field); log.info('removed existing %s',filename_field)

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

    n_psfs = len(config['bins_fwhm_centers'])*len(config['bins_ell_centers'])**2
    
    iall = 0
    for ifwhm,fwhm in enumerate(config['bins_fwhm_centers']):
        for ie1,e1 in enumerate(config['bins_ell_centers']):
            for ie2,e2 in enumerate(config['bins_ell_centers']):
                                                          
                log.debug('getting single PSF at the pixel scale of a galaxy')


                config_copy1=copy.deepcopy(config_psf)
                config_copy1['psf']['fwhm'] = fwhm
                config_copy1['psf']['ellip']['g1'] = e1
                config_copy1['psf']['ellip']['g2'] = e2
                img_gal,img_psf,_,_  = galsim.config.BuildImages(config=config_copy1,image_num=0,obj_num=0,make_psf_image=True,nimages=1)   
                img_psf = img_psf[0]
                img_psf = img_psf[galsim.BoundsI(1, orig_image_size, 1, orig_image_size)]
                pyfits.append(filename_lores,img_psf.array)
                # filename_lores = 'nbc.psf.lores.%03d.fits' % iall
                # img_psf.write(filename_lores)
                                                              
                # now the hires PSF, centers in the middle
                # log.debug('getting single PSF at high resolution')                 
                # config_copy2=copy.deepcopy(config_psf)
                # n_sub = config['upsampling']
                # n_pad = config['padding']
                # n_pix_hires = (orig_image_size + n_pad) * n_sub
                # pixel_scale_hires = float(config_copy2['image']['pixel_scale']) / float(n_sub)
                # config_copy2['image']['pixel_scale'] = pixel_scale_hires
                # config_copy2['image']['size'] = n_pix_hires
                # config_copy2['psf']['fwhm'] = fwhm
                # config_copy2['psf']['ellip']['g1'] = e1
                # config_copy2['psf']['ellip']['g2'] = e2            

                # # no pixel convolut
                # config_copy2['pix'] = {}
                # config_copy2['pix']['type'] = 'Pixel'
                # config_copy2['pix']['xw'] = orig_pixel_scale

                # img_gal,img_psf,_,_ = galsim.config.BuildImages(config=config_copy2,image_num=0,obj_num=0,make_psf_image=True,nimages=1)      
                # img_psf = img_psf[0]
                # img_psf = img_psf[galsim.BoundsI(1, int(n_pix_hires), 1, int(n_pix_hires))]             
                # pyfits.append(filename_hires,img_psf.array)
                
                log.debug('getting single PSF at high resolution')                 
                config_copy2=copy.deepcopy(config_psf)
                n_sub = config['upsampling']
                n_pad = config['padding']
                n_pix_hires = (orig_image_size + n_pad) * n_sub
                pixel_scale_hires = float(config_copy2['image']['pixel_scale']) / float(n_sub)

                img_psf=galsim.ImageD(n_pix_hires,n_pix_hires)
                psf = galsim.Kolmogorov(fwhm=fwhm)
                psf.applyShear(g1=e1,g2=e2)
                pix = galsim.Pixel(scale=orig_pixel_scale)
                psfpix = galsim.Convolve([psf,pix])
                psfpix.draw(img_psf,scale=pixel_scale_hires)

                pyfits.append(filename_hires,img_psf.array)



                # img_loaded=pyfits.getdata(filename_hires,iall)
                # pl.subplot(1,3,1)
                # pl.imshow(img_psf.array,interpolation='nearest');
                # pl.subplot(1,3,2)
                # pl.imshow(img_loaded,interpolation='nearest'); 
                # pl.subplot(1,3,3)
                # pl.imshow(img_loaded-img_psf.array,interpolation='nearest'); 
                # pl.show()

                # filename_hires = 'nbc.psf.hires.%03d.fits' % iall
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
                # filename_field = 'nbc.psf.field.%03d.fits' % iall
                # img_psf.write(filename_field)

                log.info('generated id=%3d %3d/%d psfs fwhm=%2.5f e1=% 2.5f e2=% 2.5f ifwhm=%d ie1=%d ie2=%d ' , iall , iall+1, n_psfs , fwhm, e1, e2, ifwhm, ie1,ie2)
                iall += 1



    hdus_lores=pyfits.open(filename_lores)
    hdus_hires=pyfits.open(filename_hires)
    hdus_field=pyfits.open(filename_field)

    log.info('finished writing %s with %d hdus' , filename_lores , len(hdus_lores))
    log.info('finished writing %s with %d hdus' , filename_hires , len(hdus_hires))
    log.info('finished writing %s with %d hdus' , filename_field , len(hdus_field))
    



def get_meds(noise=True):

    # get the start and end index of files

    id_first = args.first
    id_last = id_first + args.num

    for ip in range(id_first,id_last):  
        for ig,vg in enumerate(config['shear']):

            filename_cat = os.path.join(args.out_dir,'nbc.truth.%03d.g%02d.fits' % (ip,ig))
            if noise: filename_meds = os.path.join(args.out_dir,'nbc.meds.%03d.g%02d.fits' % (ip,ig))
            else: filename_meds = os.path.join(args.out_dir,'nbc.meds.%03d.g%02d.noisefree.fits' % (ip,ig))

            config_copy = copy.deepcopy(config)
            if noise==False: del(config_copy['image']['noise'])
            config_copy['input']['catalog']['file_name'] = filename_cat
            config_copy['output']['file_name'] = filename_meds

            log.info('getting %s, noise=%s' % (filename_meds , noise ))
            if args.verbosity > 2:
                galsim.config.Process(config_copy,logger=log)
            else:
                galsim.config.Process(config_copy)
            fpack(filename_meds)
    
    log.info('done all meds')

def fpack(filename):

    filename_fz = filename + '.fz'
    if os.path.isfile(filename_fz):
        os.remove(filename_fz)

    cmd=['fpack' , '-t' , '10240,1' , filename]
    try:
        subprocess.call(cmd)
    except:
        raise Exception('command failed: %s', cmd )
    os.remove(filename)
    os.rename(filename_fz,filename)
    log.debug('compressed file %s ...' % filename)


def get_psf_index(ifwhm,ie1,ie2):

    filename_key = 'psf_key.fits'
    psf_table = tabletools.loadTable(filename_key)
    key_all , key_fwhm , key_e1 , key_e2  = psf_table['id_psf'] , psf_table['id_psf_fwhm'] , psf_table['id_psf_e1'] , psf_table['id_psf_e2']

    # for ii,vv in enumerate(psf_table):
    #     print 'psf_table index=%3d ifwhm=%d ie1=%d ie2=%d' % (key_all[ii] , key_fwhm[ii], key_e1[ii], key_e2[ii])


    indices = np.ones_like(ifwhm)
    for ii,vv in enumerate(ifwhm):
        select = (ifwhm[ii]==key_fwhm) * (ie1[ii]==key_e1) * (ie2[ii]==key_e2)
        nz = np.nonzero(select)[0][0]
        indices[ii] = key_all[nz] 
        # print 'index=%3d ifwhm=%d ie1=%d ie2=%d' % (indices[ii] , ifwhm[ii], ie1[ii], ie2[ii])



    # import pdb; pdb.set_trace()

    return indices

def get_psf_key():

    n_psfs = len(config['bins_fwhm_centers'])*len(config['bins_ell_centers'])**2
    psf_key = np.zeros([n_psfs,7])
    iall = 0
    for ifwhm,fwhm in enumerate(config['bins_fwhm_centers']):
        for ie1,e1 in enumerate(config['bins_ell_centers']):
            for ie2,e2 in enumerate(config['bins_ell_centers']):
                                                          
                psf_key[iall,:] = iall,ifwhm,ie1,ie2,fwhm,e1,e2
                iall+=1

    filename_key = 'psf_key.fits'
    tabletools.saveTable(filename_key,psf_key,dtype=dtype_psfkey)
    log.info('saved %s' , filename_key )

def get_useful_cosmos_galaxies_ids():

    filename_cosmos_catalog = os.path.join(config['input']['real_catalog']['dir'],config['input']['real_catalog']['file_name'])
    filename_cosmos_catalog_fits = os.path.join(config['input']['real_catalog']['dir'],config['input']['real_catalog']['file_name']).replace('.fits','_fits.fits')
    cosmos_catalog = pyfits.getdata(filename_cosmos_catalog)
    cosmos_catalog_fits = pyfits.getdata(filename_cosmos_catalog_fits)
    n_cosmos_gals = len(cosmos_catalog)
    log.info('opened %s with %d images' , filename_cosmos_catalog, n_cosmos_gals)

    size_cut = 0.5
    select = np.array(cosmos_catalog_fits['sersicfit'][:,1]*ACS_PIXEL_SCALE > size_cut)
    cosmos_index = np.arange(len(cosmos_catalog_fits))
    import pdb; pdb.set_trace()
    ids=cosmos_index[select]
    print 'full catalog n=' , len(cosmos_catalog_fits)
    print 'catalog with size_cut<%f, n=' , size_cut, len(ids)

    filename_ids = 'cosmos_good_ids.txt'
    tabletools.saveTable(filename_ids,ids)




def get_truth_catalogs():

    filename_psf_fwhm = 'psf_fwhm_dist.txt'
    filename_psf_ell = 'psf_ell_dist.txt'
    
    bins_psf_fwhm,prob_psf_fwhm=np.loadtxt(filename_psf_fwhm).T
    bins_psf_e,prob_psf_e1,prob_psf_e2=np.loadtxt(filename_psf_ell).T
    
    n_shears = len(config['shear'])
    n_gals = int(float(config['n_gals_per_file'])) 
    
    filename_cosmos_catalog = os.path.join(config['input']['real_catalog']['dir'],config['input']['real_catalog']['file_name'])
    cosmos_catalog = pyfits.getdata(filename_cosmos_catalog)
    n_cosmos_gals = len(cosmos_catalog)
    log.info('opened %s with %d images' , filename_cosmos_catalog, n_cosmos_gals)  

    id_first = args.first
    id_last = id_first + args.num

    for ip in range(id_first,id_last):
  
        for ig,vg in enumerate(config['shear']):
            

            filename_cat = os.path.join(args.out_dir,'nbc.truth.%03d.g%02d.fits' % (ip,ig))
            
            catalog = np.zeros(n_gals,dtype=dtype_truth)
            
            ids = np.arange(n_gals)
            cosmos_ids = np.random.choice(a=n_cosmos_gals,size=n_gals)
            rotation_angle = np.random.uniform(low=0,high=2*np.pi,size=n_gals)                  
            shear_ids = np.ones_like(ids) * ig
            shear_g1  = np.ones_like(ids) * vg[0]
            shear_g2  = np.ones_like(ids) * vg[1]
            psf_fwhm_ids = np.random.choice(a=len(bins_psf_fwhm),size=n_gals,p=prob_psf_fwhm)
            psf_e1_ids = np.random.choice(a=len(bins_psf_e),size=n_gals,p=prob_psf_e1)
            psf_e2_ids = np.random.choice(a=len(bins_psf_e),size=n_gals,p=prob_psf_e2)
            psf_fwhm = bins_psf_fwhm[psf_fwhm_ids]
            psf_e1 = bins_psf_e[psf_e1_ids]
            psf_e2 = bins_psf_e[psf_e2_ids]
            psf_ids = get_psf_index(psf_fwhm_ids, psf_e1_ids, psf_e2_ids)

            catalog['id'] = ids
            catalog['id_cosmos'] = cosmos_ids
            catalog['id_shear'] = shear_ids
            catalog['id_psf'] = psf_ids
            catalog['g1_true'] = shear_g1
            catalog['g2_true'] = shear_g2
            catalog['psf_fwhm'] = psf_fwhm
            catalog['psf_e1'] = psf_e1
            catalog['psf_e2'] = psf_e2
            catalog['rotation_angle'] = rotation_angle

            tabletools.saveTable(filename_cat,catalog)

def update_truth_table(update_snr=True , update_cosmos=True , update_hsm=True, update_fwhm=True):

    log.info('getting snr, flux and fwhm for the truth table')

    noise_std = config['des_pixel_noise_sigma']

    id_first = args.first
    id_last = id_first + args.num

    psf_images = None

    filename_cosmos_catalog = os.path.join(config['input']['real_catalog']['dir'],config['input']['real_catalog']['file_name'])
    filename_cosmos_catalog_fits = os.path.join(config['input']['real_catalog']['dir'],config['input']['real_catalog']['file_name']).replace('.fits','_fits.fits')
    cosmos_catalog = pyfits.getdata(filename_cosmos_catalog)
    cosmos_catalog_fits = pyfits.getdata(filename_cosmos_catalog_fits)
    n_cosmos_gals = len(cosmos_catalog)
    log.info('opened %s with %d images' , filename_cosmos_catalog, n_cosmos_gals)

    filename_great3_info = os.path.join(config['input']['real_catalog']['dir'],'real_galaxy_selection_info.fits')
    great3_info = np.array(pyfits.getdata(filename_great3_info))


    for ip in range(id_first,id_last):

        # all_snr=[]
   
        for il,vl in enumerate(config['shear']):

            list_normsq = []

            filename_cat = os.path.join(args.out_dir,'nbc.truth.%03d.g%02d.fits' % (ip,il))
            filename_meds = os.path.join(args.out_dir,'nbc.meds.%03d.g%02d.noisefree.fits' % (ip,il))


            log.info('part %d shear %d : getting snr, flux, hsm, and fwhm' , ip, il)
            log.debug('using %s %s', filename_meds, filename_cat)

            cat = tabletools.loadTable(filename_cat)
            n_gals = len(cat)

            # assure backwards compatibility
            if 'hsm_obs_g1'             not in cat.dtype.names: cat=tabletools.appendColumn(rec=cat, name='hsm_obs_g1',         arr=np.zeros(len(cat)), dtype='f8')
            if 'hsm_obs_g2'             not in cat.dtype.names: cat=tabletools.appendColumn(rec=cat, name='hsm_obs_g2',         arr=np.zeros(len(cat)), dtype='f8')
            if 'hsm_cor_g1'             not in cat.dtype.names: cat=tabletools.appendColumn(rec=cat, name='hsm_cor_g1',         arr=np.zeros(len(cat)), dtype='f8')
            if 'hsm_cor_g2'             not in cat.dtype.names: cat=tabletools.appendColumn(rec=cat, name='hsm_cor_g2',         arr=np.zeros(len(cat)), dtype='f8')
            if 'hsm_obs_sigma'          not in cat.dtype.names: cat=tabletools.appendColumn(rec=cat, name='hsm_obs_sigma',      arr=np.zeros(len(cat)), dtype='f8')
            if 'hsm_cor_sigma'          not in cat.dtype.names: cat=tabletools.appendColumn(rec=cat, name='hsm_cor_sigma',      arr=np.zeros(len(cat)), dtype='f8')
            if 'hsm_centroid_x'         not in cat.dtype.names: cat=tabletools.appendColumn(rec=cat, name='hsm_centroid_x',     arr=np.zeros(len(cat)), dtype='f8')
            if 'hsm_centroid_y'         not in cat.dtype.names: cat=tabletools.appendColumn(rec=cat, name='hsm_centroid_y',     arr=np.zeros(len(cat)), dtype='f8')
            if 'hsm_mom_amp'            not in cat.dtype.names: cat=tabletools.appendColumn(rec=cat, name='hsm_mom_amp',        arr=np.zeros(len(cat)), dtype='f8')
            if 'fwhm'                   not in cat.dtype.names: cat=tabletools.appendColumn(rec=cat, name='fwhm',               arr=np.zeros(len(cat)), dtype='f8')
            if 'sf_i'                   not in cat.dtype.names: cat=tabletools.appendColumn(rec=cat, name='sf_i',               arr=np.zeros(len(cat)), dtype='f8')    
            if 'sf_hlr'                 not in cat.dtype.names: cat=tabletools.appendColumn(rec=cat, name='sf_hlr',             arr=np.zeros(len(cat)), dtype='f8')    
            if 'sf_sersicn'             not in cat.dtype.names: cat=tabletools.appendColumn(rec=cat, name='sf_sersicn',         arr=np.zeros(len(cat)), dtype='f8')        
            if 'sf_q'                   not in cat.dtype.names: cat=tabletools.appendColumn(rec=cat, name='sf_q',               arr=np.zeros(len(cat)), dtype='f8')    
            if 'sf_boxiness'            not in cat.dtype.names: cat=tabletools.appendColumn(rec=cat, name='sf_boxiness',        arr=np.zeros(len(cat)), dtype='f8')        
            if 'sf_phi'                 not in cat.dtype.names: cat=tabletools.appendColumn(rec=cat, name='sf_phi',             arr=np.zeros(len(cat)), dtype='f8')    
            if 'zphot'                  not in cat.dtype.names: cat=tabletools.appendColumn(rec=cat, name='zphot',              arr=np.zeros(len(cat)), dtype='f8')   
            if 'psf_fwhm_measured'      not in cat.dtype.names: cat=tabletools.appendColumn(rec=cat, name='psf_fwhm_measured',  arr=np.zeros(len(cat)), dtype='f8')   
            if 'cosmos_mag_auto'        not in cat.dtype.names: cat=tabletools.appendColumn(rec=cat, name='cosmos_mag_auto',    arr=np.zeros(len(cat)), dtype='f8')   
            if 'cosmos_flux_radius'     not in cat.dtype.names: cat=tabletools.appendColumn(rec=cat, name='cosmos_flux_radius', arr=np.zeros(len(cat)), dtype='f8')   
            if 'mean_rgpp_rp'           not in cat.dtype.names: cat=tabletools.appendColumn(rec=cat, name='mean_rgpp_rp',       arr=np.zeros(len(cat)), dtype='f8')   
            if 'to_use'                 not in cat.dtype.names: cat=tabletools.appendColumn(rec=cat, name='to_use',             arr=np.zeros(len(cat)), dtype='i4')   

            for ig in range(n_gals):

                if update_snr == True:
                    noisless_gals = meds.MEDS(filename_meds)
                    n_gals = len(noisless_gals._cat)
                    img_gal = noisless_gals.get_cutout(ig,0)
                    normsq= np.sum( img_gal.flatten()**2 )
                    snr = np.sqrt(normsq)/noise_std
                    flux = np.sum(img_gal.flatten())
                    cat[ig]['snr'] = snr
                    cat[ig]['flux'] = flux

                if update_cosmos == True:
                    current_id_cosmos = cat[ig]['id_cosmos']
                    cat[ig]['sf_i']               = cosmos_catalog_fits[current_id_cosmos]['sersicfit'][0]
                    cat[ig]['sf_hlr']             = cosmos_catalog_fits[current_id_cosmos]['sersicfit'][1]*ACS_PIXEL_SCALE
                    cat[ig]['sf_sersicn']         = cosmos_catalog_fits[current_id_cosmos]['sersicfit'][2]
                    cat[ig]['sf_q']               = cosmos_catalog_fits[current_id_cosmos]['sersicfit'][3]
                    cat[ig]['sf_boxiness']        = cosmos_catalog_fits[current_id_cosmos]['sersicfit'][4]
                    cat[ig]['sf_phi']             = cosmos_catalog_fits[current_id_cosmos]['sersicfit'][7]
                    cat[ig]['zphot']              = cosmos_catalog_fits[current_id_cosmos]['zphot']
                    cat[ig]['cosmos_mag_auto']    = cosmos_catalog_fits[current_id_cosmos]['mag_auto']
                    cat[ig]['cosmos_flux_radius'] = cosmos_catalog_fits[current_id_cosmos]['flux_radius']
                    cat[ig]['to_use']             = great3_info[current_id_cosmos]['to_use']

                if update_hsm==True:
                    img_gal = noisless_gals.get_cutout(ig,0)
                    img_psf = pyfits.getdata(os.path.join(args.out_dir,'nbc.psf.lores.fits'),cat[ig]['id_psf'])
                    gs_img_gal = image_array_to_galsim(img_gal)
                    gs_img_psf = image_array_to_galsim(img_psf)

                    try:
                        shearobj1=galsim.hsm.EstimateShear(gs_img_gal,gs_img_psf)
                        cat[ig]['hsm_cor_g1'] = shearobj1.corrected_e1 / 2.
                        cat[ig]['hsm_cor_g2'] = shearobj1.corrected_e2 / 2.
                        cat[ig]['hsm_mom_amp'] = shearobj1.moments_amp
                    except:
                        log.error('HSM failed for object ig=%d ip=%d id_cosmos=%d psf_fwhm=%2.2f' , ig, ip , cat['id_cosmos'][ig] , cat['psf_fwhm'][ig])
                        cat[ig]['hsm_cor_g1'] = -99
                        cat[ig]['hsm_cor_g1'] = -99
                        cat[ig]['hsm_mom_amp'] = -99

                    try:
                        shearobj2=galsim.hsm.FindAdaptiveMom(gs_img_gal)
                        cat[ig]['hsm_obs_g1'] = shearobj2.observed_shape.g1
                        cat[ig]['hsm_obs_g2'] = shearobj2.observed_shape.g2
                        cat[ig]['hsm_obs_sigma'] = shearobj2.moments_sigma
                        cat[ig]['hsm_cor_sigma'] = shearobj2.moments_sigma
                        cat[ig]['hsm_centroid_x'] = shearobj2.moments_centroid.x
                        cat[ig]['hsm_centroid_y'] = shearobj2.moments_centroid.y
                    except:
                        cat[ig]['hsm_obs_g1'] = -99
                        cat[ig]['hsm_obs_g2'] = -99
                        cat[ig]['hsm_obs_sigma'] = -99
                        cat[ig]['hsm_cor_sigma'] = -99
                        cat[ig]['hsm_centroid_x'] = -99
                        cat[ig]['hsm_centroid_y'] = -99


                if update_fwhm==True:

                    try:
                        noisless_gals = meds.MEDS(filename_meds)
                        img_gal = noisless_gals.get_cutout(ig,0)
                        import mathstools
                        cat[ig]['fwhm'] = mathstools.get_2D_fwhm(img_gal)
                    except:
                        log.error('getting FWHM failed for galaxy %d in %s' , ig , filename_meds )
                        cat[ig]['fwhm'] = 666

                    try:
                        # if psf_images == None: psf_images = pyfits.open(os.path.join(args.out_dir,'nbc.psf.hires.fits'))
                        # img_hires_psf = psf_images[cat[ig]['id_psf']].data
                        # cat[ig]['psf_fwhm_measured'] = mathstools.get_2D_fwhm(img_hires_psf)
                        cat[ig]['psf_fwhm_measured'] = cat[ig]['psf_fwhm'] 


                    except:
                        log.error('getting PSF FWHM failed for galaxy %d in %s' , ig , filename_meds )
                        cat[ig]['psf_fwhm_measured'] = 666

                    if (cat[ig]['fwhm'] != 666) & (cat[ig]['psf_fwhm_measured'] != 666):
                        cat[ig]['mean_rgpp_rp'] = cat[ig]['fwhm']/cat[ig]['psf_fwhm_measured']
                    else:
                        cat[ig]['mean_rgpp_rp'] = 666

                                
                if ig % 100 == 0: log.debug('getting snr, flux, hsm and fwhm of galaxy %d' , ig)


            tabletools.saveTable(filename_cat, cat)



def get_snr_from_cosmos():

    # filename_noiseless = 'nbc.meds.%03d.g%02d.fits' % (0,0)
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
    m=meds.MEDS('nbc.meds.000.g05.fits.fz')
    for i in range(100): print tru['snr'][i];pl.imshow(m.get_cutout(i,0),interpolation='nearest');pl.show()

def main():


    global log , config , args , cat

    valid_actions = ['prepare', 'generate-psf', 'generate-truth', 'generate-noiseless', 'update-truth', 'generate-noisy']
    description = 'Get input catalogs for GREAT-DES NBC2 simulation'
    parser = argparse.ArgumentParser(description=description, add_help=True)
    parser.add_argument('-v', '--verbosity', type=int, action='store', default=2, choices=(0, 1, 2, 3 ), help='integer verbosity level: min=0, max=3 [default=2]')
    parser.add_argument('-c', '--filename_config', type=str, action='store', default='nbc.yaml', help='name of the config file')
    parser.add_argument('-f', '--first', type=int, action='store', default=0, help='index of the first file to create')
    parser.add_argument('-n', '--num', type=int, action='store', default=-1, help='number of files to create, if -1 then =config[n_files]')
    parser.add_argument('-a', '--actions', nargs='+' , type=str, action='store',  help='valid actions %s' % str(valid_actions))
    parser.add_argument('-o', '--out_dir', type=str, default='./',  action='store',  help='directory which stores output truth and meds files')

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

    if args.actions==None:
        raise Exception('supply at least one of the actions -a %s' % str(valid_actions))

    for act in args.actions:
        if act not in valid_actions:
            raise Exception('%s not a valid action. Choose from %s' % (act,valid_actions))

    config['bins_fwhm_centers'] = np.linspace(config['grid_fwhm']['min'],config['grid_fwhm']['max'],config['grid_fwhm']['n_grid'])
    config['bins_ell_centers']  = np.linspace(config['grid_ell']['min'],config['grid_ell']['max'],config['grid_ell']['n_grid'])
    config['bins_snr_centers']  = GRID_SNR
    config['bins_fwhm'] = plotstools.get_bins_edges( config['bins_fwhm_centers'] )
    config['bins_ell']  = plotstools.get_bins_edges( config['bins_ell_centers'] )
    config['bins_box']  = plotstools.get_bins_edges( BOX_SIZES )
    config['bins_snr']  = plotstools.get_bins_edges( GRID_SNR )
   
    if 'prepare' in args.actions:
        # get_useful_cosmos_galaxies_ids()
        if config['population_source'] == 'flat':
            get_psf_dist();     
        # if config['population_source'] == 'des':
        #     get_params_dist_from_DES();     
        # get_noise_level()
    if 'generate-psf' in args.actions:
        get_psf_key()
        get_psf_images()
    if 'generate-truth' in args.actions:
        get_truth_catalogs()
    if 'generate-noiseless' in args.actions:
        get_meds(noise=False)
    if 'update-truth' in args.actions:
        update_truth_table(update_snr=config['update_snr'] , update_cosmos=config['update_cosmos'] , update_hsm=config['update_hsm'], update_fwhm=config['update_fwhm'])
    if 'generate-noisy' in args.actions:
        get_meds(noise=True)
    
    log.info(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))

main()

