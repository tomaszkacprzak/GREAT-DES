import os, galsim, numpy, argparse, copy, logging, sys, yaml, subprocess, time, itertools, copy, galsim.des, pyfits

dtype_table_stats = { 'names'   : ['index','n_gals','n_fail','g1','g2','size','stdv_g1','stdv_g2','stdm_g1','stdm_g2','stdv_size','stdm_size'] ,
                     'formats' : ['i8'] *3+ ['f8']*9 } 

N_GALS_COLS = 100
DEFAULT_STD_E = 0.25
N_GALS_DEBUG = 1000

def fpack(filename):

    cmd=['fpack' , '-t' , '10240,1' , filename]
    subprocess.call(cmd)
    os.remove(filename)
    logger.debug('compressing file %s ...' % filename + '.fz')

def compute_rgpp_rp(config):
    """
    @brief
    Computes rgpp_rp, i.e. pixel and psf convolved galaxy FWHM / psf FWHM
    for sersics and multisersics. Other models are not implemented in this
    function yet.
    @param result i3_results structure for either Sersics or Multisersics (must implement get_params())
    @psf i3_image image of the PSF with correct resolution (should have been used earlier for fitting)
    @options i3_options struct 
    @return fwhm_gal/fwhm_psf , fwhm_gal , fwhm_psf
    """

    upsampling = 11
    config = copy.deepcopy(config)
    config['gal']['shear']['g'] = 0.
    config['psf']['ellip']['g'] = 0.
    config['image']['draw_method'] = 'fft'
    config['image']['pixel_scale'] /= upsampling
    config['image']['size'] *= upsampling
    config['gal']['signal_to_noise'] = 10000000

    galsim.config.ProcessInput(config)
    img_gal,img_psf,_,_,_ = galsim.config.BuildSingleImage(config,image_num=0,obj_num=0,make_psf_image=True)

    fwhm_gal = get_FWHM(img_gal, upsampling=upsampling , radialProfile=False)
    fwhm_psf = get_FWHM(img_psf, upsampling=upsampling , radialProfile=False)

    return fwhm_gal/fwhm_psf , fwhm_gal , fwhm_psf

def get_centroid(image):

    x = numpy.argmax(image.sum(axis=1))
    y = numpy.argmax(image.sum(axis=0))

    return x,y

def get_FWHM(image, fwxm=0.5, upsampling=1, radialProfile=True):
    """
    Computes the FWHM of an i3 image. Per default, it computes the
    radial averaged profile for determining the FWHM. Alternatively,
    it computes the FWHMs of the profiles along the x and y axes in +
    and - direction (total of 4) and returns their average as an
    estimator of the FWHM.    
    """
    # compute weighted moments to get centroid
    # hsmresult = galsim.hsm.FindAdaptiveMom(image)

    x,y = get_centroid(image.array)

    # flip x<->y as ther is a difference in galsim indexing and python array indexing
    x0 = numpy.floor(x)
    y0 = numpy.floor(y)

    profile_x = image.array[int(x0), :]
    profile_y = image.array[:, int(y0)]

    max_val = image.array[int(x0), int(y0)]
    cut_val = max_val * fwxm

    if radialProfile:
        radii, profile = get_radialProfile(image, x0, y0, fwxm)

        diff = abs(profile - cut_val)
         
        # fhwm code from Tomek
        f1 = 0.
        f2 = 0.
        x1 = 0
        x2 = 0
      
        x1 = numpy.argmin(diff)
        f1 = profile[x1]
         
        if( f1 < cut_val ):  x2 = x1+1
        else:       x2 = x1-1
        f2 = profile[x2];
         
        a = (f1-f2)/(radii[x1] - radii[x2])
        b = f1 - a*radii[x1];
        x3 = (cut_val - b)/a;
         
        fwhm = (2.* x3) / upsampling

    else:
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
            
            x1 = numpy.argmin(diff)
            f1 = profile[x1]
         
            if( f1 < cut_val ):  x2 = x1+1
            else:       x2 = x1-1
            f2 = profile[x2];
         
            a = (f1-f2)/(x1 - x2)
            b = f1 - a*x1;
            x3 = (cut_val - b)/a;
         
            fwhms.append(2.* (dc0 + x3))

        fwhm =  numpy.mean(numpy.array(fwhms))/upsampling

    return fwhm

def my_round(x, base=1):
    return int(base * round(float(x)/base))

def copy_preloaded_config(config_processed, config_orig):

    config_new = copy.deepcopy(config_orig)
    config_new['real_catalog'] = config_processed['real_catalog']
    # config_new['n_exposures'] = config_processed['n_exposures']
    config_new['input_manager'] = config_processed['input_manager']
    config_new['real_catalog_safe'] = config_processed['real_catalog_safe']
    config_new['input']['real_catalog'] = config_processed['input']['real_catalog']
    config_new['seq_index'] = config_processed['seq_index']

    return config_new


def get_data():

    global config
    config = yaml.load(open(args.filename_config))
    filename_catalog = os.path.basename(args.filename_config.replace('.yaml','') + '.cat')
    file_catalog = open(filename_catalog,'w')

    # TODO: catalog data
    catalog_header='#index filename_meds n_gals ipsf_fwhm isnr ishear psf_fwhm snr shear1 shear2\n'
    file_catalog.write(catalog_header)

    total_n_gals = 0
    total_n_files = 0

    iall = 0

    # keep track of combination so ther is no duplications
    used_list = []

    config_copy_base = copy.deepcopy(config)
    config_copy_processed = copy.deepcopy(config)
    if not args.dry:
        galsim.config.ProcessInput(config_copy_processed)
    
    n_files =  len(config['grid']['shear'])*len(config['grid']['psf_fwhm'])*len(config['grid']['snr'])
    logger.info("will create %d files" % n_files)
    for isize,vsize in enumerate(config['grid']['psf_fwhm']):
        for isnr,vsnr in enumerate(config['grid']['snr']):
            for ishear, vshear in enumerate(config['grid']['shear']):

                # calculate the number of packages for that settings
                if args.debug:
                    n_gals,std_e = N_GALS_DEBUG,0.1
                else:
                    n_gals,std_e = get_n_gals(iall, 1)
                    
                filename_meds = os.path.basename(args.filename_config.replace('.yaml','') + '.%03d.fits' % iall)

                logger.info('%3d current params psf_fwhm=%2.2f\tsnr=%2.2f\tg1=% 2.3f\tg1=% 2.3f \t processing file %s with std_e=%2.2e and n_gals=%10d' % (
                                    iall,vsize,vsnr,vshear[0],vshear[1],filename_meds,std_e,n_gals))
                                   
                

                if args.fpack:
                    filename_meds_fz = filename_meds + '.fz'
                else:
                    filename_meds_fz = filename_meds

                # all this will be skipped in "wet" mode
                if not args.dry:
                    
                    config_use = copy_preloaded_config(config_copy_processed,config_copy_base)
                    config_use['output']['file_name'] = filename_meds
                    config_use['psf']['fwhm'] = vsize
                    config_use['gal']['signal_to_noise'] = vsnr
                    config_use['gal']['shear']['g1'] = vshear[0]
                    config_use['gal']['shear']['g2'] = vshear[1]

                    if config_use['output']['type']=='Fits':
                        config_use['image']['nx_tiles'] = n_gals/N_GALS_COLS
                        size_key = 'stamp_size'
                    elif config_use['output']['type'] == 'des_meds':
                        config_use['output']['nobjects'] = n_gals                                  
                        size_key = 'size'
                    
                    galsim.config.Process(config_use,logger=logger_config)              
                    
                    if args.fpack:
                        # compress the meds file
                        logger.debug('created file %s, compressing...' % filename_meds)
                        fpack(filename_meds)                

                    if config_use['output']['type']=='Fits':
                        save_psf_img_tiled(config_use,filename_meds_fz)
                    elif config_use['output']['type'] == 'des_meds':
                        save_psf_img_meds(config_use,filename_meds_fz)

                # save input catalog
                if config_copy_base['output']['type']=='Fits':
                    save_tiled_input_catalog(filename_meds,n_gals)

                # write the file details in the catalog
                #index filename_meds n_gals ipsf_fwhm isnr ishear psf_fwhm snr shear1 shear2 
                line_fmt = '%d\t%s\t%d\t'+ '%d\t'*3 +'% 2.8f\t'*4 +'\n' 
                line = line_fmt % ( iall,filename_meds_fz,n_gals,isize,isnr,ishear,vsize,vsnr,vshear[0],vshear[1] )
                file_catalog.write(line)
                file_catalog.flush()

                iall+=1
                total_n_files += 1
                total_n_gals += n_gals
            
    logger.info('total_n_gals=%d' % total_n_gals)
    logger.info('total_n_files=%d' % total_n_files)
    file_catalog.close()
    logger.info('saved %s' % filename_catalog)

def save_tiled_input_catalog(filename_meds,n_gals):

    filename_cat = filename_meds + '.txt'
    nx_tiles = n_gals / N_GALS_COLS 
    ny_tiles = N_GALS_COLS

    file_cat = open(filename_cat,'w')
    header = '# id x y mag(dummy) size(dummy)'
    iall = 0
    for ix in range(nx_tiles):
        for iy in range(ny_tiles):
            x = ix*config['cutout_size'] + config['cutout_size']/2.
            y = iy*config['cutout_size'] + config['cutout_size']/2.
            line = '%d\t%2.4f\t%2.4f\t0\t0\n' % (iall,x,y)
            file_cat.write(line)
            iall+=1
    file_cat.close()
    logger.info('saved %s' % filename_cat)



def get_n_gals(ident,order):

    filename_stats = os.path.basename('%s.cat.hsm.stats.stde.cat' % (args.filename_config.replace('.yaml','')))
    stats = numpy.loadtxt(filename_stats,dtype=dtype_table_stats)
    std_e = stats[ident]['stdv_g1']

    
    # if HSM had all errors, then use default
    if stats[ident]['stdv_g1'] > 0.01:
        std_e = stats[ident]['stdv_g1']
    else:
        std_e = DEFAULT_STD_E

    shears =  numpy.array([x[0] for x in config['grid']['shear']])

    n_shears = len(shears)
    var_e = std_e**2;
    var_s = numpy.mean(shears**2) - numpy.mean(shears)**2
    var_m_target = args.m_accuracy**2;

    n_gals_target = var_e/(var_s*var_m_target) 
    n_gals_per_shear = n_gals_target/n_shears

    # n_gals = (std_e / STD_DE_TARGET)**2
    n_gals = my_round(n_gals_per_shear,config['gal']['num']*1000)

    return n_gals , std_e

def save_psf_img_meds(config,filename_meds):


    config_copy=config
    hdulist=pyfits.HDUList()
    filename_psf = '%s.psf' % (filename_meds)

    orig_pixel_scale = config_copy['pixel_scale']
    orig_image_size = config_copy['cutout_size']

    logger.debug('getting single PSF at the pixel scale of a galaxy')
      
    img_gal,img_psf,_,_  = galsim.config.BuildImages(config=config,image_num=0,obj_num=0,make_psf_image=True,nimages=1,logger=logger_config)   
    img_psf[0].write(hdu_list=hdulist)
    
    # now the hires PSF, centered in the middle
    logger.debug('getting single PSF at high resolution')
    
    n_sub = 5
    n_pad = 4
    n_pix_hires = (orig_image_size + n_pad) * n_sub
    pixel_scale_hires = float(config['image']['pixel_scale']) / float(n_sub)
    config_copy['image']['pixel_scale'] = pixel_scale_hires
    config_copy['image']['size'] = n_pix_hires
    img_gal,img_psf,_,_ = galsim.config.BuildImages(config=config_copy,image_num=0,obj_num=0,make_psf_image=True,nimages=1,logger=logger_config)      
    img_psf[0].write(hdu_list=hdulist)

    # now field
    logger.debug('getting low res PSF in a field')
    config_copy['image']['pixel_scale'] = orig_pixel_scale
    config_copy['image']['stamp_size'] = orig_image_size
    config_copy['image']['type'] = 'Tiled'
    config_copy['image']['nx_tiles'] = 10
    config_copy['image']['ny_tiles'] = 10

    # This is the size of the postage stamps.
    config_copy['image']['stamp_xsize'] = orig_image_size
    config_copy['image']['stamp_ysize'] = orig_image_size

    if 'size' in config_copy['image']:
        del(config_copy['image']['size'])

def save_psf_img_tiled(config,filename_meds):


    config_copy=config
    hdulist=pyfits.HDUList()
    filename_psf = '%s.psf' % (filename_meds)

    orig_pixel_scale = config_copy['pixel_scale']
    orig_image_size = config_copy['cutout_size']

    config_copy['image']['nx_tiles'] = 2
    config_copy['image']['ny_tiles'] = 2
    config_copy['image']['nproc'] = 1

    logger.debug('getting single PSF at the pixel scale of a galaxy')
      
    img_gal,img_psf,_,_  = galsim.config.BuildImages(config=config,image_num=0,obj_num=0,make_psf_image=True,nimages=1,logger=logger_config)   
    img_psf = img_psf[0]
    img_psf = img_psf[galsim.BoundsI(1, orig_image_size, 1, orig_image_size)]
    img_psf.write(hdu_list=hdulist)
    
    # now the hires PSF, centered in the middle
    logger.debug('getting single PSF at high resolution')
    
    n_sub = 5
    n_pad = 4
    n_pix_hires = (orig_image_size + n_pad) * n_sub
    pixel_scale_hires = float(config['image']['pixel_scale']) / float(n_sub)
    config_copy['image']['pixel_scale'] = pixel_scale_hires
    config_copy['image']['stamp_size'] = n_pix_hires
    img_gal,img_psf,_,_ = galsim.config.BuildImages(config=config_copy,image_num=0,obj_num=0,make_psf_image=True,nimages=1,logger=logger_config)      
    img_psf = img_psf[0]
    img_psf = img_psf[galsim.BoundsI(1, int(n_pix_hires), 1, int(n_pix_hires))]
    img_psf.write(hdu_list=hdulist)

    # now field
    logger.debug('getting low res PSF in a field')
    config_copy['image']['pixel_scale'] = orig_pixel_scale
    config_copy['image']['stamp_size'] = orig_image_size
    config_copy['image']['type'] = 'Tiled'
    config_copy['image']['nx_tiles'] = 10
    config_copy['image']['ny_tiles'] = 10

    # This is the size of the postage stamps.
    config_copy['image']['stamp_xsize'] = orig_image_size
    config_copy['image']['stamp_ysize'] = orig_image_size

    if 'size' in config_copy['image']:
        del(config_copy['image']['size'])


    img_gal,img_psf,_,_ = galsim.config.BuildImage(config=config_copy,image_num=0,obj_num=0,make_psf_image=True,logger=logger_config)    
    img_psf.write(hdu_list=hdulist)

    hdulist.writeto(filename_psf,clobber=True)
    logger.info('saved %s' % filename_psf)


def main():

    global logger , config , args , logger_config

    description = 'Plots for the results of nmb_main. To use in ipython, create a variable global results_array, global truth_array to speed up loading'
    parser = argparse.ArgumentParser(description=description, add_help=True)
    parser.add_argument('-v', '--verbosity', type=int, action='store', default=2, choices=(0, 1, 2, 3, 4 ), help='integer verbosity level: min=0, max=4 [default=2]')
    parser.add_argument('-c', '--filename_config', default='sha1.yaml',type=str, action='store', help='name of the yaml config file')
    parser.add_argument('-m', '--m_accuracy', default=0.01,type=float, action='store', help='desired sigma_m')
    parser.add_argument('-d', '--dry', default=False,  action='store_true', help='Dry run, do not generate data')
    parser.add_argument('--debug', default=False, action='store_true', help='debug mode, runs on only a subset of galaxies')
    parser.add_argument('--fpack', default=False, action='store_true', help='compress images')
    
    args = parser.parse_args()
    # Parse the integer verbosity level from the command line args into a logging_level string
    logging_levels = { 0: logging.CRITICAL, 
                       1: logging.WARNING,
                       2: logging.INFO,
                       3: logging.DEBUG,
                       4: logging.DEBUG }
    logging_level = logging_levels[args.verbosity]
    logging.basicConfig(format="%(message)s", level=logging_level, stream=sys.stdout)
    logger = logging.getLogger("generate_catalog.py") 
    logger.setLevel(logging_level)

    if logging_level == logging.DEBUG:
        logging_level_config = logging.INFO
        if args.verbosity == 4:
            logging_level_config = logging.DEBUG
    else:
        logging_level_config = logging.WARNING

    logger_config = logging.getLogger("config") 
    logger_config.setLevel(logging_level_config)

    logger.info(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))

    if args.debug: 
        logger.critical('running in DEBUG MODE on %d galaxies' % N_GALS_DEBUG)

    STD_M_TARGET = args.m_accuracy

    get_data()
        
    logger.info(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))    

if __name__ == "__main__":

    main()