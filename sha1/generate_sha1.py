import os, galsim, numpy, argparse, copy, logging, sys, yaml, subprocess, time, itertools, copy

dtype_table_stats = { 'names'   : ['n_gals','n_fail','g1','g2','size','stdv_g1','stdv_g2','stdm_g1','stdm_g2','stdv_size','stdm_size'] ,
                     'formats' : ['i8'] *2+ ['f8']*9 } 

STD_DE_TARGET = 0.001
DEFAULT_STD_E = 0.25

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
      
        import pdb; pdb.set_trace()
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

def get_catalog_O2():

    global config
    config = yaml.load(open(args.filename_config))
    filename_catalog = os.path.basename(args.filename_config.replace('.yaml','') + '-O2.cat')
    file_catalog = open(filename_catalog,'w')
    catalog_header='#index filename_meds n_gals\
                    ihlr isnr iellip iangdeg insersic imoffat_beta imoffat_g imoffat_fwhm \
                    hlr snr ellip angdeg nsersic moffat_beta moffat_g moffat_fwhm fwhm_obj_over_fwhm_psf fwhm_obj fwhm_psf\n'
    file_catalog.write(catalog_header)

    total_n_gals = 0
    total_n_files = 0

    iall = 0

    order2_params = config['order2']['deviations'].keys()

    # keep track of combination so ther is no duplications

    for ihlr_snr,vhlr_snr in enumerate(config['order2']['hlr_snr']):

        used_list = []

        vhlr = vhlr_snr[0]
        vsnr = vhlr_snr[1]

        logger.debug('order 1 param hlr %f' % vhlr)
        logger.debug('order 1 param snr %f' % vsnr)


        # we run full grid only for start and end of the param space for size and snr

        for ip, vp in enumerate(order2_params):

            logger.debug('order 2 param %s' % vp)
            grid_param = config['order2']['deviations'][vp]['list']


            for ig, vg in enumerate(grid_param):                
                
                # assign the fiducial
                iellip = config['order2']['deviations']['ellip']['fid']
                iangdeg = config['order2']['deviations']['angdeg']['fid']
                insersic = config['order2']['deviations']['nsersic']['fid']
                imoffat_beta = config['order2']['deviations']['moffat_beta']['fid']
                imoffat_g = config['order2']['deviations']['moffat_g']['fid']
                imoffat_fwhm = config['order2']['deviations']['moffat_fwhm']['fid']

                vellip = config['order2']['deviations']['ellip']['list'][iellip]
                vangdeg = config['order2']['deviations']['angdeg']['list'][iangdeg]
                vnsersic = config['order2']['deviations']['nsersic']['list'][insersic]
                vmoffat_beta = config['order2']['deviations']['moffat_beta']['list'][imoffat_beta]
                vmoffat_g = config['order2']['deviations']['moffat_g']['list'][imoffat_g]
                vmoffat_fwhm = config['order2']['deviations']['moffat_fwhm']['list'][imoffat_fwhm]


                # assign the deviation - modify one already written by fiducial, deviations index starts from one
                exec('v%s=vg' % vp)
                exec('i%s=ig' % vp)
                logger.debug('current param %s %d %f' % (vp,ig,eval('v%s' % vp)))

                # avoid duplications
                current_list = [iellip,iangdeg,insersic,imoffat_beta,imoffat_g,imoffat_fwhm]
                if current_list in used_list:
                    logger.info('skipping this combination as already present')
                    continue

                iall += 1

                # calculate the number of packages for that settings
                if args.debug:
                    n_gals,std_e = 100,0.1
                else:
                    n_gals,std_e = get_n_gals(iall-1, 2)                    
                                                                                       
                # filename_meds = 'sha1-hlr%1d-snr%1d-ell%1d-ang%1d-ser%1d-psfb%1d-psfg%1d-psfs%1d.meds.fits' % ( 
                                    # ihlr,isnr,iellip,iangdeg,insersic,imoffat_beta,imoffat_g,imoffat_fwhm)
                    
                filename_meds = os.path.basename(args.filename_config.replace('.yaml','') + '.O2-%03d.meds.fits' % iall)

                logger.info('%3d current param %s=%2.2e\t processing file %s with std_e=%2.2e and n_gals=%10d' % (
                                    iall,vp,vg,filename_meds,std_e,n_gals))

                                   
                config_copy = copy.deepcopy(config)
                config_copy['gal']['half_light_radius'] = vhlr 
                config_copy['gal']['signal_to_noise'] = vsnr
                config_copy['gal']['shear']['g'] = vellip
                config_copy['gal']['shear']['beta'] = '%f degrees' % vangdeg
                config_copy['gal']['n'] = vnsersic
                config_copy['psf']['beta'] = vmoffat_beta
                config_copy['psf']['fwhm'] = vmoffat_fwhm
                config_copy['psf']['ellip']['g'] = vmoffat_g
                
                config_copy['output']['file_name'] = filename_meds
                config_copy['output']['nobjects'] = n_gals                                  

                filename_meds_fz = filename_meds + '.fz'

                if not args.dry:

                    galsim.config.Process(config_copy,logger=logger_config)              

                    # compress the meds file
                    fpack(filename_meds)

                    save_psf_img(config_copy,filename_meds_fz)

                fwhm_obj_over_fwhm_psf,fwhm_obj,fwhm_psf = compute_rgpp_rp(config_copy)              

                # write the file details in the catalog
                line_fmt = '%d\t%s\t%d\t'+ '%d\t'*8 +'% .8f'*11 +'\n' 

                line = line_fmt % ( iall,filename_meds_fz,n_gals,
                                    ihlr_snr,ihlr_snr,iellip,iangdeg,insersic,imoffat_beta,imoffat_g,imoffat_fwhm,
                                    vhlr,vsnr,vellip,vangdeg,vnsersic,vmoffat_beta,vmoffat_g,vmoffat_fwhm,fwhm_obj_over_fwhm_psf,fwhm_obj,fwhm_psf)
                file_catalog.write(line)
                file_catalog.flush()

                total_n_files += 1
                total_n_gals += n_gals
                used_list.append(current_list)

    logger.info('total_n_gals=%d' % total_n_gals)
    logger.info('total_n_files=%d' % total_n_files)
    file_catalog.close()


def get_catalog_O1():

    global config
    config = yaml.load(open(args.filename_config))
    filename_catalog = os.path.basename(args.filename_config.replace('.yaml','') + '-O1.cat')
    file_catalog = open(filename_catalog,'w')
    catalog_header='#index filename_meds n_gals\
                    ihlr isnr iellip iangdeg insersic imoffat_beta imoffat_g imoffat_fwhm \
                    hlr snr ellip angdeg nsersic moffat_beta moffat_g moffat_fwhm fwhm_obj_over_fwhm_psf fwhm_obj fwhm_psf\n'
    file_catalog.write(catalog_header)

    total_n_gals = 0
    total_n_files = 0

    iall = 0

    order2_params = config['order2'].keys()

    # keep track of combination so ther is no duplications
    used_list = []

    # assign the fiducial
    iellip = config['order2']['deviations']['ellip']['fid']
    iangdeg = config['order2']['deviations']['angdeg']['fid']
    insersic = config['order2']['deviations']['nsersic']['fid']
    imoffat_beta = config['order2']['deviations']['moffat_beta']['fid']
    imoffat_g = config['order2']['deviations']['moffat_g']['fid']
    imoffat_fwhm = config['order2']['deviations']['moffat_fwhm']['fid']

    vellip = config['order2']['deviations']['ellip']['list'][iellip]
    vangdeg = config['order2']['deviations']['angdeg']['list'][iangdeg]
    vnsersic = config['order2']['deviations']['nsersic']['list'][insersic]
    vmoffat_beta = config['order2']['deviations']['moffat_beta']['list'][imoffat_beta]
    vmoffat_g = config['order2']['deviations']['moffat_g']['list'][imoffat_g]
    vmoffat_fwhm = config['order2']['deviations']['moffat_fwhm']['list'][imoffat_fwhm]


    for ihlr,vhlr in enumerate(config['order1']['hlr']):
        for isnr,vsnr in enumerate(config['order1']['snr']):

            iall+=1

            logger.debug('order 1 param hlr %f' % vhlr)
            logger.debug('order 1 param snr %f' % vsnr)

            # calculate the number of packages for that settings
            if args.debug:
                n_gals,std_e = 100,0.1
            else:
                n_gals,std_e = get_n_gals(iall-1, 1)
                
                
            filename_meds = os.path.basename(args.filename_config.replace('.yaml','') + '.O1-%03d.meds.fits' % iall)

            logger.info('%3d current params hlr=%2.2f\tsnr=%2.2f\t processing file %s with std_e=%2.2e and n_gals=%10d' % (
                                iall,vhlr,vsnr,filename_meds,std_e,n_gals))
                               
            config_copy = copy.deepcopy(config)
            config_copy['gal']['half_light_radius'] = vhlr 
            config_copy['gal']['signal_to_noise'] = vsnr
            config_copy['gal']['shear']['g'] = vellip
            config_copy['gal']['shear']['beta'] = '%f degrees' % vangdeg
            config_copy['gal']['n'] = vnsersic
            config_copy['psf']['beta'] = vmoffat_beta
            config_copy['psf']['fwhm'] = vmoffat_fwhm
            config_copy['psf']['ellip']['g'] = vmoffat_g

            config_copy['output']['file_name'] = filename_meds
            config_copy['output']['nobjects'] = n_gals                                  

            filename_meds_fz = filename_meds + '.fz'

            if not args.dry:

                galsim.config.Process(config_copy,logger=logger_config)              

                # compress the meds file
                logger.debug('created file %s, compressing...' % filename_meds)
                fpack(filename_meds)

                save_psf_img(config_copy,filename_meds_fz)

            fwhm_obj_over_fwhm_psf,fwhm_obj,fwhm_psf = compute_rgpp_rp(config_copy)              

            # write the file details in the catalog
            line_fmt = '%d\t%s\t%d\t'+ '%d\t'*8 +'% .8f'*11 +'\n' 

            line = line_fmt % ( iall,filename_meds_fz,n_gals,
                                ihlr,isnr,iellip,iangdeg,insersic,imoffat_beta,imoffat_g,imoffat_fwhm,
                                vhlr,vsnr,vellip,vangdeg,vnsersic,vmoffat_beta,vmoffat_g,vmoffat_fwhm,fwhm_obj_over_fwhm_psf,fwhm_obj,fwhm_psf)
            file_catalog.write(line)
            file_catalog.flush()


            total_n_files += 1
            total_n_gals += n_gals
            
    logger.info('total_n_gals=%d' % total_n_gals)
    logger.info('total_n_files=%d' % total_n_files)
    file_catalog.close()


def get_n_gals(ident,order):

    if order == 1:
        # sha1-O2.cat.hsm.stats.stde.cat
        filename_stats = os.path.basename('%s-O1.cat.hsm.stats.stde.cat' % (args.filename_config.replace('.yaml','')))
        stats = numpy.loadtxt(filename_stats,dtype=dtype_table_stats)
        std_e = stats[ident]['stdv_g1']

    if order == 2:
        # sha1-O2.cat.hsm.stats.stde.cat
        filename_stats = os.path.basename('%s-O2.cat.hsm.stats.stde.cat' % (args.filename_config.replace('.yaml','')))
        stats = numpy.loadtxt(filename_stats,dtype=dtype_table_stats)

    # if HSM had all errors, then use default
    if stats[ident]['stdv_g1'] > 0.01:
        std_e = stats[ident]['stdv_g1']
    else:
        std_e = DEFAULT_STD_E

    n_gals = (std_e / STD_DE_TARGET)**2
    n_gals = my_round(n_gals,1000)

    return n_gals , std_e

def save_psf_img(config,filename_meds):

    logger.debug('getting single PSF at the pixel scale of a galaxy')
    config_copy = copy.deepcopy(config)
    galsim.config.ProcessInput(config_copy)
    img_gal,img_psf,_,_,_  = galsim.config.BuildImage(config=config_copy,image_num=0,obj_num=0,make_psf_image=True)

    filename_psf = '%s.psf.single.%dx%d.fits' % (filename_meds,config['image']['size'],config['image']['size'])
    img_psf.write(filename_psf)
    fpack(filename_psf)

    # now the hires PSF, centered in the middle
    logger.debug('getting single PSF at high resolution')
    config_copy = copy.deepcopy(config)

    n_sub = 5
    n_pad = 4
    n_pix_hires = (config['image']['size'] + n_pad) * n_sub
    pixel_scale_hires = float(config['image']['pixel_scale']) / float(n_sub)
    config_copy['image']['pixel_scale'] = pixel_scale_hires
    config_copy['image']['size'] = n_pix_hires

    img_gal,img_psf,_,_,_ = galsim.config.BuildImage(config=config_copy,image_num=0,obj_num=0,make_psf_image=True)    

    filename_psf_hires = '%s.psf.single.%dx%d.fits' % (filename_meds,n_pix_hires,n_pix_hires)
    img_psf.write(filename_psf_hires)
    fpack(filename_psf_hires)
    # logger.info('saved %s' % filename_psf_hires)


    # now field
    logger.debug('getting low res PSF in a field')
    config_copy = copy.deepcopy(config)
    config_copy['image']['type'] = 'Tiled'
    config_copy['image']['nx_tiles'] = 10
    config_copy['image']['ny_tiles'] = 10

    # This is the size of the postage stamps.
    config_copy['image']['stamp_xsize'] = config_copy['image']['size']
    config_copy['image']['stamp_ysize'] = config_copy['image']['size']
    del(config_copy['image']['size'])

    galsim.config.ProcessInput(config_copy)
    img_gal,img_psf,_,_ = galsim.config.BuildImage(config=config_copy,image_num=0,obj_num=0,make_psf_image=True)    

    filename_psf_field = '%s.psf.field.%dx%d.fits' % (filename_meds,config['image']['size'],config['image']['size'])
    img_psf.write(filename_psf_field)
    fpack(filename_psf_field)


def main():

    global logger , config , args , logger_config

    description = 'Plots for the results of nmb_main. To use in ipython, create a variable global results_array, global truth_array to speed up loading'
    parser = argparse.ArgumentParser(description=description, add_help=True)
    parser.add_argument('-v', '--verbosity', type=int, action='store', default=2, choices=(0, 1, 2, 3 ), help='integer verbosity level: min=0, max=3 [default=2]')
    parser.add_argument('-c', '--filename_config', default='sha1.yaml',type=str, action='store', help='name of the yaml config file')
    parser.add_argument('-d', '--dry', default=False,  action='store_true', help='Dry run, do not generate data')
    parser.add_argument('-r', '--redo', default=False, action='store_true', help='produce the files even if the exist')
    parser.add_argument('--debug', default=False, action='store_true', help='debug mode, runs on only a subset of galaxies')
    parser.add_argument('--o1', default=False, action='store_true', help='run order 1 set')
    parser.add_argument('--o2', default=False, action='store_true', help='run order 2 set')
    parser.add_argument('--o3', default=False, action='store_true', help='run nbc')

    args = parser.parse_args()
    # Parse the integer verbosity level from the command line args into a logging_level string
    logging_levels = { 0: logging.CRITICAL, 
                       1: logging.WARNING,
                       2: logging.INFO,
                       3: logging.DEBUG }
    logging_level = logging_levels[args.verbosity]
    logging.basicConfig(format="%(message)s", level=logging_level, stream=sys.stdout)
    logger = logging.getLogger("generate_catalog.py") 
    logger.setLevel(logging_level)

    if logging_level == logging.DEBUG:
        logging_level_config = logging.DEBUG
    else:
        logging_level_config = logging.WARNING

    logger_config = logging.getLogger("config") 
    logger_config.setLevel(logging_level_config)

    logger.info(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))

    if args.debug: 
        logger.critical('running in DEBUG MODE on 1000 galaxies')

    if args.o1:
        logger.info('------------- generating order 1 -------------')
        get_catalog_O1()
    if args.o2:
        logger.info('------------- generating order 2 ------------- ')
        get_catalog_O2()
    if args.o3:
        logger.info('------------- noise bias calibration not ready yet ------------- ')
        
    

if __name__ == "__main__":

    main()