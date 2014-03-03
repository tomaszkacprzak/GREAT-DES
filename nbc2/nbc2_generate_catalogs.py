import os, sys, logging, yaml, argparse, time, meds, pyfits, plotstools, tabletools, nbc2_dtypes, galsim
import pylab as pl
import numpy as np
from astropy.stats import sigma_clip

logging_level = logging.INFO
log = logging.getLogger("nbc2_generate_catalogs") 
log.setLevel(logging_level)  
log_formatter = logging.Formatter("%(asctime)s %(name)s %(levelname)s   %(message)s", "%Y-%m-%d %H:%M:%S")
stream_handler = logging.StreamHandler(sys.stdout)
stream_handler.setFormatter(log_formatter)
log.addHandler(stream_handler)

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
    n_bins=1000
    bins_fwhm=np.linspace(0.5,3,n_bins)
    bins_ell=np.linspace(-0.2,0.2,n_bins)
    bins_snr=np.linspace(0,100,n_bins)

    hist_fwhm_all = np.zeros(n_bins-1)
    hist_e1_all   = np.zeros(n_bins-1)
    hist_e2_all   = np.zeros(n_bins-1)
    hist_snr_all  = np.zeros(n_bins-1)

    for fi , fv in enumerate(files_im3[:50]):
        try:
            cat=pyfits.getdata(fv)
        except:
            print 'opening failed' , fv
            continue

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

        list_fwhm = np.array(list_fwhm)*0.27
        list_e1 = np.array(list_e1)
        list_e2 = np.array(list_e2)
        
        select= (np.array(list_fwhm) > 0) * (np.abs(np.array(list_e1) + 1j*np.array(list_e2)) < 1) 
        list_fwhm = list_fwhm[select]
        list_e1 = list_e1[select]
        list_e2 = list_e2[select]
        
        select= (np.array(list_snr)>0)
        list_snr = list_snr[select]

        hist_fwhm,bins_fwhm=pl.histogram(list_fwhm,bins=bins_fwhm)
        hist_e1,bins_e1=pl.histogram(list_e1,bins=bins_ell)
        hist_e2,bins_e2=pl.histogram(list_e2,bins=bins_ell)
        hist_snr,bins_snr=pl.histogram(list_snr,bins=bins_snr)

        hist_fwhm_all += hist_fwhm
        hist_e1_all   += hist_e1
        hist_e2_all   += hist_e2
        hist_snr_all  += hist_snr

        log.info('%3d %50s %8d' , fi, fv, len(list_fwhm))

        del(cat)


   
    pl.figure()
    pl.clf()

    bins_fwhm_centered = plotstools.get_bins_centers(bins_fwhm)
    pl.plot(bins_fwhm_centered,hist_fwhm_all,'+-')
    pl.xlabel('exposure fwhm [arsec]')
    filename_fig = 'psf_fwhm_dist.png'
    pl.savefig(filename_fig)
    pl.close()
    print 'saved' , filename_fig

    bins_e1_centered = plotstools.get_bins_centers(bins_e1)
    pl.plot(bins_e1_centered,hist_e1_all,'+-')
    pl.xlabel('psf e1')
    filename_fig = 'psf_e1_dist.png'
    pl.savefig(filename_fig)
    pl.close()
    print 'saved' , filename_fig

    bins_e2_centered = plotstools.get_bins_centers(bins_e2)
    pl.plot(bins_e2_centered,hist_e2_all,'+-')
    pl.xlabel('psf e2')
    filename_fig = 'psf_e2_dist.png'
    pl.savefig(filename_fig)
    pl.close()
    print 'saved' , filename_fig

    bins_snr_centered = plotstools.get_bins_centers(bins_snr)
    pl.plot(bins_snr_centered,hist_snr_all,'+-')
    pl.xlabel('snr')
    filename_fig = 'snr_dist.png'
    pl.savefig(filename_fig)
    pl.close()
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

def get_catalogs():

    filename_psf_fwhm = 'psf_fwhm_dist.txt'
    filename_psf_ell = 'psf_ell_dist.txt'
    filename_snr = 'snr_dist.txt'

    bins_psf_fwhm,prob_psf_fwhm=np.loadtxt(filename_psf_fwhm).T
    bins_psf_e,prob_psf_e1,prob_psf_e2=np.loadtxt(filename_psf_ell).T
    bins_snr,prob_snr=np.loadtxt(filename_snr).T

    n_shears = len(config['shear'])
    n_gals = int(float(config['n_gals_per_file'])) 
    n_pairs = n_gals/2 

    filename_cosmos_catalog = os.path.join(config['input']['real_catalog']['dir'],config['input']['real_catalog']['file_name'])
    cosmos_catalog = pyfits.getdata(filename_cosmos_catalog)
    n_cosmos_gals = len(cosmos_catalog)
    log.info('opened %s with %d images' , filename_cosmos_catalog, n_cosmos_gals)

    for ip in range(config['n_files']):

        filename_cat = 'nbc2_truth.%03d.fits' % ip

        shear_ids = np.random.choice(n_shears,size=n_pairs)
        shear_g1 = np.array(config['shear'])[shear_ids,0]
        shear_g2 = np.array(config['shear'])[shear_ids,1]
        psf_fwhm = np.random.choice(a=bins_psf_fwhm,size=n_pairs,p=prob_psf_fwhm)
        psf_e1 = np.random.choice(a=bins_psf_e,size=n_pairs,p=prob_psf_e1)
        psf_e2 = np.random.choice(a=bins_psf_e,size=n_pairs,p=prob_psf_e2)
        obj_snr = np.random.choice(a=bins_snr,size=n_pairs,p=prob_snr)
        cosmos_ids = np.random.choice(n_cosmos_gals,size=n_pairs)
        rotation_angle = np.random.uniform(low=0,high=2*np.pi,size=n_pairs)

        # create pairs
        shear_ids = np.concatenate([shear_ids[:,None] , shear_ids[:,None]] , axis=1) ; shear_ids = shear_ids.flatten()[:,None]
        shear_g1 = np.concatenate([shear_g1[:,None] , shear_g1[:,None]] , axis=1) ; shear_g1 = shear_g1.flatten()[:,None]
        shear_g2 = np.concatenate([shear_g2[:,None] , shear_g2[:,None]] , axis=1) ; shear_g2 = shear_g2.flatten()[:,None]
        psf_fwhm = np.concatenate([psf_fwhm[:,None] , psf_fwhm[:,None]] , axis=1) ; psf_fwhm = psf_fwhm.flatten()[:,None]
        psf_e1 = np.concatenate([psf_e1[:,None] , psf_e1[:,None]] , axis=1) ; psf_e1 = psf_e1.flatten()[:,None]
        psf_e2 = np.concatenate([psf_e2[:,None] , psf_e2[:,None]] , axis=1) ; psf_e2 = psf_e2.flatten()[:,None]
        obj_snr = np.concatenate([obj_snr[:,None] , obj_snr[:,None]] , axis=1) ; obj_snr = obj_snr.flatten()[:,None]
        cosmos_ids = np.concatenate([cosmos_ids[:,None] , cosmos_ids[:,None]] , axis=1) ; cosmos_ids = cosmos_ids.flatten()[:,None]
        # rotate one of the pair nby 90 deg
        rotation_angle = np.concatenate([rotation_angle[:,None] , rotation_angle[:,None] + np.pi/2.] , axis=1) ; rotation_angle = rotation_angle.flatten()[:,None]
        ids = np.arange(n_gals)[:,None]

        # 'names' : [ 'id' ,  'id_cosmos' , 'id_shear' , 'g1_true' , 'g2_true' ,  'snr' ,  'psf_fwhm' , 'psf_e1' , 'psf_e2' , 'rotation_angle']  ,
        catalog = np.concatenate([ids,cosmos_ids,shear_ids,shear_g1,shear_g2,obj_snr,psf_fwhm,psf_e1,psf_e2,rotation_angle],axis=1)
        catalog = tabletools.array2recarray(catalog,dtype=nbc2_dtypes.dtype_truth)

        tabletools.saveTable(filename_cat,catalog)


def main():


    global log , config , args , cat

    description = 'Get input catalogs for GREAT-DES NBC2 simulation'
    parser = argparse.ArgumentParser(description=description, add_help=True)
    parser.add_argument('-v', '--verbosity', type=int, action='store', default=2, choices=(0, 1, 2, 3 ), help='integer verbosity level: min=0, max=3 [default=2]')
    parser.add_argument('-c', '--filename_config', type=str, action='store', default='nbc2.yaml', help='name of the config file')

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
    # get_noise_level() 16.7552914481
    # get_psf_snr_dist()
    get_catalogs()

    log.info(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))

main()