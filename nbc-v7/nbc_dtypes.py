dtype_truth = { 
        'names' : [ 'id' ,  'id_cosmos' , 'id_shear' , 'id_psf' ,  'g1_true' , 'g2_true' ,  'flux' , 'snr' ,  'fwhm', 'psf_fwhm' , 'psf_e1' , 'psf_e2' , 'rotation_angle'
                        ,'hsm_obs_g1','hsm_obs_g2','hsm_cor_g1','hsm_cor_g2','hsm_obs_sigma','hsm_cor_sigma','hsm_centroid_x','hsm_centroid_y','hsm_mom_amp' ,
                         'psf_fwhm_measured' , 'cosmos_mag_auto'  , 'cosmos_flux_radius' 
                           , 'sf_i'         
                           , 'sf_hlr'       
                           , 'sf_sersicn'   
                           , 'sf_q'         
                           , 'sf_boxiness'  
                           , 'sf_phi'
                           , 'zphot'
                           , 'psf_file' ] ,    
        'formats' : ['i8']*4 + ['f8']*28 + ['a2048'] }

dtype_psfkey = {
        'names' : [ 'id_psf' ,  'id_psf_fwhm' , 'id_psf_e1' , 'id_psf_e2' ,  'psf_fwhm' , 'psf_e1' , 'psf_e2' ]  ,
        'formats' : ['i8']*4 + ['f4']*3 }

dtype_stats = { 'names' : ['ig' , 'tru_g1' , 'tru_g2' ,'n_gals','n_fail', 'est_g1','est_g2','est_size','est_stdv_g1','est_stdv_g2','est_stdm_g1','est_stdm_g2','est_stdv_size','est_stdm_size'],
                'formats' : ['i4']*1 + ['f8']*2 + ['i4']*2 + ['f8']*9}

dtype_bias = { 'names'   : ['n_fail', 'm1', 'm2', 'c1', 'c2', 'm1_std', 'm2_std', 'c1_std', 'c2_std', 'g1_std' , 'g2_std' ],
                     'formats' : ['i8']*1 + ['f8']*10 } 

