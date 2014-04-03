dtype_truth = { 
        'names' : [ 'id' ,  'id_cosmos' , 'id_shear' , 'id_psf' ,  'g1_true' , 'g2_true' ,  'flux' , 'snr' ,  'psf_fwhm' , 'psf_e1' , 'psf_e2' , 'rotation_angle']  ,
        'formats' : ['i8']*4 + ['f4']*8 }

dtype_psfkey = {
        'names' : [ 'id_psf' ,  'id_psf_fwhm' , 'id_psf_e1' , 'id_psf_e2' ,  'psf_fwhm' , 'psf_e1' , 'psf_e2' ]  ,
        'formats' : ['i8']*4 + ['f4']*3 }

