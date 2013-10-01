import meds
import galsim
import numpy
import galsim.des

DES_PIXEL_SCALE = 0.27
UNIT_PIXEL_SCALE = 1.

# ( u ) = ( dudrow  dudcol ) ( row )
# ( v )   ( dvdrow  dvdcol ) ( col )

class TestGal(object):

    self.Q

    def __init__(self,circular_gso):

        self.gso = circular_gso;
        self.Q = numpy.array([0,0;0,0])


    def apply_transform(self,wcstrans):

        # x'*x --> x'*A'*A*x
        self.Q = numpy.dot(wcstrans.Jxy.T , self.Q )
        self.Q = numpy.dot(self.Q,wcstrans)

    def get_elliptical_object(self,wcstrans):

        Qxx = self.Q[0,0]
        Qxy = self.Q[1,0]
        Qyy = self.Q[1,1]

        g1 = (Qxx - Qyy)/(Qxx+Qyy+2*numpy.sqrt(Qxx*Qyy - Qxy**2))
        g2 = (2*1j*Qxy)/(Qxx+Qyy+2*numpy.sqrt(Qxx*Qyy - Qxy**2))

        gso_copy = self.gso.copy()
        gso_copy.applyShear(g1=g1,g2=g2)

        return gso_copy()


class WCSTransform(object):
    """
    Very simple class which holds a WCS transformation, including a Jacobian and a shifted centroid.
    The (u,v) coordinate plane is the local tangent plane at the location of the galaxy with 
    North = +v, West = +u. The 'u' coordinace scales as cos(DEC) * dRA.
    Here, 'row' is a row in the image; 'col' is a column in the image. 
    The row/col notation is the same as arrays indexed in C and python, arr[row,col].

    Available fields:
    self.dudrow  -  element 1,1 of the Jacobian matrix
    self.dudcol  -  element 1,2 of the Jacobian matrix
    self.dvdrow  -  element 2,1 of the Jacobian matrix
    self.dvdcol  -  element 2,2 of the Jacobian matrix
    self.row0    -  horizontal position of the centroid as given by SExtractor, in pixel coordinates in the SE cutout
    self.col0    -  vertical position of the centroid as given by SExtractor, in pixel coordinates in the SE cutoust
    """



    def __init__(self , stamp_size, **kwargs):

        self.stamp_size = stamp_size

        coords_colrow = ['dudrow', 'dudcol', 'dvdrow', 'dvdcol', 'row0', 'col0']
        coords_xy = ['x0','y0','dudy','dudx','dvdy','dvdx']
        if all([x in coords_colrow for x in kwargs]):

            self.dudrow =      kwargs['dudrow']
            self.dudcol =      kwargs['dudcol']
            self.dvdrow =      kwargs['dvdrow']
            self.dvdcol =      kwargs['dvdcol']
            self.row0 =        kwargs['row0']
            self.col0 =        kwargs['col0']

            self.y0   =  stamp_size[1] - self.col0
            self.x0   =  stamp_size[0] - self.row0
            self.dudx = -self.dudrow
            self.dudy = -self.dudcol
            self.dvdx = -self.dvdrow
            self.dvdy = -self.dvdcol

        elif all([x in coords_xy for x in kwargs]):

            self.x0   = kwargs['x0']
            self.y0   = kwargs['y0']
            self.dudx = kwargs['dudx']
            self.dudy = kwargs['dudy']
            self.dvdx = kwargs['dvdx']
            self.dvdy = kwargs['dvdy']         
            # self.col0   =  stamp_size[1]-self.y0
            # self.row0   =  stamp_size[0]-self.x0    
            self.col0   =  self.y0
            self.row0   =  self.x0    
            self.dudrow = -self.dudx   
            self.dudcol = -self.dudy  
            self.dvdrow = -self.dvdx  
            self.dvdcol = -self.dvdy  


        else:
            raise ValueError('follwing keys are valid : \n' + str(coords_colrow) + '\n' + str(coords_xy) + '\n' + str(kwargs) + '\n' + str([x in coords_xy for x in kwargs]))

        self.Jcolrow = numpy.array(
            [ [self.dudrow , self.dudcol] ,
              [self.dvdrow , self.dvdcol] ] )

        self.Jxy = numpy.array(
            [ [self.dudx , self.dudy] ,
              [self.dvdx , self.dvdy] ] )

        self.invJxy = numpy.linalg.inv(self.Jxy)
        self.invJcolrow = numpy.linalg.inv(self.Jcolrow)


    def check_zero(self,f,eps=1e-15):

        if abs(f) < eps:
            return 0.0
        else:
            return f

    def get_galsim_shift(self):
        """
        Get galsim shift in arcsec, which corresponds to this transformation
        """
        
        # galsim shift 
        dx = self.col0 
        dy = self.row0 

        return dx,dy


    def get_galsim_shears(self):
        # get galsim shear to apply to object which will correspond to the WCS transformation

        # find the major and minor axis of the ellipse
        # U,S,V = numpy.linalg.svd(self.invJcolrow)
        U,S,V = numpy.linalg.svd(self.Jxy)
        VT = V.T;
        print 'Jxy'
        print self.Jxy
        print 'U'
        print U
        print 'S'
        print S
        print 'V'
        print V
        # get shear
        g = (S[0]-S[1])/(S[0]+S[1])
        # get rotation matrix1
        theta1 = numpy.angle(VT[0,0]+1j*VT[1,0]);
        theta2 = numpy.angle(U[0,0]+1j*U[1,0]);
        
        return g,theta1,theta2  
        

    def get_galsim_dilation(self):

        # dilation will correspond to the sqrt of the determinant of inverse transform
        return numpy.sqrt(numpy.abs(numpy.linalg.det(self.Jxy))) 

    def apply_wcs_transformations(self,gsobj):
        # apply shift, shear and dilation to a GSObject

        dx,dy = self.get_galsim_shift()
        g,theta1,theta2 = self.get_galsim_shears()
        dr = self.get_galsim_dilation()
        
        print 'shift : dx=%2.4f dy=  %2.4f' % (dx,dy)
        print 'shear : g =%2.4f theta1=%2.4f (%2.4f deg) theta2=%2.4f (%2.4f deg)' % (g,theta1,theta1*180./numpy.pi,theta2,theta2*180./numpy.pi)
        print 'scale : r =%2.4f' % (dr)

        gsobj.applyRotation(theta=theta1*galsim.radians)
        gsobj.applyShear(g=g,beta=0*galsim.radians)
        gsobj.applyRotation(theta=theta2*galsim.radians)
        gsobj.applyDilation(dr)
        gsobj.applyShift(dx=dx,dy=dy)

        return gsobj

    def show(self):

        print '(x,y) coordinates'
        print 'x0   % 2.4f ' % self.x0
        print 'y0   % 2.4f ' % self.y0
        print '[ dudx dudy     ] = [ % 2.4f % 2.4f ] ' % (self.dudx, self.dudy)
        print '[ dvdx dvdy     ] = [ % 2.4f % 2.4f ] ' % (self.dvdx, self.dvdy)

        print '(col,row) coordinates'
        print 'col0 % 2.4f ' % self.col0
        print 'row0 % 2.4f ' % self.row0
        print '[ dudrow dudcol ] = [ % 2.4f % 2.4f ] ' % (self.dudrow, self.dvdrow)
        print '[ dvdrow dvdcol ] = [ % 2.4f % 2.4f ] ' % (self.dudcol, self.dvdcol)

    def get_invese(self):

        # self.Jxy = numpy.array(
        #     [ [self.dudx , self.dudy] ,
        #       [self.dvdx , self.dvdy] ] )

        x0 = -self.x0
        y0 = -self.y0
        dudx = self.invJxy[0,0] 
        dudy = self.invJxy[0,1]
        dvdx = self.invJxy[1,0] 
        dvdy = self.invJxy[1,1] 
    
        return WCSTransform(self.stamp_size, x0 = x0, y0 = y0, dudx = dudx, dudy = dudy, dvdx = dvdx, dvdy = dvdy)
    

def get_image(gsgal,gspsf=None,n_pix=100,colrow=False):
    """
    Get GalSim image of GSObject.
    If PSF was supplied then use it, otherwise use only pixel kernel.
    @param gsgal GSObject to be Image
    @param gspsf GSObject of the PSF
    @param n_pix pixel size of the image
    @param center if true, then the galaxy pre-shift position will be at the center of postage 
            stamp, if not then use offset = -trueCenter() while drawing, 
    """

    gsp = galsim.GSParams()
    gsp.maximum_fft_size=20000

    img = galsim.ImageD(n_pix,n_pix)
    pix = galsim.Pixel(UNIT_PIXEL_SCALE)
    if gspsf==None:
        gsobj = galsim.Convolve([gsgal,pix],gsparams=gsp)
    else:
        gsobj = galsim.Convolve([gsgal,pix,gspsf],gsparams=gsp)

    if rowcol:
        gsobj.draw(img,dx= UNIT_PIXEL_SCALE,offset=(-img.bounds.trueCenter()))
    else:
        gsobj.draw(img,dx= UNIT_PIXEL_SCALE)

    return img



def get_multiexp_object():
    """
    Save MultiExposureObject to a meds file. 
    First create the coadd with pixel size = 0.27, using the galaxy parameters in the config file.
    The image of the coadd will look the same as the DES coadd, in DS9.
    Then produce the SE images using the WCS transforms supplied in the config file.
    In DS9 they will look the same as DES SE images.
    """

    # lists init

    truth = config['truth']
    transforms = config['transforms']

    list_images = []
    list_wcstrans = []
    list_gso = []

    gal_g1 = truth['gal_g1']  
    gal_g2 = truth['gal_g2']  
    gal_r = truth['gal_r'] 
    gal_x = truth['gal_x']
    gal_y = truth['gal_y']
    n_pix = truth['n_pix']  
    
    # first the coadd in sky coordinates - no PSF on this one

    gal_gs = galsim.Exponential(half_light_radius = gal_r)
    gal_gs.applyShear(g1=gal_g1,g2=gal_g2)
    gal_gs.applyShift(dx=gal_x,dy=gal_y)
    img_coadd = get_image(gal_gs,n_pix=n_pix)
    img_coadd.write('coadd.fits')

    wcstr = WCSTransform(stamp_size = (n_pix,n_pix), dudy=1, dvdx=1, dudx=0.0, dvdy=0.0, x0=0, y0=0) # here using arcsec for row0, col0#

    list_images += [img_coadd]
    list_wcstrans += [wcstr]
    list_gso += [gal_gs]

    psf = galsim.Moffat(beta=truth['psf_beta'],fwhm=truth['psf_fwhm']) 

    for itr,tr in enumerate(transforms):
        
        # then the exposures 

        # exposure - change scale

        print 'getting exposure %d ---------------------------- ' % itr

        if 'dudy' in tr:
            wcstr = WCSTransform(stamp_size = (n_pix,n_pix), dudy=tr['dudy'], dvdx=tr['dvdx'], dudx=tr['dudx'], dvdy=tr['dvdy'], x0=tr['x0'], y0=tr['y0']) 
        else:
            wcstr = WCSTransform(stamp_size = (n_pix,n_pix), dudrow=tr['dudrow'], dvdcol=tr['dvdcol'], dudcol=tr['dudcol'], dvdrow=tr['dvdrow'], row0=tr['row0'], col0=tr['col0'])
        
        logger.info('------------------------------------- fwd -------------------------------------')
        wcstr.show()

        gal_gs_copy = gal_gs.copy()
        gal_gs_copy = wcstr.apply_wcs_transformations(gal_gs_copy)
        img = get_image(gal_gs_copy,n_pix=n_pix)
        filename_fits = 'SE-%03d.fits' % itr
        img.write(filename_fits)

        list_wcstrans += [wcstr]
        list_images += [img]
        list_gso += [gal_gs_copy]

        inv_wcs = wcstr.get_invese()
        logger.info('------------------------------------- inv -------------------------------------')
        inv_wcs.show()
        svd_use='u'
        inv_gso = inv_wcs.apply_wcs_transformations(gal_gs_copy.copy())
        inv_img = get_image(inv_gso,n_pix=n_pix)

        # numpy.testing.assert_array_almost_equal(img.array,coadd.array)
        import pylab
        pylab.subplot(1,4,1)
        pylab.imshow(img_coadd.array)
        pylab.colorbar()
        pylab.subplot(1,4,2)
        pylab.imshow(inv_img.array)
        pylab.colorbar()
        pylab.subplot(1,4,3)
        pylab.imshow(inv_img.array - img_coadd.array)
        pylab.colorbar()
        pylab.subplot(1,4,4)
        pylab.imshow(img.array)
        pylab.show()




    # create MultiExposureObject

    meo = galsim.des.MultiExposureObject(list_images, wcstrans=list_wcstrans)

    galsim.des.write_meds(args.filename_output,[meo])
    logger.info('saved %s ' % args.filename_output)
    return list_images,list_wcstrans,list_gso


def save_psf_images():

    # this script gets the PSF for shear test 1
    # two resolutons : same as galaxy, saved in psf.fits, upsampled and padded, saved in psf.hires.fits

    import galsim
    import yaml
    import numpy

    # PSF at resolution of the galaxy
    logger.info('getting single PSF at the pixel scale of a galaxy')

    n_pix = config['truth']['n_pix']
    pixel_scale = UNIT_PIXEL_SCALE

    psf_type = galsim.Moffat
    psf_fwhm = config['truth']['psf_fwhm']
    psf_beta = config['truth']['psf_beta']

    psf = psf_type(fwhm=psf_fwhm,beta=psf_beta)
    pix = galsim.Pixel(xw=pixel_scale)
    img = galsim.ImageD(n_pix,n_pix)
    final = galsim.Convolve([psf,pix])
    dx = (numpy.random.random() - 0.5)*0.5*pixel_scale
    dy = (numpy.random.random() - 0.5)*0.5*pixel_scale
    final.applyShift(dx=dx,dy=dy)
    final.draw(img,dx=pixel_scale)

    filename_psf = 'psf.single.%dx%d.fits' % (config['truth']['n_pix'],config['truth']['n_pix'])
    img.write(filename_psf)
    print 'saved %s' % filename_psf

    logger.info('getting single PSF at high resolution')
    # now the hires PSF, centered in the middle

    n_sub = 5
    n_pad = 4
    n_pix_hires = (config['truth']['n_pix'] + n_pad) * n_sub
    pixel_scale_hires = float(UNIT_PIXEL_SCALE) / float(n_sub)

    psf = psf_type(fwhm=psf_fwhm,beta=psf_beta)
    pix = galsim.Pixel(xw=pixel_scale)
    img = galsim.ImageD(n_pix_hires,n_pix_hires)
    dx = 0.5*pixel_scale
    dy = 0.5*pixel_scale
    final = galsim.Convolve([psf,pix])
    final.applyShift(dx=dx,dy=dy)
    final.draw(img,dx=pixel_scale_hires)

    filename_psf_hires = 'psf.single.%dx%d.fits' % (n_pix_hires,n_pix_hires)
    img.write(filename_psf_hires)
    print 'saved %s' % filename_psf_hires

    logger.info('getting low res PSF in a field')

    # now create the PSF - fields for the shapelet pipeline
    ny_tiles = nx_tiles = 30

    psf_image = galsim.ImageF(n_pix * nx_tiles-1 , n_pix * ny_tiles-1)
    psf_image.setScale(pixel_scale)
    random_seed=123123

    shift_radius = pixel_scale *0.5     
    shift_radius_sq = shift_radius**2

    k=0;
    for iy in range(ny_tiles):
        for ix in range(nx_tiles):
            ud = galsim.UniformDeviate(random_seed+k)

            # Apply a random shift_radius:
            rsq = 2 * shift_radius_sq
            while (rsq > shift_radius_sq):
                dx = (2*ud()-1) * shift_radius
                dy = (2*ud()-1) * shift_radius
                rsq = dx**2 + dy**2

            this_psf = final.createShifted(dx,dy)

            b = galsim.BoundsI(ix*n_pix+1 , (ix+1)*n_pix-1, iy*n_pix+1 , (iy+1)*n_pix-1)
            sub_psf_image = psf_image[b]
            this_psf.draw(sub_psf_image,dx=pixel_scale)
            k+=1;

    # save the tiled PSF image
    filename_psf_field = 'psf.field.%dx%d.fits' % (config['truth']['n_pix'],config['truth']['n_pix'])
    psf_image.write(filename_psf_field)
    print 'saved %s' % filename_psf_field

    
    

def test_inverse(image,gso,wcstrans):
    """
    @brief Check if we can reconstruct coadd using inverse transforms from the single exposures.
    @param list_gso list of gsobjects, first is coadd, and rest is SE.
    @param list_wcstrans list of transformations that produced SE from coadd. First corresponds to coadd, so it doesn't count.
    """

    logger.info('testing inverse')

    coadd = list_images[0]

    n_pix = list_images[0].bounds.xmax

    list_se_wcstrans = list_wcstrans[1:]
    list_se_gso = list_gso[1:]
    list_se_images = list_images[1:]

    for itr,tr in enumerate(list_se_wcstrans):

        inv_transform = tr.get_invese()
        inv_transform.show()

        gal_gs_copy = list_se_gso[itr].copy()
        gal_gs_copy = inv_transform.apply_wcs_transformations(gal_gs_copy)
        img = get_image(gal_gs_copy,n_pix=n_pix)
        # filename_fits = 'SE-%03d-inv.fits' % itr 
        # img.write(filename_fits)

        # numpy.testing.assert_array_almost_equal(img.array,coadd.array)
        import pylab
        pylab.subplot(1,3,1)
        pylab.imshow(img.array)
        pylab.colorbar()
        pylab.subplot(1,3,2)
        pylab.imshow(coadd.array)
        pylab.colorbar()
        pylab.subplot(1,3,3)
        pylab.imshow(img.array - coadd.array)
        pylab.colorbar()
        pylab.show()



####################################################################################################

def main():

    import yaml
    import argparse
    import logging
    import sys

    global logger , config , args

    description = 'Save MultiExposureObject to a meds file.  \
    First create the coadd with pixel size = 0.27, using the galaxy parameters in the config file. \
    The image of the coadd will look the same as the DES coadd, in DS9. \
    Then produce the SE images using the WCS transforms supplied in the config file. \
    In DS9 they will look the same as DES SE images. '
    parser = argparse.ArgumentParser(description=description, add_help=True)
    parser.add_argument('-v', '--verbosity', type=int, action='store', default=2, choices=(0, 1, 2, 3 ), help='integer verbosity level: min=0, max=3 [default=2]')
    parser.add_argument('-o', '--filename_output', default='wcs_test_card_01.meds.fits',type=str, action='store', help='name of the output meds file')
    parser.add_argument('-c', '--filename_config', default='transforms_list.yaml',type=str, action='store', help='name of the yaml config file')
    
    args = parser.parse_args()
    # Parse the integer verbosity level from the command line args into a logging_level string
    logging_levels = { 0: logging.CRITICAL, 
                       1: logging.WARNING,
                       2: logging.INFO,
                       3: logging.DEBUG }
    logging_level = logging_levels[args.verbosity]
    logging.basicConfig(format="%(message)s", level=logging_level, stream=sys.stdout)
    logger = logging.getLogger("WCS_TEST") 
    logger.setLevel(logging_level)


    config = yaml.load(open(args.filename_config))

    list_images,list_wcstrans,list_gso = get_multiexp_object()
    # save_psf_images()


main()