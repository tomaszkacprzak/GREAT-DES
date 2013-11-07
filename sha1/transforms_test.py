import meds 
import galsim
import numpy
import galsim.des
import copy
import pyfits

DES_PIXEL_SCALE = 0.27
UNIT_PIXEL_SCALE = 1.

# ( u ) = ( dudcol dudrow ) ( col )
# ( v )   ( dvdcol dvdrow ) ( row )


class TestGal(object):
    """
    Class TestGal is a container for a circular GSObject and its covariance_matrix.

    It can be initialised:
    gal_gso = galsim.Exponential(half_light_radius = gal_r)
    test_gal = TestGal(gal_gso)
    test_gal.set_shear(g1=gal_g1,g2=gal_g2)
    test_gal.set_shift(dx=gal_x,dy=gal_y)
    img_coadd = test_gal.get_image(n_pix=n_pix)

    Then apply coordinate transformation:
    test_gal_copy.apply_transform(wcstr)
    img_se = test_gal_copy.get_image(n_pix=n_pix)

    Show current state:
    test_gal_copy_inv.show()
    """

    def __init__(self,circular_gso):
        """
        @brief Object constructor
        @param circular_gso GalSim GSObject, circular, with no shear applied
        """

        # use circular_gso
        self.gso = circular_gso;
        # use identity covariance matrix
        self.Q = numpy.array([[1,0],[0,1]],dtype=numpy.float64)
        # use origin of coordinate system
        self.xy = numpy.array([0,0],dtype=numpy.float64)

    def show(self):
        """
        @brief print covariance and center
        """

        print 'C  = [% 2.4f % 2.4f]' % (self.Q[0,0],self.Q[0,1])
        print '     [% 2.4f % 2.4f]' % (self.Q[1,0],self.Q[1,1])
        print 'xy = [% 2.4f % 2.4f]' % (self.xy[0],self.xy[1])
        
    def apply_transform(self,wcstrans):
        """
        @brief Apply WCS transform to the object
        @param wcstrans galsim.des.WCSTransform obejct to apply, going from (col,row) -> (u,v) coords
                in the transform, (col0,row0) is the position in (col,row) of the current (0,0) in (u,v)
        """

        # get the transformation matrix
        J = numpy.array(
            [ [wcstrans.dudcol , wcstrans.dudrow] ,
              [wcstrans.dvdcol , wcstrans.dvdrow] ] )

        # inverse to go from (u,v) -> (col,row)
        invJ = numpy.linalg.inv(J)

        # initialise new Q
        Q = self.Q

        # x'*x --> x'*A'*A*x
        Q = numpy.dot(invJ.T , Q )
        Q = numpy.dot(Q , invJ)
        # 
        self.xy[0] = wcstrans.col0
        self.xy[1] = wcstrans.row0

        logger.info('applying transform')
        logger.info('[% 2.4f,% 2.4f] --> [% 2.4f % 2.4f]'%( self.Q[0,0],self.Q[0,1], Q[0,0],Q[0,1]))
        logger.info('[% 2.4f,% 2.4f] --> [% 2.4f % 2.4f]'%( self.Q[1,0],self.Q[1,1], Q[1,0],Q[1,1]))
        self.Q = Q

    def get_g_from_Q(self,Q):
        """
        @brief calculate galsim shear from the covariance matrix Q
        @params Q covariance matrix
        """

        # unpack
        Qxx = Q[0,0]
        Qxy = Q[1,0]
        Qyy = Q[1,1]

        # apply equations
        div = Qxx+Qyy+2.*numpy.sqrt(Qxx*Qyy - Qxy**2)
        g1 = (Qxx - Qyy)/div
        g2 = (2*Qxy)/div

        return g1,g2

    def get_Q_from_g(self,e1,e2):
        """
        @brief calculate covariance matrix of elliptical object with elliticity e1,e2
        @param e1 elliticity e1
        @param e2 elliticity e2
        """

        # get precision matrix
        e_sq = (e1*e1+e2*e2);
        p1 = (1.+e_sq-2.*e1)/(1.-e_sq)
        p2 = (1.+e_sq+2.*e1)/(1.-e_sq)
        p3 = -2.*e2/(1.-e_sq)
        precision_matrix = numpy.array([[p1,p3],[p3,p2]])

        # calculate covariance matrix
        covariance_matrix = numpy.linalg.inv(precision_matrix)
        m1 = covariance_matrix[0,0]
        m2 = covariance_matrix[1,1]
        m3 = covariance_matrix[1,0]

        # do a check of consistency, go back to g1,g2        
        g = self.get_g_from_Q(covariance_matrix)

        print 'Q from g'
        print 'g = [% 2.4f % 2.4f]' % (e1,e2)
        
        print 'C = [% 2.4f % 2.4f]' % (m1,m3)
        print '    [% 2.4f % 2.4f]' % (m3,m2)
        print '      det = %2.4f'   % (numpy.linalg.det(covariance_matrix))

        print 'P = [% 2.4f % 2.4f]' % (p1,p3)
        print '    [% 2.4f % 2.4f]' % (p3,p2)

        print 'g from Q - consistency check'
        print 'g = [% 2.4f % 2.4f]' % (g[0],g[1])

        return covariance_matrix


    def get_r_from_Q(self,Q):
        """
        @brief Get the galsim dilation from covariance matrix
        @param Q covariance matrix
        """

        # get the determinant
        detQ = numpy.linalg.det(Q)
        # dilation is a 1/4 power of det
        r = detQ**(0.25)
        print 'det(Q) = [%2.4f]' % (detQ)
        print 'r      = [%2.4f]' % (r)
        return r

    def get_gso(self):
        """
        @brief get a galsim GSObject which has current covariance and shift
        @return new GSObject 
        """

        # get a copy of GSObject
        gso_copy = self.gso.copy()

        # get galsim transformations
        g1,g2=self.get_g_from_Q(self.Q)
        r = self.get_r_from_Q(self.Q)

        # apply transformations
        gso_copy.applyShear(g1=g1,g2=g2)
        gso_copy.applyDilation(scale=r)
        gso_copy.applyShift(dx=self.xy[0],dy=self.xy[1])

        logger.info('applying shear    (% 2.4f,% 2.4f)' % (g1,g2))
        logger.info('applying dilation (% 2.4f)' % r)
        logger.info('applying shift    (% 2.4f,% 2.4f)' % (self.xy[0],self.xy[1]))

        return gso_copy

    def set_shear(self,g1,g2):
        """
        @brief set shear of the object to g1,g2
        @param g1 shear to apply
        @param g2 shear to apply
        """

        # get covariance
        self.Q = self.get_Q_from_g(g1,g2)
        

    def set_shift(self,dx,dy):
        """
        @brief set shift of the object to dx,dy
        @param dx shift to apply
        @param dy shift to apply
        """

        self.xy[0]=dx
        self.xy[1]=dy

    def set_dilation(self,scale):
        """
        @brief set dilation of the object by setting current covariance unit determinant 
        and multiplying it by scale
        @param scale to set
        """

        self.Q /= get_r_from_Q()**2
        self.Q *= scale

    def copy(self):
        """
        @brief copy operator
        @return copied objects
        """

        clone = TestGal(self.gso)
        clone.Q = self.Q.copy()
        clone.xy = self.xy.copy()
        
        return clone

    
    def get_image(self,gspsf=None,n_pix=100,coordinate_origin='center'):
        """
        Get GalSim image of test gal.
        If PSF was supplied then use it, otherwise use only pixel kernel.
        @param gspsf GSObject of the PSF
        @param coordinate_origin {center,bottom left}, if center then the galaxy pre-shift position 
            will be at the center of postage stamp, if 'bottom left' then use offset = -trueCenter()
            while drawing, 
        @return galsim.Image of the object
        """

        # get a new GSO
        gsgal = self.get_gso()

        # get GSParams to allow large images
        gsp = galsim.GSParams()
        gsp.maximum_fft_size=20000

        # get image and pixel
        img = galsim.ImageD(n_pix,n_pix)
        pix = galsim.Pixel(UNIT_PIXEL_SCALE)

        # convolve with PSF
        if gspsf==None:
            gsobj = galsim.Convolve([gsgal,pix],gsparams=gsp)
        else:
            gsobj = galsim.Convolve([gsgal,pix,gspsf],gsparams=gsp)

        # set the center
        if coordinate_origin == 'center':
            gsobj.draw(img,dx= UNIT_PIXEL_SCALE)
        elif coordinate_origin == 'bottom left':
            # offset = -img.bounds.trueCenter() + (-1,-1)
            x0 = - float(n_pix)/2.
            y0 = - float(n_pix)/2.
            # gsobj.draw(img,dx= UNIT_PIXEL_SCALE,offset=(-img.bounds.trueCenter()))
            gsobj.draw(img,dx= UNIT_PIXEL_SCALE,offset=galsim.PositionD(x0,y0))
            # gsobj.draw(img,dx= UNIT_PIXEL_SCALE)
            # print img.bounds.trueCenter()
        else:
            raise ValueError('coordinate_origin should be {center,bottom left}')

        return img
  
def show_wcstr(self):
    """
    @brief Print galsim.des.wcstransform to screen
    @param self wcstransform to print
    """

    print 'col0 % 2.4f ' % self.col0
    print 'row0 % 2.4f ' % self.row0
    print '[ dudcol dudrow  ] = [ % 2.4f % 2.4f ] ' % (self.dudcol,self.dudrow)
    print '[ dvdcol dvdrow  ] = [ % 2.4f % 2.4f ] ' % (self.dvdcol,self.dvdrow)

def get_invese_wcstr(wcstrans,col0,row0):
    """
    @brief get get_invese_wcstr transform of a galsim.des.wcstransform
    @param wcstrans transform to invert
    @param col0 x position in new coordinate system
    @param row0 y position in new coordinate system
    @return transformed wcstransform
    """

    # get covariance 
    J = numpy.array(
        [ [ wcstrans.dudcol, wcstrans.dudrow ] ,
          [ wcstrans.dvdcol, wcstrans.dvdrow ] ] )

    # invert covariance
    invJ = numpy.linalg.inv(J)

    # assign
    row0 = row0
    col0 = col0
    dudcol = invJ[0,0]
    dudrow = invJ[0,1] 
    dvdcol = invJ[1,0] 
    dvdrow = invJ[1,1] 

    return galsim.des.WCSTransform(dudcol=dudcol, dvdrow=dvdrow, dudrow=dudrow, dvdcol=dvdcol, row0=row0, col0=col0)
   

def se_to_coadd(n_pix,wcstr):
    """
    @brief transform to image coordinate system as seen in the coadd    
                    (see des_image_coords_convention.png)
    @param n_pix number of pixel in the image
    @param wcstrans wstransform from (u,v) -> (col0,row0) in SE
    """

    row0 = float(n_pix) - wcstr.col0
    col0 = float(n_pix) - wcstr.row0

    # row0 = wcstr.row0
    # col0 = wcstr.col0

    dudcol = -wcstr.dvdcol
    dudrow = -wcstr.dvdrow
    dvdcol = -wcstr.dudcol
    dvdrow = -wcstr.dudrow

    return galsim.des.WCSTransform(dudcol=dudcol, dvdrow=dvdrow, dudrow=dudrow, dvdcol=dvdcol, row0=row0, col0=col0)


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

    list_images   = []
    list_wcstrans = []
    list_test_gso = []
    list_test_gal = []

    # get the initial parameters
    gal_g1 = truth['gal_g1']  
    gal_g2 = truth['gal_g2']  
    gal_r  = truth['gal_r'] 
    gal_x  = truth['gal_x']
    gal_y  = truth['gal_y']
    n_pix  = truth['n_pix']  
    
    # first the coadd in sky coordinates - no PSF on this one

    gal_gso = galsim.Exponential(half_light_radius = gal_r)
    test_gal = TestGal(gal_gso)
    test_gal.set_shear(g1=gal_g1,g2=gal_g2)
    test_gal.set_shift(dx=gal_x,dy=gal_y)
    flux = test_gal.gso.getFlux()
    test_gal.gso.scaleFlux(flux)

    # get the image
    img_uv = test_gal.get_image(n_pix=n_pix)
    img_uv.write('coadd.fits')

    # dummpy WCSTrans for coadd
    wcstr = galsim.des.WCSTransform(dudcol=1, dvdrow=1, dudrow=0.0, dvdcol=0.0, row0=0, col0=0) # here using arcsec for row0, col0#

    list_images   += [img_uv]
    list_wcstrans += [wcstr]
    list_test_gal += [test_gal]

    # get a simle circular PSF
    gso_psf = galsim.Moffat(beta=truth['psf_beta'],fwhm=truth['psf_fwhm']) 

    # loop over exposures
    for itr,tr in enumerate(transforms):
        
        # then the exposures 

        # exposure - change scale

        print 'getting exposure %d ---------------------------- ' % itr

        # get the WCSTransform
        wcstr = galsim.des.WCSTransform(dudrow=tr['dudrow'], dvdcol=tr['dvdcol'], dudcol=tr['dudcol'], dvdrow=tr['dvdrow'], row0=tr['row0'], col0=tr['col0'])
        
        logger.info('------------------------------------- fwd -------------------------------------')
       
        logger.info('---------------------- reg ----------------------')
        # get the exposure - copy the coadd galaxy
        test_gal_copy = test_gal.copy()
        # apply the transform to the copy
        test_gal_copy.apply_transform(wcstr)
        # get image
        img = test_gal_copy.get_image(n_pix=n_pix,coordinate_origin='bottom left',gspsf = gso_psf)
        # show the WCSTransform
        show_wcstr(wcstr)
        # save the image to fits file
        filename_fits = 'SE-%03d.fits' % itr
        img.write(filename_fits)

        # add to lists
        list_wcstrans += [wcstr]
        list_images   += [img]
        list_test_gal += [test_gal_copy]

        # coaddvew - get image of how they would like in the coadd coordinate system (see des_image_coords_convention.png)
        logger.info('---------------------- coaddview, no PSF ----------------------')
        # copy test_gal
        test_gal_copy_coaddview = test_gal.copy()
        # make a transformation from (u,v) to coadd (x,y)
        wcstr_coaddview = se_to_coadd(n_pix,wcstr)
        # print the transform we got
        show_wcstr(wcstr_coaddview)
        # apply transform
        test_gal_copy_coaddview.apply_transform(wcstr_coaddview)
        # get image
        img = test_gal_copy_coaddview.get_image(n_pix=n_pix,coordinate_origin='bottom left')
        # save file
        filename_fits = 'SE-coaddview-%03d.fits' % itr
        img.write(filename_fits)

        logger.info('------------------------------------- inv -------------------------------------')
        # now get the inverse transform, to get back to (u,v) -- sanity check
        inv_wcs = get_invese_wcstr( wcstr , gal_x , gal_y )
        # print transform
        show_wcstr(inv_wcs)
        # copy test gal
        test_gal_copy_inv = test_gal_copy.copy()

        logger.info('orig Q')
        test_gal.show()
        logger.info('fwd Q')
        test_gal_copy_inv.show()
        logger.info('inv Q')

        # apply inverse transform
        test_gal_copy_inv.apply_transform(inv_wcs)
        # see what we got
        test_gal_copy_inv.show()

        # get inverse image
        inv_img = test_gal_copy_inv.get_image(n_pix=n_pix,coordinate_origin='center')
        # check if uv and inverse are equal
        numpy.testing.assert_array_almost_equal(inv_img.array,img_uv.array)

        # diagnose errors
        # import pylab
        # pylab.subplot(1,4,1)
        # pylab.imshow(img_uv.array)
        # pylab.colorbar()
        # pylab.subplot(1,4,2)
        # pylab.imshow(img.array)
        # pylab.colorbar()
        # pylab.subplot(1,4,3)
        # pylab.imshow(inv_img.array)
        # pylab.colorbar()
        # pylab.subplot(1,4,4)
        # pylab.imshow(inv_img.array - img_uv.array)
        # pylab.colorbar()
        # pylab.show()



    # create MultiExposureObject and save to meds

    meo = galsim.des.MultiExposureObject(list_images, wcstrans=list_wcstrans)
    galsim.des.write_meds(args.filename_output,[meo])
    logger.info('saved %s ' % args.filename_output)
    return list_images,list_wcstrans,list_test_gal

def consistency_check():

    m=meds.MEDS(args.filename_output)
    logger.info('opened %s with %d exposures' % (args.filename_output,m.get_cat()['ncutout'][0]))

    import pylab as pl

    for ie in range(1,m.get_cat()['ncutout'][0]):

        filename_fits = 'SE-%03d.fits' % (ie-1)

        img_fits = pyfits.getdata(filename_fits)
        img_fits2 = galsim.fits.read(filename_fits)
        img_meds = m.get_cutout(0,ie)

        logger.info('exposure %d filename %s' % (ie, filename_fits))
        numpy.testing.assert_array_almost_equal(img_fits,img_meds)

    logger.info('all asserts successful')


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

    get_multiexp_object()
    # save_psf_images()
    consistency_check()


main()