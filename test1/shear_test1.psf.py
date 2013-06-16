# this script gets the PSF for shear test 1
# two resolutons : same as galaxy, saved in psf.fits, upsampled and padded, saved in psf.hires.fits




import galsim
import yaml
import numpy

filename_yaml = 'shear_test1.meds.yaml'

config = yaml.load(open(filename_yaml))

# PSF at resolution of the galaxy

n_pix = config['image']['size']
pixel_scale = config['image']['pixel_scale']

psf_type = eval('galsim.%s' % config['psf']['type'])
psf_fwhm = config['psf']['fwhm']

psf = psf_type(fwhm=psf_fwhm)
pix = galsim.Pixel(xw=pixel_scale)
img = galsim.ImageD(n_pix,n_pix)
final = galsim.Convolve([psf,pix])
dx = (numpy.random.random() - 0.5)*0.5*pixel_scale
dy = (numpy.random.random() - 0.5)*0.5*pixel_scale
final.applyShift(dx=dx,dy=dy)
final.draw(img,dx=pixel_scale)

filename_psf = 'psf.fits'
img.write(filename_psf)
print 'saved %s' % filename_psf

# now the hires PSF, centered in the middle

n_sub = 5
n_pad = 4
n_pix_hires = (config['image']['size'] + n_pad) * n_sub
pixel_scale_hires = float(config['image']['pixel_scale']) / float(n_sub)

psf = psf_type(fwhm=psf_fwhm)
pix = galsim.Pixel(xw=pixel_scale)
img = galsim.ImageD(n_pix_hires,n_pix_hires)
final = galsim.Convolve([psf,pix])
final.draw(img,dx=pixel_scale_hires)

filename_psf_hires = 'psf.hires.fits'
img.write(filename_psf_hires)
print 'saved %s' % filename_psf_hires

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
filename_psf_field = 'psf.field.fits'
psf_image.write(filename_psf_field)
print 'saved %s' % filename_psf_field


