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

