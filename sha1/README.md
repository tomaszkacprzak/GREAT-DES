SHA1 test 

(1) Goals of the test:
	Test the recovery of elliptical isophote galaxy from multiple exposure data, transformed by 
	a coordinate system transformation.

(2) Current status:
	First test implemented.
	Test transforms has noiseless images of Exponential galaxies (Sersic index 1).
	True parameters are in file transforms_list.yaml, including true shear.
	MEDS file corresponding to this file is 'wcs_test_card_01.meds.fits'.
	First image in the MEDS file corresponds to the galaxy in (u,v) coordinates, no PSF.
	Rest of images is SE data, convolved with circullar Moffat PSF.
	Images of PSF are available:
	psf.single.128x128.fits - single PSF at native resolution, 128x128 pixels
	psf.single.660x660.fits - single PSF at high resolution, with upsampling 5 and padding pixels (mainly for im3shape)
	psf.field.128x128.fits  - field of PSFs, dithered

(2) How to generate test sha1
	python transforms_test.py -c transforms_list.yaml
