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

(2) How to generate WCS test card
	python transforms_test.py -c transforms_list.yaml
	
(3) How to generate SHA1

	for debug mode with 1000 galaxies per file
	python generate_sha1.py -c sha1.yaml --debug	

	testing:
	rm *fits*
	python $CODE/GREAT-DES/sha1/generate_sha1.py -c sha1-test.yaml --debug
	python $CODE/GREAT-DES/sha1/run_hsm_sha1.py -i sha1-test-O1.cat
	python $CODE/GREAT-DES/sha1/run_hsm_sha1.py -i sha1-test-O2.cat
	python $CODE/GREAT-DES/sha1/plots_sha1.py -c sha1-test.yaml -i sha1-test-O1.cat
	python $CODE/GREAT-DES/sha1/plots_sha1.py -c sha1-test.yaml -i sha1-test-O2.cat


	now estimate stde by running only 1000 galaxies (use debug flag)
	python generate_sha1.py -c sha1.yaml --debug
	python run_hsm_sha1.py -i sha1-O1.cat
	python run_hsm_sha1.py -i sha1-O2.cat
	python plots_sha1.py -c sha1.yaml -i sha1-O1.cat
	python plots_sha1.py -c sha1.yaml -i sha1-O2.cat

	save the obtained stats files, so that we can use them to get the stde
	cp sha1-O1.cat.hsm.stats.cat sha1-O1.cat.hsm.stats.stde.cat
	cp sha1-O2.cat.hsm.stats.cat sha1-O2.cat.hsm.stats.stde.cat

	now create the main sample
	python generate_sha1.py -c sha1-test.yaml 
