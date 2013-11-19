SHA1 test 

(1) Goals of the test:

	There are three goals of the test:
	1.1	 Test the recovery of elliptical isophote galaxy from multiple exposure data, transformed by 
	   	 a coordinate system transformation to test the WCS convention.
	1.2	 Test the limits of SNR and galaxy size where the code performs 
	1.3	 For this limits, check the dependence of biases as a function of various parameters

(2) Current status:

	1.1 Implemented
	1.2 Implemented
	1.3 In progress

(3) Details of tests

	1.1 WCS transforms test

	Test transforms has noiseless images of Exponential galaxies (Sersic index 1).
	True parameters are in file transforms_list.yaml, including true shear.
	MEDS file corresponding to this file is 'wcs_test_card_01.meds.fits'.
	First image in the MEDS file corresponds to the galaxy in (u,v) coordinates, no PSF.
	Rest of images is SE data, convolved with circullar Moffat PSF.
	Images of PSF are available:
	psf.single.128x128.fits - single PSF at native resolution, 128x128 pixels
	psf.single.660x660.fits - single PSF at high resolution, with upsampling 5 and padding pixels (mainly for im3shape)
	psf.field.128x128.fits  - field of PSFs, dithered

	To create the test, run
	transforms_test.py -c transforms_list.yaml

	1.2 SNR vs size cut

	Produce galaxies with various SNR and size.
	Use command:
	python generate_sha1.py -c sha1-O1.yaml --o1
	Program saves a set of meds files, one meds file per snr-size combination.
	Additionally saves catalog sha1-O1.cat, which is a truth table with a row for each file.
	Run code on the data produced, for example KSB:
	python run_hsm_sha1.py -i sha1-O1.cat
	This saves a catalog with results for each meds file.
	Then make plots given all the results are present in current dir:
	python plots_sha1.py -c sha1-O1.yaml -i sha1-O1.cat --o1

	Commands I used to generate the data:

	testing:
	rm *fits*
	python $CODE/GREAT-DES/sha1/generate_sha1.py -c sha1-test.yaml --debug
	python $CODE/GREAT-DES/sha1/run_hsm_sha1.py -i sha1-test-O1.cat
	python $CODE/GREAT-DES/sha1/plots_sha1.py -c sha1-test.yaml -i sha1-test-O1.cat
	
	now estimate stde by running only 100 galaxies (use debug flag)
	rm *fz 
	python $CODE/GREAT-DES/sha1/generate_sha1.py -c sha1-O1.yaml --debug --o1
	python $CODE/GREAT-DES/sha1/run_hsm_sha1.py -i sha1-O1.cat
	python $CODE/GREAT-DES/sha1/plots_sha1.py -c sha1-O1.yaml -i sha1-O1.cat
	
	save the obtained stats files, so that we can use them to get the stde
	cp sha1-O1.cat.hsm.stats.cat sha1-O1.cat.hsm.stats.stde.cat
	
	now create the main sample
	python $CODE/GREAT-DES/sha1/generate_sha1.py -c sha1-O1.yaml --o1
	python $CODE/GREAT-DES/sha1/run_hsm_sha1.py -i sha1-O1.cat
	python $CODE/GREAT-DES/sha1/plots_sha1.py -c sha1-O1.yaml -i sha1-O1.cat

(4) PSF information

	PSF is constant for every galaxy in a meds file.
	For every meds file, there are three PSF files: 
	meds_filename.psf.single.64x64.fits.fz    -- single PSF postage stamp
	meds_filename.psf.field.64x64.fits.fz     -- field of 10x10 PSFs, each dithered within the central pixel
	meds_filename.psf.single.340x340.fits.fz  -- high resolution PSF, for Im3shape, upsampling=5, padding=4
