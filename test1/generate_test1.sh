DESDATA=.
meds_version=005
coadd_run=run1 
galsim_yaml shear_test1.meds.yaml
python shear_test1.psf.py
# profile04ellip00hlr04-test1-meds-005.fits.fz
DIRPATH=./meds/$meds_version/$coadd_run/
mkdir -p $DIRPATH
mv *fits.fz $DIRPATH
mv *psf*.fits $DIRPATH