DESDATA=.
meds_version=005
coadd_run=run1 
galsim_yaml shear_test1.meds.yaml
# test1-profile04ellip00hlr04-meds-005.fits.fz
DIRPATH=./meds/$meds_version/$coadd_run/
mkdir -p $DIRPATH
mv *fits.fz $DIRPATH