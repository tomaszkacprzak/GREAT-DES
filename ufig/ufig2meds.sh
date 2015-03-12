# ufig2meds
# convert ufig and sextractor output to meds files

# inputs
coadd_file=DES0441-4414_ufig.fits
coadd_bkg_file=DES0441-4414_ufig_bg.fits
coadd_seg_file=DES0441-4414_ufig_seg.fits
coadd_srclist=DES0441-4414_ufig2.srccat
cutout_file=DES0441-4414_ufig.meds2.fits
sexcat_file=DES0441-4414_ufig.fixed.sexcat
coaddseg_file=DES0441-4414_ufig_seg.fits
magzp_ref=30.5650332328

# outputs
cat_file=DES0441-4414_ufig.medscat
coadd_fixed_file=DES0441-4414_ufig.fixed.fits

# run commands
cmd="/cluster/home02/phys/tomaszk/local/bin/make-meds-input $sexcat_file 32 128" 
echo $cmd; 
$cmd > $cat_file

# coadd
COADD_HDU=1; # that is for image 
COADD_SEG_HDU=2;
COADD_WT_HDU=3;
# single exp
SE_HDU=$COADD_HDU; # that is for image
SE_WT_HDU=$COADD_WT_HDU;
SE_BADPIX_HDU=4;
# other
SKY_HDU=1; # in the bkg file
SEG_HDU=1; # in the seg file

coadd_hdu=$COADD_HDU coadd_wt_hdu=$COADD_WT_HDU coadd_seg_hdu=$COADD_SEG_HDU se_hdu=$SE_HDU se_wt_hdu=$SE_WT_HDU se_badpix_hdu=$SE_BADPIX_HDU sky_hdu=$SKY_HDU seg_hdu=$SEG_HDU

# python fix_hdus.py data/ufig_bcc_0.fits data/ufig_bcc_0_bkg.fits data/ufig_bcc_0.fixed.fits 
cmd="python /cluster/home02/phys/tomaszk/code/GREAT-DES/ufig/fix_hdus.py $coadd_file $coadd_seg_file $coadd_fixed_file"
echo $cmd; 
$cmd

cmd="/cluster/home02/phys/tomaszk/local/deswl_shapelets/bin/make-cutouts coadd_file=$coadd_fixed_file cat_file=$cat_file coadd_srclist=$coadd_srclist cutout_file=$cutout_file coaddseg_file=$coadd_fixed_file magzp_ref=$magzp_ref sky_hdu=$SKY_HDU coadd_hdu=$COADD_HDU coadd_wt_hdu=$COADD_WT_HDU coadd_seg_hdu=$COADD_SEG_HDU se_hdu=$SE_HDU se_wt_hdu=$SE_WT_HDU se_badpix_hdu=$SE_BADPIX_HDU sky_hdu=$SKY_HDU seg_hdu=$SEG_HDU" 
echo $cmd; $cmd

