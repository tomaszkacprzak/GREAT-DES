# ufig2meds
# convert ufig and sextractor output to meds files

data_dir=data
meds_dir=meds
make_cutouts=/home/tomek/Work/code/desdm/wl/trunk/src/bin/make-cutouts
coadd_file=$data_dir/ufig_bcc_0.fixed.fits
cat_file=$data_dir/ufig_bcc_0.meds.cat
coadd_srclist=$data_dir/ufig_bcc_0.srccat
cutout_file=$meds_dir/ufig_bcc_0.meds.fits
sexcat_file=$data_dir/ufig_bcc_0.sexcat
coaddseg_file=$data_dir/ufig_bcc_0_seg.fits
magzp_ref=31.9

make-meds-input $sexcat_file 32 128 > $cat_file

# coadd
COADD_HDU=2; # that is for image 
COADD_WT_HDU=3;
COADD_SEG_HDU=1;
# single exp
SE_HDU=$COADD_HDU; # that is for image
SE_WT_HDU=$COADD_WT_HDU;
SE_BADPIX_HDU=4;
# other
SKY_HDU=1; # in the bkg file
SEG_HDU=1; # in the seg file

coadd_hdu=$COADD_HDU coadd_wt_hdu=$COADD_WT_HDU coadd_seg_hdu=$COADD_SEG_HDU se_hdu=$SE_HDU se_wt_hdu=$SE_WT_HDU se_badpix_hdu=$SE_BADPIX_HDU sky_hdu=$SKY_HDU seg_hdu=$SEG_HDU

python fix_hdus.py data/ufig_bcc_0.fits data/ufig_bcc_0_bkg.fits data/ufig_bcc_0.fixed.fits 

cmd="/home/tomek/Work/code/desdm/wl/trunk/src/bin/make-cutouts coadd_file=$coadd_file cat_file=$cat_file coadd_srclist=$coadd_srclist cutout_file=$cutout_file coaddseg_file=$coaddseg_file magzp_ref=$magzp_ref sky_hdu=$SKY_HDU coadd_hdu=$COADD_HDU coadd_wt_hdu=$COADD_WT_HDU coadd_seg_hdu=$COADD_SEG_HDU se_hdu=$SE_HDU se_wt_hdu=$SE_WT_HDU se_badpix_hdu=$SE_BADPIX_HDU sky_hdu=$SKY_HDU seg_hdu=$SEG_HDU" 
echo $cmd
$cmd

