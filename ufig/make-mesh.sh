# ufig2meds
# convert ufig and sextractor output to meds files

display_usage() { 
	echo -e "\nUsage:\n$0 [tilename magzp_ref] \n" 
	echo "example: sh ufig2meds.sh DES0441-4414/DES0441-4414_r 30.5650332328"
	} 
# if less than two arguments supplied, display usage 
	if [  $# -ne 2 ] 
	then 
		display_usage
		exit 1
	fi 
 
# check whether user had supplied -h or --help . If yes display usage 
	if [[ ( $# == "--help") ||  $# == "-h" ]] 
	then 
		display_usage
		exit 0
	fi


TILENAME=$1
SUFF_BG=_bg.fits
SUFF_SEG=_seg.fits
SUFF_SRC=.srccat
SUFF_SEX=.sexcat
SUFF_MEDSCAT=.medscat
SUFF_CUTOUT=.meds.fits
SUFF_FIX=.fix.fits

# inputs
coadd_file=$TILENAME.fits
coadd_bkg_file=$TILENAME$SUFF_BG
coadd_seg_file=$TILENAME$SUFF_SEG
sexcat_file=$TILENAME$SUFF_SEX
magzp_ref=$2

# outputs
coadd_srclist=$TILENAME$SUFF_SRC
cat_file=$TILENAME$SUFF_MEDSCAT
cutout_file=$TILENAME$SUFF_CUTOUT
coadd_fixed_file=$TILENAME$SUFF_FIX

# create source cat
cmd="echo 0 0 $coadd_fixed_file $coadd_bkg_file $coadd_seg_file $magzp_ref" 
echo $cmd
$cmd > $coadd_srclist

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

cmd="python /cluster/home02/phys/tomaszk/code/GREAT-DES/ufig/fix_hdus.py $coadd_file $coadd_seg_file $coadd_fixed_file"
echo $cmd; 
$cmd


cmd="/cluster/home02/phys/tomaszk/local/deswl_shapelets/bin/make-cutouts coadd_file=$coadd_fixed_file cat_file=$cat_file coadd_srclist=$coadd_srclist cutout_file=$cutout_file coaddseg_file=$coadd_fixed_file magzp_ref=$magzp_ref sky_hdu=$SKY_HDU coadd_hdu=$COADD_HDU coadd_wt_hdu=$COADD_WT_HDU coadd_seg_hdu=$COADD_SEG_HDU se_hdu=$SE_HDU se_wt_hdu=$SE_WT_HDU se_badpix_hdu=$SE_BADPIX_HDU sky_hdu=$SKY_HDU seg_hdu=$SEG_HDU" 
echo $cmd; $cmd

# coadd_hdu=$COADD_HDU coadd_wt_hdu=$COADD_WT_HDU coadd_seg_hdu=$COADD_SEG_HDU se_hdu=$SE_HDU se_wt_hdu=$SE_WT_HDU se_badpix_hdu=$SE_BADPIX_HDU sky_hdu=$SKY_HDU seg_hdu=$SEG_HDU
# coadd_file=$coadd_fixed_file image, weight
# coaddseg_file=$coadd_fixed_file seg
# sky_hdu=$SKY_HDU coadd_hdu=$COADD_HDU coadd_wt_hdu=$COADD_WT_HDU coadd_seg_hdu=$COADD_SEG_HDU se_hdu=$SE_HDU se_wt_hdu=$SE_WT_HDU se_badpix_hdu=$SE_BADPIX_HDU sky_hdu=$SKY_HDU seg_hdu=$SEG_HDU" 