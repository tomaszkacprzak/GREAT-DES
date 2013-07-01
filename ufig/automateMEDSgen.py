import os
import sys

#Directories (coadds, MEDS, SEs)
data_dir = '/Users/cbruderer/TEMP/'
meds_dir = '/Users/cbruderer/TEMP/'

#Directory make-cutout
make_cutouts = '/Users/cbruderer/DES/devel/wl/trunk/src/bin/make-cutouts'

#HDUS
#coadd
COADD_HDU = '2' # that is for image
COADD_WT_HDU = '3'
COADD_SEG_HDU = '1'
#single exp
SE_HDU = COADD_HDU # that is for image
SE_WT_HDU = COADD_WT_HDU
SE_BADPIX_HDU = '4'
#other
SKY_HDU = '1' # in the bkg file
SEG_HDU = '1' # in the seg file


#Indices of tile to create MEDS file
startindex = 1
endindex = 1

#Template name of BCC tiles (name = template1 + index + ending)
template1 = 'ufig_bcc_'

#File name of dummy source list
fname_dummy = data_dir + 'coadd_srclist.srccat'

# for j in range(startindex,endindex+1):
# 	template2 = template1 + str(j)
# 	print ''
# 	print template2

# 	#Different files
# 	coadd_file = template2 + '.fits'
# 	coadd_cat = template2 + '.sexcat'
# 	coadd_fixed = template2 + '.fixed.fits'
# 	coadd_seg = template2 + '_seg.fits'
# 	coadd_bkg = template2 + '_bkg.fits'

# 	se_fixed = coadd_fixed
# 	se_seg = coadd_seg #for now same as coadd
# 	se_bkg = coadd_bkg #for now same as coadd

# 	magzp_ref = str(31.9) #for now fixed for all the BCC images

# 	#Target files
# 	cat_file = template2 + '.meds.cat'
# 	cutout_file = template2 + '.meds.cat'

# 	#Write source list in dummy srclist-file
# 	srclist = open(fname_dummy,'w')
# 	templine = data_dir+se_fixed + ' ' + data_dir+se_bkg + ' ' + data_dir+se_seg + ' ' + magzp_ref + '\n'
# 	srclist.write(templine)
# 	srclist.close()

# 	#Create catalogs out of sexcat
# 	os.system('make-meds-input ' + data_dir+coadd_cat + ' 32 128 > ' + meds_dir+cat_file)

# 	#Fix HDU (add weight, create bad pixel map (both trivial atm))
# 	os.system('python fix_hdus.py '+ data_dir+coadd_file +' '+ data_dir+coadd_bkg +' '+ data_dir+coadd_fixed)

# 	#Create MEDS
# 	hdus = 'sky_hdu='+SKY_HDU+' '+'coadd_hdu='+COADD_HDU+' '+'coadd_wt_hdu='+COADD_WT_HDU+' '+'coadd_seg_hdu='+COADD_SEG_HDU+' '+'se_hdu='+SE_HDU+' '+'se_wt_hdu='+SE_WT_HDU+' '+'se_badpix_hdu='+SE_BADPIX_HDU+' '+'sky_hdu='+SKY_HDU+' '+'seg_hdu='+SEG_HDU
# 	files = 'coadd_file='+data_dir+coadd_file+' '+'cat_file='+meds_dir+cat_file+' '+'coadd_srclist='+fname_dummy+' '+'cutout_file='+meds_dir+cutout_file+' '+'coaddseg_file='+data_dir+coadd_seg+' '+'magzp_ref='+magzp_ref
# 	cmd = make_cutouts + ' ' + files + ' ' + hdus
# 	print cmd
# 	#sys.exit()
# 	os.system(cmd)


for j in range(startindex,endindex+1):
	template2 = 'ufig_bcc_1'
	print ''
	print template2

	#Different files
	coadd_file = template2 + '.fits'
	coadd_cat = template2 + '.sexcat'
	coadd_fixed = template2 + '.fixed.fits'
	coadd_seg = template2 + '_seg.fits'
	coadd_bkg = template2 + '_bkg.fits'

	se_fixed = coadd_fixed
	se_seg = coadd_seg #for now same as coadd
	se_bkg = coadd_bkg #for now same as coadd

	magzp_ref = 31.9 #for now fixed for all the BCC images

	#Target files
	cat_file = template2 + '.meds.cat'
	cutout_file = template2 + '.meds.cat'

	#Write source list in dummy srclist-file
	#srclist = open(fname_dummy,'w')
	#templine = data_dir+se_fixed + ' ' + data_dir+se_bkg + ' ' + data_dir+se_seg + ' ' + str(magzp_ref) + '\n'
	#srclist.write(templine)
	#srclist.close()

	#Create catalogs out of sexcat
	os.system('make-meds-input ' + data_dir+coadd_cat + ' 32 128 > ' + meds_dir+cat_file)

	#Fix HDU (add weight, create bad pixel map (both trivial atm))
	os.system('python fix_hdus.py '+ data_dir+coadd_file +' '+ data_dir+coadd_bkg +' '+ data_dir+coadd_fixed)

	#Create MEDS
	hdus = 'sky_hdu='+SKY_HDU+' '+'coadd_hdu='+COADD_HDU+' '+'coadd_wt_hdu='+COADD_WT_HDU+' '+'coadd_seg_hdu='+COADD_SEG_HDU+' '+'se_hdu='+SE_HDU+' '+'se_wt_hdu='+SE_WT_HDU+' '+'se_badpix_hdu='+SE_BADPIX_HDU+' '+'sky_hdu='+SKY_HDU+' '+'seg_hdu='+SEG_HDU
	files = 'coadd_file='+data_dir+coadd_file+' '+'cat_file='+meds_dir+cat_file+' '+'coadd_srclist='+fname_dummy+' '+'cutout_file='+meds_dir+cutout_file+' '+'coaddseg_file='+data_dir+coadd_seg+' '+'magzp_ref='+str(magzp_ref)
	cmd = make_cutouts + ' ' + files + ' ' + hdus
	print cmd
	#sys.exit()
	os.system(cmd)



# coadd=$data_dir/ufig_bcc_0.fits
# coadd_file=$data_dir/ufig_bcc_0.fixed.fits
# cat_file=$data_dir/ufig_bcc_0.meds.cat

# coadd_srclist=$data_dir/ufig_bcc_0.srccat

# cutout_file=$meds_dir/ufig_bcc_0.meds.fits
# sexcat_file=$data_dir/ufig_bcc_0.sexcat
# coaddseg_file=$data_dir/ufig_bcc_0_seg.fits
# se_bkg_file=$data_dir/ufig_bcc_0_bkg.fits
# magzp_ref=31.9

# echo $coadd_file $se_bkg_file $coaddseg_file $magzp_ref >  $coadd_srclist

# make-meds-input $sexcat_file 32 128 > $cat_file

# #coadd_hdu=$COADD_HDU coadd_wt_hdu=$COADD_WT_HDU coadd_seg_hdu=$COADD_SEG_HDU se_hdu=$SE_HDU se_wt_hdu=$SE_WT_HDU se_badpix_hdu=$SE_BADPIX_HDU sky_hdu=$SKY_HDU seg_hdu=$SEG_HDU

# python fix_hdus.py $data_dir/ufig_bcc_0.fits $data_dir/ufig_bcc_0_bkg.fits $data_dir/ufig_bcc_0.fixed.fits 

# cmd="$make_cutouts coadd_file=$coadd_file cat_file=$cat_file coadd_srclist=$coadd_srclist cutout_file=$cutout_file coaddseg_file=$coaddseg_file magzp_ref=$magzp_ref sky_hdu=$SKY_HDU coadd_hdu=$COADD_HDU coadd_wt_hdu=$COADD_WT_HDU coadd_seg_hdu=$COADD_SEG_HDU se_hdu=$SE_HDU se_wt_hdu=$SE_WT_HDU se_badpix_hdu=$SE_BADPIX_HDU sky_hdu=$SKY_HDU seg_hdu=$SEG_HDU" 
# echo $cmd
# $cmd