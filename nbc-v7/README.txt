Calibration commands:

# after getting the results, run this command to produce table of snr,size,m,c etc.
python ~/code/GREAT-DES/nbc-v7/nbc_v7.py -c nbc-v7.yaml -m im3shape -a get_calibration -n 600

# then use the table to make a model of m,v as a function of snr, size
python ~/code/GREAT-DES/nbc-v7/nbc_v7.py -c nbc-v7.yaml -m im3shape -a get_bias_model 

# then apply calibration to DES data
python ~/code/GREAT-DES/nbc-v7/nbc_v7.py -c nbc-v7.yaml -m im3shape -a apply_calibration_des 

# then apply calibration to simulated data
python ~/code/GREAT-DES/nbc-v7/nbc_v7.py -c nbc-v7.yaml -m im3shape -a apply_calibration_sim

# check if calibration works: plot m, alpha vs redshift in COSMOS catalogs
python ~/code/GREAT-DES/nbc-v7/nbc_v7.py -c nbc-v7.yaml -m im3shape -a get_bias_vs_redshift
python ~/code/GREAT-DES/nbc-v7/nbc_v7.py -c nbc-v7.yaml -m im3shape -a plot_bias_vs_redshift






Other commands which work (or used to work):

python ~/code/GREAT-DES/nbc-v7/nbc_v7.py -c nbc-v7.yaml -m im3shape -a get_distributions -n 200
python ~/code/GREAT-DES/nbc-v7/nbc_v7.py -c nbc-v7.yaml -m im3shape -a get_calibration -n 200
python ~/code/GREAT-DES/nbc-v7/nbc_v7.py -c nbc-v7.yaml -m im3shape -a get_bias_model 
python ~/code/GREAT-DES/nbc-v7/nbc_v7.py -c nbc-v7.yaml -m im3shape -a apply_calibration_des 
python ~/code/GREAT-DES/nbc-v7/nbc_v7.py -c nbc-v7.yaml -m im3shape -a get_calibration --use_calibration -n 200
python ~/code/GREAT-DES/nbc-v7/nbc_v7.py -c nbc-v7.yaml -m im3shape -a plot_bias_vs_redshift
python ~/code/GREAT-DES/nbc-v7/nbc_v7.py -c nbc-ngmix009.yaml -m ngmix -a plot_bias_vs_redshift

ipy ~/code/GREAT-DES/nbc-v7/nbc_v7.py -c nbc-v7.yaml -a get_distributions
ipy ~/code/GREAT-DES/nbc-v7/nbc_v7.py -c nbc-v7.yaml -a get_PSF_leakage -n 200 -v3

# [ucabtok@login08 141016_nbc_v7]$ qsub gen.legion.missing.82.sh
# [ucabtok@login08 141016_nbc_v7]$ qsub gen.legion.missing.83.sh
# [ucabtok@login08 141016_nbc_v7]$ qsub gen.legion.missing.101.sh

# source ~/source_ucl_des_shear.sh
# python -m py3shape.analyze_meds /home/ucabtok/Scratch/141016_nbc_v7/data/nbc2.meds.000.g00.fits /home/ucabtok/ucl_des_shear//utils/des.ini /home/ucabtok/Scratch/141016_nbc_v7/outputs/nbc2.meds.000.g00.fits-0 --file_selected_objects /home/ucabtok/Scratch/141016_nbc_v7/dummy.txt --first 0 --count 2 -v3 --save-source-paths=/home/ucabtok/Scratch/141016_nbc_v7/outputs/nbc2.meds.000.g00.txt.sources --no-clobber --psf_input_method=great_des
# don't use nearest_pixel_mask

# test post processing
# python -m py3shape.analyze_meds /home/ucabtok//projects/141016_nbc_v7/data/nbc2.meds.000.g00.fits /home/ucabtok/ucl_des_shear//utils/des.ini /home/ucabtok//projects/141016_nbc_v7/outputs/nbc2.meds.000.g00.fits-0 --file_selected_objects /home/ucabtok/Scratch/141016_nbc_v7/dummy.txt --first 0 --count 3  -v3 --save-source-paths=/home/ucabtok//projects/141016_nbc_v7/outputs/nbc2.meds.000.g00.txt.sources --no-clobber --psf_input_method=great_des
# python -m py3shape.analyze_meds /home/ucabtok//projects/141016_nbc_v7/data/nbc2.meds.000.g00.fits /home/ucabtok/ucl_des_shear//utils/des.ini /home/ucabtok//projects/141016_nbc_v7/outputs/nbc2.meds.000.g00.fits-1000 --file_selected_objects /home/ucabtok/Scratch/141016_nbc_v7/dummy.txt --first 0 --count 3  -v3 --save-source-paths=/home/ucabtok//projects/141016_nbc_v7/outputs/nbc2.meds.000.g00.txt.sources --no-clobber --psf_input_method=great_des

# mkdir outputs
# mkdir wcs
# mkdir main_cats
# mkdir epoch_cats
# echo "/home/ucabtok//projects/141016_nbc_v7/data/nbc2.meds.000.g00.fits" > meds_list.test.txt
# python ~/ucl_des_shear/des_post/create_wcs_cat.py meds_list.test.txt
# export UCL_DES_SHEAR=~/code/ucl_des_shear/
# export PYTHONPATH=$PYTHONPATH:/Users/tomek/code/ucl_des_shear:/Users/tomek/code/ucl_des_shear/im3shape/py3shape
# python -m des_post.postprocess meds_list.test.txt --no-blinding -v --mode GREATDES
# python -m des_post.postprocess meds_list.txt --no-blinding -v --mode GREATDES

# post process at plempy
export UCL_DES_SHEAR=/home/ast/refreg/temp/tomek/code/ucl_des_shear
export PYTHONPATH=$PYTHONPATH:/home/ast/refreg/temp/tomek/code/ucl_des_shear/im3shape/py3shape/:/home/ast/refreg/temp/tomek/code/ucl_des_shear/
# python -m des_post.postprocess meds_list.txt --no-blinding -v --mode GREATDES

python ~/code/GREAT-DES/nbc-v7/nbc_v7.py -c nbc-v7.yaml -m im3shape -a get_calibration
python ~/code/GREAT-DES/nbc-v7/nbc_v7.py -c nbc-v7.yaml -m im3shape -a get_bias_model -n 126
python ~/code/GREAT-DES/nbc-v7/nbc_v7.py -c nbc-v7.yaml -m im3shape -a apply_calibration_des
python ~/code/ucl_des_shear/des_post/add_weights.py -c weights.des.yaml -n440 -v2 -a apply_weights

python ~/code/GREAT-DES/nbc-v7/nbc_v7.py -c nbc-v7.yaml -m im3shape -a apply_calibration_sim
python ~/code/ucl_des_shear/des_post/add_weights.py -c weights.sim.yaml -n5000 -v2 -a apply_weights --clobber

python ~/code/GREAT-DES/nbc-v7/nbc_v7.py -c nbc-v7.yaml -m im3shape -a get_bias_vs_redshift
python ~/code/GREAT-DES/nbc-v7/nbc_v7.py -c nbc-v7.yaml -m im3shape -a plot_bias_vs_redshift

python ~/code/GREAT-DES/nbc-v7/nbc_v7.py -c nbc-ngmix009.yaml -m ngmix -a get_bias_vs_redshift
python ~/code/GREAT-DES/nbc-v7/nbc_v7.py -c nbc-ngmix009.yaml -m ngmix -a plot_bias_vs_redshift

python ~/code/GREAT-DES/nbc-v7/nbc_v7.py -c nbc-v7.yaml -m im3shape -a get_jacknife_regions

minion$ python ~/code/GREAT-DES/nbc-v7/nbc_v7.py -c nbc-ngmix009-002.yaml -m ngmix -a get_meane_vs_size

