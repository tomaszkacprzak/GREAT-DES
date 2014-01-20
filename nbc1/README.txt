rm *fits*

# generate small amount of data to get the standard deviaiton of ellipticity estimators
python $CODE/GREAT-DES/nbc1/nbc1_generate_data.py -c $CODE/GREAT-DES/nbc1/nbc1.yaml --debug  -v2
# run hsm on the small amount of data
python $CODE/GREAT-DES/nbc1/nbc1_run_hsm.py -i nbc1.cat
# analyse these results
python $CODE/GREAT-DES/nbc1/nbc1_plots.py -i nbc1.cat -c $CODE/GREAT-DES/nbc1/nbc1.yaml
# copy the resultus to use it later for calculations of number of needed galaxies using ellipticity standard deviation 
cp nbc1.cat nbc1.stde.cat
cp nbc1.cat.hsm.stats.cat nbc1.cat.hsm.stats.stde.cat
# now generate the main data sample
python $CODE/GREAT-DES/nbc1/nbc1_generate_data.py -c $CODE/GREAT-DES/nbc1/nbc1.yaml  -v2
# use results from im3shape
python $CODE/GREAT-DES/nbc1/nbc1_plots.py -i nbc1.tiled.cat -c nbc1.tiled.yaml --method_id im3
# 20140114 use run-tiled-001, im3shape 'cleaned' catalogs from Michael Hirsch
python $CODE/GREAT-DES/nbc1/nbc1_plots.py -i nbc1.tiled.cat -c nbc1.tiled.yaml --method_id im3_cleaned
# then run the calibration 
python $CODE/GREAT-DES/nbc1/nbc1_calibrate.py -i run-tiled-001/nbc1.tiled.cat -c run-tiled-001/nbc1.tiled.yaml --method_id im3_cleaned