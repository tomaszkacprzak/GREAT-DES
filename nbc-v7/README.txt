python ~/code/GREAT-DES/nbc-v7/nbc_v7.py -c nbc-v7.yaml -m im3shape -a get_distributions -n 200
python ~/code/GREAT-DES/nbc-v7/nbc_v7.py -c nbc-v7.yaml -m im3shape -a get_calibration -n 200
python ~/code/GREAT-DES/nbc-v7/nbc_v7.py -c nbc-v7.yaml -m im3shape -a get_bias_model 
python ~/code/GREAT-DES/nbc-v7/nbc_v7.py -c nbc-v7.yaml -m im3shape -a apply_calibration_des 
python ~/code/GREAT-DES/nbc-v7/nbc_v7.py -c nbc-v7.yaml -m im3shape -a get_calibration --use_calibration -n 200
python ~/code/GREAT-DES/nbc-v7/nbc_v7.py -c nbc-v7.yaml -m im3shape -a plot_bias_vs_redshift
python ~/code/GREAT-DES/nbc-v7/nbc_v7.py -c nbc-ngmix009.yaml -m ngmix -a plot_bias_vs_redshift

# bias table used for v7.2
bias_table.v7.2.fits