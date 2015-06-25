n_first=000
n_meds=600
n_shears=8
n_gals=1000
n_split=10
dir_home='/cluster/home02/phys/tomaszk/'
dir_wdir='/cluster/scratch_xl/shareholder/refregier/gamperl/bruclaud/tomaszk/150104_nbc_v8/001-greatdes005/'


templ="python -m py3shape.analyze_meds %s %s/code/ucl_des_shear/utils/des.ini %s --file_selected_objects %s/dummy.txt --first %d --count %d  --save-source-paths=%s --psf_input_method=great_des \n"

filename_cmds = 'commands.sh'
f = open(filename_cmds,'w')

n_last = n_first+n_meds
for im in range(n_first,n_last):
	for ig in range(n_shears):
		for ip in range(n_split):
			first = ip*n_gals
			filename_meds = '%s/data/nbc.meds.%03d.g%02d.fits' % (dir_wdir,im,ig)
			filename_out = 'nbc.meds.%03d.g%02d.fits-%d'% (im,ig,first)
			filename_sources = 'nbc.meds.%03d.g%02d.txt.sources' % (im,ig)
			line = templ % (filename_meds,dir_home,filename_out,dir_wdir,first,n_gals,filename_sources)
			f.write(line)

f.close()
print 'wrote', filename_cmds

