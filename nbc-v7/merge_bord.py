from astropy.table import Table, vstack, Column
# astropy FTW
import numpy as np
import glob, os, pyfits, time

import os 
import matplotlib as mpl
if 'DISPLAY' not in os.environ:
    mpl.use('agg')
    print 'using backend ' , mpl.get_backend()
import sys, logging, yaml, argparse, time, pyfits, warnings, subprocess
import pylab as pl
import numpy as np
from nbc_dtypes import *
warnings.simplefilter('once')

logging_level = logging.INFO
log = logging.getLogger("merge_bord") 
log.setLevel(logging_level)  
log_formatter = logging.Formatter("%(asctime)s %(name)s %(levelname)s   %(message)s", "%Y-%m-%d %H:%M:%S")
stream_handler = logging.StreamHandler(sys.stdout)
stream_handler.setFormatter(log_formatter)
log.addHandler(stream_handler)
log.propagate = False


def add_col(rec, name, arr, dtype=None):
    import numpy
    arr = numpy.asarray(arr)
    if dtype is None:
        dtype = arr.dtype
    newdtype = numpy.dtype(rec.dtype.descr + [(name, dtype)])
    newrec = numpy.empty(rec.shape, dtype=newdtype)
    for field in rec.dtype.fields:
        newrec[field] = rec[field]
    newrec[name] = arr
    return newrec

 
def get_filelists(b_path,d_path):

    b_files = glob.glob(b_path)
    d_files = glob.glob(d_path)
    b_files.sort()
    d_files.sort()

    b_files_base = [ os.path.basename(f) for f in b_files]
    d_files_base = [ os.path.basename(f) for f in d_files]

    print len(b_files), len(d_files)

    b_files_dir = os.path.dirname(b_files[0])
    d_files_dir = os.path.dirname(d_files[0])

    intersec = set.intersection(set(b_files_base),set(d_files_base))

    b_files_out = []
    d_files_out = []
    for fi in intersec:
        b_files_out.append( os.path.join(b_files_dir,fi) )
        d_files_out.append( os.path.join(d_files_dir,fi) )

    # Strip out the one tile in bulge but not disc
    # bulge_files = [f for f in bulge_files if "DES1000+0209" not in f]

    b_files_out.sort()
    d_files_out.sort()

    return b_files_out, d_files_out

def merge_files(bulge_files,disc_files,truth_files=None,tag='test',greatdes=False):

    n_total_sva1 = 0
    n_b_sva1 = 0
    n_d_sva1 = 0

    cols = ['coadd_objects_id','e1','e2','mean_psf_e1_sky','mean_psf_e2_sky','ra','dec','ra_as','dec_as','info_flag','error_flag','snr','mean_rgpp_rp','radius','disc_flux','bulge_flux','n_exposure','stamp_size','mean_mask_fraction','is_bulge','is_disc']
    cols_tru = ['psf_e1','psf_e2','id_shear','cosmos_mag_auto','g1_true','g2_true','sf_hlr']

    ifiles = 0

    list_cats = []

    print 'n_disc_files',len(disc_files)
    print 'n_bulge_files',len(bulge_files)

    bulge_files.sort()
    disc_files.sort()

    for bulge_filename, disc_filename in zip(bulge_files, disc_files):

        if not os.path.split(disc_filename)[1]==os.path.split(bulge_filename)[1]:
            import pdb; pdb.set_trace()

        filename_new = os.path.join(args.dir_out,os.path.basename(bulge_filename))
        if os.path.isfile(filename_new) and (args.clobber.lower()=='skip'):
            print('File exists, skipping: %s' % filename_new)
            continue


        truth_filename = disc_filename.replace('meds','truth').replace('results_disc','data')

        cat_b = np.array(pyfits.getdata(bulge_filename))
        cat_d = np.array(pyfits.getdata(disc_filename))
        if truth_files!=None: cat_t = np.array(Table.read(truth_filename))
        print '---------', ifiles, time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
        print bulge_filename,len(cat_b)
        print disc_filename,len(cat_d)
        if truth_files!=None: print truth_filename,len(cat_t)
        
        
        # need to match coadd_object_ids here
        intersect = np.intersect1d(cat_b['coadd_objects_id'],cat_d['coadd_objects_id'])
        
        select_b = np.in1d(cat_b['coadd_objects_id'],intersect)
        select_d = np.in1d(cat_d['coadd_objects_id'],intersect)
        
        ndiscard_bulge =  select_b.size - select_b.sum()
        ndiscard_disc =   select_d.size - select_d.sum()
        cat_b = cat_b[select_b]
        cat_d = cat_d[select_d]
        if truth_files!=None: cat_t = cat_t[cat_d['coadd_objects_id']]


        print 'len intersect b=',len(cat_b),'len intersect d=',len(cat_d),'ndiscard_disc d=',ndiscard_disc, 'ndiscard_b=', ndiscard_bulge
        
        if not (cat_b['coadd_objects_id'] == cat_d['coadd_objects_id']).all():
            import pdb; pdb.set_trace()


        bulge_good = (cat_b['info_flag']==0) & (cat_b['error_flag']==0)
        disc_good = (cat_d['info_flag']==0) & (cat_d['error_flag']==0)
    
        bulge_bad = ~bulge_good
        disc_bad = ~disc_good

        bulge_better = (cat_b['likelihood']-cat_d['likelihood']) > 0
        disc_better = (cat_b['likelihood']-cat_d['likelihood']) < 0
        
        cat_final = cat_d.copy()

        # do a pre-assignment
        cat_final[bulge_better] = cat_b[bulge_better]
        cat_final[disc_better]  = cat_d[disc_better]

        
        # SO in this bit we allow the info_flag to override
        # the likelihood in deciding which is better?
        # do we want that?
        
        select_b = (bulge_good&disc_bad) | (bulge_good&disc_good&bulge_better)
        select_d = (disc_good&bulge_bad) | (bulge_good&disc_good&disc_better)
        
        # cat_final[select_b] = cat_b[select_b]
        # cat_final[select_d] = cat_d[select_d]

        # cat_final.add_column(Column(name='is_bulge', data=select_b.astype(int)))
        # cat_final.add_col(Column(name='is_disc', data=select_d.astype(int)))
    
        cat_final = add_col(cat_final,'is_bulge', bulge_better,dtype=np.bool)
        cat_final = add_col(cat_final,'is_disc', disc_better,dtype=np.bool)

        if truth_files!=None: 
            for ct in cols_tru: 
                cat_final = add_col(cat_final,ct, cat_t[ct])


        # finally apply the info and error flags 
        # select_flags = (cat_b['info_flag']==0)&(cat_final['error_flag']==0)
        # cat_final = cat_final[select_flags]

        
        n_total = len(cat_final)
        n_common = np.sum(bulge_good& disc_good)
        n_common_b = np.sum(bulge_good& disc_good&bulge_better)
        n_common_d = np.sum(bulge_good& disc_good&disc_better)
        n_only_b = np.sum(bulge_good&(~ disc_good))
        n_only_d = np.sum(disc_good&(~bulge_good))  
        n_check = n_only_b + n_only_d +  n_common_b + n_common_d        
        
        print 'n_total =%d n_common =%d n_common_b =%d n_common_d =%d n_only_b =%d n_only_d=%d n_check=%d' % (n_total,n_common,n_common_b,n_common_d,n_only_b,n_only_d,n_check)
                
        # filename_new = bulge_filename.replace('bulge','bord')
        if (args.clobber.lower()=='true') or ( (~os.path.isfile(filename_new)) and (args.clobber.lower()=='skip')) or ( (~os.path.isfile(filename_new)) and (args.clobber.lower()=='false')):
            pyfits.writeto(filename_new,cat_final, clobber=True)
            print 'wrote',filename_new
        elif os.path.isfile(filename_new) and (args.clobber.lower()=='false'):
            raise('File exists: %s' % filename_new)
        elif os.path.isfile(filename_new) and (args.clobber.lower()=='skip'):
            print('File exists, skipping: %s' % filename_new)
        
        n_total_sva1 += n_total
        n_b_sva1 += n_common_b+n_only_b
        n_d_sva1 += n_common_d+n_only_d

        if truth_files!=None:
            part_id = int(os.path.basename(truth_filename)[10:13])
            shear_id = int(os.path.basename(truth_filename)[15:17])
            prefix = 10000000*part_id + 100000*shear_id
            cat_final['coadd_objects_id'] += prefix
            list_cats.append(cat_final[cols])

        ifiles+=1

    print 'n_total_sva1=%d n_b_sva1=%d n_d_sva1=%d' % (n_total_sva1,n_b_sva1,n_d_sva1)

    # cat_final_conc = np.concatenate(list_cats)
    # print len(cat_final_conc)
    # filename_full = '%sbord.fits' % tag

    # import pdb; pdb.set_trace()
        




def main():

    global args

    description = 'Get input catalogs for GREAT-DES NBC2 simulation'
    parser = argparse.ArgumentParser(description=description, add_help=True)
    parser.add_argument('-v', '--verbosity', type=int, action='store', default=2, choices=(0, 1, 2, 3 ), help='integer verbosity level: min=0, max=3 [default=2]')
    parser.add_argument('-do', '--dir_out', type=str,    action='store',  help='directory which stores output')
    parser.add_argument('-fb', '--files_bulge', type=str,  action='store',  help='wildcard to select bulge files ex /path/to/dir/DES*' )
    parser.add_argument('-fd', '--files_disc', type=str,   action='store',  help='wildcard to select disc  files ex /path/to/dir/DES*' )
    parser.add_argument( '--tag', type=str,   action='store',  help='not sure what this is for' )
    parser.add_argument( '--greatdes', action='store_true',  help='GREAT-DES mode (not using some bits in info flags' )
    parser.add_argument( '--clobber', action='store', default='False', type=str, help='GREAT-DES mode (not using some bits in info flags' )

    args = parser.parse_args()
    logging_levels = { 0: logging.CRITICAL, 
                       1: logging.WARNING,
                       2: logging.INFO,
                       3: logging.DEBUG }
    logging_level = logging_levels[args.verbosity]
    log.setLevel(logging_level)  
    
    log.info(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))

    # tag = 'v9test'
    # bulge_path = '/Users/tomek/data/DES/im3shape-v9/im3shape_v9/bulge/main/DES*'
    # disc_path = '/Users/tomek/data/DES/im3shape-v9/im3shape_v9/disc/main/DES*'
    tag = 'nbc006'
    bulge_path = args.files_bulge
    disc_path  = args.files_disc

    bulge_files,disc_files = get_filelists(bulge_path,disc_path)
    merge_files(bulge_files,disc_files,truth_files=None,tag=tag,greatdes=args.greatdes)

main()
