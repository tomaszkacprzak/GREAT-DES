import numpy as np; import pylab as pl; import tktools as tt;
import  sys, logging, yaml, argparse, time, copy, itertools, tktools, warnings, os, fitsio, pyfits;
warnings.simplefilter("once")
sys.path.append('/home/tomek/code/tktools');
sys.path.append('/Users/tomek/code/tktools');
logging_level = logging.INFO; logger = logging.getLogger("nbc-v7-errors"); logger.setLevel(logging_level)  
log_formatter = logging.Formatter("%(asctime)s %(name)s %(levelname)s   %(message)s", "%Y-%m-%d %H:%M:%S")
stream_handler = logging.StreamHandler(sys.stdout); stream_handler.setFormatter(log_formatter)
if logger.handlers == [] : logger.addHandler(stream_handler); logger.propagate = False
import nbc_v7_select, nbc_v7_stats
import plotstools, fitting


def load_selection():

        res_sim = tt.load( '%s/res_sim.fits' % args.output_dir)
        res_des = tt.load( '%s/res_des.fits' % args.output_dir)
        res_tru = tt.load( '%s/res_tru.fits' % args.output_dir)

        return res_sim, res_tru, res_des


def get_bias_for_each_galaxy():

        res_sim, res_tru, res_des = load_selection()

        logger.info('using final selection')
        logger.info(config['selection_string_final_sim'])
        cat_res = res_sim
        cat_tru = res_tru
        exec config['selection_string_final_sim']
        select_sim = select.copy()
        res_sim = res_sim[select_sim]
        res_tru = res_tru[select_sim]

        logger.info('getting unique cosmos ids')
        unique_cosmos_ids = np.unique(res_tru['id_cosmos'])

        list_bias = []
        n_skip = 0

        nbc_v7_stats.logger.setLevel(logging.WARNING)

        filename_str = None  # kill plots in get_mc

        for iu, id_cosmos in enumerate(unique_cosmos_ids):
            
            select = res_tru['id_cosmos']==id_cosmos
            res_sim_select = res_sim[select]
            res_tru_select = res_tru[select]
            n_gals_sim = len(res_sim_select)

            std_e = np.std(res_sim_select['e1'],ddof=1)
            bulge_fraction = np.sum(res_sim_select['disc_flux']==0)/float(len(res_sim_select))
            median_snr = np.median(res_sim_select['snr'])
            median_rgpp_rp = np.median(res_sim_select['mean_rgpp_rp'])
            median_radius = np.median(res_sim_select['radius'])

            if n_gals_sim < 50:
                mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2,mms,std_mms,ccs,std_ccs,mm1s,std_mm1s,cc1s,std_cc1s,mm2s,std_mm2s,cc2s,std_cc2s = (0,)*36
                n_skip += n_gals_sim 
            else:
                mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,pmm1,std_pmm1,pmm2,std_pmm2,mean_e1,mean_e2,stdm_e1,stdm_e2,mms,std_mms,ccs,std_ccs,mm1s,std_mm1s,cc1s,std_cc1s,mm2s,std_mm2s,cc2s,std_cc2s=nbc_v7_stats.get_mc(res_sim_select,res_tru_select,None,use_calibration=args.use_calibration,use_weights=args.use_weights,filename_str=filename_str,correct_selection_bias=config['correct_selection_bias'])
          
            logger.info('------------------------------------- id=%5d id_cosmos=%6d n_gals=%4d m=%2.3f +/- %2.3f bulge_fraction=%2.2f median_snr=%2.2f median_rgpp_rp=%2.2f' % (iu,id_cosmos,n_gals_sim,mm,std_mm,bulge_fraction,median_snr,median_rgpp_rp))

            list_bias.append( [iu,id_cosmos,n_gals_sim,std_e,mm,std_mm,cc,std_cc,mm1,std_mm1,mm2,std_mm2,cc1,std_cc1,cc1,std_cc2,pmm,std_pmm,pcc,std_pcc,bulge_fraction,mms,std_mms,ccs,std_ccs,mm1s,std_mm1s,cc1s,std_cc1s,mm2s,std_mm2s,cc2s,std_cc2s,median_snr,median_rgpp_rp,median_radius] )

        arr_bias = tktools.arr2rec(np.array(list_bias),dtype={'names': ["index","id_cosmos","n_gals_sim","std_e","mm","std_mm","cc","std_cc","mm1","std_mm1","mm2","std_mm2","cc1","std_cc1","cc2","std_cc2","pmm","std_pmm","pcc","std_pcc","bulge_fraction","mms","std_mms","ccs","std_ccs","mm1s","std_mm1s","cc1s","std_cc1s","mm2s","std_mm2s","cc2s","std_cc2s","median_snr","median_rgpp_rp","median_radius"], 'formats': ['i4']*3 + ['f8']*33 })

        filename_table = os.path.join(args.output_dir,'cosmos_gals_bias.fits')
        tt.save(filename_table,arr_bias,clobber=True)

def plot_bias_stats():

        filename_table = os.path.join(args.output_dir,'cosmos_gals_bias.fits')
        arr_bias = tt.load(filename_table)        

        select = np.isfinite(arr_bias['mm'])
        arr_bias_select = arr_bias[select]
        x_bootstrap= arr_bias_select['mm']; n_bootstrap=100; sigma_bootstrap = np.std(np.mean(x_bootstrap[np.reshape(np.random.choice(len(x_bootstrap),len(x_bootstrap)*n_bootstrap,replace=True),[len(x_bootstrap),n_bootstrap])],axis=0),ddof=1)
        mean_m = np.mean(arr_bias_select['mm'])
        median_m = np.median(arr_bias_select['mm'])
        logger.info('mean m=%2.4f +/- %2.4f median_m=%2.4f' % (mean_m,sigma_bootstrap,median_m))

        import pdb; pdb.set_trace()


def get_histmatrix():

        res_sim, res_tru, res_des = load_selection()

        logger.info('using final selection')
        logger.info(config['selection_string_final_sim'])
        cat_res = res_sim
        cat_tru = res_tru
        exec config['selection_string_final_sim']
        select_sim = select.copy()
        res_sim = res_sim[select_sim]
        res_tru = res_tru[select_sim]

        logger.info('getting unique cosmos ids')
        unique_cosmos_ids = np.unique(res_tru['id_cosmos'])
        n_unique_cosmos_ids = len(unique_cosmos_ids)

        list_bias = []
        n_skip = 0

        nbc_v7_stats.logger.setLevel(logging.WARNING)

        filename_str = None  # kill plots in get_mc

        bins_e1 = np.linspace(-1,1,100)
        bins_e2 = np.linspace(-1,1,100)
        bins_sn = np.logspace(1,2,100)
        bins_sz = np.linspace(1.15,3,100)

        n_cols = bins_e1-1+bins_e2-1+bins_sn-1+bins_sz

        histmatrix = np.zeros(n_unique_cosmos_ids,n_cols)

        for iu, id_cosmos in enumerate(unique_cosmos_ids[:100]):
            
            select = res_tru['id_cosmos']==id_cosmos
            res_sim_select = res_sim[select]
            res_tru_select = res_tru[select]
            n_gals_sim = len(res_sim_select)

            std_e = np.std(res_sim_select['e1'],ddof=1)
            bulge_fraction = np.sum(res_sim_select['disc_flux']==0)/float(len(res_sim_select))
            median_snr = np.median(res_sim_select['snr'])
            median_rgpp_rp = np.median(res_sim_select['mean_rgpp_rp'])
            median_radius = np.median(res_sim_select['radius'])

            great_des_e1 = res_sim_select['e1'] - res_tru_select['g1_true']
            great_des_e1 = res_sim_select['e2'] - res_tru_select['g2_true']

            hist_e1,_,_ = np.histogram(great_des_e1, bins=bins_e2)
            hist_e2,_,_ = np.histogram(great_des_e2, bins=bins_e2)
            hist_sn,_,_ = pl.histogram(res_sim_select['snr'] ,bins=bins_sn) 
            hist_sz,_,_ = pl.histogram(res_sim_select['mean_rgpp_rp'] ,bins=bins_sz) 

            row = np.concatenate([hist_e1,hist_e2,hist_sn,hist_sz])
            histmatrix[iu,:] = row
            import pdb; pdb.set_trace()
                          
            logger.info('------------------------------------- id=%5d id_cosmos=%6d n_gals=%4d' % (iu,id_cosmos,n_gals_sim))

        
        filename_table = os.path.join(args.output_dir,'cosmos_gals_histmatrix.cpickle')
        tt.save(filename_table,histmatrix,clobber=True)




def main():


    global logger , config , args

    description = 'Get statistics and plot results of noise bias calibration runs'
    parser = argparse.ArgumentParser(description=description, add_help=True)
    parser.add_argument('-v', '--verbosity', type=int, action='store', default=2, choices=(0, 1, 2, 3 ), help='integer verbosity level: min=0, max=3 [default=2]')
    parser.add_argument('-c', '--filename_config', default='sva1.yaml',type=str, action='store', help='name of the yaml config file')
    parser.add_argument('--use_calibration', default=False, action='store_true', help='if to apply calibration columns')
    parser.add_argument('--use_weights', default=False, action='store_true', help='if to apply weights')
    parser.add_argument('--fig_format', default='.png', type=str, action='store', help='format of the figure files')
    parser.add_argument('-o','--output_dir', default='.', type=str, action='store', help='dir to store results')
    parser.add_argument('-m', '--method', default='im3shape',type=str, action='store', help='im3shape or ngmix')

    args = parser.parse_args()
    logging_levels = { 0: logging.CRITICAL,1: logging.WARNING,2: logging.INFO,3: logging.DEBUG }
    logging_level = logging_levels[args.verbosity]; logger.setLevel(logging_level)  
    logger.info(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))

    config = yaml.load(open(args.filename_config))
    nbc_v7_select.config=config; nbc_v7_select.args = args; 
    nbc_v7_stats.config=config; nbc_v7_stats.args = args; 

    # get_bias_for_each_galaxy()
    # plot_bias_stats()
    get_histmatrix()

if __name__=='__main__':
    main()
