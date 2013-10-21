import numpy, galsim, sys, logging, yaml, argparse, time, pylab, copy, itertools

dtype_table_results =  { 'names'   : ['index','g1','g2','size','x0','y0'] ,
                         'formats' : ['i8'] + ['f8']*5 } 

dtype_table_cat =  { 'names'   : ['index', 'filename_meds', 'n_gals' , 'ihlr', 'isnr', 'iellip', 'iangdeg', 'insersic', 'imoffat_beta', 'imoffat_g', 'imoffat_fwhm', 'hlr', 'snr', 'ellip', 'angdeg', 'nsersic', 'moffat_beta', 'moffat_g', 'moffat_fwhm'],
                     'formats' : ['i8'] + ['a40'] + ['i8']*9 + ['f8']*11 } 

dtype_table_stats = { 'names'   : ['n_gals','n_fail','g1','g2','size','stdv_g1','stdv_g2','stdm_g1','stdm_g2','stdv_size','stdm_size'] ,
                     'formats' : ['i8'] *2+ ['f8']*9 } 


def write_stats(file_stats,stats):

    stats_line_fmt = '%d\t'*2 + '% .8f\t'*9 + '\n'
    line = stats_line_fmt % (
            stats['n_gals'],
            stats['n_fail'],
            stats['est_g1'],
            stats['est_g2'],
            stats['est_size'],
            stats['est_stdv_g1'],
            stats['est_stdv_g2'],
            stats['est_stdm_g1'],
            stats['est_stdm_g2'],
            stats['est_stdv_size'],
            stats['est_stdm_size']
        )
    file_stats.write(line)


def get_shear(filename_result):

    res = numpy.loadtxt(filename_result,dtype=dtype_table_results)
    n_res = len(res)
    # select_successful = res['flag']==0
    select_successful = numpy.abs(res['g1'] + 1j*res['g2'])<1
    n_success = sum(select_successful)
    n_fail = n_res - n_success
    res = res[select_successful]

    stats = {}
    
    if n_fail == n_res:
        stats['est_g1']  = 0.
        stats['est_g2'] = 0.
        stats['est_size'] = 1.
        stats['est_stdv_g1'] = 0.1
        stats['est_stdv_g2'] = 0.1
        stats['est_stdm_g1'] = 0.1
        stats['est_stdm_g2'] = 0.1
        stats['est_stdv_size'] = 0.1
        stats['est_stdm_size'] = 0.1
        stats['n_fail'] = n_res 
        stats['n_gals'] = n_res
        logger.error('file %s - all measurements are error' % filename_result)
    else:
        stats['est_g1'] = numpy.mean(res['g1'])
        stats['est_g2'] = numpy.mean(res['g2'])
        stats['est_size'] = numpy.mean(res['size'])
        stats['est_stdv_g1'] = numpy.std(res['g1'],ddof=1)
        stats['est_stdv_g2'] = numpy.std(res['g2'],ddof=1)
        stats['est_stdm_g1'] = numpy.std(res['g1'],ddof=1)/numpy.sqrt(n_success)
        stats['est_stdm_g2'] = numpy.std(res['g2'],ddof=1)/numpy.sqrt(n_success)
        stats['est_stdv_size'] = numpy.std(res['size'],ddof=1)
        stats['est_stdm_size'] = numpy.std(res['size'],ddof=1)/numpy.sqrt(n_success)
        stats['n_fail'] = n_fail
        stats['n_gals'] = n_res

    # save histograms

    n_bins = 100
    bins_e = numpy.linspace(-1,1,n_bins)
    pylab.hist(res['g1'],bins=bins_e,histtype='step',color='r',label='g1')
    pylab.hist(res['g2'],bins=bins_e,histtype='step',color='b',label='g2')
    pylab.legend()
    pylab.title('%s ellipticity' % filename_result)
    pylab.xlabel('ellip [g]')
    pylab.ylabel('histogram')
    filename_fig = 'fig.ellip.%s.png' % filename_result
    pylab.savefig(filename_fig)
    logger.info('saved %s' % filename_fig)
    pylab.close()

    return stats

def get_stats():

    truth_cat = numpy.loadtxt(args.filename_input,dtype=dtype_table_cat)

    filename_stats = '%s.%s.stats.cat' % (args.filename_input,args.method_id)
    file_stats = open(filename_stats,'w')
    logger.debug('opened file %s' % filename_stats)

    stats_header = '# n_gals n_fail g1 g2 size stdv_g1 stdv_g2 stdm_g1 stdm_g2 stdv_size stdm_size\n'
    file_stats.write(stats_header)

    for it,vt in enumerate(truth_cat):

        filename_result = '%s.%s.cat' % (vt['filename_meds'],args.method_id)
        logger.debug('getting stats for file %s' % filename_result)
        stats = get_shear(filename_result)
        write_stats(file_stats,stats)

    file_stats.close()
    logger.info('saved file %s' % filename_stats)

def plotsO1():

    filename_stats = '%s.%s.stats.cat' % (args.filename_input,args.method_id)
    truth_cat = numpy.loadtxt(args.filename_input,dtype=dtype_table_cat)
    stats_cat = numpy.loadtxt(filename_stats,dtype=dtype_table_stats)
    logger.debug('opened %s with %d rows' % (args.filename_input,len(truth_cat)))
    logger.debug('opened %s with %d rows' % (filename_stats,len(stats_cat)))

    
    # order 1 plots
    for ihlr,vhlr in enumerate(config['order1']['hlr']):


        select = truth_cat['ihlr'] == ihlr
        snr = truth_cat[select]['snr']
        angrad = truth_cat[select]['angdeg']/180.*numpy.pi
        e_true = truth_cat[select]['ellip'] * numpy.exp(1j*angrad*2.)   # angle of galaxy orientation, not complex shear angle
        bias_g1 = (stats_cat[select]['g1'] - e_true.real)/e_true.real
        bias_g2 = (stats_cat[select]['g2'] - e_true.imag)/e_true.imag

        print e_true, bias_g1, bias_g2

        g1_err  = stats_cat[select]['stdm_g1']/e_true.real
        g2_err  = stats_cat[select]['stdm_g2']/e_true.imag

        pylab.errorbar(snr,bias_g1,yerr=g1_err,fmt='r+')
        pylab.errorbar(snr,bias_g2,yerr=g2_err,fmt='rx')
        pylab.xlim([0.8*min(snr),1.1*max(snr)])
        pylab.xlabel('SNR')
        pylab.ylabel('de/e')
        title_str = 'hlr=%2.2f' % vhlr
        pylab.title(title_str)
        filename_fig = 'fig.%s.hlr%d.eps' % (filename_stats,ihlr)
        pylab.savefig(filename_fig)
        logger.info('saved %s' % filename_fig)
        pylab.close()

def plotsO2():

    filename_stats = '%s.%s.stats.cat' % (args.filename_input,args.method_id)
    truth_cat = numpy.loadtxt(args.filename_input,dtype=dtype_table_cat)
    stats_cat = numpy.loadtxt(filename_stats,dtype=dtype_table_stats)
    logger.debug('opened %s with %d rows' % (args.filename_input,len(truth_cat)))
    logger.debug('opened %s with %d rows' % (filename_stats,len(stats_cat)))
    params = config['order2']['deviations'].keys()
    
    # order 1 plots
    for ihlr_snr,vhlr_snr in enumerate(config['order2']['hlr_snr']):

        for iparam,vparam in enumerate(params):

            select = truth_cat['ihlr'] == ihlr_snr  
            hlr,snr = vhlr_snr[0],vhlr_snr[1]

            current_truth_cat = truth_cat[select]
            current_stats_cat = stats_cat[select]
            
            fid_params = copy.deepcopy(params)
            fid_params.remove(vparam)

            select = numpy.ones(len(current_truth_cat),dtype=bool)
            # [ print ( current_truth_cat['i'+p] == config['order2']['deviations'][vparam]['fid']  ) for p in fid_params ]
            for p in fid_params:
                select *= (current_truth_cat['i'+p] == config['order2']['deviations'][p]['fid'])

            current_stats_cat = current_stats_cat[select]
            current_truth_cat = current_truth_cat[select]

            angrad = current_truth_cat['angdeg']/180.*numpy.pi
            e_true = current_truth_cat['ellip'] * numpy.exp(1j*angrad*2.)   # angle of galaxy orientation, not complex shear angle
            bias_g1 = (current_stats_cat['g1'] - e_true.real)/e_true.real
            bias_g2 = (current_stats_cat['g2'] - e_true.imag)/e_true.imag
            g1_err  = current_stats_cat['stdm_g1']/e_true.real
            g2_err  = current_stats_cat['stdm_g2']/e_true.imag

            # bias_g1 = (current_stats_cat['g1'] - e_true.real)
            # bias_g2 = (current_stats_cat['g2'] - e_true.imag)
            # g1_err  = current_stats_cat['stdm_g1']
            # g2_err  = current_stats_cat['stdm_g2']

            param_grid = current_truth_cat[vparam]

            pylab.errorbar(param_grid,bias_g1,yerr=g1_err,fmt='b+')
            pylab.errorbar(param_grid,bias_g2,yerr=g2_err,fmt='bx')
            pylab.xlim([0.8*min(param_grid),1.1*max(param_grid)])
            pylab.xlabel(vparam)
            pylab.ylabel('de/e')
            title_str = 'hlr=%2.2f snr=%2.2f' % (hlr,snr)
            pylab.title(title_str)

            filename_fig = 'fig.%s.%d.%s.eps' % (filename_stats,ihlr_snr,vparam)
            pylab.savefig(filename_fig)
            logger.info('saved %s' % filename_fig)
            pylab.close()




def main():

    import numpy, galsim, sys, logging, yaml, argparse, time

    global logger , config , args

    description = 'Run on sha1 data with HSM'
    parser = argparse.ArgumentParser(description=description, add_help=True)
    parser.add_argument('-v', '--verbosity', type=int, action='store', default=2, choices=(0, 1, 2, 3 ), help='integer verbosity level: min=0, max=3 [default=2]')
    parser.add_argument('-i', '--filename_input', default='sha1-O1.cat',type=str, action='store', help='name of the output catalog')
    parser.add_argument('-c', '--filename_config', default='sha1.yaml',type=str, action='store', help='name of the yaml config file')
    parser.add_argument('-m', '--method_id', default='hsm',type=str, action='store', help='name of the yaml config file')
    
    args = parser.parse_args()
    # Parse the integer verbosity level from the command line args into a logging_level string
    logging_levels = { 0: logging.CRITICAL, 
                       1: logging.WARNING,
                       2: logging.INFO,
                       3: logging.DEBUG }
    logging_level = logging_levels[args.verbosity]
    logging.basicConfig(format="%(message)s", level=logging_level, stream=sys.stdout)
    logger = logging.getLogger("run_hsm_sha1") 
    logger.setLevel(logging_level)

    logger.info(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))

    config = yaml.load(open(args.filename_config))


    get_stats()
    if 'O1' in args.filename_input:
        plotsO1()
    if 'O2' in args.filename_input:
        plotsO2()

    logger.info(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))


main()