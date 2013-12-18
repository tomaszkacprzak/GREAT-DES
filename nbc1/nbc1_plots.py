import numpy, galsim, sys, logging, yaml, argparse, time, pylab, copy, itertools


dtype_table_cat     =  { 'names'   : ['index', 'filename_meds', 'n_gals' , 'ipsf_fwhm', 'isnr', 'ig', 'psf_fwhm', 'snr', 'g1' , 'g2'],
                         'formats' : ['i8'] + ['a40'] + ['i8']*3 + ['f8']*5 } 

dtype_table_bias     =  { 'names'   : ['index', 'ipsf_fwhm', 'isnr' ,  'nfail', 'psf_fwhm', 'snr', 'm1', 'm2', 'c1', 'c2', 'm1_std', 'm2_std', 'c1_std', 'c2_std', 'g_std' ],
                         'formats' : ['i8']*4 + ['f8']*11 } 

dtype_table_stats = { 'names'   : ['index','n_gals','n_fail','g1','g2','size','stdv_g1','stdv_g2','stdm_g1','stdm_g2','stdv_size','stdm_size'] ,
                     'formats' : ['i8'] *3+ ['f8']*9 } 

dtype_table_results =  { 'names'   : ['index','g1','g2','size','x0','y0'] ,
                         'formats' : ['i8'] + ['f8']*5 } 

dtype_table_results_im3shape = { 
        'names' : [ 'ID' ,  'catalogx' , 'catalogy' , 'Likelihood' , 'x0' , 'y0' , 'g1' , 'g2' , 'radius' , 'bulge_A' , 'disc_A' , 'flux_ratio' , 'size' , 'snr' , 'min_residuals' , 'max_residuals' , 'model_min' , 'model_max' , 'psf_e1' , 'psf_e2' , 'psf_fwhm' , 'psf_beta' , 'Rgpp_Rp' , 'levmar_e0' , 'levmar_e' , 'levmar_J' , 'levmar_Dp' , 'levmar_mu' , 'levmar_De' , 'levmar_nit' , 'levmar_reason' , 'levmar_neval' , 'levmar_njac' , 'levmar_nlin' , 'time_taken']  ,
        'formats' : ['i8'] + ['f4']*28 + ['i8']*5 + ['f4']*1 }
            

dtype_results = { 'hsm' : dtype_table_results , 'im3' : dtype_table_results_im3shape }

req1_dg = 0.003
req2_dg = 0.02
req3_dg = 0.1
req1_nfail = 0.01
req2_nfail = 0.05

def get_line_fit(x,y,sig):
        """
        @brief get linear least squares fit with uncertainity estimates
        y(X) = b*X + a
        see numerical recipies 15.2.9
        @param X    function arguments 
        @param y    function values
        @param sig  function values standard deviations   
        @return a - additive 
        @return b - multiplicative
        @return C - covariance matrix
        """
        
        import numpy

        invsig2 = sig**-2;
        
        S  = numpy.sum(invsig2)
        Sx = numpy.inner(x,invsig2)
        Sy = numpy.inner(y,invsig2)
        Sxx = numpy.inner(invsig2*x,x)
        Sxy = numpy.inner(invsig2*x,y)

        D = S*Sxx - Sx**2
        a = (Sxx*Sy - Sx*Sxy)/D
        b = (S*Sxy  - Sx*Sy)/D
        
        Cab = numpy.zeros((2,2))
        Cab[0,0] = Sxx/D
        Cab[1,1] = S/D
        Cab[1,0] = -Sx/D
        Cab[0,1] = Cab[1,0]
        
        return a,b,Cab


def surf_interp(arg1,arg2,value):

    grid_x, grid_y = numpy.mgrid[min(arg1):max(arg1):100j,min(arg2):max(arg2):100j]
    array_hlr  = numpy.array(arg1,ndmin=2).T
    array_snr  = numpy.array(arg2,ndmin=2).T
    points = numpy.concatenate((array_hlr,array_snr),axis=1)
    import scipy.interpolate
    grid = scipy.interpolate.griddata(points, value, (grid_x, grid_y), method='cubic')
    return grid.T


def plot_add_req(mult=1.):

    corner = pylab.xlim()[0]
    length = abs(pylab.xlim()[1]) + abs(pylab.xlim()[0])
    pylab.gca().add_patch(pylab.Rectangle(  (corner,-mult*req3_dg), length , 2*mult*req3_dg , facecolor = '0.9' , edgecolor='k' ))
    pylab.gca().add_patch(pylab.Rectangle(  (corner,-mult*req2_dg), length , 2*mult*req2_dg , facecolor = '0.8' , edgecolor='k' ))
    pylab.gca().add_patch(pylab.Rectangle(  (corner,-mult*req1_dg), length , 2*mult*req1_dg , facecolor = '0.7' , edgecolor='k' ))

def plot_add_nfail_req():

    corner = pylab.xlim()[0]
    length = abs(pylab.xlim()[1]) + abs(pylab.xlim()[0])
    pylab.gca().add_patch(pylab.Rectangle(  (corner,-req2_nfail), length , 2*req2_nfail , facecolor = '0.9' , edgecolor='k' ))
    pylab.gca().add_patch(pylab.Rectangle(  (corner,-req1_nfail), length , 2*req1_nfail , facecolor = '0.7' , edgecolor='k' ))

def write_stats(file_stats,stats):

    stats_line_fmt = '%d\t'*3 + '% .8f\t'*9 + '\n'
    line = stats_line_fmt % (
            stats['index'],
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

    use_dtype = dtype_results[args.method_id]
    res = numpy.loadtxt(filename_result,dtype=use_dtype)
    n_res = len(res)
    select_successful = numpy.abs(res['g1'] + 1j*res['g2'])<1
    n_success = sum(select_successful)
    n_fail = n_res - n_success
    res = res[select_successful]

    stats = {}
    
    if n_fail == n_res:
        stats['est_g1']  = 0.01
        stats['est_g2'] = 0.01
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
    try:
        pylab.hist(res['g1'],bins=bins_e,histtype='step',color='r',label='g1')
        pylab.hist(res['g2'],bins=bins_e,histtype='step',color='b',label='g2')
    except Exception,errmsg:
        logger.error('hist failed')
        pylab.hist(numpy.random.rand(1000)*2-1,bins=bins_e,histtype='step',color='r',label='g1')
        pylab.hist(numpy.random.rand(1000)*2-1,bins=bins_e,histtype='step',color='b',label='g2')
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

    stats_header = '# index n_gals n_fail g1 g2 size stdv_g1 stdv_g2 stdm_g1 stdm_g2 stdv_size stdm_size\n'
    file_stats.write(stats_header)

    for it,vt in enumerate(truth_cat):

        filename_result = '%s.%s.cat' % (vt['filename_meds'],args.method_id)
        logger.debug('getting stats for file %s' % filename_result)
        stats = get_shear(filename_result)
        stats['index'] = it
        write_stats(file_stats,stats)

    file_stats.close()
    logger.info('saved file %s' % filename_stats)

def write_bias_results(file_bias,bias_result=None,header=False):

    if header:
        header = '#  index ipsf_fwhm isnr n_fail psf_fwhm snr m1 m2 c1 c2 m1_std m2_std c1_std c2_std g_std \n'
        file_bias.write(header)
    else:
        fmt = '%d\t'*4 + '% 2.8f\t'*11 + '\n'
        line = fmt % (
            bias_result['index'] ,
            bias_result['ipsf_fwhm'],
            bias_result['isnr'] ,
            bias_result['n_fail'],
            bias_result['psf_fwhm'],
            bias_result['snr'],
            bias_result['m1'],
            bias_result['m2'],
            bias_result['c1'],
            bias_result['c2'],
            bias_result['m1_std'],
            bias_result['m2_std'],
            bias_result['c1_std'],
            bias_result['c2_std'],
            bias_result['g_std']
        )

        file_bias.write(line)


def get_mc():

    filename_stats = '%s.%s.stats.cat' % (args.filename_input,args.method_id)
    truth_cat = numpy.loadtxt(args.filename_input,dtype=dtype_table_cat)
    stats_cat = numpy.loadtxt(filename_stats,dtype=dtype_table_stats)
    n_gals_total = sum(truth_cat['n_gals'])
    logger.info('opened %s with %d rows and %d galaxies total' % (args.filename_input,len(truth_cat),n_gals_total))
    logger.info('opened %s with %d rows' % (filename_stats,len(stats_cat)))

    filename_bias = '%s.%s.bias.cat' % (args.filename_input,args.method_id)
    file_bias = open(filename_bias,'w')
    write_bias_results(file_bias,header=True)

    iall = 0
    
    for ipsf_fwhm,vpsf_fwhm in enumerate(config['grid']['psf_fwhm']):

        for isnr,vsnr in enumerate(config['grid']['snr']):

            select = numpy.array(truth_cat['ipsf_fwhm'] == ipsf_fwhm) * numpy.array(truth_cat['isnr'] == isnr)
            current_stats_cat = stats_cat[select]
            current_truth_cat = truth_cat[select]

            g1_true = current_truth_cat['g1']
            g2_true = current_truth_cat['g2']
            g1_bias = (current_stats_cat['g1'] - g1_true)
            g2_bias = (current_stats_cat['g2'] - g2_true)
            g1_err  = stats_cat[select]['stdm_g1']
            g2_err  = stats_cat[select]['stdm_g2']

            [c1,m1,C1cm] = get_line_fit(g1_true,g1_bias,g1_err)
            [c2,m2,C2cm] = get_line_fit(g2_true,g2_bias,g2_err)
            m1_std = numpy.sqrt(C1cm[1,1])
            m2_std = numpy.sqrt(C2cm[1,1])
            c1_std = numpy.sqrt(C1cm[0,0])
            c2_std = numpy.sqrt(C2cm[0,0])

            g_std = numpy.mean(current_stats_cat['stdv_g1'])
            n_fail = float(sum(current_stats_cat['n_fail'])) / float(sum(current_stats_cat['n_gals']))

            pylab.figure()
            pylab.errorbar(g1_true,g1_bias,yerr=g1_err,fmt='b.',label='m1=%2.3f +/- %2.3f' % (m1,m1_std))
            pylab.plot(g1_true,g1_true*m1 + c1,'b-')
            pylab.errorbar(g2_true,g2_bias,yerr=g2_err,fmt='r.',label='m1=%2.3f +/- %2.3f' % (m2,m2_std))
            pylab.plot(g2_true,g2_true*m2 + c2,'r-')
            pylab.title('snr=%2.2f psf_fwhm=%2.2f' % (vsnr,vpsf_fwhm))
            pylab.legend()
            filename_fig = 'fig.bias.snr%03d.psf%03d.png' % (isnr,ipsf_fwhm)
            pylab.savefig(filename_fig)
            pylab.close()
            logger.info('saved %s' % filename_fig)

            bias_result = {}
            bias_result['index']  = iall
            bias_result['ipsf_fwhm'] = ipsf_fwhm
            bias_result['isnr']  = isnr
            bias_result['psf_fwhm'] = vpsf_fwhm
            bias_result['snr'] = vsnr
            bias_result['m1'] = m1
            bias_result['m2'] = m2
            bias_result['c1'] = c1
            bias_result['c2'] = c2
            bias_result['m1_std'] = m1_std
            bias_result['m2_std'] = m2_std
            bias_result['c1_std'] = c1_std
            bias_result['c2_std'] = c2_std
            bias_result['g_std'] = g_std
            bias_result['n_fail'] = n_fail
            logger.info('psf_fwhm=%2.2f\tsnr=%2.2f\tm=% 2.2e\tc=% 2.2f\tm_err=%2.2e\tc_err=%2.2e\tg_std=%2.2f\tn_fail=%2.2f' % (vpsf_fwhm, vsnr,m1,c1,m1_std,c1_std,g_std,n_fail))                

            write_bias_results(file_bias, bias_result=bias_result)

            iall+=1

def plot_mc():

        filename_bias = '%s.%s.bias.cat' % (args.filename_input,args.method_id)
        results_mc = numpy.loadtxt(filename_bias,dtype=dtype_table_bias)  



        # snr = config['grid']['snr']
        # for ig in range(len(snr)):
        #     logger.info('%d e_tru=(% 0.3f,% 0.3f) , e_est=(% 0.3f,% 0.3f) , bias_g1=(% 0.3f,% 0.3f)' % (
        #         ig, e_true.real, e_true.imag, 
        #         current_stats_cat['g1'][ig], current_stats_cat['g2'][ig], 
        #         bias_g1[ig], bias_g2[ig]))           

        psf_fwhm_colors=['r','g','b']
        
        for ipsf_fwhm,vpsf_fwhm in enumerate(config['grid']['psf_fwhm']):
            # for isnr,vsnr in enumerate(config['grid']['snr']):

            select = (results_mc['ipsf_fwhm'] == ipsf_fwhm)
            results_select = results_mc[select]
            sort = numpy.argsort(results_select['snr'])
            results_select=results_select[sort] 

            pylab.figure(1)
            pylab.errorbar(results_select['snr'],results_select['m1'],yerr=results_select['m1_std'],fmt='r+--')
            pylab.errorbar(results_select['snr'],results_select['m2'],yerr=results_select['m2_std'],fmt='rx:')
            xlim_add = (max(results_select['snr']) - min(results_select['snr']))*0.1
            pylab.xlim([min(results_select['snr'])-xlim_add,max(results_select['snr'])+xlim_add])
            plot_add_req()
            pylab.xlabel('SNR')
            pylab.ylabel('multiplicative bias m')
            pylab.xticks(results_select['snr'])
            title_str = 'PSF FWHM=%2.2f ' % (vpsf_fwhm)
            pylab.title(title_str)
            pylab.legend(['%s g1' % args.method_id,'%s g2' % args.method_id],loc='lower left',ncol=2,mode='expand')
            filename_fig = 'fig.%s.psf_fwhm%d.png' % (filename_bias,ipsf_fwhm)
            pylab.savefig(filename_fig)
            logger.info('saved %s' % filename_fig)
            pylab.close()

            pylab.figure(2)
            pylab.errorbar(results_select['snr'],results_select['m1'],yerr=results_select['m1_std'],fmt=psf_fwhm_colors[ipsf_fwhm]+'--',label='psf_fwhm=%2.2f m1' % vpsf_fwhm)
            pylab.errorbar(results_select['snr'],results_select['m2'],yerr=results_select['m2_std'],fmt=psf_fwhm_colors[ipsf_fwhm]+':',label='psf_fwhm=%2.2f m2' % vpsf_fwhm)

        pylab.figure(2)
        xlim_add = (max(results_select['snr']) - min(results_select['snr']))*0.1
        pylab.xlim([min(results_select['snr'])-xlim_add,max(results_select['snr'])+xlim_add])
        plot_add_req()
        pylab.xlabel('SNR')
        pylab.ylabel('multiplicative bias m')
        pylab.xticks(results_select['snr'])
        pylab.legend(loc='lower left',ncol=2,mode='expand')
        filename_fig = 'fig.%s.png' % (filename_bias)
        pylab.savefig(filename_fig)
        logger.info('saved %s' % filename_fig)
        pylab.close()

        # pylab.figure()
        # pylab.plot(snr[sort],n_fail[sort],'ro-')
        # xlim_add = (max(snr) - min(snr))*0.1
        # pylab.xlim([min(snr)-xlim_add,max(snr)+xlim_add])
        # plot_add_nfail_req()
        # pylab.xlabel('SNR')
        # pylab.ylabel('n_fail [%]')
        # pylab.xticks(snr)
        # title_str = 'hlr=%2.2f' % vhlr
        # pylab.title(title_str)
        # filename_fig = 'fig.%s.hlr%d.nfail.png' % (filename_stats,ihlr)
        # pylab.savefig(filename_fig)
        # logger.info('saved %s' % filename_fig)
        # pylab.close()

    # # plot the surface

    # angrad = float(config['gal']['shear']['beta'].split()[0])/180.*numpy.pi
    # e_true = config['gal']['shear']['g'] * numpy.exp(1j*angrad*2.)   # angle of galaxy orientation, not complex shear angle

    # n_grid = 100
    # pylab.figure()
    # bias = (stats_cat['g1'] - e_true.real)/e_true.real
    # grid=surf_interp(truth_cat['hlr'],truth_cat['snr'],bias)
    # pylab.imshow(grid, origin='lower',extent=(min(truth_cat['hlr']),max(truth_cat['hlr']),min(truth_cat['snr']),max(truth_cat['snr'])), aspect='auto')
    # pylab.ylabel('SNR')
    # pylab.xlabel('hlr [arcsec]')
    # title_str = '%s' % args.method_id
    # pylab.title(title_str)
    # filename_fig = 'fig.%s.hlr_snr.%s.png' % (filename_stats,args.method_id)
    # pylab.colorbar()
    # pylab.savefig(filename_fig)
    # logger.info('saved %s' % filename_fig)
    # pylab.close()


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


    # get_stats()
    # get_mc()
    plot_mc()

    logger.info(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))


main()