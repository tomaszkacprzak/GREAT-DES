import argparse
import logging
import sys
import numpy
import meds
import pylab
import pyfits

n_bins = 100


def saveHistograms():

    medsdata = meds.MEDS(args.filename_meds)
    results = numpy.loadtxt(args.filename_results)

    logger.info('opened meds file %s with %d galaxies' % (args.filename_meds,medsdata.size))
    logger.info('opened results file %s with %d galaxies' % (args.filename_meds,len(results)))

    e1 = results[:,3]
    e2 = results[:,4]

    pylab.hist(e1,bins=n_bins)
    pylab.show()


def saveResidualPlots():

    medsdata = meds.MEDS(args.filename_meds)
    results = numpy.loadtxt(args.filename_results)

    logger.info('opened meds file %s with %d galaxies' % (args.filename_meds,medsdata.size))
    logger.info('opened results file %s with %d galaxies' % (args.filename_meds,len(results)))

    for i in range(args.n_residuals):

        # load image, weight
        filename_image  = 'image.%d.fits' % i
        filename_model  = 'model.%d.fits' % i
        filename_weight = 'weight.%d.fits' % i
        gal_image = pyfits.getdata(filename_image)
        gal_model = pyfits.getdata(filename_model)
        gal_weight = pyfits.getdata(filename_weight)
        gal_residual = gal_model - gal_image

        subplot_nx = 2
        subplot_ny = 2
        subplot_ni  = 0

        pylab.figure()

        subplot_ni+=1; pylab.subplot(subplot_nx,subplot_ny,subplot_ni)
        pylab.imshow(gal_image,interpolation='nearest')


        subplot_ni+=1; pylab.subplot(subplot_nx,subplot_ny,subplot_ni)
        pylab.imshow(gal_model,interpolation='nearest')

        subplot_ni+=1; pylab.subplot(subplot_nx,subplot_ny,subplot_ni)
        pylab.imshow(gal_weight,interpolation='nearest')

        subplot_ni+=1; pylab.subplot(subplot_nx,subplot_ny,subplot_ni)
        pylab.imshow(gal_residual,interpolation='nearest')

        filename_figure = 'figs/fig.%03d.png' % i
        pylab.savefig(filename_figure)
        logger.info('saved figure %s' % filename_figure)

        pylab.close()





def main():

    global logger , config , args

    description = 'make residual plots'
    parser = argparse.ArgumentParser(description=description, add_help=True)
    parser.add_argument('-v', '--verbosity', type=int, action='store', default=2, choices=(0, 1, 2, 3 ), help='integer verbosity level: min=0, max=3 [default=2]')
    parser.add_argument('-r', '--filename_results', default='results.cat',type=str, action='store', help='name of the im3shape output catalog')
    parser.add_argument('-m', '--filename_meds', default='meds.fits',type=str, action='store', help='name of the meds file with the data')
    parser.add_argument('-n', '--n_residuals',default=100,type=int, action='store', help='number of images to produce')

    # parser.add_argument('-d', '--dry', default=False,  action='store_true', help='Dry run, dont generate data')

    args = parser.parse_args()
    # Parse the integer verbosity level from the command line args into a logging_level string
    logging_levels = { 0: logging.CRITICAL, 
                       1: logging.WARNING,
                       2: logging.INFO,
                       3: logging.DEBUG }
    logging_level = logging_levels[args.verbosity]
    logging.basicConfig(format="%(message)s", level=logging_level, stream=sys.stdout)
    logger = logging.getLogger("results_residuals") 
    logger.setLevel(logging_level)

    saveResidualPlots()


main()






