import numpy
import argparse
import logging
import sys

bin = '/project/projectdirs/des/wl/desdata/users/cogs/ucl_des_shear_dev/utils/launch_im3shape_GREATDES2.sh'


def genCommands():

    filename_commands = 'commands.%d.%d.sh' % (args.start,args.count)
    file_commands = open(filename_commands,'w')

    filename_list = 'test2.filelist.cat'
    file_list = open(filename_list)

    n_meds = 0
    for f in file_list:

        f= f.replace('\n','')
        file_in =   args.data_dir + '/' + f
        file_out = "%s/%s.%d.%d.im3.cat" % (args.output_dir,f,args.start,args.count)
        cmd = "%s %s %d %d %s\n" % (bin, file_in, args.start, args.count, file_out)

        file_commands.write(cmd)

        n_meds += 1

    file_commands.close()
    logger.info('created %s' % filename_commands)

    time_per_galaxy = 5.   # sec
    number_galaxies = float(n_meds * args.count)
    n_nodes = float(args.n_nodes)

    wallclock_hours = time_per_galaxy * number_galaxies / 60. / 60. / n_nodes

    print 'set wallclock to %f hours with %d nodes' % (wallclock_hours , n_nodes)






        # /project/projectdirs/des/wl/desdata/users/cogs/ucl_des_shear_dev/utils/launch_im3shape_GREATDES2.sh test2.hlr00.snr00.shear00.part00.meds.fits.fz 700  800  /global/homes/t/tomaszk/projects/130821_minions_test/output/test2.hlr00.snr00.shear00.part00.meds.fits.fz.7.cat

def main():

    global logger , config , args , logger_config

    description = 'Plots for the results of nmb_main. To use in ipython, create a variable global results_array, global truth_array to speed up loading'
    parser = argparse.ArgumentParser(description=description, add_help=True)
    parser.add_argument('-v', '--verbosity', type=int, action='store', default=2, choices=(0, 1, 2, 3 ), help='integer verbosity level: min=0, max=3 [default=2]')
    parser.add_argument('-s', '--start', type=int, default=0, action='store', help='integer starting galaxy index')
    parser.add_argument('-c', '--count', type=int, default=10, action='store', help='integer number of galaxies to process')
    parser.add_argument('-d', '--data_dir', default='/global/project/projectdirs/des/wl/desdata/users/cogs/GREAT-DES/test2/meds/',type=str, action='store', help='name of the dir for output catalog')
    parser.add_argument('-o', '--output_dir', default='/global/project/projectdirs/des/wl/desdata/users/cogs/GREAT-DES/test2/im3shape',type=str, action='store', help='name of the dir for output catalog')
    parser.add_argument('-n', '--n_nodes', default=192,type=int, action='store', help='number of nodes, used for calculating wallclock time')

    
    args = parser.parse_args()
    # Parse the integer verbosity level from the command line args into a logging_level string
    logging_levels = { 0: logging.CRITICAL, 
                       1: logging.WARNING,
                       2: logging.INFO,
                       3: logging.DEBUG }
    logging_level = logging_levels[args.verbosity]
    logging.basicConfig(format="%(message)s", level=logging_level, stream=sys.stdout)
    logger = logging.getLogger("generate_catalog.py") 
    logger.setLevel(logging_level)

    genCommands()



main()