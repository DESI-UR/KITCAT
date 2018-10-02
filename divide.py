""" Script for submitting jobs to calculate DD(s), DR(s), and RR(s) """

# Standard Python modules
import time
import argparse

# Python modules
import numpy

# User-defined modules
import lib.io as lio
import lib.helper as lhelper
import lib.analysis as lanalysis

if __name__ == "__main__":
    """ Main"""
    print('')

    def parse_command_line():
        parser = argparse.ArgumentParser(description='preprocess data')
        parser.add_argument('--prefix', '-p',
                            help    = 'output prefix.',
                            dest    = 'prefix',
                            type    = str)
        parser.add_argument('--ijob', '-i',
                            help    = 'job index',
                            default = 0,
                            dest    = 'ijob',
                            type    = int)
        parser.add_argument('--njob', '-n',
                            help    = 'total number of z-slice',
                            default = 1,
                            dest    = 'njob',
                            type    = int)
        parser.add_argument('--time', '-t',
                            help='save runtime',
                            action='store_true',
                            default=False)
        params = parser.parse_args()
        return params

    args = parse_command_line()

    # set job helper
    job_helper = lhelper.JobHelper(args.njob)
    job_helper.set_current_job(args.ijob)

    # set timer
    time_dd = 0
    time_dr = 0
    time_rr = 0

    # load preprocess data
    preprocess_params = lio.load('%s_preprocess.pkl' % args.prefix)
    rr_tree = preprocess_params['rr']['tree']
    rr_catalog = preprocess_params['rr']['catalog']
    dd_tree = preprocess_params['dd']['tree']
    dd_catalog = preprocess_params['dd']['catalog']
    dr_tree = preprocess_params['dr']['tree']
    dr_catalog = preprocess_params['dr']['catalog']
    bins = preprocess_params['bins']

    # calculate f(theta)
    start_time = time.time()
    ftheta = lanalysis.get_ftheta(
        catalog = rr_catalog,
        tree    = rr_tree,
        theta_max   = bins.max('theta'),
        theta_nbins = bins.num_bins('theta'),
        job_helper  = job_helper,
        same        = True)
    time_rr = time.time()-start_time
    print("--- %f seconds ---" % time_rr)


    # save("%s_divide_%03d.pkl" % (args.prefix, args.ijob), correlation)
