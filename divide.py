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
                            help    = 'save runtime',
                            action  = 'store_true',
                            default = False)
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
    rr_params = preprocess_params['rr']
    dd_params = preprocess_params['dd']
    dr_params = preprocess_params['dr']
    bins = preprocess_params['bins']
    cosmos_list = preprocess_params['cosmos_list']
    helper = preprocess_params['helper']

    # calculate f(theta)
    print('')
    start_time = time.time()
    ftheta = None
    ftheta = lanalysis.get_ftheta(
        catalog     = rr_params['catalog'],
        tree        = rr_params['tree'],
        theta_max   = bins.max('theta'),
        theta_nbins = bins.num_bins('theta'),
        job_helper  = job_helper,
        same        = True)
    time_rr = time.time()-start_time
    print("--- %f seconds ---" % time_rr)

    # calculate ztheta
    print('')
    start_time = time.time()
    ztheta = None
    ztheta = lanalysis.get_ztheta(
        tree_catalog    = dr_params['tree_catalog'],
        pair_catalog    = dr_params['pair_catalog'],
        tree            = dr_params['tree'],
        z_min           = bins.min('z'),
        z_max           = bins.max('z'),
        z_nbins         = bins.num_bins('z'),
        theta_max       = bins.max('theta'),
        theta_nbins     = bins.num_bins('theta'),
        job_helper      = job_helper,)
    time_dr = time.time()-start_time
    print("--- %f seconds ---" % time_dr)

    # calculate zztheta
    print('')
    start_time = time.time()
    zztheta = None
    zztheta = lanalysis.get_zztheta(
        catalog         = dd_params['catalog'],
        tree            = dd_params['tree'],
        z_min           = bins.min('z'),
        z_max           = bins.max('z'),
        z_nbins         = bins.num_bins('z'),
        theta_max       = bins.max('theta'),
        theta_nbins     = bins.num_bins('theta'),
        job_helper      = job_helper,
        same            = True,)
    time_dd = time.time()-start_time
    print("--- %f seconds ---" % time_dd)

    # store to helper object
    helper.ftheta = ftheta
    helper.ztheta = ztheta
    helper.zztheta = zztheta
    if args.ijob == 0:
        helper.cosmos_list = cosmos_list
        helper.bins = bins
    lio.save("%s_divide_%03d-%03d.pkl" % (args.prefix, args.ijob, args.njob),
             helper)
