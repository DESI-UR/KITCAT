""" Script for submitting jobs to calculate DD(s), DR(s), and RR(s) """

# Standard Python modules
import time
import argparse

# Python modules
import numpy

# User-defined modules
from lib.myio import save, load
from lib.helper import JobHelper


def main():
    """ Main
    Cmd arguments are job number, total number of jobs, configuration file,
    and output prefix
      + --ijob: Job number must be an integer from 0 to total_job-1.
      + --njob: Total number of jobs that will be submitted.
      + --config: Setting for correlation function. See below.
      + --prefix: Prefix of the output files (can include folder path, folder
      must exist).
    If job number is 0, will also save comoving distribution P(r) and
    normalization factor."""
    print("DIVIDE module")

    # Read in cmd argument
    parser = argparse.ArgumentParser(
        description='Divide calculation into multiple parts.')
    parser.add_argument('-p', '-P', '--prefix', type=str, help='Output prefix.')
    parser.add_argument('-i', '-I', '--ijob', type=int, default=0,
                        help='Index of job. From 0 to N-1.')
    parser.add_argument('-n', '-N', '--njob', type=int, default=1,
                        help='Total number of jobs.')
    parser.add_argument('-t', '-T', '--time', action='store_true',
                        default=False, help='Enable save runtime.')
    parser.add_argument('--version', action='version', version='KITCAT 1.10')
    args = parser.parse_args()

    # Set job helper
    job_helper = JobHelper(args.njob)
    job_helper.set_current_job(args.ijob)

    # Set timer
    timer = {'dd': None, 'dr': None, 'rr': None}

    # Load pickle file
    loader = next(load('{}_preprocess.pkl'.format(args.prefix)))

    # Set correlation helper
    correlation = loader['helper']

    # Parallelize DD if only one cosmology model is given
    if len(correlation.models) == 1:
        # Keeping track of time
        start_time = time.time()

        # Import data and calculate DD(s)
        data = loader['dd']
        tree, catalog = data.build_balltree(
            'euclidean', model=correlation.models[0], return_catalog=True)
        correlation.set_dd(tree, catalog, job_helper)

        # Save and print out time
        timer['dd'] = time.time()-start_time
        print("--- %f seconds ---" % (timer['dd']))

    # Calculate f(theta)
    # Keeping track of time
    start_time = time.time()

    # Import data and calculate f(theta)
    tree = loader['rr']['tree']
    catalog = loader['rr']['catalog']
    correlation.set_theta_distr(tree, catalog, job_helper)

    # Save and print out time
    timer['rr'] = time.time()-start_time
    print("--- %f seconds ---" % (timer['rr']))

    # Calculate g(theta, rz)
    # Keeping track of time
    start_time = time.time()

    # Import data and calculate g(theta, rz)
    tree = loader['dr']['tree']
    data_catalog = loader['dr']['data_catalog']
    angular_catalog = loader['dr']['angular_catalog']
    mode = loader['dr']['mode']
    correlation.set_z_theta_distr(tree, data_catalog, angular_catalog,
                                  mode, job_helper)

    # Save and print out time
    timer['dr'] = time.time()-start_time
    print("--- %f seconds ---" % (timer['dr']))

    # Save object and timer
    if args.time:
        runtime = [timer['rr'], timer['dr'], timer['dd']]
        numpy.savetxt("%s_timer.txt" % args.prefix, runtime,
                      header='RR(s)  DR(s)  DD(s)')

    save("%s_divide_%03d.pkl" % (args.prefix, args.ijob), correlation)


if __name__ == "__main__":
    main()
