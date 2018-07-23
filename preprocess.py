""" Module to create preprocessed catalog object"""

# Standard Python module
import os
import configparser
import argparse

# User-defined module
from lib.myio import save, read_cosmology
from lib.catalog import GalaxyCatalog
from lib.helper import CorrelationHelper
from lib.bins import Bins

def main():
    """ Convert .fits catalogs and other binning, cosmology information into
    data structures. """

    # Read in cmd argument
    parser = argparse.ArgumentParser(
        description='Preprocess galaxy and random catalogs.')
    parser.add_argument('-c', '-C', '--config', type=str,
                        help='Path to configuration file.')
    parser.add_argument(
        '-p', '-P', '--prefix', type=str,
        help='Output prefix. Output is saved as PREFIX_preprocess.pkl')
    parser.add_argument('-i', '-I', '--islice', type=int, default=0,
                        help='Index of Z-slice. From 0 to N-1.')
    parser.add_argument('-n', '-N', '--nslice', type=int, default=1,
                        help='Total number of Z-slices.')
    parser.add_argument('-a', '-A', '--auto', action='store_true',
                        default=False, help='Set automatic binning')
    parser.add_argument('-b', '-B', '--binwidth', type=float, default=2.00,
                        help='Separation binwidth. Disable if auto is False.')
    parser.add_argument('--version', action='version', version='KITCAT 1.10')
    args = parser.parse_args()

    print('PREPROCESS module')

    # Check if arguments are valid
    if args.islice < 0 or args.islice >= args.nslice:
        raise ValueError('islice must be at least 0 and less than nslice.')

    # Read the provided configuration file
    print('- Reading configuration file from %s' % args.config)
    config = configparser.RawConfigParser(os.environ)
    config.read(args.config)

    # Setting up cosmology models
    print(' - Setting up cosmology models')
    models = read_cosmology(config['COSMOLOGY'])

    # Setting up binning
    print(' - Setting up binning')
    nbins_params = None if args.auto else config['NBINS']
    auto = args.binwidth if args.auto else None
    bins = Bins(config['LIMIT'], models,
                nbins=nbins_params,
                islice=args.islice, nslice=args.nslice,
                auto=auto)


    # Initialize catalog and save dictionary
    print('- Initialize catalog')
    gal_data = GalaxyCatalog(config['GAL_DATA'], bins.limit) # data catalog
    gal_rand = GalaxyCatalog(config['GAL_RAND'], bins.limit) # random catalog
    gal_rand = gal_rand.to_rand(bins.limit, bins.nbins)
    if config['VOID_DATA'] is not None and config['VOID_RAND'] is not None:
        void_data = GalaxyCatalog(config['VOID_DATA'], bins.limit)
        void_rand = GalaxyCatalog(config['VOID_RAND'], bins.limit)
        void_rand = void_rand.to_rand(bins.limit, bins.nbins)
    gal_save_dict = {'dd': None, 'dr': None, 'rr': None, 'helper': None}


    # Create a kd-tree for DD calculation and pickle
    print('- Setting up DD')
    gal_save_dict['dd'] = gal_data

    # Create a balltree for f(theta), RR calculation and pickle
    print('- Setting up RR')
    tree, catalog = gal_rand.build_balltree(return_catalog=True)
    gal_save_dict['rr'] = {'tree': tree, 'catalog': catalog}

    # Create a balltree for g(theta, r), DR calculation and pickle
    # Runtime is O(N*logM); M = size of tree; N = size of catalog
    # Create tree from catalog with smaller size
    print('- Setting up DR')
    mode = 'angular_tree'
    if rand.angular_distr.shape[0] <= data.ntotal:
        # Tree from Data, catalog from Random
        tree = data.build_balltree(metric='haversine', return_catalog=False)
        mode = 'data_tree'
    save_dict['dr'] = {'tree': tree,
                       'data_catalog': data.get_catalog(),
                       'angular_catalog': catalog,
                       'mode': mode}

    # Save helper object
    helper = CorrelationHelper(bins, models)
    helper.set_z_distr(gal_rand.z_distr)
    helper.set_norm(norm_dd=gal_data.norm(),
                    norm_rr=gal_rand.norm(),
                    norm_dr=gal_rand.norm(data))
    save_dict['helper'] = helper

    # Save
    save('{}_preprocess.pkl'.format(args.prefix), save_dict)


if __name__ == "__main__":
    main()
