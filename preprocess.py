""" Module to create preprocessed catalog object"""

# Standard Python module
import argparse

import numpy as np

# User-defined module
import lib.io as lio
import lib.catalog as lcatalog
import lib.helper as lhelper
import lib.bins as lbins
import lib.cosmology as lcosmology


if __name__ == '__main__':
    """ preprocess data: convert catalogs, binning, cosmology into a data structure
    """
    def parse_command_line():
        parser = argparse.ArgumentParser(description='preprocess data')
        parser.add_argument('--config', '-c',
                            help    = 'config file to read',
                            dest    = 'config',
                            type    = str)
        parser.add_argument('--prefix', '-p',
                            help    = 'output prefix',
                            dest    = 'prefix',
                            type    = str)
        parser.add_argument('--islice', '-i',
                            help    = 'z-slice index',
                            default = 0,
                            dest    = 'islice',
                            type    = int)
        parser.add_argument('--nslice', '-n',
                            help    = 'total number of z-slice',
                            default = 1,
                            dest    = 'nslice',
                            type    = int)
        parser.add_argument('--auto_binning', '-a',
                            help    = 'enable automatic bins',
                            action  = 'store_true',
                            dest    = 'auto_binning',
                            default = False)
        parser.add_argument('--binwidth_distance', '-b',
                            help    = 'binw of s if auto is enabled',
                            default = 2.00,
                            dest    = 'binw_s',
                            type    = float)
        params = parser.parse_args()
        return params

    args = parse_command_line()

    # check if arguments are valid
    if args.islice < 0 or args.islice >= args.nslice:
        raise ValueError('islice must be at least 0 and less than nslice.')

    # set up cosmology
    cosmos_params = lio.parse_config(args.config, 'COSMOLOGY')
    cosmos_list = []
    print('')
    print('setting up cosmology')
    print('- number of cosmology models: %d' %cosmos_params['n_cosmos'])
    for i in range(cosmos_params['n_cosmos']):
        hubble0 = cosmos_params['hubble0'][i]
        omega_m0 = cosmos_params['omega_m0'][i]
        omega_de0 = cosmos_params['omega_de0'][i]
        cosmos_list.append(lcosmology.Cosmology(hubble0   = hubble0,
                                                omega_m0  = omega_m0,
                                                omega_de0 = omega_de0))

    # set up binning
    print('')
    print('setting up binning')
    messg = 'auto' if args.auto_binning else 'manual'
    print('- binning mode: %s' %messg)

    limit_params = lio.parse_config(args.config, 'LIMIT')
    nbins_params = lio.parse_config(args.config, 'NBINS')
    bins = lbins.Bins(
        limit_params = limit_params,
        nbins_params = nbins_params,
        min_cosmo    = lcosmology.min_cosmo(cosmos_list),
        max_cosmo    = lcosmology.max_cosmo(cosmos_list),
        islice       = args.islice,
        nslice       = args.nslice,
        auto_binning = args.auto_binning,
        binw_s       = args.binw_s)

    # initialize catalog
    print('')
    print('initialize catalog')
    same = True
    d1_params = lio.parse_config(args.config, 'GALAXY_1')
    d2_params = lio.parse_config(args.config, 'GALAXY_2')
    r1_params = lio.parse_config(args.config, 'RANDOM_1')
    r2_params = lio.parse_config(args.config, 'RANDOM_2')

    d1 = lcatalog.GalaxyCatalog(d1_params, bins.limit)
    d2 = lcatalog.GalaxyCatalog(d2_params, bins.limit)
    r1 = lcatalog.GalaxyCatalog(r1_params, bins.limit)
    r2 = lcatalog.GalaxyCatalog(r2_params, bins.limit)

    # calculate normalization constant
    norm_dd = lcatalog.get_norm(d1, d2, same=same)
    norm_rr = lcatalog.get_norm(r1, r2, same=same)
    norm_d1r2 = lcatalog.get_norm(d1, r2, same=False)
    norm_d2r1 = norm_d1r2 if same else lcatalog.get_norm(d2, r1, same=False)

    print('- normalize factor:')
    print(' +   norm_dd: %.4e, %.4e' % norm_dd)
    print(' +   norm_rr: %.4e, %.4e' % norm_rr)
    print(' + norm_d1r2: %.4e, %.4e' % norm_d1r2)
    print(' + norm_d2r1: %.4e, %.4e' % norm_d2r1)

    # convert to random catalog
    r1 = r1.to_rand(
        z_min       = bins.min('z'),
        z_max       = bins.max('z'),
        z_nbins     = bins.num_bins('z'),
        ra_min      = bins.min('ra'),
        ra_max      = bins.max('ra'),
        ra_nbins    = bins.num_bins('ra'),
        dec_min     = bins.min('dec'),
        dec_max     = bins.max('dec'),
        dec_nbins   = bins.num_bins('dec'))
    r2 = r2.to_rand(
        z_min       = bins.min('z'),
        z_max       = bins.max('z'),
        z_nbins     = bins.num_bins('z'),
        ra_min      = bins.min('ra'),
        ra_max      = bins.max('ra'),
        ra_nbins    = bins.num_bins('ra'),
        dec_min     = bins.min('dec'),
        dec_max     = bins.max('dec'),
        dec_nbins   = bins.num_bins('dec'))

    # set up catalog and tree for RR(s)
    print('')
    print('setting up for RR(s)')
    rr_tree = r1.build_tree()
    if r1.ngals < r2.ngals:
        rr_pair_catalog = r1.get_catalog()
        rr_tree_catalog = r2.get_catalog()
    else:
        rr_pair_catalog = r2.get_catalog()
        rr_tree_catalog = r1.get_catalog()


    # set up catalog and tree for DD(s)
    print('')
    print('setting up for DD(s)')
    dd_tree = d1.build_tree(metric='haversine')
    if d1.ngals < d2.ngals:
        dd_pair_catalog = d1.get_catalog()
        dd_tree_catalog = d2.get_catalog()
    else:
        dd_pair_catalog = d2.get_catalog()
        dd_pair_catalog = d1.get_catalog()

    # set up catalog and tree for DR(s)
    print('')
    print('setting up for DR(s)')
    d1r2_tree = r2.build_tree()
    d1r2_pair_catalog = d1.get_catalog()
    d1r2_tree_catalog = r2.get_catalog()
    d2r1_tree = r1.build_tree()
    d2r1_pair_catalog = d2.get_catalog()
    d2r1_tree_catalog = r1.get_catalog()

    # set up helper object
    print('')
    print('setting up helper object')
    helper = lhelper.CorrelationHelper()
    helper.z1_distr = r1.z_distr
    helper.z2_distr = r2.z_distr
    helper.norm_dd = np.array(norm_dd)
    helper.norm_rr = np.array(norm_rr)
    helper.norm_d1r2 = np.array(norm_d1r2)
    helper.norm_d2r1 = np.array(norm_d2r1)

    # save into pickle file
    save_params = {
        'rr': {'tree': rr_tree,
               'tree_catalog': rr_tree_catalog,
               'pair_catalog': rr_pair_catalog},
        'dd': {'tree':dd_tree,
               'tree_catalog': dd_tree_catalog,
               'pair_catalog': dd_pair_catalog},
        'd1r2': {'tree': d1r2_tree,
               'tree_catalog': d1r2_tree_catalog,
               'pair_catalog': d1r2_pair_catalog,},
        'd2r1': {'tree': d2r1_tree,
               'tree_catalog': d2r1_tree_catalog,
               'pair_catalog': d2r1_pair_catalog,},
        'cosmos_list': cosmos_list,
        'bins': bins,
        'helper': helper,
        'same': same
    }
    lio.save('%s_preprocess.pkl' % args.prefix, save_params)

    print('')
