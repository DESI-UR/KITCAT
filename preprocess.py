""" Module to create preprocessed catalog object"""

# Standard Python module
import argparse

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
    data_params = lio.parse_config(args.config, 'GALAXY')
    rand_params = lio.parse_config(args.config, 'RANDOM')

    data = lcatalog.GalaxyCatalog(data_params, bins.limit)
    rand = lcatalog.GalaxyCatalog(rand_params, bins.limit)
    rand = rand.to_rand(
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
    rr_tree = rand.build_tree()
    rr_catalog = rand.get_catalog()

    # set up catalog and tree for DD(s)
    print('')
    print('setting up for DD(s)')
    dd_tree = data.build_tree(metric='haversine')
    dd_catalog = data.get_catalog()

    # set up catalog and tree for DR(s)
    print('')
    print('setting up for DR(s)')
    dr_tree = None
    dr_catalog = None

    # save into pickle file
    save_params = {
        'rr': {'tree': rr_tree, 'catalog': rr_catalog},
        'dd': {'tree': dd_tree, 'catalog': dd_catalog},
        'dr': {'tree': dr_tree, 'catalog': dr_catalog},
        'cosmos_list': cosmos_list,
        'bins': bins,
    }
    lio.save('%s_preprocess.pkl' % args.prefix, save_params)