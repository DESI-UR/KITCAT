""" Module with class for binnings """

import numpy as np

class Bins(object):
    """ Class to handle uniform binnings """
    def __init__(
        self,
        limit_params = None,
        nbins_params = None,
        min_cosmo    = None,
        max_cosmo    = None,
        islice       = 0,
        nslice       = 1,):

        # set up bin limit
        min_max = limit_params.copy()
        min_max['s_min'] = 0

        # convert to radian
        if limit_params['unit'] == 'deg':
            min_max['ra_max'] = np.deg2rad(min_max['ra_max'])
            min_max['ra_min'] = np.deg2rad(min_max['ra_min'])
            min_max['dec_max'] = np.deg2rad(min_max['dec_max'])
            min_max['dec_min'] = np.deg2rad(min_max['dec_min'])

        # for angular separation
        r_min = max_cosmo.z2r(min_max['z_min'])
        theta_max = np.arccos(1. - min_max['s_max']**2/(2 * r_min**2))

        # for redshift slice
        diff = (min_max['z_max'] - min_max['z_min']) / nslice
        z_min = min_max['z_min'] + diff * islice
        z_max = z_min + diff + max_cosmo.dels_to_delz(min_max['s_max'], z_min)
        z_max = min(z_max, min_max['z_max'])

        # create class dictionary
        self.limit = {}
        self.limit['s'] = (0., min_max['s_max'])
        self.limit['dec'] = (min_max['dec_min'], min_max['dec_max'])
        self.limit['ra'] = (min_max['ra_min'], min_max['ra_max'])
        self.limit['z'] = (z_min, z_max)
        self.limit['theta'] = (0., theta_max)

        # set up number of bins
        self.nbins = {}

        auto_binning = nbins_params['auto']
        if not auto_binning:
            self.nbins['s'] = nbins_params['s']
            self.nbins['dec'] = nbins_params['dec']
            self.nbins['ra'] = nbins_params['ra']
            self.nbins['z'] = nbins_params['z']
            self.nbins['theta'] = nbins_params['theta']
        else:
            self.nbins['s'] = nbins_params['s']
            self._set_auto_nbins(min_cosmo)

        # Print out number of bins
        self.print_info()

    def __eq__(self, other):
        """ Comparing one bins with other """
        for key, val in self.limit.items():
            if val != other.limit[key]:
                return False
        for key, val in self.nbins.items():
            if val != other.nbins[key]:
                return False
        return True

    def _set_auto_nbins(self, cosmo):
        """ Set number of bins based on binwidths """

        # Separation
        binw_s = (self.max('s') - self.min('s')) /self.num_bins('s')

        # Redshift/Comoving distance distribution
        binw_r = binw_s/2.
        self.nbins['z'], binw_r = self._find_nbins('r', binw_r, cosmo)

        # Angular variable
        binw_angl = binw_r/cosmo.z2r(self.max('z'))
        for key in ['dec', 'ra', 'theta']:
            self.nbins[key], _ = self._find_nbins(key, binw_angl)

    def _find_nbins(self, key, binw, model=None):
        """ Return the number of bins given key and binwidth """
        if key == 'r' and model is not None:
            low, high = model.z2r([self.min('z'), self.max('z')])
        else:
            low = self.min(key)
            high = self.max(key)

        nbins = int(np.ceil((high-low)/binw))
        binw = (high-low)/nbins

        return nbins, binw

    def find_zslice(self, index, total, model):
        """ Find zslice """
        diff = (self.max('z') - self.min('z')) / total
        z_low = index * diff + self.min('z')
        z_high = z_low + diff + model.dels_to_delz(self.max('s'))
        return z_low, z_high

    def min(self, key):
        """ Return binning lower bound """
        return self.limit[key][0]

    def max(self, key):
        """ Return binning upper bound """
        return self.limit[key][1]

    def num_bins(self, key):
        """ Return number of bins """
        return self.nbins[key]

    def binw(self, key):
        """ Return binwidth """
        return (self.max(key) - self.min(key)) / self.num_bins(key)

    def bins(self, key, cosmo=None):
        """ Return uniform binning """
        bins_min = self.min(key)
        bins_max = self.max(key)
        nbins = self.num_bins(key)

        # Uniformly over 'r'
        if key == 'z' and cosmo is not None:
            bins_min = cosmo.z2r(bins_min)
            bins_max = cosmo.z2r(bins_max)

        return np.linspace(bins_min, bins_max, nbins + 1)

    def print_info(self):
        """ Print out binning information """
        print('- bin range: ')
        for key in sorted(self.limit.keys()):
            print(' + %6s: %.5f, %.5f' % (key, self.min(key), self.max(key)))

        print('- nbins:')
        for key in sorted(self.nbins.keys()):
            print(' + %6s: %4d' % (key, self.num_bins(key)))
