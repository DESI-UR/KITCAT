""" Module to handle galaxy survey catalogs """

# Python modules
import numpy as np
from astropy.table import Table
from sklearn.neighbors import BallTree, KDTree

def hist2point(hist, bins_x, bins_y, exclude_zeros=True):
    """ Convert 2D histogram into a set of weighted data points.
        Use bincenter as coordinate.
        Inputs:
        + hist: ndarray
            2D histogram
        + bins_x: ndarray
            Binedges along the  x-axis
        + bins_y: ndarray
            Binedges along the y-axis
        + exclude_zeros: bool (default=True)
            If True, return non-zero weight points only.
        Outputs:
        + catalog: ndarrays
            Catalog of weighted data points. Format: [X, Y, Weight]"""

    # Get bins center and create a grid
    center_x = 0.5*(bins_x[:-1]+bins_x[1:])
    center_y = 0.5*(bins_y[:-1]+bins_y[1:])
    grid_x, grid_y = np.meshgrid(center_x, center_y)
    grid_x = grid_x.flatten()
    grid_y = grid_y.flatten()

    # Convert the grid into a catalog
    hist = hist.T.flatten()
    catalog = np.array([grid_x, grid_y, hist]).T

    # Return catalog
    if exclude_zeros:
        zeros = (hist == 0)
        return catalog[np.logical_not(zeros)]
    return catalog

class GalaxyCatalog(object):
    """ Class to handle galaxy catalogs. """

    def __init__(self, reader, limit):
        """ Initialize galaxy catalog.
        Inputs:
        + reader: dict
            Dictionary with path and properties of catalog file.
            Keys: path, dec, ra, z. Values format: str
        + limit: dict
            Dictionary with boundaries (inclusive) of catalog.
            If None, use full catalog.
            Keys: dec, ra, z. Values format (min, max) """
        self.ntotal = 0
        self.catalog = None
        self.set_catalog(reader)
        if limit is not None:
            self.set_limit(limit)

    def set_catalog(self, reader):
        """ Read in galaxy catalog from FITS file. Convert redshift to
        comoving distance.
        Inputs:
        + reader: dict
            Dictionary with path and properties of catalog file.
            Keys: path, dec, ra, z. Values format: str """

        # Read in configuration file
        print(" - Import catalog from: {}".format(reader['path']))

        # Read in table
        table = Table.read(reader['path'])
        dec = np.deg2rad(table[reader['dec']].data)
        ra = np.deg2rad(table[reader['ra']].data)
        z = table[reader['z']].data
        try:
            w = table[reader['weight']]
        except KeyError:
            w_fkp = table[reader['weight_fkp']].data
            w_noz = table[reader['weight_noz']].data
            w_cp = table[reader['weight_cp']].data
            w_sdc = table[reader['weight_sdc']].data
            w = w_sdc*w_fkp*(w_noz+w_cp-1)

        self.catalog = np.array([dec, ra, z, w]).T
        self.ntotal = self.catalog.shape[0]

    def get_catalog(self, model=None):
        """ Get catalog. If cosmology is given, convert redshift to comoving
        distance.
        Inputs:
        + model: cosmology.Cosmology (default=None)
            Cosmology model to convert redshift to comoving distance.
        Outputs:
        + catalog: ndarray
            Data catalog. """
        # If cosmology is given, convert redshift to comoving distance
        catalog = np.copy(self.catalog)
        if model is not None:
            catalog[:, 2] = model.z2r(catalog[:, 2])
        return catalog

    def set_limit(self, limit):
        """ Set boundaries of catalog
        Inputs:
        + limit: dict
            Dictionary with boundaries (inclusive) of catalog
            Keys: dec, ra, z. Values format: (min, max) """

        # Extract min and max values from dictionary
        min_dec, max_dec = limit['dec']
        min_ra, max_ra = limit['ra']
        min_z, max_z = limit['z']

        # Set boundaries
        dec, ra, z = self.catalog[:, :3].T
        mask = ((min_dec <= dec) & (dec <= max_dec) &
                (min_ra <= ra) & (ra <= max_ra) &
                (min_z <= z) & (z <= max_z))

        self.catalog = self.catalog[mask]
        self.ntotal = self.catalog.shape[0]

    def to_rand(self, limit, nbins):
        """ Convert into RandomCatalog and return.
        Inputs:
        + limit: list, tuple, ndarray, dict
            Bins range in order 'dec', 'ra', 'z'.
            If dict, use values from keys 'dec', 'ra', 'z'
        + nbins: list, tuple, ndarray, dict
            Number of bins in order 'dec', 'ra', 'z'.
            If dict, use values from keys 'dec', 'ra', 'z'
        Outputs:
        + rand: RandomCatalog
            Distribution catalog """

        # Initialize bins range and number of bins
        num_bins = nbins
        bins_range = limit
        if isinstance(bins_range, dict):
            bins_range = (limit['dec'], limit['ra'], limit['z'])
        if isinstance(num_bins, dict):
            num_bins = (nbins['dec'], nbins['ra'], nbins['z'])

        # Calculate comoving distribution
        z_distr_w, bins_z = self.redshift_distr(
            bins_range[2], num_bins[2], weighted=True, normed=True)
        z_distr_uw, _ = self.redshift_distr(
            bins_range[2], num_bins[2], weighted=False, normed=True)

        # Calculate angular distribution as a set of weighted data points
        angular_distr, bins_dec, bins_ra = self.angular_distr(
            bins_range[:2], num_bins[:2], weighted=False, normed=False)

        # Set up DistrCatalog attributes
        rand = RandomCatalog()
        rand.z_distr = np.array([z_distr_w, z_distr_uw])
        rand.angular_distr = hist2point(angular_distr, bins_dec, bins_ra)
        rand.bins_z = bins_z
        rand.bins_dec = bins_dec
        rand.bins_ra = bins_ra
        rand.norm_vars['ntotal'] = self.ntotal
        rand.norm_vars['sum_w'] = np.sum(self.catalog[:, 3])
        rand.norm_vars['sum_w2'] = np.sum(self.catalog[:, 3]**2)

        return rand

    def redshift_distr(self, limit, nbins, model=None,
                       weighted=False, normed=False):
        """ Calculate redshift distribution.
        If cosmology is given, convert redshift to comoving distance.
        Inputs:
        + limit: list, tuple, ndarray, dict
            Binning range for redshift histogram.
            If dict, use value of key 'z'.
        + nbins: int
            Number of bins for comoving histogram.
        + model: cosmology.Cosmology (default=None)
            Cosmology model to convert comoving to redshift.
        + weighted: bool (default=False)
            If True, return weighted histogram.
            Else return unweighted histogram.
        + normed: bool (default=False)
            If True, normalized by number of galaxies.
        Outputs:
        + z_distr: ndarray
            Redshift (comoving) distribution.
        + bins_z: ndarray
            Redshift (comoving) binedges. """
        bins_range = limit['z'] if isinstance(limit, dict) else limit
        weights = self.catalog[:, 3] if weighted else None
        z = self.catalog[:, 2]

        # If given cosmological model, calculate comoving distribution
        if model is not None:
            bins_range = model.z2r(bins_range)
            z = model.z2r(z)

        z_distr, bins_z = np.histogram(z, bins=nbins, range=bins_range,
                                       weights=weights)

        # Normalized by number of galaxies
        if normed:
            z_distr = z_distr/self.ntotal

        return z_distr, bins_z

    def angular_distr(self, limit, nbins, weighted=False, normed=False):
        """ Calculate angular distribution.
        Inputs:
        + limit: list, tupple, ndarray, dict
            Bins range in order 'dec', 'ra'.
            If dict, use value from keys 'dec' and 'ra'
        + nbins: list, tupple, ndarray, dict
            Number of bins in order 'dec', 'ra'.
            If dict, use value from keys 'dec' and 'ra'
        + weighted: bool (default=False)
            If True, return weighted histogram. Else return unweighted.
        + normed: bool (default=False)
            If True, normalized by number of galaxies.
        Outputs:
        + angular_distr: ndarray
            Angular distribution
        + bins_dec: ndarray
            Angular bins for declination
        + bins_ra: ndarray
            Angular bins for right ascension """

        # Initialize bins range and number of bins
        num_bins = nbins
        bins_range = limit
        if isinstance(bins_range, dict):
            bins_range = (limit['dec'], limit['ra'])
        if isinstance(num_bins, dict):
            num_bins = (nbins['dec'], nbins['ra'])

        # Calculate distribution
        weights = self.catalog[:, 3] if weighted else None
        angular_distr, bins_dec, bins_ra = np.histogram2d(
            self.catalog[:, 0], self.catalog[:, 1],
            bins=num_bins, range=bins_range, weights=weights)

        # Normalized by number of galaxies
        if normed:
            angular_distr = 1. * angular_distr / self.catalog.shape[0]

        return angular_distr, bins_dec, bins_ra

    def norm(self):
        """ Return unweighted and weighted normalization factor
        for pairwise separation distribution
        Unweighted equation:
        - norm = 0.5(ntotal^2-ntotal); ntotal is the size of catalog
        Weighted equation:
        - norm = 0.5(sum_w^2-sum_w2)
        where sum_w and sum_w2 are the sum of weights and weights squared
        Outputs:
        + w_norm: float
            Weighted normalization factor
        + uw_norm: float
            Unweighted normalization factor """

        sum_w = np.sum(self.catalog[:, 3])
        sum_w2 = np.sum(self.catalog[:, 3]**2)

        w_norm = 0.5 * (sum_w**2 - sum_w2)
        uw_norm = 0.5 * (self.ntotal**2 - self.ntotal)

        return w_norm, uw_norm

    def to_cartesian(self, model):
        """ Return galaxy catalog in Cartesian coordinates
        Inputs:
        + cosmo: cosmology.Cosmology (default=None)
            Cosmology model to convert redshift to comoving.
        Outputs:
        + catalog: np.ndarray
            Galaxy catalog in Cartesian coordinate system. """
        dec, ra, z = self.catalog[:, :3].T
        r = model.z2r(z)
        catalog = np.array([r * np.cos(dec) * np.cos(ra),
                            r * np.cos(dec) * np.sin(ra),
                            r * np.sin(dec),
                            self.catalog[:, 3]]).T
        return catalog

    def build_balltree(self, metric, model=None, return_catalog=False, leaf=40):
        """ Build a balltree from catalog. If metric is 'euclidean',
        cosmology is required.
        Inputs:
        + metric: str
            Metric must be either 'haversine' or 'euclidean'.
            If metric is 'haversine', build a balltree from DEC and RA
            coordinates of galaxies.
            If metric is 'euclidean', build a 3-dimensional kd-tree
        + return_catalog: bool (default=False)
            If True, return the catalog in balltree
        + cosmo: cosmology.Cosmology (default=None)
            Cosmology model to convert redshift to comoving.
        + leaf: int (default=40)
            Number of points at which KD-tree switches to brute-force. A leaf
            node is guaranteed to satisfy leaf_size <= n_points <= 2*leaf_size,
            except in the case that n_samples < leaf_size.
            More details can be found at sklearn.neightbors.BallTree. """
        if metric == 'euclidean':
            if model is None:
                raise TypeError('Cosmology must be given if metric ' +
                                'is "euclidean".')
            # Convert Celestial coordinate into Cartesian coordinate
            catalog = self.to_cartesian(model)
            tree = KDTree(catalog[:, :3], leaf_size=leaf, metric='euclidean')
        elif metric == 'haversine':
            catalog = self.catalog[:, :-2]
            tree = BallTree(catalog, leaf_size=leaf, metric=metric)
        else:
            raise ValueError('Metric must be either "haversine" ' +
                             'or "euclidean".')

        print(" - Creating BallTree with metric %s " % metric)

        # Return KD-tree and the catalog
        if return_catalog:
            return tree, catalog
        return tree


class RandomCatalog(object):
    """ Class to handle random catalog. Random catalog has the angular and
    redshif (comoving) distribution, but not the coordinates of each galaxy. """

    def __init__(self):
        """ Initialize angular, redshift (comoving) distribution,
        and normalization variables """
        self.z_distr = None
        self.angular_distr = None
        self.bins_z = None
        self.bins_ra = None
        self.bins_dec = None
        self.norm_vars = {'ntotal': None, 'sum_w': None, 'sum_w2': None}

    def norm(self, data_catalog=None):
        """ Return unweighted and weighted normalization factor
        for pairwise separation distribution
        Unweighted equation:
        - norm = 0.5(ntotal^2-ntotal); ntotal is the size of catalog
        Weighted equation:
        - norm = 0.5(sum_w^2-sum_w2);
        where sum_w and sum_w2 are the sum of weights and weights squared
        Inputs:
        + data_catalog: DataCatalog (default=None)
            If None, calculate the normalization factor for itself (i.e. RR).
            Otherwise, calculate the normalization factor for
            correlation distribution with the input catalog (i.e. DR).
        Outputs:
        + w_norm: float
            Weighted normalization factor
        + uw_norm: float
            Unweighted normalization factor """

        if data_catalog is not None:
            w_norm = (np.sum(data_catalog.catalog[:, 3]) *
                      self.norm_vars['sum_w'])
            uw_norm = data_catalog.ntotal * self.norm_vars['ntotal']
            return w_norm, uw_norm

        w_norm = 0.5 * (self.norm_vars['sum_w']**2 - self.norm_vars['sum_w2'])
        uw_norm = 0.5 * (self.norm_vars['ntotal']**2 - self.norm_vars['ntotal'])
        return w_norm, uw_norm

    def build_balltree(self, return_catalog=False, leaf=40):
        """ Build a balltree using DEC and RA from angular distributions.
        Metric: haversine.
        + return_catalog: bool (default=False)
            If True, return the angular distribution catalog.
        + leaf: int (default=40)
            Number of points at which KD-tree switches to brute-force. A leaf
            node is guaranteed to satisfy leaf_size <= n_points <= 2*leaf_size,
            except in the case that n_samples < leaf_size.
            More details can be found at sklearn.neightbors.BallTree."""
        print("  - Creating BallTree")
        balltree = BallTree(self.angular_distr[:, :2], leaf_size=leaf,
                            metric='haversine')
        if return_catalog:
            return balltree, self.angular_distr
        return balltree
