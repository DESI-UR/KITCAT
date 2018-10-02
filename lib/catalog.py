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

    def __init__(self, catalog_params, limit_params):
        """ initialize galaxy catalog

        Parameters:
        -----------
        catalog_params: dict
            path and parameters of .fits
        limit_params: dict """

        # import catalog from .fits file
        print('- import catalog from %s' %catalog_params['path'])

        table = Table.read(catalog_params['path'])
        dec = np.deg2rad(table[catalog_params['dec']])
        ra = np.deg2rad(table[catalog_params['ra']])
        z = table[catalog_params['z']]
        try:
            w = table[catalog_params['weight']]
        except KeyError:
            w_fkp = table[catalog_params['weight_fkp']]
            w_noz = table[catalog_params['weight_noz']]
            w_cp = table[catalog_params['weight_cp']]
            w_sdc = table[catalog_params['weight_sdc']]
            w = w_sdc*w_fkp*(w_noz+w_cp-1)
        self.catalog = np.array([dec, ra, z, w]).T

        # apply limit cut
        min_dec, max_dec = limit_params['dec']
        min_ra, max_ra = limit_params['ra']
        min_z, max_z = limit_params['z']

        # Set boundaries
        mask = ((min_dec <= dec) & (dec <= max_dec) &
                (min_ra <= ra) & (ra <= max_ra) &
                (min_z <= z) & (z <= max_z))

        self.catalog = self.catalog[mask]
        self.ngals = self.catalog.shape[0]

    def get_catalog(self, cosmo=None):
        """ return catalog. convert z to d if cosmology is given """
        catalog = np.copy(self.catalog)
        if cosmo is not None:
            catalog[:, 2] = cosmo.z2r(catalog[:, 2])
        return catalog

    def norm(self):
        """ return weighted and unweighted normalization factor for DD """

        w = np.copy(self.catalog[:, 3])
        sum_w = np.sum(w)
        sum_w2 = np.sum(w**2)
        w_norm = 0.5 * (sum_w**2 - sum_w2)
        uw_norm = 0.5 * (self.ngals**2 - self.ngals)

        return w_norm, uw_norm

    def to_cartesian(self, cosmo):
        """ return galaxy catalog in Cartesian coordinates"""

        dec, ra, z, w = np.copy(self.catalog).T
        r = cosmo.z2r(z)
        catalog = np.array([r * np.cos(dec) * np.cos(ra),
                            r * np.cos(dec) * np.sin(ra),
                            r * np.sin(dec),
                            w]).T
        return catalog

    def to_rand(
        self,
        z_min       = None,
        z_max       = None,
        z_nbins     = None,
        ra_min      = None,
        ra_max      = None,
        ra_nbins    = None,
        dec_min     = None,
        dec_max     = None,
        dec_nbins   = None,
        ):

        """ convert into RandomCatalog

        Parameters:
        -----------
        z_min, z_max, ra_min, ra_max, dec_min, dec_max: float
        z_nbins, ra_nbins, dec_nbins: int

        Returns:
        --------
        rand: RandomCatalog """

        # calculate z distr (weighted and unweighted)
        z_distr_w, bins_z = self.redshift_distr(
            z_min       = z_min,
            z_max       = z_max,
            z_nbins     = z_nbins,
            weighted    = True,
            normed      = True)
        z_distr_uw, _ = self.redshift_distr(
            z_min       = z_min,
            z_max       = z_max,
            z_nbins     = z_nbins,
            weighted    = False,
            normed      = True)

        # calculate angular distribution
        angular_distr, bins_dec, bins_ra = self.angular_distr(
            ra_min      = ra_min,
            ra_max      = ra_max,
            ra_nbins    = ra_nbins,
            dec_min     = dec_min,
            dec_max     = dec_max,
            dec_nbins   = dec_nbins,
            weighted    = False,
            normed      = False)

        # Set up DistrCatalog attributes
        w = self.catalog[:, 3]

        rand = RandomCatalog()
        rand.ngals = self.ngals
        rand.sum_w = np.sum(w)
        rand.sum_w2 = np.sum(w**2)
        rand.z_distr = np.array([z_distr_w, z_distr_uw])
        rand.angular_distr = hist2point(angular_distr, bins_dec, bins_ra)
        rand.bins_z = bins_z
        rand.bins_dec = bins_dec
        rand.bins_ra = bins_ra

        return rand

    def redshift_distr(
        self,
        z_min    = None,
        z_max    = None,
        z_nbins  = None,
        cosmo    = None,
        weighted = False,
        normed   = False,
        ):
        """ calculate z distribution. convert z to d if cosmology is given

        Parameters:
        -----------
        z_min, z_max: float
        nbins: int
        cosmo:
        weighted: bool
            if true, return weighted histogram
        normed: bool
            if true, normalized histogram by number of data

        Returns:
        --------
        z_distr: array of shape (nbins, )
        bins_z: array of shape (nbins+1, )
            binsedge of distribution """

        w = self.catalog[:, 3] if weighted else None
        z = np.copy(self.catalog[:, 2])
        norm = 1.*self.ngals if normed else 1.

        # convert z to d
        if cosmo is not None:
            z_max = cosmo.z2r(z_max)
            z_min = cosmo.z2r(z_min)
            z = model.z2r(z)

        z_distr, bins_z = np.histogram(z,
                                       bins     = z_nbins,
                                       range    = (z_min, z_max),
                                       weights  = w)
        z_distr = z_distr/norm

        return z_distr, bins_z

    def angular_distr(
        self,
        ra_min      = None,
        ra_max      = None,
        ra_nbins    = None,
        dec_min     = None,
        dec_max     = None,
        dec_nbins   = None,
        weighted    = False,
        normed      = False):
        """ calculate angular distribution.

        Parameters:
        -----------
        ra_min, ra_max, dec_min, dec_max: float
        ra_nbins, dec_nbins: int
        weighted: bool
            if true, return weighted histogram
        normed: bool
            if true, normalized histogram by number of data

        Returns:
        --------
        angular_distr: array of shape (dec_nbins, ra_nbins)
        bins_dec: array of shape (dec_nbins, )
            binsedge of distribution
        bins_dec: array of shape (ra_nbins, )
            binsedge of distribution """

        dec, ra, _, w = np.copy(self.catalog).T
        w = w if weighted else None
        norm = 1.*self.ngals if normed else 1.

        angular_distr, bins_dec, bins_ra = np.histogram2d(
            dec, ra,
            bins        = (dec_nbins, ra_nbins),
            range       = ([dec_min, dec_max], [ra_min, ra_max]),
            weights     = w)
        angular_distr = angular_distr/norm

        return angular_distr, bins_dec, bins_ra

    def build_tree(self, metric, cosmo=None, return_catalog=False, leaf=40):
        """ a balltree from catalog. If metric is 'euclidean', cosmology is required.

        Parameters:
        ----------
        metric: str
            Metric must be either 'haversine' or 'euclidean'.
            If metric is 'haversine', build a tree from dec and ra.
            If metric is 'euclidean', build a tree from x, y, z.
        cosmo: cosmology.Cosmology (default=None)
            Cosmology model to convert redshift to comoving.
        leaf: int (default=40)
            Number of points at which KD-tree switches to brute-force. A leaf
            node is guaranteed to satisfy leaf_size <= n_points <= 2*leaf_size,
            except in the case that n_samples < leaf_size.
            More details can be found at sklearn.neightbors.BallTree. """

        if metric == 'euclidean':
            # convert Celestial coordinate into Cartesian coordinate
            catalog = self.to_cartesian(cosmo)
            tree = KDTree(catalog[:, :3], leaf_size=leaf, metric=metric)

        elif metric == 'haversine':
            catalog = self.catalog[:, :2]
            tree = BallTree(catalog, leaf_size=leaf, metric=metric)

        else:
            raise ValueError('metric must be "haversin" or "euclidean".')
        print("- building tree: %s " % metric)

        # return KD-tree and the catalog
        if return_catalog:
            return tree, catalog
        return tree


class RandomCatalog(object):
    """ Class to handle random catalog. Random catalog has the angular and
    redshif (comoving) distribution, but not the coordinates of each galaxy. """

    def __init__(self):
        """ initialize angular, z distribution, and norm variables """

        self.ngals = 0
        self.sum_w = 0
        self.sum_w2 =  0

        self.z_distr = None
        self.angular_distr = None
        self.bins_z = None
        self.bins_ra = None
        self.bins_dec = None


    def norm(self, data_catalog=None):
        """ return unweighted and weighted normalization factor """

        if data_catalog is not None:
            w_norm = np.sum(data_catalog.catalog[:, 3]) * sum_W
            uw_norm = data_catalog.ntotal * ngals
            return w_norm, uw_norm

        w_norm = 0.5 * (self.norm_vars['sum_w']**2 - self.norm_vars['sum_w2'])
        uw_norm = 0.5 * (self.norm_vars['ngals']**2 - self.norm_vars['ngals'])
        return w_norm, uw_norm

    def get_catalog(self, cosmo=None):
        """ return angular distribution """
        return np.copy(self.angular_distr)

    def build_tree(self, leaf=40):
        """ build a balltree from angular distributions using haversine.

        Parameters:
        -----------
        leaf: int (default=40)
            Number of points at which KD-tree switches to brute-force. A leaf
            node is guaranteed to satisfy leaf_size <= n_points <= 2*leaf_size,
            except in the case that n_samples < leaf_size.
            More details can be found at sklearn.neightbors.BallTree.

        Returns:
        --------
        balltree """

        print("- building tree: haversine")
        balltree = BallTree(self.angular_distr[:, :2],
                            leaf_size=leaf,
                            metric='haversine')
        return balltree
