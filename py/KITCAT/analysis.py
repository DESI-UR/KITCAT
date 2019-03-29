import numpy as np
from KITCAT.helper import JobHelper

def get_dd(
    catalog, tree,
    s_max       = 200.,
    s_nbins     = 50.,
    job_helper  = None,
    same        = False,
    checkpoint  = 10000
    ):
    """ calculate weighted and unweighted DD(s) of catalog and tree

    Parameters:
    -----------
    catalog:
    tree: kd-tree
    s_max: float
    s_nbins: int
    job_helper:
    same: bool
        set True if the tree is built from catalog .
        if True, will not apply double counting correction.

    Returns:
    --------
    dd: array of shape (2, s_nbins) """

    dd = np.zeros((2, s_nbins))

    # if job_helper is None, assume one job
    if job_helper is None:
        job_helper = JobHelper(1)
        job_helper.set_current_job(0, verbose=False)
    start, end = job_helper.get_index_range(catalog.shape[0])
    n = end - start - 1

    print('')
    print('calculate DD(s) from index %d to %d' % (start, end - 1))

    for i, pt in enumerate(catalog[start:end]):

        # print out checkpoint
        if i % checkpoint is 0:
            print('- index: %d/%d' % (i, n))

        # start query
        index, s = tree.query_radius(pt[:3].reshape(1, -1),
                                     r=s_max,
                                     return_distance=True)
        index = index[0]
        s = s[0]

        # fill weighted distribution
        # w =  w1 * w2
        w = catalog[:, 3][index]*pt[3]
        hist, _ = np.histogram(s, bins=s_nbins, range=(0., s_max), weights=w)
        dd[0] += hist

        # fill unweighted distribution
        hist, _ = np.histogram(s, bins=num_bins_s, range=(0., s_max))
        dd[1] += hist

    # double counting correction
    if same:
        dd = dd / 2.

    return dd

def get_ftheta(
    pair_catalog, tree_catalog, tree,
    theta_max   = 0.18,
    theta_nbins = 100,
    job_helper  = None,
    same        = False,
    checkpoint  = 10000
    ):
    """ calculate f(theta) of catalog and tree.

    Parameters:
    -----------
    catalog:
    tree: kd-tree
    theta_max: float
    theta_nbins: int
    job_helper:
    same: bool
        set True if the tree is built from catalog .
        if True, will not apply double counting correction.

    Returns:
    --------
    ftheta: array of shape (theta_nbins,) """

    ftheta = np.zeros(theta_nbins)

    # if job_helper is None, assume one job
    if job_helper is None:
        job_helper = JobHelper(1)
        job_helper.set_current_job(0, verbose=False)
    start, end = job_helper.get_index_range(pair_catalog.shape[0])
    n = end - start - 1

    print('')
    print('calculate f(theta) from index %d to %d' % (start, end - 1))

    for i, pt in enumerate(pair_catalog[start:end]):
        if i % checkpoint is 0:
            print('- index: %d/%d' % (i, n))
        index, theta = tree.query_radius(pt[:2].reshape(1, -1),
                                         r=theta_max,
                                         return_distance=True)
        index = index[0]
        theta = theta[0]

        # w = w1 * w2
        w = pt[2] * tree_catalog[:, 2][index]
        hist, _ = np.histogram(theta,
                               bins     = theta_nbins,
                               range    = (0., theta_max),
                               weights  = w)
        ftheta += hist

    if same:
        # Correction for double counting
        ftheta = ftheta / 2.

    return ftheta

def get_ztheta(
    pair_catalog, tree_catalog, tree,
    z_min       = 0.4,
    z_max       = 0.7,
    z_nbins     = 600,
    theta_max   = 0.18,
    theta_nbins = 100,
    job_helper  = None,
    checkpoint  = 10000
    ):
    """ calculate f(theta) of catalog and tree.

    Parameters:
    -----------
    pair_catalog:
    tree_catalog:
    tree: kd-tree
    z_min: float
    z_max: float
    z_nbins: int
    theta_max: float
    theta_nbins: int
    job_helper:

    Returns:
    --------
    ztheta: array of shape (theta_nbins,) """

    ztheta = np.zeros((2, theta_nbins, z_nbins))

    # if job_helper is None, assume one job
    if job_helper is None:
        job_helper = JobHelper(1)
        job_helper.set_current_job(0, verbose=False)
    start, end = job_helper.get_index_range(pair_catalog.shape[0])
    n = end - start - 1

    # set up binning parameters
    nbins = (theta_nbins, z_nbins)
    bins_range = ((0., theta_max), (z_min, z_max))

    print("calculate ztheta from index %d to %d" % (start, end - 1))

    for i, pt in enumerate(pair_catalog[start:end]):
        if i % checkpoint is 0:
            print('- index: %d/%d' % (i, n))

        index, theta = tree.query_radius(pt[:2].reshape(1, -1),
                                         r=theta_max,
                                         return_distance=True)
        index = index[0]
        theta = theta[0]
        iz = int(z_nbins * (pt[2]-z_min)/(z_max - z_min))

        # fill unweighted histogram
        w = tree_catalog[:, 2][index]
        hist, _ = np.histogram(theta,
                               bins     = theta_nbins,
                               range    = (0., theta_max),
                               weights  = w)
        ztheta[1][:, iz] += hist

        # fill weighted histogram
        w *= pt[3]
        hist, _ = np.histogram(theta,
                               bins     = theta_nbins,
                               range    = (0., theta_max),
                               weights  = w)
        ztheta[0][:, iz] += hist

    return ztheta


def get_zztheta(
    pair_catalog, tree_catalog, tree,
    z_min       = 0.4,
    z_max       = 0.7,
    z_nbins     = 600,
    theta_max   = 0.18,
    theta_nbins = 100,
    job_helper  = None,
    same        = False,
    checkpoint  = 10000
    ):
    """ calculate f(theta) of catalog and tree.

    Parameters:
    -----------
    galaxy_catalog:
    random_catalog:
    tree: kd-tree
    z_min: float
    z_max: float
    z_nbins: int
    theta_max: float
    theta_nbins: int
    job_helper:
    same: bool
        set True if the tree is built from catalog .
        if True, will not apply double counting correction.

    Returns:
    --------
    zztheta: array of shape (2, theta_nbins, z_nbins, z_nbins) """

    zztheta = np.zeros((2, theta_nbins, z_nbins, z_nbins))

    # if job_helper is None, assume one job
    if job_helper is None:
        job_helper = JobHelper(1)
        job_helper.set_current_job(0, verbose=False)
    start, end = job_helper.get_index_range(pair_catalog.shape[0])
    n = end - start - 1

    nbins = (theta_nbins, z_nbins)
    bins_range = ((0., theta_max), (z_min, z_max))

    print('')
    print('calculate zztheta from index %d to %d' % (start, end - 1))

    for i, pt in enumerate(pair_catalog[start:end]):
        if i % checkpoint is 0:
            print('- index: %d/%d' % (i, n))
        index, theta = tree.query_radius(pt[:2].reshape(1, -1),
                                         r=theta_max,
                                         return_distance=True)
        index = index[0]
        theta = theta[0]
        iz = int(z_nbins * (pt[2]-z_min)/(z_max - z_min))
        z = tree_catalog[:, 2][index]

        # fill weighted histogram
        # w = w1 * w2
        w = pt[3] * tree_catalog[:, 3][index]
        hist, _, _  = np.histogram2d(theta, z,
                                     bins     = nbins,
                                     range    = bins_range,
                                     weights  = w)
        zztheta[0][:, :, iz] += hist

        # fill unweighted histogram
        hist, _, _  = np.histogram2d(theta, z,
                                     bins     = nbins,
                                     range    = bins_range)
        zztheta[1][:, :, iz] += hist

    # double counting correction
    if same:
        zztheta = zztheta / 2.
    return zztheta
