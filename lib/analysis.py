

import numpy as np
from lib.helper import JobHelper


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
    catalog: array if shape (s_nbins, 3)
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
    n = start - end + 1

    print('')
    print('calculate DD(s) from index %d to %d' % (start, end - 1))

    w_catalog = catalog[:, 3]
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
        w = w_catalog[index]*pt[3]
        hist, _ = np.histogram(s, bins=s_nbins, range=(0., s_max), weights=w)
        dd[0] += hist

        # fill unweighted distribution
        hist, _ = np.histogram(s, bins=num_bins_s, range=(0., s_max))
        dd[1] += hist

    # double counting correction
    if same:
        dd[0][0] -= np.sum(catalog[start:end, 3]**2)
        dd[1][0] -= end - start
        dd = dd / 2.

    return dd

def get_ftheta(
    catalog, tree,
    theta_max   = 0.18,
    theta_nbins = 100,
    job_helper  = None,
    same        = False,
    checkpoint  = 10000
    ):
    """ calculate f(theta) of catalog and tree.

    Parameters:
    -----------
    catalog: array if shape (n_gal, 3)
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
    start, end = job_helper.get_index_range(catalog.shape[0])
    n = start - end + 1

    print('')
    print('calculate f(theta) from index %d to %d' % (start, end - 1))

    catalog_w = catalog[:, 2]
    for i, pt in enumerate(catalog[start:end]):
        if i % 10000 is 0:
            print('- index: %d/%d' % (i, n))
        index, theta = tree.query_radius(pt[:2].reshape(1, -1),
                                         r=theta_max,
                                         return_distance=True)
        index = index[0]
        theta = theta[0]

        # w = w1 * w2
        w = pt[2] * catalog_w[index]
        hist, _ = np.histogram(theta, bins=theta_nbins, range=(0., theta_max), weights=w)
        ftheta += hist

    if same:
        # Correction for double counting
        ftheta = ftheta / 2.

    return ftheta