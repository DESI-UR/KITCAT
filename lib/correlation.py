""" Module to handle correlation function """

# Python modules
import numpy as np

def tpcf(rr, dd, d1r2, d2r1, norm_rr, norm_dd, norm_d1r2, norm_d2r1):
    """ Calculate 1d or 2d tpcf

    Arguments:
    ----------
    rr, dd, d1r2, d2r1: array of shape (2, n, 1) or (2, n, n)
        Unormalized distribution
    norm_rr, norm_dd, norm_d1r2, norm_d2r1: (2, )
        Normalization factor

    Returns:
    --------
    xi: array of shape (2, n, 1) or (2, n, n)
        two-point correlation function
    xi_err: array of shape (2, n, 1) or (2, n, n)
        statistical error of the two-point correlation
    """

    # normalize
    norm_rr = norm_rr.reshape(2, 1, 1)
    norm_dd = norm_dd.reshape(2, 1, 1)
    norm_d1r2 = norm_d1r2.reshape(2, 1, 1)
    norm_d2r1 = norm_d2r1.reshape(2, 1, 1)

    rr /= norm_rr
    dd /= norm_dd
    d1r2 /= norm_d1r2
    d2r1 /= norm_d2r1

    # calculate error
    rr_err = get_error(rr) / norm_rr
    dd_err = get_error(dd) / norm_dd
    d1r2_err = get_error(d1r2) / norm_d1r2
    d2r1_err = get_error(d2r1) / norm_d2r1

    # calculate tpcf
    xi = dd - d1r2 -d2r1 + rr
    xi = np.where(rr != 0, xi/rr, 0)
    xi_err = dd_err/rr

    return xi, xi_err

def get_error(dist):
    """ calculate error for rr, dr, dd """
    err = np.zeros_like(dist)
    err[1] = np.sqrt(dist[1])
    err[0] = np.where(err[1] > 0, dist[0]/err[1], 0)
    return err
