#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Script for combining job results and calculate DD(s), DR(s), and RR(s) """

# Standard Python module
import argparse
import glob

# User-defined module
import lib.io as lio
import lib.correlation as lcorrelation

if __name__ == "__main__":
    """ integration """

    print('')

    # read in helper
    print('reading file')
    fname_list = sorted(glob.glob("test_divide_*.pkl"))
    helper = None
    for i, fname in enumerate(fname_list):
        print("- %s" % fname)
        if i != 0:
            helper.add(lio.load(fname))
            continue
        helper = lio.load(fname)
        import numpy as np
        helper.ztheta_d1r2 = np.copy(helper.ztheta_d1r2)


    cosmo = helper.cosmos_list[0]
    bins = helper.bins

    theta = bins.bins('theta')
    theta = 0.5*(theta[:-1]+theta[1:])
    r = cosmo.z2r(bins.bins('z'))
    r = 0.5 * (r[:-1] + r[1:])
    s = np.linspace(0, 200, 51)

    ftheta = helper.ftheta
    ztheta = helper.ztheta_d1r2
    zztheta = helper.zztheta
    z_distr = helper.z1_distr

    rr = np.zeros((2, 50, 50))
    dd = np.zeros((2, 50, 50))
    dr = np.zeros((2, 50, 50))

    # calculate rr
    # loop over cosmology models
    print('')
    print('calculate rr')

    for i in range(2):

        # calculate 2-d weight matrix
        w = ftheta[:, None]*z_distr[i][None, :]
        w = w.ravel()

        sin_2 = np.sin(theta/2.)
        cos_2 = np.cos(theta/2.)

        # Calculate RR(sigma, pi)
        for j, pt_r in enumerate(r):
            # calculate distance
            sigma = (pt_r + r[None, :])*sin_2[:, None]
            sigma = sigma.ravel()

            pi = np.abs(pt_r - r[None, :])*cos_2[:, None]
            pi = pi.ravel()

            hist, _, _ = np.histogram2d(sigma, pi,
                                        bins    = (s, s),
                                        weights = z_distr[i][j]*w)
            rr[i] += hist

    # calculate dr
    print('')
    print('calculate dr')

    for i in range(2):

        w = ztheta[i].ravel()

        sin_2 = np.sin(theta/2.)
        cos_2 = np.cos(theta/2.)

        # Calculate DR(s)
        for j, pt_r in enumerate(r):
            # calculate distance
            sigma = (pt_r + r[None, :])*sin_2[:, None]
            sigma = sigma.ravel()

            pi = np.abs(pt_r - r[None, :])*cos_2[:, None]
            pi = pi.ravel()

            hist, _, _ = np.histogram2d(sigma, pi,
                                        bins    = (s, s),
                                        weights = w)
            dr[i] += hist

    # calculate dd
    print('')
    print('calculate dd')

    for i in range(2):
        for j, pt_theta in enumerate(theta):

            sin_2 = np.sin(pt_theta/2)
            cos_2 = np.cos(pt_theta/2)

            for k, pt_r in enumerate(r):
                w = zztheta[i, j, k]

                # calculate distance
                sigma = (pt_r + r)*sin_2
                pi = np.abs(pt_r - r)*cos_2

                hist, _, _ = np.histogram2d(sigma, pi,
                                            bins    = (s, s),
                                            weights = w)
                dd[i] += hist

    import seaborn as sns
    import matplotlib.pyplot as plt

    i = 0
    rr = rr[0]
    dd = dd[0]
    dr = dr[0]
    rr = rr/helper.norm_rr.reshape(-1, 1, 1)
    dd = dd/helper.norm_dd.reshape(-1, 1, 1)
    dr = dr/helper.norm_d1r2.reshape(-1, 1, 1)
    rr = rr[i]
    dd = dd[i]
    dr = dr[i]

    tpcf = dd - 2*dr + rr
    tpcf = np.where(rr != 0, tpcf/rr, 0)

    sns.heatmap(tpcf)
    plt.show()
