#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Script for combining job results and calculate DD(s), DR(s), and RR(s) """

# Standard Python module
import argparse
import glob

# User-defined module
import numpy as np
import lib.io as lio
import lib.correlation as lcorrelation

if __name__ == "__main__":
    """ integration """

    print('')

    # read in helper
    print('reading file')
    fname_list = sorted(glob.glob("output/test/test_divide_*.pkl"))
    helper = None
    for i, fname in enumerate(fname_list):
        print("- %s" % fname)
        if i != 0:
            helper.add(lio.load(fname))
            continue
        helper = lio.load(fname)
        helper.ztheta_d1r2 = np.copy(helper.ztheta_d1r2)


    cosmo = helper.cosmos_list[0]
    bins = helper.bins

    theta = bins.bins('theta')
    theta = 0.5*(theta[:-1]+theta[1:])
    sin_2 = np.sin(theta/2.)
    cos_2 = np.cos(theta/2.)
    r = cosmo.z2r(bins.bins('z'))
    r = 0.5 * (r[:-1] + r[1:])
    s = np.linspace(0, 200, 51)

    ftheta = helper.ftheta
    ztheta = helper.ztheta_d1r2[0]
    zztheta = helper.zztheta[0]
    z_distr = helper.z1_distr[0]

    # calculate rr
    # create the weight matrix
    rr = np.zeros((50, 50))
    w = ftheta[:, None]*z_distr[None, :]
    for k, pt_r in enumerate(r):
        sigma = sin_2[:, None]*(pt_r + r[None, :])
        pi = cos_2[:, None]*np.abs(pt_r - r[None, :])
        hist, _, _ = np.histogram2d(sigma.ravel(), pi.ravel(),
                                    bins=(s,s),
                                    weights=w.ravel()*z_distr[k])
        rr += hist

    # calculate dr
    dr = np.zeros((50, 50))
    for k, pt_r in enumerate(r):
        sigma = sin_2[:, None]*(pt_r + r[None, :])
        pi = cos_2[:, None]*np.abs(pt_r - r[None, :])
        hist, _, _ = np.histogram2d(sigma.ravel(), pi.ravel(),
                                    bins=(s,s),
                                    weights=ztheta.ravel()*z_distr[k])
        dr += hist

    dd = np.zeros((50, 50))
    for k, pt_r in enumerate(r):
        sigma = sin_2[:, None]*(pt_r + r[None, :])
        pi = cos_2[:, None]*np.abs(pt_r - r[None, :])
        hist, _, _ = np.histogram2d(sigma.ravel(), pi.ravel(),
                                    bins=(s,s),
                                    weights=zztheta[:, :, k].ravel())
        dd += hist

    lio.save('test2d.pkl', {'rr': rr, 'dr': dr, 'dd': dd})
