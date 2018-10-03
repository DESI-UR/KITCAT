""" Module to handle correlation function """

# Python modules
import numpy as np

class Correlation(object):
    """ class to handle correlation function """

    def __init__(
        self,
        rr          = None,
        dd          = None,
        d1r2        = None,
        d2r1        = None,
        bins        = None,
        norm_rr     = 1.,
        norm_dd     = 1.,
        norm_d1r2   = 1.,
        norm_d2r1   = 1.):
        """ constructor """

        # store parameters
        self.rr = rr
        self.dd = dd
        self.d1r2 = d1r2
        self.d2r1 = d2r1
        self.bins = bins
        self.norm_rr = norm_rr.reshape(-1, 1)
        self.norm_dd = norm_dd.reshape(-1, 1)
        self.norm_d1r2 = norm_d1r2.reshape(-1, 1)
        self.norm_d2r1 = norm_d2r1.reshape(-1, 1)

        # normalize
        self.rr /= self.norm_rr
        self.dd /= self.norm_dd
        self.d1r2 /= self.norm_d1r2
        self.d2r1 /= self.norm_d2r1

        # calculate error
        self.rr_err = self.get_error(rr) # / self.norm_rr
        self.dd_err = self.get_error(dd) # / self.norm_dd
        self.d1r2_err = self.get_error(d1r2) # / self.norm_dr
        self.d2r1_err = self.get_error(d2r1) # / self.norm_dr

        # calculate tpcf
        s = 0.5*(self.bins[1:] + self.bins[:-1])
        self.tpcf = self.dd - self.d1r2 - self.d2r1 + self.rr
        self.tpcf = np.where(self.rr != 0, self.tpcf / self.rr, 0)
        self.tpcf_err = self.dd_err/self.rr

        # calculate tpcf s^2
        self.tpcfss = self.tpcf * s**2
        self.tpcfss_err = self.tpcf_err * s**2

    def get_error(self, distr):
        """ calculate error for rr, dr, dd """
        err = np.zeros_like(distr)
        err[1] = np.sqrt(distr[1])
        err[0] = np.where(err[1] > 0, distr[0]/err[1], 0)
        return err
