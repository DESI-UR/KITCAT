""" Module to handle correlation function """

# Python modules
import numpy as np

class Correlation(object):
    """ class to handle correlation function """

    def __init__(self, rr, dr, dd, bins, norm_rr, norm_dr, norm_dd):
        """ constructor """

        # store parameters
        self.rr = rr
        self.dr = dr
        self.dd = dd
        self.bins = bins
        self.norm_rr = norm_rr.reshape(-1, 1)
        self.norm_dr = norm_dr.reshape(-1, 1)
        self.norm_dd = norm_dd.reshape(-1, 1)

        # normalize
        self.rr /= self.norm_rr
        self.dr /= self.norm_dr
        self.dd /= self.norm_dd

        # calculate error
        self.rr_err = self.get_error(rr) / self.norm_rr
        self.dr_err = self.get_error(dr) / self.norm_dr
        self.dd_err = self.get_error(dd) / self.norm_dd

        # calculate tpcf
        s = 0.5*(self.bins[1:] + self.bins[:-1])
        self.tpcf = self.dd - 2*self.dr + self.rr
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
