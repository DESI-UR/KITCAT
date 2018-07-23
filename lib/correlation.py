""" Module to handle correlation function """

# Python modules
import numpy as np

class Correlation(object):
    """ Class to handle correlation function """
    def __init__(self, data_data, data_rand, rand_rand, bins, norm):
        """ Constructor sets up normalized DD, DR, RR from unnormalized"""
        self.bins = bins

        # Initialize separation distribution and its statistical errors
        self.w_distr = {}
        self.u_distr = {}

        # Set statistical errors
        keys = ('dd', 'dr', 'rr')
        for i, distr in enumerate((data_data, data_rand, rand_rand)):
            distr_err = self.distr_error(distr[0], distr[1])
            self.w_distr[keys[i]] = [distr[0], distr_err[0]]
            self.u_distr[keys[i]] = [distr[1], distr_err[1]]

        # Normalize and set separation distribution
        for i, distr in enumerate((self.w_distr, self.u_distr)):
            for key in keys:
                distr[key][0] = distr[key][0]/norm[key][i]
                distr[key][1] = distr[key][1]/norm[key][i] # propagate error

    def get_distr(self, key, weighted=True):
        """ Get distribution """
        if weighted:
            return self.w_distr[key]
        return self.u_distr[key]

    def distr_error(self, w_distr, uw_distr):
        """ Get statistical bin error of separation distribution
        Inputs:
        + w_dist: ndarray
            Weighted distribution
        + uw_dist: ndarray
            Unweighted distribution
        Outputs:
        + error: ndarray
            Bin errors """
        uw_err = np.sqrt(uw_distr)
        w_err = np.where(uw_err != 0, w_distr/uw_err, 0.)
        return w_err, uw_err

    def tpcf(self, weighted=True):
        """ Calculate two-point correlation and its error
        Outputs:
        + tpcf: ndarray
            Two-point correlation function
            Equation: (DD-2DR+RR)/RR
        + tpcf_err ndarray
            Statistical error of two-point correlation function """
        distr = self.w_distr if weighted else self.u_distr

        # Calculate tpcf
        # tpcf = (dd-2dr+rr)/rr
        tpcf = distr['dd'][0]-2*distr['dr'][0]+distr['rr'][0]
        tpcf = np.where(distr['rr'][0] != 0, tpcf/distr['rr'][0], 0.)

        # Calculate error of tpcf
        # err^2 = err_dd^2/rr^2
        tpcf_err = distr['dd'][1]/distr['rr'][0]
        # tpcf_err += np.where(
            # distr['rr'] != 0,
            # (distr['rr'][1]*(2*distr['dr'][0]-distr['dd'][0]))**2/distr['rr'][0]**2, 0)
        # tpcf_err = np.where(distr['rr'][0], tpcf_err/distr['rr'][0]**2, 0.)
        # tpcf_err = np.sqrt(tpcf_err)
        return tpcf, tpcf_err

    def tpcfss(self, weighted=True):
        """ Calculate two-point correlation s^2 and its error
        Outputs:
        + tpcfss: ndarray
            Two-point correlation function s^2
            Equation: s^2*(DD-2DR+RR)/RR """
        s = (self.bins[1:]+self.bins[:-1])/2.
        tpcf, tpcf_err = self.tpcf(weighted)
        return tpcf*s**2, tpcf_err*s**2
