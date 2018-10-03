""" Module with helper class """

import numpy as np

def distance(theta, r1, r2):
    """ Calculate distance between two points at radius r1, r2 separated by
    angle theta """
    return np.sqrt(r1**2 + r2**2 - 2 * r1 * r2 * np.cos(theta))

class JobHelper(object):
    """ Class to handle multiprocess job """
    def __init__(self, total_jobs):
        """ Constructor """
        if total_jobs <= 0:
            raise ValueError('total jobs must be at least 1')
        self.total_jobs = total_jobs
        self.current_job = 0

    def increment(self, verbose=True):
        """ Increment current job index by 1 """
        if self.current_job != self.total_jobs - 1:
            self.current_job += 1
        else:
            print('already at the last job index')

        if verbose:
            print("job %d/%d" % (self.current_job+1, self.total_jobs))

    def set_current_job(self, current_job, verbose=True):
        """ Set current job index to input """
        if current_job < 0 or current_job >= self.total_jobs:
            raise ValueError('job must be at least 0 and less than total job')
        self.current_job = current_job

        if verbose:
            print("job %d/%d" % (self.current_job+1, self.total_jobs))

    def get_index_range(self, size):
        """ Calculate the start and end indices given job size
        Inputs:
        + size: int
            Size of job.
        Outputs:
        + job_range: tuple
            Return the start and end indices """
        job_index = np.floor(np.linspace(0, size, self.total_jobs + 1))
        job_index = job_index.astype(int)
        job_range = (job_index[self.current_job],
                     job_index[self.current_job + 1])
        return job_range

class CorrelationHelper(object):
    """ Class to handle multiprocess correlation function calculation """

    def __init__(self):
        """ constructor """
        self.norm_dd = 1.
        self.norm_dr = 1.
        self.norm_rr = 1.
        self.zztheta = None
        self.ftheta = None
        self.ztheta = None
        self.z_distr = None
        self.cosmos_list = None
        self.bins = None

    def add(self, other):
        self.ftheta  += other.ftheta
        self.ztheta  += other.ztheta
        self.zztheta += other.zztheta

    def get_rr(self, nloops=1):

        n_models = len(self.cosmos_list)
        rr = np.zeros((n_models, 2, self.bins.num_bins('s')))

        theta = self.bins.bins('theta')
        theta = 0.5*(theta[:-1] + theta[1:])

        s = self.bins.bins('s')

        # calculate rr
        # loop over cosmology models
        for i, cosmo in enumerate(self.cosmos_list):
            r = cosmo.z2r(self.bins.bins('z'))
            r = 0.5 * (r[:-1] + r[1:])

            if nloops == 0:
                # no explicit for-loop to calculate distance
                # fast at the expense of memory
                # one for-loop is needed for weighting

                for j in range(2):

                    # calculate 3-d weight matrix
                    w = (self.z_distr[j][None, :, None] *
                         self.z_distr[j][None, None, :] *
                         self.ftheta[:, None, None])

                    # calculate 3-d distance matrix
                    dist = distance(theta = theta[:, None, None],
                                    r1    = r[None, :, None],
                                    r2    = r[None, None, :])

                    # calculate RR(s)
                    rr[i][j], _ = np.histogram(dist,
                                               bins     = s,
                                               weights  = w)

            elif nloops == 1:

                for j in range(2):

                    # calculate 2-d weight matrix
                    w = self.ftheta[:, None] * self.z_distr[j][None, :]

                    # Build a 2-d distance matrices
                    temp = np.cos(theta[:, None]) * r[None, :]
                    r_sq = r[None, :]**2

                    # Calculate RR(s)
                    for k, pt_r in enumerate(r):
                        dist = np.sqrt(pt_r**2 + r_sq - 2*pt_r*temp)
                        hist, _ = np.histogram(dist,
                                               bins    = s,
                                               weights = self.z_distr[j][k]*w)
                        rr[i][j] += hist
            else:
                raise ValueError('Invalid input for num_loops')

        rr = np.squeeze(rr)
        return rr


    def get_dr(self, nloops=1):

        n_models = len(self.cosmos_list)
        dr = np.zeros((n_models, 2, self.bins.num_bins('s')))

        theta = self.bins.bins('theta')
        theta = 0.5*(theta[:-1] + theta[1:])

        s = self.bins.bins('s')

        # calculate dr
        # loop over cosmology models
        for i, cosmo in enumerate(self.cosmos_list):
            r = cosmo.z2r(self.bins.bins('z'))
            r = 0.5 * (r[:-1] + r[1:])

            if nloops == 0:
                # no explicit for-loop to calculate distance
                # fast at the expense of memory
                # one for-loop is needed for weighting

                for j in range(2):
                    # Calculate 3-dimensional weights matrix
                    w = self.z_distr[j][None, None, :]*self.ztheta[j][:, :, None]

                    # calculate 3-d distance matrix
                    dist = distance(theta = theta[:, None, None],
                                    r1    = r[None, :, None],
                                    r2    = r[None, None, :])

                    # calculate RR(s)
                    dr[i][j], _ = np.histogram(dist,
                                               bins     = s,
                                               weights  = w)

            elif nloops == 1:
                for j in range(2):

                    w = self.ztheta[j]

                    # Build a 2-d distance matrices
                    temp = np.cos(theta[:, None]) * r[None, :]
                    r_sq = r[None, :]**2

                    # Calculate RR(s)
                    for k, pt_r in enumerate(r):
                        dist = np.sqrt(pt_r**2 + r_sq - 2*pt_r*temp)
                        hist, _ = np.histogram(dist,
                                               bins    = s,
                                               weights = self.z_distr[j][k]*w)
                        dr[i][j] += hist
            else:
                raise ValueError('Invalid input for num_loops')

        dr = np.squeeze(dr)
        return dr

    def get_dd(self):

        n_models = len(self.cosmos_list)
        dd = np.zeros((n_models, 2, self.bins.num_bins('s')))

        theta = self.bins.bins('theta')
        theta = 0.5*(theta[:-1] + theta[1:])

        s = self.bins.bins('s')

        for i, cosmo in enumerate(self.cosmos_list):

            # convert z to r
            r = cosmo.z2r(self.bins.bins('z'))
            r = 0.5 * (r[:-1] + r[1:])

            # calculate weighted and unweighted
            for j in range(2):
                for k, pt_theta in enumerate(theta):
                    for l, pt_r in enumerate(r):
                        w = self.zztheta[j, k, l]
                        dist = np.sqrt(pt_r**2 + r**2 - 2*r*pt_r * np.cos(pt_theta))
                        hist, _ = np.histogram(dist,
                                               bins    = s,
                                               weights = w)
                        dd[i][j] += hist

        dd = np.squeeze(dd)
        return dd
