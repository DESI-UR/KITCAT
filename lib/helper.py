""" Module with helper class """

import numpy as np

class JobHelper(object):
    """ class to handle multiprocess job """

    def __init__(self, total_jobs):
        """ constructor """
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
        """ calculate the start and end indices given job size """
        job_index = np.floor(np.linspace(0, size, self.total_jobs + 1))
        job_index = job_index.astype(int)
        job_range = (job_index[self.current_job],
                     job_index[self.current_job + 1])
        return job_range

class CorrelationHelper(object):
    """ class to handle multiprocess correlation function calculation """

    def __init__(self):
        """ constructor """

        # normalization
        self.norm_dd = 1.
        self.norm_rr = 1.
        self.norm_d1r2 = 1.
        self.norm_d2r1 = 1.

        # distribution
        self.zztheta = None
        self.ftheta = None
        self.ztheta_d1r2 = None
        self.ztheta_d2r1 = None
        self.z1_distr = None
        self.z2_distr = None

        # misc
        self.cosmos_list = None
        self.bins = None

    def add(self, other):
        self.zztheta += other.zztheta
        self.ftheta  += other.ftheta
        self.ztheta_d1r2 += other.ztheta_d1r2
        self.ztheta_d2r1 += other.ztheta_d2r1

    def get_rr(self):

        n_models = len(self.cosmos_list)
        rr = np.zeros((n_models, 2, self.bins.num_bins('s')))

        theta = self.bins.bins('theta')
        theta = 0.5*(theta[:-1] + theta[1:])

        s = self.bins.bins('s')

        # calculate rr
        # loop over cosmology models
        for i, cosmo in enumerate(self.cosmos_list):
            cosmo_params = list(cosmo.params.values())
            print('- h0, om0, ode0: %s' % cosmo_params)
            r = cosmo.z2r(self.bins.bins('z'))
            r = 0.5 * (r[:-1] + r[1:])

            for j in range(2):

                # calculate 2-d weight matrix
                w = self.ftheta[:, None] * self.z1_distr[j][None, :]

                # Build a 2-d distance matrices
                temp = np.cos(theta[:, None]) * r[None, :]
                r_sq = r[None, :]**2

                # Calculate RR(s)
                for k, pt_r in enumerate(r):
                    dist = np.sqrt(pt_r**2 + r_sq - 2*pt_r*temp)
                    hist, _ = np.histogram(dist,
                                           bins    = s,
                                           weights = self.z2_distr[j][k]*w)
                    rr[i][j] += hist

        return rr

    def get_dr(self, mode='r1'):

        if mode == 'r1':
            ztheta = self.ztheta_d2r1
            z_distr = self.z1_distr
        elif mode == 'r2':
            ztheta = self.ztheta_d1r2
            z_distr = self.z2_distr

        n_models = len(self.cosmos_list)
        dr = np.zeros((n_models, 2, self.bins.num_bins('s')))

        theta = self.bins.bins('theta')
        theta = 0.5*(theta[:-1] + theta[1:])

        s = self.bins.bins('s')

        # calculate dr
        # loop over cosmology models
        for i, cosmo in enumerate(self.cosmos_list):
            cosmo_params = list(cosmo.params.values())
            print('- h0, om0, ode0: %s' % cosmo_params)

            r = cosmo.z2r(self.bins.bins('z'))
            r = 0.5 * (r[:-1] + r[1:])

            for j in range(2):

                w = ztheta[j]

                # Build a 2-d distance matrix
                temp = np.cos(theta[:, None]) * r[None, :]
                r_sq = r[None, :]**2

                # Calculate RR(s)
                for k, pt_r in enumerate(r):
                    dist = np.sqrt(pt_r**2 + r_sq - 2*pt_r*temp)
                    hist, _ = np.histogram(dist,
                                           bins    = s,
                                           weights = z_distr[j][k]*w)
                    dr[i][j] += hist

        return dr

    def get_dd(self):

        n_models = len(self.cosmos_list)
        dd = np.zeros((n_models, 2, self.bins.num_bins('s')))

        theta = self.bins.bins('theta')
        theta = 0.5*(theta[:-1] + theta[1:])

        s = self.bins.bins('s')

        for i, cosmo in enumerate(self.cosmos_list):
            cosmo_params = list(cosmo.params.values())
            print('- h0, om0, ode0: %s' % cosmo_params)

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

        return dd
