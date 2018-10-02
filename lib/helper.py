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
            raise ValueError('total jobs must be at least 1.')
        self.total_jobs = total_jobs
        self.current_job = 0

    def increment(self, verbose=True):
        """ Increment current job index by 1 """
        if self.current_job != self.total_jobs - 1:
            self.current_job += 1
        else:
            print('already at the last job index.')

        if verbose:
            print("job %d/%d." % (self.current_job, self.total_jobs))

    def set_current_job(self, current_job, verbose=True):
        """ Set current job index to input """
        if current_job < 0 or current_job >= self.total_jobs:
            raise ValueError('job must be at least 0 and less than total job.')
        self.current_job = current_job

        if verbose:
            print("job %d/%d." % (self.current_job, self.total_jobs))

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
    def __init__(self, bins, models):
        """ Constructor sets up attributes """

        # Initialize binnings
        self.bins = bins

        # Initalize norm factor
        self.norm = {"dd": None, "dr": None, "rr": None}

        # Initialize histogram
        self.data_data = np.zeros((2, bins.num_bins('s')))
        self.theta_distr = np.zeros(bins.num_bins('theta'))
        self.z_theta_distr = np.zeros((2, bins.num_bins('theta'),
                                       bins.num_bins('z')))
        self.z_distr = None

        # Initalize cosmology models
        self.models = models

    def set_z_distr(self, z_distr):
        """ Setting up P(r) """
        self.z_distr = z_distr

    def set_norm(self, norm_dd, norm_dr, norm_rr):
        """ Setting up normalization """
        self.norm["dd"] = norm_dd
        self.norm["dr"] = norm_dr
        self.norm["rr"] = norm_rr

    def set_dd(self, tree, catalog, job_helper=None, same=True):
        """ Calculate DD(s) using a modified nearest-neighbors kd-tree search.
        Metric: Euclidean
        Inputs:
        + tree: sklearn.neighbors.BallTree or sklearn.neighbors.KDTree
            KD-tree built from data catalog
        + catalog: np.ndarray
            Data catalog in KD-tree (with weights)
        + job_helper: helper.JobHelper (default=None)
            Job manager to handle multiprocess indices.
            If None, assume one job.
        + same: bool (default=True)
            If True, catalog and tree.catalog are the same. Account for
            double counting.
        """

        # Reset previous value of dd
        self.data_data = np.zeros((2, self.bins.num_bins('s')))

        # If job_helper is None, assume one job
        if job_helper is None:
            job_helper = JobHelper(1)
            job_helper.set_current_job(0)

        # Get start and end indices
        start, end = job_helper.get_index_range(catalog.shape[0])

        s_max = self.bins.max('s')
        num_bins_s = self.bins.num_bins('s')
        print("- Calculate DD from index %d to %d" % (start, end - 1))
        for i, point in enumerate(catalog[start:end]):
            if i % 10000 is 0:
                print(' - Index: %d' % i)
            index, s = tree.query_radius(point[: 3].reshape(1, -1), r=s_max,
                                         return_distance=True)

            # Fill weighted distribution
            # weight is the product of the weights of two points
            weights = catalog[:, 3][index[0]]*point[3]
            hist, _ = np.histogram(s[0], bins=num_bins_s, range=(0., s_max),
                                   weights=weights)
            self.data_data[0] += hist

            # Fill unweighted distribution
            hist, _ = np.histogram(s[0], bins=num_bins_s, range=(0., s_max))
            self.data_data[1] += hist

        if same:
            # Correction for double counting in the first bin from pairing a galaxy
            # with itself
            self.data_data[0][0] -= np.sum(catalog[start:end, 3]**2)
            self.data_data[1][0] -= end - start

            # Correction for double counting
            self.data_data = self.data_data / 2.

    def set_theta_distr(self, tree, catalog, job_helper=None, same=True):
        """ Calculate f(theta) using a modified nearest-neighbors search
        BallTree algorithm. Metric = 'haversine'.
        Inputs:
        + tree: sklearn.neighbors.BallTree
            Balltree built from data catalog
        + catalog: np.ndarray
            Angular catalog in Balltree (with weights)
        + job_helper: helper.JobHelper (default=None)
            Job manager to handle multiprocess indices.
            If None, assume one job.
        + same: bool (default=True)
            If True, catalog and tree.catalog are the same. Account for
            double counting."""

        # Reset previous value of f(theta)
        self.theta_distr = np.zeros(self.bins.num_bins('theta'))

        # If job_helper is None, assume one job
        if job_helper is None:
            job_helper = JobHelper(1)
            job_helper.set_current_job(0)

        # Get start and end indices
        start, end = job_helper.get_index_range(catalog.shape[0])

        theta_max = self.bins.max('theta')
        nbins_theta = self.bins.num_bins('theta')
        print("- Construct f(theta) from index %d to %d" % (start, end - 1))
        for i, point in enumerate(catalog[start:end]):
            if i % 10000 is 0:
                print('  - Index: %d' % i)
            index, theta = tree.query_radius(point[:2].reshape(1, -1),
                                             r=theta_max,
                                             return_distance=True)
            # weight is the product of the weights of each point
            weights = point[2] * catalog[:, 2][index[0]]
            hist, _ = np.histogram(theta[0], bins=nbins_theta,
                                   range=(0., theta_max), weights=weights)
            self.theta_distr += hist

        if same:
            # Correction for double counting
            self.theta_distr = self.theta_distr / 2.

    def set_z_theta_distr(self, tree, data_catalog, angular_catalog, mode,
                          job_helper=None):
        """ Calculate g(theta, z) using a modified nearest-neighbors BallTree
        search. Metric = 'haversine'.
        NOTE: assume uniformed comoving bins
        Inputs:
        + tree: sklearn.neighbors.BallTree
            Balltree built from DEC and RA coordinates
        + data_catalog: np.ndarray
            Catalog from galaxy data (with weights)
        + angular_catalog: np.ndarray
            Angular catalog from random data (with weights)
        + job_helper: helper.JobHelper (default=None)
            Job manager to handle multiprocess indices.
            If None, assume one job. """

        # Reset previous values of g(theta, r)
        self.z_theta_distr = np.zeros((2, self.bins.num_bins('theta'),
                                       self.bins.num_bins('z')))

        # Get start and end indices
        if mode == "angular_tree":
            start, end = job_helper.get_index_range(data_catalog.shape[0])
        elif mode == "data_tree":
            start, end = job_helper.get_index_range(angular_catalog.shape[0])

        # Initialize some binning variables
        num_bins = (self.bins.num_bins('theta'), self.bins.num_bins('z'))
        theta_max = self.bins.max('theta')
        z_min = self.bins.min('z')
        z_max = self.bins.max('z')
        limit = ((0., theta_max), (z_min, z_max))

        print("- Construct angular-redshift from index %d to %d" % (start,
                                                                    end - 1))
        if mode == "angular_tree":
            for i, point in enumerate(data_catalog[start:end]):
                if i % 10000 is 0:
                    print(' - Index: %d' % i)
                index, theta = tree.query_radius(
                    point[:2].reshape(1, -1), r=theta_max, return_distance=True)
                z = np.repeat(point[2], index[0].size)

                # Fill unweighted histogram
                weights = angular_catalog[:, 2][index[0]]
                hist, _, _ = np.histogram2d(
                    theta[0], z, bins=num_bins, range=limit, weights=weights)
                self.z_theta_distr[1] += hist

                # Fill weighted histogram
                # weight is the product of the weight of the data point
                # and the weight of the angular point.
                weights = weights * point[3]
                hist, _, _ = np.histogram2d(
                    theta[0], z, bins=num_bins, range=limit, weights=weights)
                self.z_theta_distr[0] += hist

        elif mode == "data_tree":
            for i, point in enumerate(angular_catalog[start:end]):
                if i % 10000 is 0:
                    print(' - Index: %d' % i)
                index, theta = tree.query_radius(
                    point[:2].reshape(1, -1), r=theta_max, return_distance=True)
                z = data_catalog[:, 2][index[0]]

                # Fill weighted histogram
                # weight is the product of the weight of the data point
                # and the weight of the angular point
                weights = point[2] * data_catalog[:, 3][index[0]]
                hist, _, _ = np.histogram2d(
                    theta[0], z, bins=num_bins, range=limit, weights=weights)
                self.z_theta_distr[0] += hist

                # Fill unweighted histogram
                # weight is the weight of the angular point
                hist, _, _ = np.histogram2d(
                    theta[0], z, bins=num_bins, range=limit)
                self.z_theta_distr[1] += hist*point[2]

    def add(self, other):
        """ Add another part """
        self.data_data = self.data_data + other.data_data
        self.theta_distr = self.theta_distr + other.theta_distr
        self.z_theta_distr = self.z_theta_distr + other.z_theta_distr

    def get_dd(self):
        """ Return DD(s) """
        print("- Calculate DD(s)")
        return self.data_data

    def get_rr(self, num_loops=1):
        """ Calculate and return RR(s) """
        if self.z_distr is None:
            return RuntimeError("Redshift/comoving distribution is None.")

        print("- Calculate RR(s)")

        # Initialize separation distribution and binning
        n_models = len(self.models)
        rand_rand = np.zeros((n_models, 2, self.bins.num_bins('s')))

        # Angular separation bins
        bins_theta = self.bins.bins('theta')
        bins_theta = 0.5 * (bins_theta[:-1] + bins_theta[1:])

        # Calculate DD with weight
        print("  - Number of loops: %d" % num_loops)
        if num_loops == 0:
            # No explicit for loop in distance calculation
            # Fast calculation at the expense of memory
            # For loop over weights included due to memory
            for i in range(2):
                # Calculate 3-dimensional weights matrix
                weights = (self.z_distr[i][None, :, None] *
                           self.z_distr[i][None, None, :] *
                           self.theta_distr[:, None, None])

                # Loop over cosmology models
                for j, model in enumerate(self.models):
                    # Apply cosmology to calculate comoving bins
                    bins_r = model.z2r(self.bins.bins('z'))
                    bins_r = 0.5 * (bins_r[:-1] + bins_r[1:])

                    # Calculate 3-dimensional separation matrix
                    dist = distance(bins_theta[:, None, None],
                                    bins_r[None, :, None],
                                    bins_r[None, None, :])

                    # Calculate RR
                    rand_rand[j][i], _ = np.histogram(
                        dist, bins=self.bins.bins('s'), weights=weights)

        elif num_loops == 1:
            for i in range(2):
                if i == 0:
                    print(' - Weighted')
                else:
                    print(' - Unweighted')
                # Caculate 2-d weight matrix
                weights = (self.theta_distr[:, None] * self.z_distr[i][None, :])

                for j, model in enumerate(self.models):
                    # Apply cosmology to calculate comoving bins
                    bins_r = model.z2r(self.bins.bins('z'))
                    bins_r = 0.5 * (bins_r[:-1] + bins_r[1:])

                    # Build a 2-d distance matrices
                    temp = np.cos(bins_theta[:, None]) * bins_r[None, :]
                    bins_r_sq = bins_r[None, :]**2

                    # Calculate RR
                    for k, r in enumerate(bins_r):
                        if k % 100 is 0:
                            print('  - Index: %d' % k)
                        dist = np.sqrt(r**2 + bins_r_sq - 2*r*temp)
                        hist, _ = np.histogram(
                            dist, bins=self.bins.bins('s'),
                            weights=self.z_distr[i][k] * weights)
                        rand_rand[j][i] += hist
        else:
            raise ValueError('Invalid input for num_loops')

        rand_rand = np.squeeze(rand_rand)
        return rand_rand

    def get_dr(self, num_loops=1):
        """ Calculate and return DR(s) """
        if self.z_distr is None:
            raise RuntimeError("Comoving distribution is None.")

        print("- Calculate DR(s)")

        # Initialize separation distribution and binning
        n_models = len(self.models)
        data_rand = np.zeros((n_models, 2, self.bins.num_bins('s')))

        # Angular separation bins
        bins_theta = self.bins.bins('theta')
        bins_theta = 0.5 * (bins_theta[:-1] + bins_theta[1:])

        # Calculate DR weighted and unweighted
        print(" - Number of loops: %d" % num_loops)
        if num_loops == 0:
            # No explicit for loop in distance calculation
            # Fast calculation at the expense of memory
            # For loop over weights included due to memory
            for i in range(2):
                # Calculate 3-dimensional weights matrix
                weights = (self.z_distr[i][None, None, :] *
                           self.z_theta_distr[i][:, :, None])

                # Loop over cosmology
                for j, model in enumerate(self.models):
                    # Apply cosmology to calculate comoving bins
                    bins_r = model.z2r(self.bins.bins('z'))
                    bins_r = 0.5 * (bins_r[:-1]  + bins_r[1:])

                    # Calculate 3-dimensional separation matrix
                    dist = distance(bins_theta[:, None, None],
                                    bins_r[None, :, None],
                                    bins_r[None, None, :])

                    # Calculate DR
                    data_rand[j][i], _ = np.histogram(
                        dist, bins=self.bins.bins('s'), weights=weights)
        elif num_loops == 1:
            for i in range(2):
                if i == 0:
                    print(' - Weighted')
                else:
                    print(' - Unweighted')
                # Caculate 2-d weight matrix
                weights = self.z_theta_distr[i]

                for j, model in enumerate(self.models):
                    # Apply cosmology to calculate comoving bins
                    bins_r = model.z2r(self.bins.bins('z'))
                    bins_r = 0.5 * (bins_r[:-1] + bins_r[1:])

                    # Build a 2-d distance matrices
                    temp = np.cos(bins_theta[:, None]) * bins_r[None, :]
                    bins_r_sq = bins_r[None, :]**2

                    # Caculate DR
                    for k, r in enumerate(bins_r):
                        if k % 100 is 0:
                            print('  - Index: %d' % k)
                        dist = np.sqrt(r**2 + bins_r_sq - 2*r*temp)
                        hist, _ = np.histogram(
                            dist, bins=self.bins.bins('s'),
                            weights=self.z_distr[i][k] * weights)
                        data_rand[j][i] += hist

        data_rand = np.squeeze(data_rand)
        return data_rand
