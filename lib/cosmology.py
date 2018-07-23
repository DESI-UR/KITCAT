""" Convert redshift into comoving distance given cosmology

AstroPy has internally defined cosmologies based on recent CMB missions
(see http://docs.astropy.org/en/stable/cosmology/#built-in-cosmologies).
We'll use those parameters here to define comoving distance, but with
two changes:

1) Define H0 = 100h km/s/Mpc; astropy default uses the measured value of H0.
2) Allow for the possibile definition of nonzero curvature (not by default).
"""

# Python modules
import numpy
from astropy import cosmology
from scipy import interpolate


class Cosmology():
    """ Class to manage cosmological parameters and convert redshift z into
        comoving distance using linear interpolation technique. """
    def __init__(self, cosmo=None):
        """ Initialize a cosmological model. Note that H0 = 100 h km/s/Mpc
        is used. For now the cosmological parameters measured by Planck
        (P.A.R. Ade et al., Paper XIII, A&A 594:A13, 2016) are used.
        """
        self.model = None
        self.comoving_table = None
        self.set_model(cosmo)
        self._z2r = interpolate.pchip(*self.comoving_table.T)
        self._r2z = interpolate.pchip(*self.comoving_table[:, ::-1].T)

    def _set_comoving_table(self):
        """ Set redshift-comoving table """
        # Default parameters
        z_min = 0.
        z_max = 3.0
        step = 0.00005

        # Set table
        n =  int(numpy.ceil((z_max-z_min)/step))
        z = numpy.linspace(z_min, z_max, n)
        r = self.model.comoving_distance(z)
        self.comoving_table = numpy.array([z, r]).T

    def set_model(self, cosmo=None):
        """ Read cosmologies from configuration file and reset table
        Inputs:
        + cosmo: dictionary
            Dictionary with comsology paramaters
            Key: hubble0, omega_m0, omega_b0, omega_de0, temp_cmb, nu_eff
        """
        # Set astropy cosmology model
        if cosmo is None:
            self.model = cosmology.LambdaCDM(H0=100, Om0=0.307, Ode0=0.693)
        else:
            params = {}
            for key, val in cosmo.items():
                if key in ['hubble0', 'omega_m0', 'omega_de0']:
                    params[key] = float(val)
            self.model = cosmology.LambdaCDM(
                H0=params['hubble0'], Om0=params['omega_m0'], Ode0=params['omega_de0'],)

        # Set up redshift-comoving table
        self._set_comoving_table()

    def z2r(self, z):
        """Convert redshift to comoving distance by linear interpolating from table.
        Inputs:
        + z: list, tuple, ndarray, float
            Redshift 0 < z < 3.0.
        Outputs:
        + r: list, tuple, ndarray, float
            Return comoving distance given set cosmology. """

        # Convert list and tuple into numpy array
        if isinstance(z, (list, tuple)):
            z = numpy.asarray(z)

        # Check if input exceeds limit of table
        if numpy.all(z < 0) or numpy.all(z > 3.0):
            raise ValueError('Redshift must be between 0 and 3.0')

        r = self._z2r(z)
        if isinstance(z, numpy.ndarray):
            return r
        # Avoid 0-dimensional array
        return float(r)

    def r2z(self, r):
        """ Convert comoving distance to redshift by linear interpolating from table.
        Inputs:
        + r: list, tuple, ndarray, float
            Comoving distance within limit of table
        Outputs:
        + z: list, tuple, ndarray, float
            Return redshift given set cosmology. """
        # Convert list and tuple into numpy array
        if isinstance(r, (list, tuple)):
            r = numpy.asarray(r)

        # Check if input exceeds limit of table
        r_min = self.comoving_table[:, 1][0]
        r_max = self.comoving_table[:, 1][-1]
        if numpy.all(r < r_min) or numpy.all(r > r_max):
            raise ValueError('Comoving distance exceeds limit of table.')

        z = self._r2z(r)
        if isinstance(r, (list, tuple, numpy.ndarray)):
            return z
        # Avoid 0-dimensional array
        return float(z)

    def dels_to_delz(self, delta_s, z_min):
        """ Compute delta z given delta s and the lower redshift limit. """
        return self.r2z(self.z2r(z_min)+delta_s)-z_min

    @staticmethod
    def min_cosmo(cosmo_list, return_index=False):
        """ Find cosmology with the minimum value of Omega_m0 """
        min_value = 100000 # sentinel value
        min_cosmo = None
        index = None
        for i, cosmo in enumerate(cosmo_list):
            if cosmo.model.Om0 <= min_value:
                min_value = cosmo.model.Om0
                min_cosmo = cosmo
                index = i
        if return_index:
            return min_cosmo, index
        return min_cosmo

    @staticmethod
    def max_cosmo(cosmo_list, return_index=False):
        """ Find cosmology with the maximum value of Omega_m0 """
        max_value = -100000 # sentinel value
        max_cosmo = None
        index = None
        for i, cosmo in enumerate(cosmo_list):
            if cosmo.model.Om0 >= max_value:
                max_value = cosmo.model.Om0
                max_cosmo = cosmo
                index = i
        if return_index:
            return max_cosmo, index
        return max_cosmo
