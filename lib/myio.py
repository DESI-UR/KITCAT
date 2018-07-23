""" Modules for handling simple I/O functions"""

import pickle
from lib.cosmology import Cosmology

def load(fname):
    """ Load pickle """
    with open(fname, "rb") as f:
        while True:
            try:
                yield pickle.load(f)
            except EOFError:
                break

def save(fname, *save_list):
    """ Pickle a list of objects """
    pickle_out = open(fname, 'wb')
    for save_object in save_list:
        pickle.dump(save_object, pickle_out, protocol=-1)
    pickle_out.close()

def read_cosmology(params):
    """ Read and set up cosmology models from a configuration section """
    # Temporary dictionary to store parameters
    models_params = {'hubble0': [], 'omega_m0': [], 'omega_de0': []}
    n = -1  # keep track of the number of models
    for key, val in params.items():
        if key in models_params.keys():
            # Convert to float
            pars = [float(p) for p in val.split(', ')]
            models_params[key] = pars

            # Check number of models
            n = max(len(pars), 0)
            if n != len(pars):
                raise ValueError('Number of parameters inconsistent.')

    # Initialize cosmology models
    print('  - Number of models: %d' %n)
    models = []
    for i in range(n):
        temp = {}
        for key, val in models_params.items():
            temp[key] = val[i]
        models.append(Cosmology(temp))

    return models
