""" Modules for handling simple I/O functions"""

import os
import pickle
import configparser

from lib.cosmology import Cosmology

def parse_config(config_file, section):
    """ parse config file into a dictionary """
    parser = configparser.RawConfigParser(os.environ)
    parser.read(config_file)

    params_dict = {}

    if section in ["GALAXY_1", "GALAXY_2", "RANDOM_1", "RANDOM_2"]:
        params_dict['path'] = parser.get(section, 'path')
        params_dict['ra'] = parser.get(section, 'ra')
        params_dict['dec'] = parser.get(section, 'dec')
        params_dict['z'] = parser.get(section, 'z')
        try:
            params_dict['weight'] = parser.get(section, 'weight')
        except:
            params_dict['weight_fkp'] = parser.get(section, 'weight_fkp')
            params_dict['weight_noz'] = parser.get(section, 'weight_noz')
            params_dict['weight_sdc'] = parser.get(section, 'weight_sdc')
            params_dict['weight_cp'] = parser.get(section, 'weight_cp')
    elif section == "NBINS":
        params_dict['ra'] = parser.getint(section, 'ra')
        params_dict['dec'] = parser.getint(section, 'dec')
        params_dict['theta'] = parser.getint(section, 'theta')
        params_dict['z'] = parser.getint(section, 'z')
        params_dict['s'] = parser.getint(section, 's')
    elif section == "LIMIT":
        params_dict['unit'] = parser.get(section, 'unit')
        params_dict['ra_max'] = parser.getfloat(section, 'ra_max')
        params_dict['dec_max'] = parser.getfloat(section, 'dec_max')
        params_dict['ra_min'] = parser.getfloat(section, 'ra_min')
        params_dict['dec_min'] = parser.getfloat(section, 'dec_min')
        params_dict['z_max'] = parser.getfloat(section, 'z_max')
        params_dict['z_min'] = parser.getfloat(section, 'z_min')
        params_dict['s_max'] = parser.getfloat(section, 's_max')
    elif section == "COSMOLOGY":
        def getlist(parser, section, key, dtype):
            val = parser.get(section, key)
            if ',' in val:
                val = [v.strip() for v in val.split(',')]
                if dtype == str:
                    return val
                return list(map(dtype, val))
            return [dtype(val)]
        params_dict['hubble0'] = getlist(parser, section, 'hubble0', float)
        params_dict['omega_m0'] = getlist(parser, section, 'omega_m0', float)
        params_dict['omega_de0'] = getlist(parser, section, 'omega_de0', float)
        params_dict['n_cosmos'] = len(params_dict['hubble0'])

    return params_dict

def load(fname):
    """ load pickle """
    with open(fname,'rb') as f:
        save_object = pickle.load(f)
    return save_object


def save(fname, save_object):
    """ save pickle """
    with open(fname, 'wb') as f:
        pickle.dump(save_object, f, protocol=-1)


