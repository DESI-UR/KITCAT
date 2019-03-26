#!/usr/bin/env python

# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function

#
# Standard imports
#
import glob
import os
import sys
import pip

#
# setuptools' sdist command ignores MANIFEST.in
#
from distutils.command.sdist import sdist as DistutilsSdist
from setuptools import setup, find_packages
from setuptools.command.install import install as InstallCommand
from py.KITCAT import versioning as ver

class Install(InstallCommand):
    """ Customized setuptools install command which uses pip."""
    def run(self, *args, **kwargs):
        print("called function with {}".format(args))
        try:
            # for older versions of pip
            from pip._internal import main
        except:
            # for newer versions of pip
            from pip import main
        main(['install', '.'])
        InstallCommand.run(self, *args, **kwargs)

def pip_install(package):
    if hasattr(pip, 'main'):
        pip.main(['install', package])
    else:
        pip._internal.main(['install', package])
        
#
# Begin setup
#
setup_keywords = dict()
#
# THESE SETTINGS NEED TO BE CHANGED FOR EVERY PRODUCT.
#
setup_keywords['name'] = 'KITCAT'
setup_keywords['description'] = 'Kd-tree Implementation for Two-point Correlation AlgoriThm'
setup_keywords['author'] = 'Tri Nguyen, Tolga Yapici'
setup_keywords['author_email'] = 'tnguy51@u.rochester.edu, tyapici@ur.rochester.edu'
setup_keywords['license'] = 'BSD'
setup_keywords['url'] = 'https://github.com/DESI-UR/KITCAT'
setup_keywords['version'] = ver.get_version(out_type='string')

print("Version is set to {}".format(ver.get_version(out_type='string')))

# this is a dirty way of installing the dependencies
# but when I let everything handled by pip/setuptools,
# the wheel generation spawns thousands of processes and kills
# the computer
print("Installing necessary packages with pip")
fid = open('requirements.txt')
for line in fid:
    print(line)
    pip_install(line)

#
# END OF SETTINGS THAT NEED TO BE CHANGED.
#

# Set other keywords for the setup function.  These are automated, & should
# be left alone unless you are an expert.
#
# Treat everything in bin/ except *.rst as a script to be installed.
#
if os.path.isdir('bin'):
    setup_keywords['scripts'] = [fname for fname in glob.glob(os.path.join('bin', '*'))
        if not os.path.basename(fname).endswith('.rst')]

setup_keywords['provides'] = [setup_keywords['name']]
#setup_keywords['setup_requires'] = [']
#setup_keywords['install_requires'] = ['numpy>=1.13.1', 'healpy>=1.11.0', 'configparser>=3.5', \
#                                      'astropy>=1.2.1', 'scipy>=0.19.1', 'matplotlib>=2.0.0', \
#                                      'scikit-learn>=0.18.1', 'pip>19.0.0']
setup_keywords['zip_safe'] = False
setup_keywords['use_2to3'] = True
setup_keywords['packages'] = find_packages('py')
setup_keywords['package_dir'] = {'': 'py'}
#setup_keywords['cmdclass'] = {'install': Install,}

# Add internal data directories.
#
#setup_keywords['package_data'] = {'KITCAT': ['data/*',]}

print(setup_keywords)

# Run setup command.
#
setup(**setup_keywords)
