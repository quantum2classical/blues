from __future__ import absolute_import, division, print_function
from os.path import join as pjoin

# Format expected by setup.py and doc/source/conf.py: string of form "X.Y.Z"
_version_major = 0
_version_minor = 0
_version_micro = 1  # use '' for first of series, number for 1 and above
_version_extra = 'dev1'
# _version_extra = ''  # Uncomment this for full releases

# Construct full version string from these.
_ver = [_version_major, _version_minor]
if _version_micro:
    _ver.append(_version_micro)
if _version_extra:
    _ver.append(_version_extra)

__version__ = '.'.join(map(str, _ver))

CLASSIFIERS = ["Development Status :: 3 - Alpha",
               "Environment :: Console",
               "Intended Audience :: Science/Research",
               "License :: OSI Approved :: MIT License",
               "Operating System :: OS Independent",
               "Programming Language :: Python",
               "Topic :: Scientific/Engineering"]

# Description should be a one-liner:
description = "BLUES: Bond Locator Utilizing Electronic Structure"
# Long description will go up on the pypi page
long_description = """

BLUES
========
BLUES is a pipeline written in Python that runs JANPA with given MOLDEN files
and outputs an accurate InChI or SMILES string.

License
=======
``blues`` is licensed under the terms of the MIT license. See the file
"LICENSE" for information on the history of this software, terms & conditions
for usage, and a DISCLAIMER OF ALL WARRANTIES.

All trademarks referenced herein are property of their respective holders.

Copyright (c) 2018--, Luke D Gibson, The University of Washington
Department of Chemical Engineering.
"""

NAME = "blues"
MAINTAINER = "Luke D Gibson"
MAINTAINER_EMAIL = "ldgibson@uw.edu"
DESCRIPTION = description
LONG_DESCRIPTION = long_description
URL = "http://github.com/quantum2classical/blues"
DOWNLOAD_URL = ""
LICENSE = "MIT"
AUTHOR = "Luke D Gibson"
AUTHOR_EMAIL = "ldgibson@uw.edu"
PLATFORMS = "OS Independent"
MAJOR = _version_major
MINOR = _version_minor
MICRO = _version_micro
VERSION = __version__
PACKAGE_DATA = {'blues': [pjoin('janpa', '*.jar')]}
REQUIRES = ['numpy']
