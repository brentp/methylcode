from setuptools import setup
from distutils.extension import Extension
import numpy as np

import os.path as op
import sys
sys.path.insert(0, op.dirname(__file__))

import methylcoder.version as V

ext_modules = [ Extension("methylcoder/cbowtie",
                sources=["methylcoder/cbowtie.c"],
                include_dirs=[np.get_include(), "methylcoder"]) ]

setup(
    license="BSD",
    name = "methylcoder",
    version = V.version,
    ext_modules = ext_modules,
    packages=['methylcoder'],
    zip_safe=False,
    requires=['numpy', 'pyfasta'],
    test_suite="nose.collector",
    entry_points = {
        'console_scripts': ['methylcoder = methylcoder:main']
    }
)
