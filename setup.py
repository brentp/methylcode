from setuptools import setup
from distutils.extension import Extension
import numpy as np

#import methylcoder

ext_modules = [ Extension("methylcoder/cbowtie",
                sources=["methylcoder/cbowtie.c"], 
                include_dirs=[np.get_include(), "methylcode"]) ]

setup(
    license="BSD",
    name = "methylcoder",
    #version = methylcoder.__version__,
    ext_modules = ext_modules,
    packages=['methylcoder'],
    zip_safe=False,
    install_requires=['numpy', 'pyfasta'],
    entry_points = {
        'console_scripts': ['methylcoder = methylcoder:main']
    }
)
