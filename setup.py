import sys
from distutils.core import setup
from setuptools import find_packages

# TODO add loggers and warnings
# TODO lazy module importing (https://github.com/bwesterb/py-demandimport)

if sys.version_info < (3, 6):
    print("SciLK requires Python >= 3.6")
    sys.exit(1)

# from Cython.Build import cythonize
#
# os.environ['CFLAGS'] = '-O3 -Wall'

setup(
    name="biomisc",
    packages=find_packages(),
    scripts=['primercut.py'],
    install_requires=[
        'fn',
        'biopython',
        'multipledispatch',
        'numba',
        'numpy',
        'hypothesis',
        'click',
        'pandas',
        'regex'
    ]
)
