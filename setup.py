"""
Installation script for GSForge.

This script checks if the ``$GSFORGE_MINIMAL`` environment variable exists,
if so the requirements installed are a reduced set, so as to make containers
made this way smaller.

This means that each workflow script made in this way should have its own
requirements.txt file containing any packages used (if you wish to make a
docker image for your script).
"""

import os
from setuptools import setup, find_packages

# These are the packages used by the classes in ''GSForge.models''.
requirements = """
numpy
pandas
xarray
param
""".split()

viz_requirements = """
Boruta
bokeh
click
datashader
h5py
holoviews
jupyter
matplotlib
methodtools
netcdf4
panel
scikit-learn
scipy
seaborn
statsmodels
tqdm
umap_learn
""".split()


if "GSFORGE_MINIMAL" not in os.environ:
    requirements += viz_requirements

# requirements = """\
# Boruta
# bokeh
# click
# datashader
# h5py
# holoviews
# jupyter
# matplotlib
# methodtools
# netcdf4
# numpy
# pandas
# panel
# param
# scikit-learn
# scipy
# seaborn
# statsmodels
# tqdm
# umap_learn
# xarray""".split()

setup(
    name='GSForge',
    version='0.2',
    packages=find_packages(),
    url='https://systemsgenetics.github.io/GSForge/',
    license='LICENSE.txt',
    author='Tyler Biggs',
    author_email='tyler.biggs@wsu.edu',
    description='Feature (gene) selection package for gene expression data.',
    python_requires='>3.6',
    install_requires=requirements,
)
