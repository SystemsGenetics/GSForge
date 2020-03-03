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
core_requirements = """
numpy
pandas
xarray
param
""".split()

doc_requirements = """
nbsite
nbsphinx
selenium
sphinx_ioam_theme
""".split()

full_requirements = """
Boruta
bokeh
click
datashader
h5py
holoviews
jupyter
matplotlib
graphviz
methodtools
netcdf4
panel
scikit-learn
scipy
seaborn
statsmodels
umap_learn
""".split()


installation_modes = {
    "normal": core_requirements + full_requirements,
    "docs": core_requirements + doc_requirements + full_requirements,
    "container": core_requirements
}

requirements = installation_modes[os.environ.get("GSFORGE_INSTALL_MODE", "normal")]

setup(
    name='GSForge',
    version='0.7.1',
    packages=find_packages(),
    url='https://systemsgenetics.github.io/GSForge/',
    license='LICENSE.txt',
    author='Tyler Biggs',
    author_email='tyler.biggs@wsu.edu',
    description='Feature (gene) selection package for gene expression data.',
    python_requires='>=3.8',
    install_requires=requirements,
)
