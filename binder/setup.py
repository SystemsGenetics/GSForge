from setuptools import setup, find_packages

setup(
    name='GSForge',
    version='0.2',
    packages=find_packages(),
    url='https://systemsgenetics.github.io/GSForge/',
    license='LICENSE.txt',
    author='Tyler Biggs',
    author_email='tyler.biggs@wsu.edu',
    description='Feature (gene) selection package for gene expression data.',
    install_requires="""\
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
numpy
pandas
panel
param
scikit-learn
scipy
seaborn
statsmodels
tqdm
umap_learn
xarray""".split(),
    #
    # ['boruta', 'umap_learn', 'h5py', 'netcdf4',
    #                   'xarray', 'pandas', 'numpy', 'lightgbm',
    #                   'param', 'scipy', 'scikit-learn', 'click',
    #                   'methodtools', 'datashader', 'seaborn'],
)

