from setuptools import setup, find_packages

setup(
    name='GSForge',
    version='0.2',
    packages=find_packages(),
    url='',
    license='',
    author='Tyler Biggs',
    author_email='tyler.biggs@wsu.edu',
    description='Feature (gene) selection package for gene expression data.',
    install_requires=['boruta', 'umap_learn', 'h5py', 'netcdf4',
                      'xarray', 'pandas', 'numpy', 'lightgbm',
                      'param', 'scipy', 'scikit-learn', 'click',
                      'joypy', 'methodtools', 'dask[dataframe]', 'cloudpickle'],
)

