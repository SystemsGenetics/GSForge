"""
This module contains functions for preparing data for transfer to and from the R programming language.

This functionality is powered by the [`rpy2`](https://rpy2.readthedocs.io/en/version_2.8.x/) library.

This means:
+ The count and label dataframes will have the same indexes.
+ The shape will be transposed (genes as rows and samples as columns).
"""
import xarray as xr

__all__ = [
    "Py_counts_to_R",
    "Py_labels_to_R",
    "R_counts_to_Py_counts",
]


def Py_counts_to_R(counts: xr.DataArray):
    """Prepare a count `xarray.DataArray` as a `pandas.DataFrame` for transfer to the R programming language.

    Removes dashes from index names, as R does something(?) with them.
    """
    # Get a copy and stack the data so it is no longer flat.
    count_df = counts.transpose().copy().to_dataframe().unstack()
    # Since we are stacking a named array, the outer level index will remain.
    # Remove this index level via:
    count_df = count_df.droplevel(level=0, axis=1)
    return count_df


def Py_labels_to_R(labels: xr.DataArray):
    """
    :param labels:
    :return:
    """
    label_df = labels.copy().to_dataframe()
    return label_df


def R_counts_to_Py_counts(r_count_array, original_count_array):
    values = r_count_array.transpose()
    return xr.DataArray(values, coords=original_count_array.coords)
