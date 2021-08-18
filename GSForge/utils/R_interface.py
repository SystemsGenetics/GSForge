"""
This module contains functions for preparing data for transfer to and from the R programming language.

Such functionality is powered by the `rpy2 <https://rpy2.readthedocs.io/en/version_2.8.x/>`_ library.

No conversion function is needed for preparing labels, as one can use the builtin ``pandas.DataFrame``
function ``to_dataframe``.::

    label_df = labels.to_dataframe()

"""

import xarray as xr

__all__ = [
    "Py_counts_to_R",
    "R_counts_to_Py",
    "Py_labels_to_R",
]


def Py_counts_to_R(counts: xr.DataArray):
    """
    Prepare a count ``xarray.DataArray`` as a ``pandas.DataFrame`` for transfer to the R programming language.

    This function transposes the data to have genes as rows and samples as columns.
    It then converts to a ``pandas.DataFrame`` and removes extraneous index levels.

    The inverse of this function is ``R_counts_to_Py_counts``.

    :param counts:
        An ``xr.DataArray`` count matrix.

    :returns count_df:
        A ``pandas.DataFrame`` ready to transfer to an R environment.
    """
    # Get a copy and stack the data so it is no longer flat.
    count_df = counts.transpose().copy().to_dataframe().unstack()
    # Since we are stacking a named array, the outer level index will remain.
    # Remove this index level via:
    count_df = count_df.droplevel(level=0, axis=1)

    # Ensure the samples are returned in their original order.
    count_df = count_df[counts[counts.dims[0]].values]
    return count_df


def R_counts_to_Py(r_count_array, original_count_array):
    """
    Prepares a ``numpy`` array (count matrix) for use in ``GSForge``.

    This function transposes the data (so that it has samples as rows and genes as columns).

    Inverts the conversion provided by ``Py_counts_to_R``.

    :param r_count_array:
        A ``numpy`` array of count values. Presumed to be oriented with genes as rows,
        and samples as columns.

    :param original_count_array:
        A copy of the original count array from which coordinates will be drawn.

    :return:
        An ``xarray.DataArray`` of the count values.
    """
    values = r_count_array.transpose()
    return xr.DataArray(values, coords=original_count_array.coords)


def Py_labels_to_R(label_ds):
    df = label_ds.to_dataframe()
    df.index.name = None
    return df
