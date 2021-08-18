import csv

import numpy as np
import pandas as pd
import xarray as xr


def sniff_csv_sep(file_path: str) -> str:
    """Determine and return the separator used in a given .csv file.

    Works by sniffing the first line of the given file.

    :param file_path:
        The file path of the .csv file.

    :return:
        Delimiter used in the given .csv file.

    """

    with open(file_path, 'r') as csv_file:
        dialect = csv.Sniffer().sniff(csv_file.readline())

    return dialect.delimiter


def infer_xarray_variables(xr_dataset: xr.Dataset,
                           quantile_size: int = 10,
                           skip=None):
    """Build a dictionary of inferred column categories.

    :param xr_dataset:
        The xarray.Dataset to be examined.
    :param quantile_size:
        Unique count of floats and integers to be considered
        as categorical values instead of continuous.
    :param skip:
        The names of variables to skip. This should include the name
        of the gene expression matrix.

    :return:
        A dictionary of {type: variable_names}. Where 'type' is the
        inferred type.
    """
    if skip is None:
        skip = []

    selected_variables = [var for var in xr_dataset.variables
                          if var not in skip
                          and len(xr_dataset[var].dims) >= 1]  # 0 dim variables are problematic.

    incomplete = [var for var in selected_variables
                  if xr_dataset[var].to_series().isnull().values.any()]

    discrete_variables = [var for var in selected_variables if
                          xr_dataset[var].dtype == "object"]

    if quantile_size is not None:
        quantile_variables = [var for var in selected_variables
                              if xr_dataset[var].to_series().nunique() < quantile_size]
        continuous_variables = [var for var in selected_variables
                                if var not in quantile_variables
                                and np.issubdtype(xr_dataset[var].dtype, np.number)]

        return {'all_labels': selected_variables, 'discrete': discrete_variables,
                'continuous': continuous_variables, 'quantile': quantile_variables,
                'incomplete': incomplete}

    else:
        continuous_variables = [var for var in selected_variables
                                if var not in discrete_variables]
        return {'all_labels': selected_variables, 'discrete': discrete_variables,
                'continuous': continuous_variables, 'incomplete': incomplete}


def infer_pandas_columns(df: pd.DataFrame,
                         quantile: bool = True,
                         quantile_size: int = 10,
                         skip=None) -> dict:
    """Build a dictionary of inferred column categories.

    These categories are: "quantile", "continuous" and an optional
    "quantile".

    :param df:
        DataFrame to infer columns from.

    :param quantile:
        Boolean control of weather to construct the 'quantile' set.

    :param quantile_size:
        The max number of unique items a column can have and still be
        assigned to the 'quantile' category.

    :return:
        A dictionary with lists of columns in each category.

    """
    if skip is None:
        skip = ["counts", "Sample"]

    selected_columns = sorted([col for col in df.columns
                               if col not in skip])

    discrete = [x for x in selected_columns if df[x].dtype == object]

    if quantile:
        quantile = [x for x in selected_columns if df[x].nunique() < quantile_size
                    if x not in discrete]
        continuous = [x for x in selected_columns if x not in discrete]
        return {'all_labels': selected_columns, 'discrete': discrete,
                'continuous': continuous, 'quantile': quantile}

    else:
        continuous = [x for x in selected_columns if x not in discrete]
        return {'all_labels': selected_columns, 'discrete': discrete,
                'continuous': continuous}


def load_count_df(count_path, sniff=True, **kwargs):
    """Load the gene expression matrix as a `pandas.DataFrame` object.

    :param count_path:
        The path to the gene expression matrix.

    :param sniff:
        Infer the separator of the file.

    :param kwargs:
        Arguments to be passed to `pandas.read_csv`.

    :return:
        The gene expression matrix as a `pandas.DataFrame`.

    """

    if sniff:
        sep = sniff_csv_sep(count_path)
        kwargs = {**{"sep": sep}, **kwargs}

    count_df = pd.read_csv(count_path, **kwargs)
    count_df = count_df.fillna(0)

    return count_df


def load_label_df(label_path: str, sniff: bool = True, **kwargs) -> pd.DataFrame:
    """Load the gene annotation data as a `pandas.DataFrame` object.

    :param label_path:
        The path to the gene annotation data.

    :param sniff:
        Infer the separator of the file.

    :param kwargs:
        Arguments to be passed to `pandas.read_csv`.

    :return:
        The gene annotation data as a `pandas.DataFrame`.
    """

    if sniff:
        sep = sniff_csv_sep(label_path)
        kwargs = {**{"sep": sep}, **kwargs}

    label_df = pd.read_csv(label_path, **kwargs)
    label_df.index.name = "Sample"

    return label_df


def xrarray_gem_from_pandas(count_df: pd.DataFrame,
                            label_df: pd.DataFrame = None,
                            transpose_count_df: bool = True) -> xr.Dataset:
    """Stitch together a gene expression and annotation DataFrames into
    a single `xarray.Dataset` object.

    :param count_df: The gene expression matrix as a `pandas.DataFrame`.
        This file is assumed to have genes as rows and samples as columns.

    :param label_df: The gene annotation data as a `pandas.DataFrame`.
        This file is assumed to have samples as rows and annotation observations
        as columns.

    :param transpose_count_df: Whether to transpose the count matrix from the typical creation
        format (genes as rows, samples as columns) to the more standard (samples as rows,
        observations as columns).

    :return: An `xarray.Dataset` containing the gene expression matrix and
        the gene annotation data.
    """
    # TODO: Add input validation.
    # Shape validation.
    # Index matching validation.
    # Index name validation.
    count_array = xr.Dataset(
        {"counts": (("Gene", "Sample"), count_df.values)},
        coords={
            "Sample": count_df.columns.values,
            "Gene": count_df.index.values
        }
    )

    if transpose_count_df is True:
        count_array = count_array.transpose()

    if label_df is None:
        return count_array

    else:
        label_df.index.name = "Sample"
        label_ds = label_df.to_xarray()
        return label_ds.merge(count_array)
