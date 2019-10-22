"""
Utility functions for GEMprospector.
"""

# Standard library imports.
import csv
import numpy as np
import json

# Data science imports.
import xarray as xr
import pandas as pd
import inspect


def download_demo_data():
    pass


def get_by_json_attr(xr_dataset, key_pair):
    main_key, json_key = key_pair
    json_str = xr_dataset.attrs[main_key]
    json_attrs = json.loads(json_str)
    return json_attrs.get(json_key)


def filter_by_json_attr(xr_dataset, key_pair, req_value):
    if get_by_json_attr(xr_dataset, key_pair) == req_value:
        return True
    else:
        return False


def params_to_json(parameterized_instance, skip: list = None):
    """Build a .json string of the string parameters from a given instance of a
    `param.Parameterized` subclass.

    :param parameterized_instance:
    :param skip:
    :return:
    """
    values = {key: value for key, value in parameterized_instance.get_param_values()
              if isinstance(value, str) and value not in skip}
    return json.dumps(values)


def clean_nans(scores):
    """
    Fixes Issue #1240: NaNs can't be properly compared, so change them to the
    smallest value of scores's dtype. -inf seems to be unreliable.
    """
    # XXX where should this function be called? fit? scoring functions
    # themselves?
    scores = scores.copy()
    scores[np.isnan(scores)] = np.finfo(scores.dtype).min
    return scores


def kwargs_overlap(paramaterized, func):
    key_set = set(inspect.signature(func).parameters.keys()
                  ).intersection(set(paramaterized.param.objects().keys()))
    return {key: getattr(paramaterized, key) for key in key_set}


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
                          if var not in skip]

    incomplete = [var for var in selected_variables
                  if xr_dataset[var].to_series().isnull().values.any()]

    discrete_variables = [var for var in selected_variables if
                          xr_dataset[var].dtype == "object"]

    if quantile_size is not None:
        quantile_variables = [var for var in selected_variables
                              if xr_dataset[var].to_series().nunique() < quantile_size]
        continuous_variables = [var for var in selected_variables
                                if var not in quantile_variables]

        return {'all_labels': selected_variables, 'discrete': discrete_variables,
                'continuous': continuous_variables, 'quantile': quantile_variables,
                'incomplete': incomplete}

    else:
        continuous_variables = [var for var in selected_variables
                                if var not in discrete_variables]
        return {'all_labels': selected_variables, 'discrete': discrete_variables,
                'continuous': continuous_variables, 'incomplete': incomplete}


def infer_columns(df: pd.DataFrame,
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
