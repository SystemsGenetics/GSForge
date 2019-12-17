import json
import pandas as pd
import xarray as xr
import inspect

__all__ = [
    "get_by_json_attr",
    "filter_by_json_attr",
    "params_to_json",
    "kwargs_overlap",
    "xrarray_gem_from_pandas",
]


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


def kwargs_overlap(paramaterized, func):
    """Gets the intersection between the parameters of a paramaterized model and
     the keyword arguments of a function, and returns the current intersection
     and their values as a dictionary."""
    key_set = set(inspect.signature(func).parameters.keys()
                  ).intersection(set(paramaterized.param.objects().keys()))
    return {key: getattr(paramaterized, key) for key in key_set}


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
