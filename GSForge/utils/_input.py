import xarray as xr
import pandas as pd
import csv

__all__ = [
    "sniff_csv_sep",
    "infer_xarray_variables",
    "infer_pandas_columns",
    "load_count_df",
    "load_label_df",
]


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
