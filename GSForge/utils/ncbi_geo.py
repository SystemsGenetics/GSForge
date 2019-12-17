"""
Convert SOFT files from the NCBI: Geo Database and convert them into xarray.Dataset objects.
"""

import warnings
from ..utils._models import xrarray_gem_from_pandas


def soft_to_pandas(soft_object):
    """Converts a file downloaded from the NCBI gene expression omnibus to a pair
    of `pandas.DataFrame` objects."""
    count_df = soft_object.table.loc[:, soft_object.columns.index]
    count_df["Gene"] = soft_object.table["ID_REF"]

    if count_df["Gene"].dtype != "object":
        count_df["Gene"] = "Gene_" + count_df["Gene"].astype("str")
    count_df = count_df.set_index("Gene")
    label_df = soft_object.columns
    return count_df, label_df


def build_geo_xarray(geo_id: str, destdir: str = "./"):
    """Construct an xarray.Dataset object from a NCBI GEO DataSet Record."""
    try:
        import GEOparse
    except ImportError:
        warnings.warn('You must have the "GEOparse" package installed.', ImportWarning)

    soft_object = GEOparse.get_GEO(geo=geo_id, destdir=destdir)
    count_df, label_df = soft_to_pandas(soft_object)
    dataset = xrarray_gem_from_pandas(count_df, label_df)
    dataset = dataset.assign_attrs(soft_object.metadata)
    return dataset
