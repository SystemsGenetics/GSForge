import functools
import json
import pathlib
from textwrap import dedent
from typing import List, TypeVar, Union

import pandas as pd
import param
import xarray as xr

from ._utils import (
    infer_xarray_variables,
    xrarray_gem_from_pandas,
    load_count_df,
    load_label_df,
)

# Declare typing hints.
# Typing hints are not used by the Python runtime in any way.
# They are used by third party tools and linters.
# In python 3.8 we can use forward references.
TypeAGEM = TypeVar("TypeAGEM", bound="AnnotatedGEM")
TypePath = Union[str, 'PathLike[Any]']


class AnnotatedGEM(param.Parameterized):
    """
    A data class for a gene expression matrix and any associated sample or
    gene annotations.

    This model holds the count expression matrix, and any associated labels
    or annotations as an ``xarray.Dataset`` object under the ``.data`` attribute.
    By default this dataset will be expected to have its indexes named "Gene"
    and "Sample", although there are parameters to override those arrays and
    index names used.

    An AnnotatedGEM object can be created with one of the class methods:

    ``from_files()``
      A helper function for loading disparate GEM and annotation files through
      pandas.read_csv().

    ``from_pandas()``
      Reads in a GEM pandas.DataFrame and an optional annotation DataFrame. These
      must share the same sample index.

    ``from_netcdf()``
      Reads in from a .nc filepath. Usually this means loading a previously
      created AnnotatedGEM.


    **Randomly generate a demo AnnotatedGEM**

    # >>> from sklearn.datasets import make_multilabel_classification
    # >>> data, labels = make_multilabel_classification()
    # >>> agem = AnnotatedGEM.from_pandas(pd.DataFrame(data), pd.DataFrame(labels), name="Generated GEM")


    # >>> agem
    # <GSForge.AnnotatedGEM>
    # Name: Generated GEM
    # Selected GEM Variable: 'counts'
    #     Gene   100
    #     Sample 100

    **View the entire gene or sample index:**

    # >>> agem.gene_index
    # <xarray.DataArray 'Gene' (Gene: 100)>...

    # >>> agem.sample_index
    # <xarray.DataArray 'Sample' (Sample: 100)>...

    # >>> agem.infer_variables()
    # {'all_labels': ...

    """

    data = param.ClassSelector(class_=xr.Dataset, allow_None=False, doc=dedent("""\
    An ``xarray.Dataset`` object that contains the Gene Expression Matrix, and any 
    needed annotations. This `xarray.Dataset` object is expected to have a count 
    array named 'counts', that has coordinates ('Gene', 'Sample')."""))

    count_array_name = param.String(default="counts", doc=dedent("""\
    This parameter controls which variable from the `xarray.Dataset` should be
    considered to be the 'count' variable.
    Consider using this if you require different index names, or wish to control 
    which count array among many should be used by default."""))

    sample_index_name = param.String(default="Sample", doc=dedent("""\
    This parameter controls which variable from the `xarray.Dataset` should be
    considered to be the 'sample' coordinate.
    Consider using this if you require different coordinate names."""))

    gene_index_name = param.String(default="Gene", doc=dedent("""\
    This parameter controls which variable from the `Xarray.Dataset` should be 
    considered to be the 'gene index' coordinate.
    Consider using this if you require different coordinate names."""))

    def __init__(self, *args, **params) -> None:
        if args:
            params = _annotated_gem_dispatch(*args, **params)
        # TODO: Ensure that a data object exists in params.
        #       Give an invalid input warning to the user.
        super().__init__(**params)

    def __repr__(self) -> str:
        """
        Display a summary of this AnnotatedGEM.
        """
        summary = [f"<GSForge.{type(self).__name__}>"]
        summary += [f"Name: {self.name}"]
        summary += [f"Selected GEM Variable: '{self.count_array_name}'"]
        summary += [f"    {self.gene_index_name}   {self.data[self.gene_index_name].shape[0]}"]
        summary += [f"    {self.sample_index_name} {self.data[self.sample_index_name].shape[0]}"]
        return "\n".join(summary)

    @property
    def gene_index(self) -> xr.DataArray:
        """
        Returns the entire gene index of this AnnotatedGEM object as an ``xarray.DataArray``.

        The actual variable or coordinate that this returns is controlled by the
        ``gene_index_name`` parameter.
        """
        return self.data[self.gene_index_name].copy(deep=True)

    @property
    def sample_index(self) -> xr.DataArray:
        """
        Returns the entire sample index of this AnnotatedGEM object as an ``xarray.DataArray``.

        The actual variable or coordinate that this returns is controlled by the
        ``sample_index_name`` parameter.
        """
        return self.data[self.sample_index_name].copy(deep=True)

    @property
    def count_array_names(self) -> List[str]:
        """
        Returns a list of all available count arrays contained within this AnnotatedGEM object.

        This is done simply by returning all data variables that have the same dimension set
        as the default count array.
        """
        default_dims = set(self.data[self.count_array_name].dims)
        return [var for var in self.data.data_vars if set(self.data[var].dims) == default_dims]

    def infer_variables(self, quantile_size: int = 10, skip: bool = None) -> dict:
        """
        Infer categories for the variables in the AnnotatedGEM's labels.

        :param quantile_size:
            The maximum number of unique elements before a variable is no
            longer considered as a `quantile-able` set of values.

        :param skip:
            The variables to be skipped.

        :return:
            A dictionary of the inferred value types.
        """
        if skip is None:
            skip = self.count_array_names + [self.sample_index_name, self.gene_index_name]

        sample_dim = {self.sample_index_name}
        gene_annots = [var for var in self.data.data_vars if set(self.data[var].dims) != sample_dim]
        if gene_annots:
            skip += gene_annots

        return infer_xarray_variables(xr_dataset=self.data, quantile_size=quantile_size, skip=skip)

    @staticmethod
    def _parse_xarray_dataset(data: xr.Dataset, **params) -> dict:
        """
        Parse arguments for AnnotatedGEM creation via an `xarray.Dataset`.

        :param data:
            An `xarray.Dataset` if this dataset has different index names than default
            (Gene, Sample, counts), be sure to explicitly set those parameters
            (`gene_index_name`, `sample_index_name`, `count_array_name`).

        :param params:
            Other parameters to set.

        :return:
            A parsed parameter dictionary.
        """
        existing_params = data.attrs.get("__GSForge.AnnotatedGEM.params")
        if existing_params:
            existing_params = json.loads(existing_params)
            params = {**existing_params, **params}
        return {"data": data, **params}

    @classmethod
    def _parse_netcdf_path(cls, netcdf_path, **params):
        """
        Parse arguments for AnnotatedGEM creation via a path to an `xarray.Dataset` saved as a .netcdf file.

        :param netcdf_path:
            A path to a `netcdf` file. If this file has different index names than default
            (Gene, Sample, counts), be sure to explicitly set those parameters
            (`gene_index_name`, `sample_index_name`, `count_array_name`).

        :param params:
            Other parameters to set.

        :return:
            A parsed parameter dictionary.
        """
        params = cls._parse_xarray_dataset(xr.open_dataset(netcdf_path), **params)
        return {"data": xr.open_dataset(netcdf_path), **params}

    @classmethod
    def from_netcdf(cls, netcdf_path: TypePath, **params) -> TypeAGEM:
        """
        Construct an `AnnotatedGEM` object from a `netcdf` (.nc) file path.

        :param netcdf_path:
            A path to a `netcdf` file. If this file has different index names than default
            (Gene, Sample, counts), be sure to explicitly set those parameters
            (`gene_index_name`, `sample_index_name`, `count_array_name`).

        :return:
            A new AnnotatedGEM.
        """
        params = cls._parse_xarray_dataset(xr.open_dataset(netcdf_path), **params)
        return cls(**params)

    @classmethod
    def _parse_pandas(cls, count_df: pd.DataFrame, label_df: pd.DataFrame = None, **params) -> dict:
        """
        Parse arguments for the construction of an `AnnotatedGEM` object from `pandas.DataFrame` objects.

        :param count_df:
            The gene expression matrix as a `pandas.DataFrame`.
            This file is assumed to have genes as rows and samples as columns.

        :param label_df: The gene annotation data as a `pandas.DataFrame`.
            This file is assumed to have samples as rows and annotation observations
            as columns.

        :return:
            A parsed parameter dictionary.
        """
        data = xrarray_gem_from_pandas(count_df=count_df, label_df=label_df)
        return {"data": data, **params}

    @classmethod
    def from_pandas(cls, count_df: pd.DataFrame, label_df: pd.DataFrame = None, **params) -> TypeAGEM:
        """
        Construct a `GEM` object from `pandas.DataFrame` objects.

        :param count_df:
            The gene expression matrix as a `pandas.DataFrame`.
            This file is assumed to have genes as rows and samples as columns.

        :param label_df:
            The gene annotation data as a `pandas.DataFrame`.
            This file is assumed to have samples as rows and annotation observations
            as columns.

        :return:
            An instance of the `AnnotatedGEM` class.
        """
        data = xrarray_gem_from_pandas(count_df=count_df, label_df=label_df)
        params = {"data": data, **params}
        return cls(**params)

    @classmethod
    def _parse_files(cls, count_path, label_path=None, count_kwargs=None, label_kwargs=None, **params) -> dict:
        if count_kwargs is None:
            count_kwargs = dict(index_col=0)

        # Expand and resolve the given paths.
        count_path = str(pathlib.Path(count_path).expanduser().resolve())

        if label_path:
            label_path = str(pathlib.Path(label_path).expanduser().resolve())

        count_df = load_count_df(count_path=count_path, **count_kwargs)

        if label_path:
            if label_kwargs is None:
                label_kwargs = dict(index_col=0)
            label_df = load_label_df(label_path=label_path, **label_kwargs)
        else:
            label_df = None

        # Check to ensure the indexes match.
        if label_path and all(label_sample not in count_df.columns.values
                              for label_sample in label_df.index.values):
            raise ValueError(dedent("""Files cannot be automatically processed, please
            load your data in to `pandas.DataFrame` objects and provide them to the
            `AnnotatedGEM.from_pandas` constructor instead."""))

        data = xrarray_gem_from_pandas(count_df, label_df)

        return {"data": data, **params}

    @classmethod
    def from_files(cls,
                   count_path: str,
                   label_path: str = None,
                   count_kwargs: dict = None,
                   label_kwargs: dict = None,
                   **params) -> TypeAGEM:
        """
        Construct a `AnnotatedGEM` object from file paths and optional parsing arguments.

        :param count_path:
            The path to the gene expression matrix.

        :param label_path:
            The path to the gene annotation data.

        :param count_kwargs:
            Arguments to be passed to `pandas.read_csv` for the count matrix.

        :param label_kwargs:
            Arguments to be passed to `pandas.read_csv` for the annotations.

        :return:
            An instance of the `AnnotatedGEM` class.
        """
        params = cls._parse_files(count_path=count_path, label_path=label_path,
                                  count_kwargs=count_kwargs, label_kwargs=label_kwargs, **params)
        return cls(**params)

    def save(self, path):
        """
        Save as a netcdf (.nc) to the file at `path`.

        :param path:
            The filepath to save to. This should use the `.nc` extension.

        :return:
            The path to which the file was saved.
        """
        if path is None:
            path = f"{self.name}.nc"

        # Save some parameters that may be unique to this instance.
        params_to_save = {key: value for key, value in self.get_param_values()
                          if isinstance(value, str)}
        params_str = json.dumps(params_to_save)
        self.data.attrs.update({"__GSForge.AnnotatedGEM.params": params_str})
        self.data.to_netcdf(path, mode="w")
        return path


# Python 3.8 will let us move this code into the class body, and add
# # add register the @classmethod functions.
@functools.singledispatch
def _annotated_gem_dispatch(*args, **params):
    raise TypeError(f"Source of type: {type(args[0])} not supported.")


_annotated_gem_dispatch.register(xr.Dataset, AnnotatedGEM._parse_xarray_dataset)
_annotated_gem_dispatch.register(str, AnnotatedGEM._parse_netcdf_path)
_annotated_gem_dispatch.register(pathlib.PosixPath, AnnotatedGEM._parse_netcdf_path)
_annotated_gem_dispatch.register(pd.DataFrame, AnnotatedGEM._parse_pandas)
