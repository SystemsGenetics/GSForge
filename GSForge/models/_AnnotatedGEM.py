from __future__ import annotations

import json
from pathlib import Path
from textwrap import dedent
from typing import List, Union, AnyStr, IO, Dict

import numpy as np
import pandas as pd
import param
import xarray as xr

from ._utils import (
    infer_xarray_variables,
    xrarray_gem_from_pandas,
    load_count_df,
    load_label_df,
)
from GSForge._singledispatchmethod import singledispatchmethod


class AnnotatedGEM(param.Parameterized):
    """
    A data class for a gene expression matrix and any associated sample or gene annotations.

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

    **Generate random a demo AnnotatedGEM:**

    .. doctest::
        :options: +SKIP

        >>> from sklearn.datasets import make_multilabel_classification
        >>> data, labels = make_multilabel_classification()
        >>> agem = AnnotatedGEM.from_pandas(pd.DataFrame(data), pd.DataFrame(labels), name='Generated GEM')
        >>> agem
        <GSForge.AnnotatedGEM>
        Name: Generated GEM
        Selected GEM Variable: ‘counts’
            Gene 100
            Sample 100

    **View the entire gene or sample index:**

    .. doctest::
        :options: +SKIP

        >>> agem.gene_index
        <xarray.DataArray ‘Gene’ (Gene: 100)>...

    .. doctest::
        :options: +SKIP

        >>> agem.sample_index
        <xarray.DataArray ‘Sample’ (Sample: 100)>...

    """

    data = param.ClassSelector(class_=xr.Dataset, allow_None=False, doc="""\
    An ``xarray.Dataset`` object that contains the Gene Expression Matrix, and any 
    needed annotations. This ``xarray.Dataset`` object is expected to have a count 
    array named 'counts', that has coordinates ('Gene', 'Sample').""")

    # TODO: Consider name change to 'selected_count_array'.
    count_array_name = param.String(default="counts", doc="""\
    This parameter controls which variable from the ``xarray.Dataset`` should be
    considered to be the 'count' variable.
    Consider using this if you require different index names, or wish to control 
    which count array among many should be used by default.""")

    sample_index_name = param.String(default="Sample", doc="""\
    This parameter controls which variable from the ``xarray.Dataset`` should be
    considered to be the 'sample' coordinate.
    Consider using this if you require different coordinate names.""")

    gene_index_name = param.String(default="Gene", doc="""\
    This parameter controls which variable from the ``xarray.Dataset`` should be 
    considered to be the 'gene index' coordinate.
    Consider using this if you require different coordinate names.""")

    @singledispatchmethod
    def __annotated_gem_dispatch(*args, **params):
        raise TypeError(f"Source of type: {type(args[0])} not supported.")

    def __init__(self, *args, **params) -> None:
        if args:
            params = self.__annotated_gem_dispatch(*args, **params)
        super().__init__(**params)

    def __repr__(self) -> str:
        """Display a summary of this AnnotatedGEM."""
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

        The variable or coordinate that this returns is controlled by the
        `gene_index_name` parameter.

        Returns
        -------
        xarray.DataArray
            The complete gene index of this AnnotatedGEM.
        """
        return self.data[self.gene_index_name].copy(deep=True)

    @property
    def sample_index(self) -> xr.DataArray:
        """
        Returns the entire sample index of this AnnotatedGEM object as an ``xarray.DataArray``.

        The actual variable or coordinate that this returns is controlled by the
        `sample_index_name` parameter.

        Returns
        -------
        xarray.DataArray
            The complete sample index of this AnnotatedGEM.
        """
        return self.data[self.sample_index_name].copy(deep=True)

    @property
    def count_array_names(self) -> List[str]:
        """
        Returns a list of all available count arrays contained within this AnnotatedGEM object.

        This is done simply by returning all data variables that have the same dimension set
        as the default count array.

        Returns
        -------
        List[str]
            A list of available count arrays in this AnnotatedGEM.
        """
        default_dims = set(self.data[self.count_array_name].dims)
        return [var for var in self.data.data_vars if set(self.data[var].dims) == default_dims]

    def infer_variables(self, quantile_size: int = 10, skip: bool = None) -> Dict[str, np.ndarray]:
        """
        Infer categories for the variables in the AnnotatedGEM's labels.

        Parameters
        ----------
        quantile_size : int
            The maximum number of unique elements before a variable is no
            longer considered as a `quantile-able` set of values.

        skip : bool
            The variables to be skipped.

        Returns
        -------
            A dictionary of the inferred value types.
        """
        if skip is None:
            skip = self.count_array_names + [self.sample_index_name, self.gene_index_name]

        sample_dim = {self.sample_index_name}
        gene_annots = [var for var in self.data.data_vars if set(self.data[var].dims) != sample_dim]
        if gene_annots:
            skip += gene_annots

        return infer_xarray_variables(xr_dataset=self.data, quantile_size=quantile_size, skip=skip)

    @__annotated_gem_dispatch.register(str)
    @__annotated_gem_dispatch.register(Path)
    @classmethod
    def _parse_netcdf_path(cls, netcdf_path: Union[str, Path, IO[AnyStr]], **params) -> dict:
        """
        Parse arguments for AnnotatedGEM creation via a path to an ``xarray.Dataset`` saved as a .netcdf file.

        Parameters
        ----------
        netcdf_path : Union[str, Path, IO[AnyStr]]
            A path to a netcdf file. If this file has different index names than default
            (Gene, Sample, counts), be sure to explicitly set those parameters
            (`gene_index_name`, `sample_index_name`, `count_array_name`).

        params : dict
            Other parameters to set.

        Returns
        -------
            A parsed parameter dictionary.
        """
        params = cls._parse_xarray_dataset(xr.open_dataset(netcdf_path), **params)
        return {"data": xr.open_dataset(netcdf_path), **params}

    @__annotated_gem_dispatch.register(xr.Dataset)
    @staticmethod
    def _parse_xarray_dataset(data: xr.Dataset, **params) -> dict:
        """
        Parse arguments for AnnotatedGEM creation via an ``xarray.Dataset``.

        Parameters
        ----------
        data : xr.Dataset
            An ``xarray.Dataset`` if this dataset has different index names than default
            (Gene, Sample, counts), be sure to explicitly set those parameters
            (`gene_index_name`, `sample_index_name`, `count_array_name`).

        params : dict
            Other parameters to set.

        Returns
        -------
            A parsed parameter dictionary.
        """
        existing_params = data.attrs.get("__GSForge.AnnotatedGEM.params")
        if existing_params:
            existing_params = json.loads(existing_params)
            params = {**existing_params, **params}
        return {"data": data, **params}

    @classmethod
    def from_netcdf(cls, netcdf_path: Union[str, Path, IO[AnyStr]], **params) -> AnnotatedGEM:
        """
        Construct an ``AnnotatedGEM`` object from a netcdf (.nc) file path.

        Parameters
        ----------
        netcdf_path : Union[str, Path, IO[AnyStr]]
            A path to a `netcdf` file. If this file has different index names than default
            (Gene, Sample, counts), be sure to explicitly set those parameters
            (`gene_index_name`, `sample_index_name`, `count_array_name`).

        Returns
        -------
        AnnotatedGEM : A new instance of the AnnotatedGEM class.
        """
        params = cls._parse_xarray_dataset(xr.open_dataset(netcdf_path), **params)
        return cls(**params)

    @__annotated_gem_dispatch.register(pd.DataFrame)
    @classmethod
    def _parse_pandas(cls, count_df: pd.DataFrame, label_df: pd.DataFrame = None, **params) -> dict:
        data = xrarray_gem_from_pandas(count_df=count_df, label_df=label_df)
        return {"data": data, **params}

    @classmethod
    def from_pandas(cls, count_df: pd.DataFrame, label_df: pd.DataFrame = None, **params) -> AnnotatedGEM:
        """
        Construct a `GEM` object from `pandas.DataFrame` objects.

        Parameters
        ----------
        count_df : pd.DataFrame
            The gene expression matrix as a `pandas.DataFrame`.
            This file is assumed to have genes as rows and samples as columns.

        label_df : pd.DataFrame
            The gene annotation data as a `pandas.DataFrame`.
            This file is assumed to have samples as rows and annotation observations
            as columns.

        Returns
        -------
        AnnotatedGEM : A new instance of the AnnotatedGEM class.
        """
        data = xrarray_gem_from_pandas(count_df=count_df, label_df=label_df)
        params = {"data": data, **params}
        return cls(**params)

    @classmethod
    def _parse_files(cls,
                     count_path: Union[str, Path, IO[AnyStr]],
                     label_path: Union[str, Path, IO[AnyStr]] = None,
                     count_kwargs: dict = None,
                     label_kwargs: dict = None, **params) -> dict:
        if count_kwargs is None:
            count_kwargs = dict(index_col=0)

        # Expand and resolve the given paths, ensures tilde '~' and environment variables are handled.
        count_path = str(Path(count_path).expanduser().resolve())

        count_df = load_count_df(count_path=count_path, **count_kwargs)

        if label_path:
            if label_kwargs is None:
                label_kwargs = dict(index_col=0)
            label_path = str(Path(label_path).expanduser().resolve())
            label_df = load_label_df(label_path=label_path, **label_kwargs)
        else:
            label_df = None

        # Check to ensure the indexes match.
        if label_path and all(label_sample not in count_df.columns.values
                              for label_sample in label_df.index.values):
            raise ValueError(dedent("""Files cannot be automatically processed, please
        load your data in to ``pandas.DataFrame`` objects and provide them to the
        ``AnnotatedGEM.from_pandas`` constructor instead."""))

        data = xrarray_gem_from_pandas(count_df, label_df)

        return {"data": data, **params}

    @classmethod
    def from_files(cls,
                   count_path: Union[str, Path, IO[AnyStr]],
                   label_path: Union[str, Path, IO[AnyStr]] = None,
                   count_kwargs: dict = None,
                   label_kwargs: dict = None,
                   **params) -> AnnotatedGEM:
        """
        Construct a ``AnnotatedGEM`` object from file paths and optional parsing arguments.

        Parameters
        ----------
        count_path : Union[str, Path, IO[AnyStr]]
            Path to the gene expression matrix.

        label_path : Union[str, Path, IO[AnyStr]]
            Path to the gene annotation data.

        count_kwargs : dict
            A dictionary of arguments to be passed to ``pandas.read_csv`` for the count matrix.

        label_kwargs : dict
            A dictionary of arguments to be passed to ``pandas.read_csv`` for the annotations.

        Returns
        -------
        AnnotatedGEM : A new instance of the AnnotatedGEM class.
        """
        params = cls._parse_files(count_path=count_path, label_path=label_path,
                                  count_kwargs=count_kwargs, label_kwargs=label_kwargs, **params)
        return cls(**params)

    def save(self, path: Union[str, Path, IO[AnyStr]]) -> str:
        """
        Save as a netcdf (.nc) to the file at ``path``.

        Parameters
        ----------
        path : Union[str, Path, IO[AnyStr]]
            The filepath to save to. This should use the ``.nc`` extension.

        Returns
        -------
        str : The path to which the file was saved.
        """
        if path is None:
            path = f"{self.name}.nc"

        # Save some parameters that may be unique to this instance.
        params_to_save = {key: value for key, value in self.get_param_values()
                          if isinstance(value, str)}
        params_str = json.dumps(params_to_save)
        self.data.attrs.update({"__GSForge.AnnotatedGEM.params": params_str})
        self.data.to_netcdf(path)
        return path
