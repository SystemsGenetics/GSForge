import pathlib
import json
import param
import functools

import pandas as pd
import xarray as xr

from textwrap import dedent

from .. import utils


# TODO: Track and infer the 'count' matrix variables, versus sample and gene annotations.
class AnnotatedGEM(param.Parameterized):
    """
    A wrapper class for a gene expression matrix and any associated phenotype
    data or labels.

    This model holds the count expression matrix, and any associated labels
    or annotations as an `xarray.Dataset` object under the `.data` attribute.
    By default this dataset will be expected to have its indexes named "Gene"
    and "Sample", although there are parameters to override those arrays and
    index names used.
    """

    data = param.ClassSelector(class_=xr.Dataset, doc=dedent("""\
    An `Xarray.Dataset` object that contains the Gene Expression Matrix, and any 
    needed annotations. This `Xarray.Dataset` object is expected to have a count 
    array named 'counts', that has coordinates ('Gene', 'Sample')."""))

    # TODO: This may not be needed, a select_count_array should be in the Interface class.
    count_array_name = param.String(default="counts", doc=dedent("""\
    This parameter controls which variable from the `Xarray.Dataset` should be
    considered to be the 'count' variable.
    Consider using this if you require different index names, or wish to control 
    which count array among many should be used by default."""))

    sample_index_name = param.String(default="Sample", doc=dedent("""\
    This parameter controls which variable from the `Xarray.Dataset` should be
    considered to be the 'sample' coordinate.
    Consider using this if you require different coordinate names."""))

    gene_index_name = param.String(default="Gene", doc=dedent("""\
    This parameter controls which variable from the `Xarray.Dataset` should be 
    considered to be the 'gene index' coordinate.
    Consider using this if you require different coordinate names."""))

    def __init__(self, *args, **params):
        if args:
            params = _annotated_gem_dispatch(*args, **params)
        super().__init__(**params)

    def __repr__(self) -> str:
        """Display a summary of this AnnotatedGEM."""
        summary = [f"<GEMprospector.{type(self).__name__}>"]
        summary += [f"Name: {self.name}"]
        summary += [f"Selected GEM Variable: '{self.count_array_name}'"]
        summary += [f"    {self.gene_index_name}   {self.data[self.gene_index_name].shape[0]}"]
        summary += [f"    {self.sample_index_name} {self.data[self.sample_index_name].shape[0]}"]
        return "\n".join(summary)

    @property
    def gene_index(self) -> xr.DataArray:
        """Returns the entire gene index of this AnnotatedGEM object as an `Xarray.DataArray`.

        The actual variable or coordinate that this returns is controlled by the
        `gene_index_name` parameter.
        """
        return self.data[self.gene_index_name].copy(deep=True)

    @property
    def sample_index(self) -> xr.DataArray:
        """Returns the entire sample index of this AnnotatedGEM object as an `Xarray.DataArray`.

        The actual variable or coordinate that this returns is controlled by the
        `sample_index_name` parameter.
        """
        return self.data[self.sample_index_name].copy(deep=True)

    @property
    def count_array_names(self) -> list:
        """Returns a list of all available count arrays contained within this AnnotatedGEM object.

        This is done simply by returning all data variables that have the same dimension set
        as the default count array.
        """
        default_dims = set(self.data[self.count_array_name].dims)
        return [var for var in self.data.data_vars if set(self.data[var].dims) == default_dims]

    # TODO: Remove from this class, move to interface.
    def zero_mask_counts(self, count_name: str = None) -> xr.DataArray:
        """Returns the count matrix as an `Xarray.DataArray` with any zero values masked as NaN.

        :param count_name: [Optional] Name of the variable to pull from the data. For use if there
            are multiple versions of the count data present (i.e. normalizations).

        :return: An `Xarray.DataArray` of counts with any zero values masked as NaN.
        """
        if count_name is None:
            count_name = self.count_array_name
        return self.data[count_name].where(self.data[count_name] > 0.0)

    # TODO: Remove from this class, move to interface.
    def zero_dropped_counts(self, count_name: str = None) -> xr.DataArray:
        """Returns the count matrix as an `Xarray.DataArray` with any zero values dropped
        along the 'gene' axis.

        :param count_name: [Optional] Name of the variable to pull from the data. For use if there
            are multiple versions of the count data present (i.e. normalizations).

        :return: An `Xarray.DataArray` of counts with any zero values dropped along the gene axis.
        """
        if count_name is None:
            count_name = self.count_array_name
        zero_masked_counts = self.zero_mask_counts(count_name=count_name)
        return zero_masked_counts.dropna(dim=self.gene_index_name)

    # TODO: Remove from this class, move to interface.
    def zero_dropped_gene_index(self, count_name: str = None) -> xr.DataArray:
        """Returns the count matrix gene index as an `Xarray.DataArray` with any genes containing
        zeros having been dropped.

        :param count_name: [Optional] Name of the variable to pull from the data. For use if there
            are multiple versions of the count data present (i.e. normalizations).

        :return: An `Xarray.DataArray` of the counts gene index with any zero values dropped along
            the gene axis.
        """
        return self.zero_dropped_counts(count_name=count_name)[self.gene_index_name]

    # TODO: Remove from this class, move to interface.
    def label_dataframe(self, labels: list = None) -> pd.DataFrame:
        """Returns the supplied list of labels as a `pandas.DataFrame` object.

        :param labels: [Optional] the list of variables to be included as columns
            in the output.

        :return: A `pandas.DataFrame` object.
        """
        if labels is None:
            labels = self.data.attrs["all_labels"]
        return self.data[labels].to_dataframe()

    # TODO: Remove from this class, move to interface?
    def infer_variables(self, quantile_size=10, skip=None) -> dict:
        """Infer categories for the variables in the AnnotatedGEM's labels.

        :param quantile_size: The maximum number of unique elements before a variable is no
            longer considered as a `quantile-able` set of values.

        :param skip: The variables to be skipped.

        :return: A dictionary of the inferred value types.
        """
        if skip is None:
            # TODO: Ensure that other potential count arrays end up here.
            skip = [self.count_array_name,
                    self.sample_index_name,
                    self.gene_index_name]

        return utils.infer_xarray_variables(xr_dataset=self.data,
                                            quantile_size=quantile_size,
                                            skip=skip)

    # def plot_label_bars(self, max_=8, labels=None):
    #     return plots.plot_label_bars(self.label_dataframe(labels), max_)

    # TODO: Consider adding an option for transposing.
    @staticmethod
    def xrarray_gem_from_pandas(count_df: pd.DataFrame,
                                label_df: pd.DataFrame = None) -> xr.Dataset:
        """Stitch together a gene expression and annotation DataFrames into
        a single `xarray.Dataset` object.

        :param count_df: The gene expression matrix as a `pandas.DataFrame`.
        :param label_df: The gene annotation data as a `pandas.DataFrame`.

        :return: An `xarray.Dataset` containing the gene expression matrix and
            the gene annotation data.
        """
        count_array = xr.Dataset(
            {"counts": (("Gene", "Sample"), count_df.values)},
            coords={
                "Sample": count_df.columns.values,
                "Gene": count_df.index.values
            }
        )

        if label_df is None:
            return count_array.transpose().to_dataset()

        else:
            label_ds = label_df.to_xarray()
            # TODO: Update `skip` parameter.
            attrs = utils.infer_xarray_variables(label_ds, skip=["counts", "Gene", "Sample"])
            label_ds = label_ds.assign_attrs(attrs)
            return label_ds.merge(count_array).transpose()

    @staticmethod
    def parse_xarray_dataset(data, **params):
        existing_params = data.attrs.get("__GEMprospector.AnnotatedGEM.params")
        if existing_params:
            existing_params = json.loads(existing_params)
            params = {**existing_params, **params}
        return {"data": data, **params}

    @classmethod
    def parse_netcdf_path(cls, netcdf_path, **params):
        params = cls.parse_xarray_dataset(xr.open_dataset(netcdf_path), **params)
        return {"data": xr.open_dataset(netcdf_path), **params}

    @classmethod
    def from_netcdf(cls, netcdf_path, **params):
        """Construct a `GEM` object from a `netcdf` file path."""
        params = cls.parse_xarray_dataset(xr.open_dataset(netcdf_path), **params)
        return cls(**params)

    @classmethod
    def parse_pandas(cls, count_df, label_df, **params):
        data = cls.xrarray_gem_from_pandas(count_df=count_df, label_df=label_df)
        return {"data": data, **params}

    @classmethod
    def from_pandas(cls,
                    count_df: pd.DataFrame,
                    label_df: pd.DataFrame = None,
                    **params):
        """Construct a `GEM` object from `pandas.DataFrame` objects.

        :param count_df: The gene expression matrix as a `pandas.DataFrame`.

        :param label_df: The gene annotation data as a `pandas.DataFrame`.

        :return: An instance of the `GEM` class.
        """
        data = cls.xrarray_gem_from_pandas(count_df=count_df, label_df=label_df)
        params = {"data": data, **params}
        instance = super().__new__(cls)
        instance.set_param(**params)
        return instance

    @classmethod
    def parse_files(cls, count_path, label_path=None, count_kwargs=None, label_kwargs=None, **params):
        if count_kwargs is None:
            count_kwargs = dict(index_col=0)

        # Expand and resolve the given paths.
        count_path = str(pathlib.Path(count_path).expanduser().resolve())

        if label_path:
            label_path = str(pathlib.Path(label_path).expanduser().resolve())

        count_df = utils.load_count_df(count_path=count_path, **count_kwargs)

        if label_path:
            if label_kwargs is None:
                label_kwargs = dict(index_col=0)
            label_df = utils.load_label_df(label_path=label_path, **label_kwargs)
        else:
            label_df = None

        # Check to ensure the indexes match.
        if label_path and all(label_sample not in count_df.columns.values
                              for label_sample in label_df.index.values):
            raise ValueError(dedent("""Files cannot be automatically processed, please
            load your data in to `pandas.DataFrame` objects and provide them to the
            `AnnotatedGEM.from_pandas` constructor instead."""))

        data = cls.xrarray_gem_from_pandas(count_df, label_df)

        return {"data": data, **params}

    @classmethod
    def from_files(cls,
                   count_path: str,
                   label_path: str = None,
                   count_kwargs: dict = None,
                   label_kwargs: dict = None,
                   **params):
        """Construct a `GEM` object from file paths and optional parsing arguments.

        :param count_path: The path to the gene expression matrix.

        :param label_path: The path to the gene annotation data.

        :param count_kwargs: Arguments to be passed to `pandas.read_csv` for the count matrix.

        :param label_kwargs: Arguments to be passed to `pandas.read_csv` for the annotations.

        :return: An instance of the `GEM` class.
        """
        params = cls.parse_files(count_path=count_path, label_path=label_path,
                                 count_kwargs=count_kwargs, label_kwargs=label_kwargs, **params)
        return cls(**params)

    def save(self, path):
        """Save as a netcdf (.nc) to the file at `path`.

        :param path: The filepath to save to. This should use the `.nc` extension.

        :return: The path to which the file was saved.
        """
        if path is None:
            path = f"{self.name}.nc"

        # Save some parameters that may be unique to this instance.
        params_to_save = {key: value for key, value in self.get_param_values()
                          if isinstance(value, str)}
        params_str = json.dumps(params_to_save)
        self.data.attrs.update({"__GEMprospector.AnnotatedGEM.params": params_str})
        self.data.to_netcdf(path, mode="w")
        return path


# Python 3.8 will let us move this code into the class body, and add
# # add register the @classmethod functions.
@functools.singledispatch
def _annotated_gem_dispatch(*args, **params):
    raise TypeError(f"Source of type: {type(args[0])} not supported.")


_annotated_gem_dispatch.register(xr.Dataset, AnnotatedGEM.parse_xarray_dataset)
_annotated_gem_dispatch.register(str, AnnotatedGEM.parse_netcdf_path)
_annotated_gem_dispatch.register(pd.DataFrame, AnnotatedGEM.parse_pandas)
