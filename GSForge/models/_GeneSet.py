import copy
import functools
import json
import os
import pathlib
from textwrap import dedent
from typing import TypeVar, Type, Union

import numpy as np
import pandas as pd
import param
import xarray as xr

# Declare typing hints.
# Typing hints are not used by the Python runtime in any way.
# They are used by third party tools and linters.
# In python 3.8 we can use forward references.
TypeGS = TypeVar("TypeGS", bound="GeneSet")
TypePath = Union[str, 'PathLike[Any]']


class GeneSet(param.Parameterized):
    """
    A data class for a the result of a gene selection or analysis.

    A GeneSet can also be a measurement or ranking of a set of genes, and this could include all of
    the 'available' genes. In such cases a boolean array 'support' indicates membership in the GeneSet.
    """

    data = param.Parameter(allow_None=False, doc=dedent("""\
    Contains a gene-index `xarray.Dataset` object, it should have
    only those genes that are considered 'within' the GeneSet
    in the index, or a boolean variable named 'support'."""))

    support_index_name = param.String(default="support", doc=dedent("""\
    This parameter controls which variable should be considered to be the
    (boolean) variable indicating membership in this GeneSet."""))

    gene_index_name = param.String(default="Gene", doc=dedent("""\
    This parameter controls which variable from the `xarray.Dataset` should be 
    considered to be the 'gene index' coordinate.
    Consider using this if you require different coordinate names."""))

    def __init__(self, *args, **params):
        if args:
            params = _geneset_dispatch(*args, **params)
        super().__init__(**params)

    def __repr__(self) -> str:
        """
        Display a summary of this GeneSet.
        """
        support_size = self.gene_support().shape[0]
        percent_support = support_size / self.data['Gene'].shape[0]
        summary = [f"<GSForge.{type(self).__name__}>"]
        summary += [f"Name: {self.name}"]
        summary += [
            f"    Supported Genes:  {support_size}, {percent_support:.2%} of {self.data[self.gene_index_name].shape[0]}"]
        return "\n".join(summary)

    @property
    def gene_index(self) -> np.ndarray:
        """
        Returns the entire gene index of this GeneSet object as an ``xarray.DataArray``.

        The variable or coordinate that this returns is controlled by the ``gene_index_name``
        parameter.

        :return:
            The entire gene index of this GeneSet as an ``xarray.DataArray``.
        """
        return self.data[self.gene_index_name].copy()

    def gene_support(self) -> np.ndarray:
        """
        Returns the list of genes 'supported in this GeneSet.

        The value that this return is (by default) controlled by the self.support_index_name
        parameter.

        :return:
            A numpy array of the genes 'supported' by this GeneSet.
        """
        if self.support_index_name in self.data:
            supported_genes = self.data[self.gene_index_name].sel(
                {self.gene_index_name: (self.data[self.support_index_name] == True)}).values.copy()
            return self.data[self.gene_index_name].sel({self.gene_index_name: supported_genes}).values.copy()
        else:
            Warning("Warning, no boolean support array found. Returning the complete "
                    "set of genes in this GeneSet.")
            return self.data[self.gene_index_name].values.copy()

    @property
    def support_exists(self) -> bool:
        """
        Returns True if a support array exists, and that it has at least one member within, returns
        False otherwise.
        """
        if self.support_index_name not in self.data:
            return False
        if self.gene_support().shape[0] > 0:
            return True
        return False

    def set_support_by_genes(self, genes: np.ndarray) -> TypeGS:
        """
        Set this GeneSet support to the given genes. This function calculates the boolean support
        array for the gene index via `np.isin(gene_index, genes)`.

        :param genes:
            An array of genes which represent the "supported" subset within the entire gene index.

        :return:
            None. This function modifies the current object in place.
        """
        gs_copy = copy.deepcopy(self)
        gs_copy.data[self.support_index_name] = ((gs_copy.gene_index_name,),
                                                 np.isin(gs_copy.gene_index, np.asarray(genes)))
        return gs_copy

    def set_support_from_boolean_array(self, boolean_array: np.ndarray) -> TypeGS:
        """
        Set this GeneSet support based on the given boolean array, which must be the same length as the existing
        gene index.

        This function calculates the boolean support array for the gene index via `np.isin(gene_index, genes)`.

        :param boolean_array:
            A boolean `numpy.ndarray` or an `xarray.DataArray`.

        :return:
            A copy of this GeneSet with the support array modified.
        """
        gs_copy = copy.deepcopy(self)
        gs_copy.data[gs_copy.support_index_name] = ((gs_copy.gene_index_name,), np.asarray(boolean_array, dtype=bool))
        return gs_copy

    # TODO: Allow this to take a user-given callable function. Document with an example.
    def get_n_top_genes(self, score_variable: str = None, mode: str = "absolute_largest",
                        within_support: bool = True, **mode_kwargs) -> np.ndarray:
        """
        Get the top n genes as an array, selected via `mode` by using `score_variable` from the GeneSet.data variables.

        **Available modes**:

        `absolute_largest`
            Returns the n absolute largest values per the `score_variable` variable.
            Useful for selecting the absolute largest log-fold-change, regardless of the sign.

        :param score_variable:
            The name of the score variable to use. If this is not given some common score variable names are
            used, and the first found is used.

        :param mode:
            The mode by which to calculate how genes are ranked.

        :param within_support:
            Whether to limit the genes returned to those within the "support" of this GeneSet. Defaults to True.

        :param mode_kwargs:
            Remaining keyword arguments are passed to the 'mode' function.

        # TODO: Comment mode functions and their arguments clearly. Include how to write such a function.

        :return:
            An array of n, top-ranked, genes.
        """

        # TODO: Consider moving these functions to utilities?
        def _abs_largest(scores: xr.DataArray, **kwargs):
            nn = kwargs["n"]
            top_idx = np.argsort(np.abs(scores.values))[::-1][:nn]
            return scores.isel({self.gene_index_name: top_idx})[self.gene_index_name].values.copy()

        def _above_threshold(scores: xr.DataArray, **kwargs):
            threshold = kwargs["threshold"]
            return scores.isel({self.gene_index_name: scores.values > threshold})[self.gene_index_name].values.copy()

        def _below_threshold(scores: xr.DataArray, **kwargs):
            threshold = kwargs["threshold"]
            return scores.isel({self.gene_index_name: scores.values < threshold})[self.gene_index_name].values.copy()

        def _above_absolute_threshold(scores: xr.DataArray, **kwargs):
            threshold = kwargs["threshold"]
            return scores.isel({self.gene_index_name: np.abs(scores.values) > threshold})[
                self.gene_index_name].values.copy()

        modes = {
            "absolute_largest": _abs_largest,
            "above_threshold": _above_threshold,
            "below_threshold": _below_threshold,
            "above_absolute_threshold": _above_absolute_threshold,
        }

        def _scan_for_score_variable(target=score_variable):
            score_variables = [
                "score",
                "logFC",  # from EdgeR.
                "log2FoldChange",  # from DESeq2.
                "pvalue",   # from DESeq2.
            ]
            if target is not None:  # An explicit variable has been given.
                return target
            else:
                for name in score_variables:
                    for var_name in list(self.data.variables.keys()):
                        if name == var_name:
                            print(f"Automatically found a score variable, {name}. Ensure this the desired variable.")
                            return name

        score_name = _scan_for_score_variable()

        # TODO: Convert this to an internal helper function.
        # Look for a variable name that ends in "_score" if none is given.
        # score_name = score_variable
        # if score_variable is None:
        #     for var_name in list(self.data.variables.keys()):
        #         # TODO: Check more common 'score' variable names. Convert to a variable somewhere...
        #         if "score" in var_name:
        #             score_name = var_name
        #             break  # Use the first match.

        if score_name is None:
            raise ValueError("A score variable could not be automatically identified, please specify one.")

        data = self.data.sel({self.gene_index_name: self.gene_support()})[score_name] \
            if within_support \
            else self.data[score_name]

        # Assign the selected mode function.
        selected_mode = modes.get(mode)
        if selected_mode is None:
            raise ValueError(f"Given mode selection: {mode} is not a valid mode. Choose from:\n{modes.keys()}")

        return selected_mode(data, **mode_kwargs)

    @staticmethod
    def _parse_xarray_dataset(data: xr.Dataset, **params) -> dict:
        """
        Parses an argument set containing an `xarray.Dataset`, combines with given params and
        returns a parameter dictionary.

        :param data:
            An `xarray.Dataset` passed by the user.

        :param params:
            Keyword-value arguments.

        :return:
            A parsed parameter dictionary.
        """
        existing_params = data.attrs.get("__GSForge.GeneSet.params")
        if existing_params:
            existing_params = json.loads(existing_params)
            params = {**existing_params, **params}
        return {"data": data, **params}

    @classmethod
    def _parse_netcdf_path(cls: Type[TypeGS], netcdf_path: TypePath, **params) -> dict:
        """
        Parses an argument set containing an `xarray.Dataset`, combines with given params and
        returns a parameter dictionary.

        :param netcdf_path:
            The path to the .netcdf file to be used.

        :param params:
            Keyword-value arguments.

        :return:
            A parsed parameter dictionary.
        """
        params = cls._parse_xarray_dataset(xr.open_dataset(netcdf_path), **params)
        return {"data": xr.open_dataset(netcdf_path), **params}

    @classmethod
    def from_xarray_dataset(cls, data: xr.Dataset, **params) -> TypeGS:
        """
        Create a GeneSet from an `xarray.Dataset`.

        :param data:
            An `xarray.Dataset` object. See the `.data` parameter of this class.

        :param params:
            Other parameters to set.

        :return:
            A new instance of a GeneSet.
        """
        params = cls._parse_xarray_dataset(data, **params)
        return cls(**params)

    @classmethod
    def from_netcdf(cls, path, **params):
        """
        Create a `GeneSet` object from a `netcdf` file path.

        :param path:
            The path to the .netcdf file to be used.

        :param params:
            Other parameters to set.

        :return:
            A new instance of a GeneSet.
        """
        path = str(pathlib.Path(path).expanduser().resolve())
        params = cls._parse_xarray_dataset(xr.open_dataset(path), **params)
        return cls(**params)

    @staticmethod
    def _parse_pandas(dataframe, genes: np.ndarray = None, attrs: dict = None, **params) -> dict:
        """
        Parse a `pandas.DataFrame` for use in a GeneSet.

        :param dataframe:
            A `pandas.DataFrame` object. Assumed to be indexed by genes names.

        :param genes:
            If you have a separate (but ordered the same!) gene array that corresponds to your data,
            it can be passed here to be set as the index appropriately.

        :param attrs:
            A dictionary of attributes to be added to the `xarray.Dataset.attrs` attribute.

        :param params:
            Other parameters to set.

        :return:
            A parsed parameter dictionary.
        """

        if genes is not None:
            dataframe["Gene"] = genes
            dataframe = dataframe.set_index("Gene")

        if dataframe.index.name is None:
            dataframe.index.name = "Gene"
        else:
            params = {"gene_index_name": dataframe.index.name, **params}

        data = dataframe.to_xarray()

        if attrs is not None:
            data = data.assign_attrs(attrs)

        return {"data": data, **params}

    @classmethod
    def from_pandas(cls, dataframe: pd.DataFrame, genes=None, attrs=None, **params):
        """
        Create a GeneSet from a `pandas.DataFrame`.

        :param dataframe:
            A `pandas.DataFrame` object. Assumed to be indexed by genes names.

        :param genes:
            If you have a separate (but ordered the same!) gene array that corresponds to your data,
            it can be passed here to be set as the index appropriately.

        :param attrs:
            A dictionary of attributes to be added to the `xarray.Dataset.attrs` attribute.

        :param params:
            Other parameters to set.

        :return:
            A new GeneSet object.
        """
        params = cls._parse_pandas(dataframe=dataframe, genes=genes, attrs=attrs, **params)
        return cls(**params)

    def to_dataframe(self, only_supported: bool = True) -> pd.DataFrame:
        """
        Convert this GeneSet.data attribute to a `pandas.DataFrame`. This restricts the data
        returned to include only those genes that are returned by `GeneSet.gene_support()`.

        :param only_supported:
            Defaults to True, set to False if you want all GeneSet data to be in the DataFrame returned.

        :return:
            A `pandas.DataFrame` of this `GeneSet.data` attribute.
        """
        selection = self.gene_support() if only_supported else self.gene_index
        return self.data.sel(Gene=selection).to_dataframe()

    @staticmethod
    def _parse_gene_array(selected_gene_array: np.ndarray,
                          complete_gene_index=None, attrs=None, **params) -> dict:
        """
        Parses arguments for a new GeneSet from an array or list of 'selected' genes. Such genes are
        assumed to be within the optionally supplied `complete_gene_index`.

        :param selected_gene_array:
            The genes 'selected' to be within the support of this GeneSet.

        :param complete_gene_index:
            Optional. The complete gene index to which those selected genes belong.

        :param attrs:
            Attributes to add to the `xarray.Dataset.attrs` attribute.

        :param params:
            Parameters to initialize the `GeneSet` instance with.

        :return:
            A parsed parameter dictionary.
        """
        coords = selected_gene_array if complete_gene_index is None else complete_gene_index

        data = xr.Dataset({"support": (["Gene"], np.isin(coords, selected_gene_array))}, coords={"Gene": coords})

        if attrs is not None:
            data = data.assign_attrs(attrs)

        return {"data": data, **params}

    @classmethod
    def from_gene_array(cls, selected_gene_array: np.ndarray,
                        complete_gene_index=None, attrs=None, **params) -> TypeGS:
        """
        Parses arguments for a new GeneSet from an array or list of 'selected' genes. Such genes are
        assumed to be within the optionally supplied `complete_gene_index`.

        :param selected_gene_array:
            The genes 'selected' to be within the support of this GeneSet.

        :param complete_gene_index:
            Optional. The complete gene index to which those selected genes belong.

        :param attrs:
            Attributes to add to the `xarray.Dataset.attrs` attribute.

        :param params:
            Parameters to initialize the `GeneSet` instance with.

        :return:
            A new GeneSet object.
        """
        params = cls._parse_gene_array(selected_gene_array=selected_gene_array,
                                       complete_gene_index=complete_gene_index,
                                       attrs=attrs, **params)
        return cls(**params)

    @staticmethod
    def _parse_GeneSets(*gene_sets, mode: str = "union", attrs=None, **params) -> dict:
        """
        Combines the GeneSet objects given using `mode` to create a single new
        GeneSet object.

        Since the complete gene index is not necessarily known, it must minimally
        be the union of all genes included in the provided gene sets.

        :param gene_sets:
            One or more `GSForge.GeneSet` objects.

        :param mode:
            Mode by which the gene_sets should be combined. Options are "union"
            or "intersection".

        :param attrs:
            Optional attributes for the combined GeneSet. These attributes are
            added to the `GeneSet.data.attrs` attribute.

        :param params:
            Keyword parameters for the GeneSet object to be initialized with.

        :return:
            A new GeneSet object that contains genes from the provided `gene_sets`.
        """
        # Get the mode function from the dictionary of options.
        modes = {"union": lambda sets: functools.reduce(np.union1d, sets),
                 "intersection": lambda sets: functools.reduce(np.intersect1d, sets)}
        mode_fn = modes[mode]

        # Construct the minimal 'complete' gene index for the new GeneSet.
        gene_index = functools.reduce(np.union1d, [gs.gene_index.values for gs in gene_sets])

        # Construct the new support array based on the mode selected.
        gene_support = mode_fn([gs.gene_support() for gs in gene_sets])

        # Construct the new xarray.Dataset object.
        data = xr.Dataset(
            {"support": (["Gene"], np.isin(gene_index, gene_support))},
            coords={"Gene": gene_index})

        if attrs is not None:
            data = data.assign_attrs(attrs)

        return {"data": data, **params}

    @classmethod
    def from_GeneSets(cls, *gene_sets, mode: str = "union", attrs=None, **params) -> TypeGS:
        """
        Create a new GeneSet by combining all the genes in the given GeneSets.

        No variables or attributes from the original GeneSets are maintained in this process.
        """
        params = cls._parse_GeneSets(*gene_sets,
                                     mode=mode,
                                     attrs=attrs,
                                     **params)
        return cls(**params)

    @classmethod
    def from_bool_array(cls, bool_array, complete_gene_index, attrs=None, **params) -> TypeGS:
        """
        Create a GeneSet object from a boolean support array. This requires a matching
        gene index array.

        :param bool_array:
            A boolean array representing support within this GeneSet.

        :param complete_gene_index:
            The complete gene index.

        :param attrs:
            Optional attributes for the combined GeneSet. These attributes are
            added to the `GeneSet.data.attrs` attribute.

        :param params:
            Keyword parameters for the `GeneSet` object to be initialized with.

        :return:
            A new `GeneSet` object.
        """
        data = xr.Dataset({"support": (["Gene"], bool_array)}, coords={"Gene": complete_gene_index})
        if attrs is not None:
            data = data.assign_attrs(attrs)
        params = {"data": data, **params}
        return cls(**params)

    # TODO: Consider overwrite protection?
    # TODO: Make output of save functions consistent.
    def save_as_netcdf(self, target_dir=None, name=None) -> str:
        """
        Save this GeneSet as a netcdf (.nc) file in the `target_dir` directory.

        The default filename will be: `{GeneSet.name}.nc`, if the GeneSet does not
        have a name, one must be provided via the `name` argument.

        :param target_dir:
            The directory to place the saved GeneSet into.

        :param name:
            The name to give the GeneSet upon saving.

        :returns output_path:
            The path to which the file was saved.
        """
        target_dir = os.getcwd() if target_dir is None else target_dir

        # Create any needed intermediate directories.
        if not os.path.exists(target_dir):
            os.makedirs(target_dir)

        if self.name is None and name is None:
            raise ValueError("The GeneSet must be named, or you must pass a name to the save function.")

        active_name = name if name is not None else self.name
        # Remove any spaces in the active name.
        active_name = active_name.replace(" ", "_")
        output_path = os.path.join(target_dir, active_name + ".nc")

        # Save some parameters that may be unique to this instance.
        params_to_save = {key: value for key, value in self.get_param_values()
                          if isinstance(value, str)}
        params_str = json.dumps(params_to_save)
        self.data.attrs.update({"__GSForge.GeneSet.params": params_str})

        self.data.to_netcdf(output_path)

        return output_path


# Python 3.8 will let us move this code into the class body, and add
# add register the @classmethod functions.
@functools.singledispatch
def _geneset_dispatch(*args, **params):
    raise TypeError(f"Source of type: {type(args[0])} not supported.")


_geneset_dispatch.register(str, GeneSet._parse_netcdf_path)
_geneset_dispatch.register(pathlib.PosixPath, GeneSet._parse_netcdf_path)
_geneset_dispatch.register(xr.Dataset, GeneSet._parse_xarray_dataset)
_geneset_dispatch.register(pd.DataFrame, GeneSet._parse_pandas)
_geneset_dispatch.register(np.ndarray, GeneSet._parse_gene_array)
_geneset_dispatch.register(list, GeneSet._parse_gene_array)
