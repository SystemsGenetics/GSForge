import json
import pathlib
import param
import os
import functools

import pandas as pd
import xarray as xr
import numpy as np

from textwrap import dedent

# Local imports.
from .. import utils


# TODO: Feature: Automatic GeneSet replicate combination?
#       + From name or some hash, and by:
#       + By membership
#       + By filter_function(variables)
# TODO: Add support inference?
#       Create the boolean support array from a variable in the given data parameter.
class GeneSet(param.Parameterized):
    """
    A `GeneSet` is at the very least a list of genes which are
    considered to be 'within' the GeneSet.

    A GeneSet can also be a measurement or ranking of a set of genes,
    and this could include all of the 'available' genes. In such cases
    a boolean array 'support' indicates membership in the GeneSet.
    """

    data = param.Parameter(allow_None=False, doc=dedent("""\
    Contains a gene-index `Xarray.Dataset` object, it should have
    only those genes that are considered 'within' the GeneSet
    in the index, or a boolean variable named 'support'."""))

    support_index_name = param.String(default="support", doc=dedent("""\
    This parameter controls which variable should be considered to be the
    (boolean) variable indicating membership in this GeneSet."""))
    
    gene_index_name = param.String(default="Gene", doc=dedent("""\
    This parameter controls which variable from the `Xarray.Dataset` should be 
    considered to be the 'gene index' coordinate.
    Consider using this if you require different coordinate names."""))

    def __init__(self, *args, **params):
        if args:
            params = _geneset_dispatch(*args, **params)
        super().__init__(**params)

    def __repr__(self):
        """Display a summary of this GeneSet."""
        support_size = self.gene_support().shape[0]
        percent_support = support_size / self.data['Gene'].shape[0]
        summary = [f"<GEMprospector.{type(self).__name__}>"]
        summary += [f"Name: {self.name}"]
        summary += [f"    Supported Genes:  {support_size}, {percent_support:.2%} of {self.data[self.gene_index_name].shape[0]}"]
        return "\n".join(summary)

    def gene_support(self) -> np.array:
        """Returns the list of genes 'supported in this GeneSet.

        The value that this return is (by default) controlled by the self.support_index_name
        parameter.

        :return: A numpy array of the genes 'supported' by this GeneSet.
        """
        # TODO: Consider implementing an override for the target variable.
        if self.support_index_name in self.data:
            supported_genes = self.data[self.gene_index_name].sel(
                {self.gene_index_name: (self.data[self.support_index_name] == True)}).values.copy()
            return self.data[self.gene_index_name].sel({self.gene_index_name: supported_genes}).values.copy()
        else:
            # TODO: Implement proper warnings.
            Warning("Warning, no boolean support array found. Returning the complete "
                    "set of genes in this GeneSet.")
            return self.data[self.gene_index_name].values.copy()

    def set_gene_support(self, genes):
        """Set this GeneSet support to the given genes."""
        self.data[self.support_index_name] = ((self.gene_index_name,), np.isin(self.data.Gene.values, genes))

    def set_boolean_support(self, variable):
        support = np.array(self.data[variable], dtype=bool)
        self.data[self.support_index_name] = ((self.gene_index_name,), support)

    # TODO: Combine with k_best_genes as a mode or option.
    def k_abs_best_genes(self, k=100, score_name=None):
        #     Fixes Issue #1240: NaNs can't be properly compared, so change them to the
        #     smallest value of scores's dtype. -inf seems to be unreliable.
        scores = self.data[score_name].values.copy()
        scores[np.isnan(scores)] = np.finfo(scores.dtype).min

        top_k_indexes = np.argsort(np.abs(scores))[-k:]
        top_k_genes = self.data.isel({self.gene_index_name: top_k_indexes})[self.gene_index_name].values.copy()
        return top_k_genes
        
    def k_best_genes(self, k=100, score_name=None) -> np.array:
        """Select the highest scoring genes from the 'score_name' variable.

        :param k: The number of genes to return.

        :param score_name: The variable name to rank genes by.

        :return: A numpy array of the top k genes based on their scores in `score_name`.
        """
        # Look for a variable name that ends in "_score" if none is given.
        if score_name is None:
            for var_name in list(self.data.variables.keys()):
                if "_scores" in var_name:
                    score_name = var_name
                    break  # Use the first match.

        scores = self.data[score_name].values.copy()
        scores = scores[np.isnan(scores)] = np.finfo(scores.dtype).min
        top_k_indexes = np.argsort(scores)[-k:]
        top_k_genes = self.data.isel({self.gene_index_name: top_k_indexes}).Gene.values.copy()
        return top_k_genes

    def q_best_genes(self, q=0.999, score_name=None) -> np.array:
        """Returns a numpy array of the q best genes based on the quantile `q`,
        and the target variable `score_name`.

        :param q: The quantile cutoff.

        :param score_name: The target variable to judge the genes by.

        :return: A numpy array of the top `q` quantile genes based on `score_name`.
        """
        if score_name is None:
            for var_name in list(self.data.variables.keys()):
                if "_scores" in var_name:
                    score_name = var_name
                    break  # Use the first match.

        quantile_score = self.data[score_name].quantile(q=q)
        quantile_selection = (self.data[score_name] >= quantile_score)
        top_q_genes = self.data.sel({"Gene": quantile_selection}).Gene.values.copy()
        return top_q_genes

    # TODO: Ensure 'attrs' get added to the dataset if they are found in params.
    @staticmethod
    def parse_xarray_dataset(data, **params):
        existing_params = data.attrs.get("__GSForge.GeneSet.params")
        if existing_params:
            existing_params = json.loads(existing_params)
            params = {**existing_params, **params}
        return {"data": data, **params}

    @classmethod
    def parse_netcdf_path(cls, netcdf_path, **params):
        params = cls.parse_xarray_dataset(xr.open_dataset(netcdf_path), **params)
        return {"data": xr.open_dataset(netcdf_path), **params}

    @classmethod
    def from_xarray_dataset(cls, data, **params):
        params = cls.parse_xarray_dataset(data, **params)
        return cls(**params)

    @classmethod
    def from_netcdf(cls, path, **params):
        """Construct a `GeneSet` object from a `netcdf` file path."""
        path = str(pathlib.Path(path).expanduser().resolve())
        params = cls.parse_xarray_dataset(xr.open_dataset(path), **params)
        return cls(**params)

    @staticmethod
    def parse_pandas(dataframe, genes=None, attrs=None, **params):
        """Parse a `pandas.DataFrame` for use in a GeneSet."""

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
        params = cls.parse_pandas(dataframe=dataframe, genes=genes, attrs=attrs, **params)
        return cls(**params)

    @staticmethod
    def parse_gene_array(selected_gene_array: np.ndarray, complete_gene_index=None, attrs=None, **params):
        coords = selected_gene_array if complete_gene_index is None else complete_gene_index

        data = xr.Dataset({"support": (["Gene"], np.isin(coords, selected_gene_array))}, coords={"Gene": coords})

        if attrs is not None:
            data = data.assign_attrs(attrs)

        return {"data": data, **params}

    @classmethod
    def from_gene_array(cls, selected_gene_array: np.ndarray, complete_gene_index=None, attrs=None, **params):
        params = cls.parse_gene_array(selected_gene_array=selected_gene_array,
                                      complete_gene_index=complete_gene_index,
                                      attrs=attrs, **params)
        return cls(**params)

    @staticmethod
    def parse_GeneSets(gene_sets, complete_gene_index, attrs=None, **params):
        union = np.array(list(set.union(*[set(lin.gene_support()) for lin in gene_sets])))
        data = xr.Dataset({"support": (["Gene"], np.isin(complete_gene_index, union))},
                          coords={"Gene": complete_gene_index})
        if attrs is not None:
            data = data.assign_attrs(attrs)
        return {"data": data, **params}

    @classmethod
    def from_GeneSets(cls, gene_sets, complete_gene_index, attrs=None, **params):
        """Create a new GeneSet by combining all the genes in the given GeneSets.
        No variables or attributes from the original GeneSets are maintained in this process."""
        params = cls.parse_GeneSets(gene_sets=gene_sets, complete_gene_index=complete_gene_index,
                                    attrs=attrs, **params)
        return cls(**params)

    @classmethod
    def from_bool_array(cls, bool_array, genes, attrs=None, **params):
        data = xr.Dataset({"support": (["Gene"], bool_array)}, coords={"Gene": genes})
        if attrs is not None:
            data = data.assign_attrs(attrs)
        params = {"data": data, **params}
        return cls(**params)

    # TODO: Consider overwrite protection?
    def save_as_netcdf(self, target_dir=None, name=None):
        """Save this GeneSet as a netcdf (.nc) file in the `target_dir` directory.

        The default filename will be: `{GeneSet.name}.nc`, if the GeneSet does not
        have a name, one must be provided via the `name` argument.

        :param target_dir: The directory to place the saved GeneSet into.

        :param name: The name to give the GeneSet upon saving.

        :returns output_path: The path to which the file was saved.
        """
        target_dir = os.getcwd() if target_dir is None else target_dir

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


_geneset_dispatch.register(str, GeneSet.parse_netcdf_path)
_geneset_dispatch.register(xr.Dataset, GeneSet.parse_xarray_dataset)
_geneset_dispatch.register(pd.DataFrame, GeneSet.parse_pandas)
_geneset_dispatch.register(np.ndarray, GeneSet.parse_gene_array)
_geneset_dispatch.register(list, GeneSet.parse_gene_array)
