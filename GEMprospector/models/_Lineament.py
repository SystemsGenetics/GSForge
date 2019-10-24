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


# TODO: Add support inference?
#       Create the boolean support array from a variable in the given data parameter.
class Lineament(param.Parameterized):
    """
    A `Lineament` is at the very least a list of genes which are
    considered to be 'within' the lineament.

    A Lineament can also be a measurement or ranking of a set of genes,
    and this could include all of the 'available' genes. In such cases
    a boolean array 'support' indicates membership in the lineament.
    """

    data = param.Parameter(allow_None=False, doc=dedent("""\
    Contains a gene-index `Xarray.Dataset` object, it should have
    only those genes that are considered 'within' the lineament
    in the index, or a boolean variable named 'support'."""))

    support_index_name = param.String(default="support", doc=dedent("""\
    This parameter controls which variable should be considered to be the
    (boolean) variable indicating membership in this Lineament."""))

    target = param.Parameter(default=None, doc=dedent("""\
    Optional. List of which target or targets this lineament describes. These values
    should reference a label category name, or a category therein.\n
    These should be provided as a list of tuples (label, category)."""))

    def __init__(self, *args, **params):
        if args:
            params = _lineament_dispatch(*args, **params)
        super().__init__(**params)

        # Infer parameters from the metadata attributes within the dataset.
        if self.target is None:
            try:
                target = self.data.attrs.get("y_variables")
                if target is not None:
                    self.set_param(target=target)
            except KeyError:
                Warning(f"Lineament {self.name} has no declared 'target', "
                        f"and one cannot be found within the dataset attrs.")

    def __repr__(self):
        """Display a summary of this Lineament."""
        support_size = self.gene_support().shape[0]
        percent_support = support_size / self.data['Gene'].shape[0]
        summary = [f"<GEMprospector.{type(self).__name__}>"]
        summary += [f"Name: {self.name}"]
        summary += [f"    Lineament Target: '{self.target}'"]
        summary += [f"    Supported Genes:  {support_size}, {percent_support:.2%} of {self.data['Gene'].shape[0]}"]
        return "\n".join(summary)

    def gene_support(self) -> np.array:
        """Returns the list of genes 'supported in this Lineament.

        The value that this return is (by default) controlled by the self.support_index_name
        parameter.

        :return: A numpy array of the genes 'supported' by this lineament.
        """
        # TODO: Consider implementing an override for the target variable.
        if self.support_index_name in self.data:
            supported_genes = self.data["Gene"].sel(
                {"Gene": (self.data[self.support_index_name] == True)}).values.copy()
            return self.data.Gene.sel(Gene=supported_genes).values.copy()
        else:
            # TODO: Implement proper warnings.
            Warning("Warning, no boolean support array found. Returning the complete "
                    "set of genes in this Lineament.")
            return self.data.Gene.values.copy()

    def set_gene_support(self, genes):
        """Set this lineaments support to the given genes."""
        self.data[self.support_index_name] = (("Gene",), np.isin(self.data.Gene.values, genes))

    def set_boolean_support(self, variable):
        support = np.array(self.data[variable], dtype=bool)
        self.data[self.support_index_name] = (("Gene",), support)

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

        scores = utils.clean_nans(self.data[score_name].values)
        top_k_indexes = np.argsort(scores)[-k:]
        top_k_genes = self.data.isel({"Gene": top_k_indexes}).Gene.values.copy()
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

    @staticmethod
    def parse_xarray_dataset(data, **params):
        existing_params = data.attrs.get("__GEMprospector.Lineament.params")
        if existing_params:
            existing_params = json.loads(existing_params)
            params = {**existing_params, **params}
        return {"data": data, **params}

    @classmethod
    def from_xarray_dataset(cls, data, **params):
        params = cls.parse_xarray_dataset(data, **params)
        return cls(**params)

    @classmethod
    def from_netcdf(cls, path, **params):
        """Construct a `Lineament` object from a `netcdf` file path."""
        path = str(pathlib.Path(path).expanduser().resolve())
        params = cls.parse_xarray_dataset(xr.open_dataset(path), **params)
        return cls(**params)

    @staticmethod
    def parse_pandas(dataframe, genes=None, attrs=None, **params):

        if genes is not None:
            dataframe["Gene"] = genes
            dataframe = dataframe.set_index("Gene")

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
    def parse_lineaments(lineaments, complete_gene_index, attrs=None, **params):
        union = np.array(list(set.union(*[set(lin.gene_support()) for lin in lineaments])))
        data = xr.Dataset({"support": (["Gene"], np.isin(complete_gene_index, union))},
                          coords={"Gene": complete_gene_index})
        if attrs is not None:
            data = data.assign_attrs(attrs)
        return {"data": data, **params}

    @classmethod
    def from_lineaments(cls, lineaments, complete_gene_index, attrs=None, **params):
        """Create a new lineament by combining all the genes in the given lineaments.
        No variables or attributes from the original lineaments are maintained in this process."""
        params = cls.parse_lineaments(lineaments=lineaments, complete_gene_index=complete_gene_index,
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
        """Save this Lineament as a netcdf (.nc) file in the `target_dir` directory.

        The default filename will be: `{lineament_name}.nc`, if the lineament does not
        have a name, one must be provided via the `name` argument.

        :param target_dir: The directory to place the saved lineament into.

        :param name: The name to give the Lineament upon saving.

        :returns output_path: The path to which the file was saved.
        """
        target_dir = os.getcwd() if target_dir is None else target_dir

        if self.name is None and name is None:
            raise ValueError("The lineament must be named, or you must pass a name to the save function.")

        active_name = name if name is not None else self.name
        # Remove any spaces in the active name.
        active_name = active_name.replace(" ", "_")
        output_path = os.path.join(target_dir, active_name + ".nc")

        # Save some parameters that may be unique to this instance.
        params_to_save = {key: value for key, value in self.get_param_values()
                          if isinstance(value, str)}
        params_str = json.dumps(params_to_save)
        self.data.attrs.update({"__GEMprospector.Lineament.params": params_str})

        self.data.to_netcdf(output_path)

        return output_path


# Python 3.8 will let us move this code into the class body, and add
# add register the @classmethod functions.
@functools.singledispatch
def _lineament_dispatch(*args, **params):
    raise TypeError(f"Source of type: {type(args[0])} not supported.")


_lineament_dispatch.register(str, Lineament.from_netcdf)
_lineament_dispatch.register(xr.Dataset, Lineament.from_xarray_dataset)
_lineament_dispatch.register(pd.DataFrame, Lineament.from_pandas)
_lineament_dispatch.register(np.ndarray, Lineament.from_gene_array)
_lineament_dispatch.register(list, Lineament.from_gene_array)
