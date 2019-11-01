"""
The GEMprospector interface model forms a basis for interacting with data
stored in the `AnnotatedGEM` and `LineamentCollection` objects.
"""

import param
from textwrap import dedent
import numpy as np
import xarray as xr
import functools

from ._AnnotatedGEM import AnnotatedGEM
from ._Lineament import Lineament
from ._LineamentCollection import LineamentCollection

from .. import utils


# TODO: Consider `gene_subset` parameter.
# TODO: Consider generalizing some LineamentCollection-specific functions
#       as found in the Lineament_Connectivity_Panel.
class Interface(param.Parameterized):
    """
    The Interface provides common API access for interacting with the `AnnotatedGEM` and
    `LineamentCollection` objects. It also accepts an `AnnotatedGEM` and a single `Lineament`
    for subset selection.
    """

    gem = param.ClassSelector(class_=AnnotatedGEM, doc=dedent("""\
    An AnnotatedGEM object."""), default=None, precedence=-1.0)

    lineament = param.ClassSelector(class_=Lineament, doc=dedent("""\
    A Lineament object which provides a gene subset."""), default=None, precedence=-1.0)

    lineament_collection = param.ClassSelector(class_=LineamentCollection, default=None,
                                               precedence=-1.0)

    lineament_keys = param.List(default=None, precedence=-1.0)

    # gene_subset = param.Parameter(default=None, precedence=-1.0, doc=dedent("""\
    # A list of genes to use in a given operation. These can be supplied
    # directly as a list of genes, or can be drawn from a given Lineament."""))

    sample_subset = param.Parameter(default=None, precedence=-1.0, doc=dedent("""\
    A list of samples to use in a given operation. These can be supplied
    directly as a list of genes, or can be drawn from a given Lineament."""))

    x_variable = param.String(default=None, precedence=-1.0, doc=dedent("""\
    The name of the count matrix used."""))

    y_variables = param.Parameter(doc=dedent("""\
    The name of the variable(s) with which the `x_variable` will be fitted.\n
    If these must be modified (binarized, normalized, etc.) this should be done
    prior to calling a GEMOperation. Passing a list of variables causes the
    interface to return an `Xarray.Dataset` object when `Interface.y_data` is
    called."""), precedence=-1.0, default=None)

    count_mask = param.ObjectSelector(doc=dedent("""\
    The type of mask to use for the count matrix.
        + 'complete' returns the entire count matrix as numbers.
        + 'masked' returns the entire count matrix with zero or missing as NaN values.
        + 'dropped' returns the count matrix without genes that have zero or missing values.
    """), default='complete', objects=["complete", "masked", "dropped"], precedence=-1.0)

    target_mask = param.ObjectSelector(doc=dedent("""\
    The type of mask to use for the target array.
        + 'complete' returns the entire target array.
        + 'masked' returns the entire target array with zero or missing as NaN values.
        + 'dropped' returns the target array without samples that have zero or missing values.
    """), default='complete', objects=["complete", "dropped"], precedence=-1.0)

    def __init__(self, *args, **params):
        if args:
            params = _interface_dispatch(*args, **params)

        super().__init__(**params)

    #         if self.x_variable is None:
    #             self.set_param(x_variable=self.gem.count_array_name)

    @staticmethod
    def parse_annotated_gem(annotated_gem, *args, **params):
        params = {"gem": annotated_gem,
                  "x_variable": annotated_gem.count_array_name,
                  **params}
        if args:
            params = _interface_dispatch(*args, **params)
        return params

    @staticmethod
    def parse_lineament_collection(lineament_collection, *args, **params):
        if lineament_collection.gem is not None:
            params = {"gem": lineament_collection.gem,
                      "x_variable": lineament_collection.gem.count_array_name,
                      **params}
        params = {"lineament_collection": lineament_collection, **params}

        if args:
            params = _interface_dispatch(*args, **params)

        return params

    @staticmethod
    def parse_lineament(lineament, *args, **params):
        params = {"lineament": lineament, **params}

        if args:
            params = _interface_dispatch(*args, **params)

        return params

    def get_gene_index(self, count_variable=None) -> np.array:
        """Get the currently selected gene index as a numpy array.

        :param count_variable: The variable to be retrieved.

        :return: A numpy array of the currently selected genes.
        """

        # Check if an explicit variable has been requested. Otherwise use the default x_variable.
        if count_variable is None:
            count_variable = self.x_variable

        if self.lineament_keys is not None:
            support = np.array(list(self.lineament_collection.union(self.lineament_keys)))
            counts = self.gem.data.sel({self.gem.gene_index_name: support})[count_variable]
        else:
            counts = self.gem.data[count_variable]

        if self.count_mask == "complete":
            pass
        elif self.count_mask == "masked":
            counts = counts.where(counts > 0.0)
        elif self.count_mask == "dropped":
            counts = counts.where(counts > 0.0).dropna(dim=self.gem.gene_index_name)

        return counts[self.gem.gene_index_name].values.copy()

    @property
    def gene_index_name(self):
        return self.gem.gene_index_name

    @property
    def sample_index_name(self):
        return self.gem.sample_index_name

    def get_sample_index(self) -> np.array:
        """Get the currently selected sample index as a numpy array.

        :param sample_variable: The variable to be retrieved, this should probably not be changed.

        :return: A numpy array of the currently selected samples.
        """

        if self.sample_subset is not None:
            subset = self.gem.data.sel({self.gem.sample_index_name: self.sample_subset})
        else:
            subset = self.gem.data

        if self.target_mask == "complete":
            pass
        elif self.target_mask == "dropped":
            if self.y_variables is not None:
                subset = subset[self.y_variables].dropna(dim=self.gem.sample_index_name)

            subset = subset.dropna(dim=self.gem.sample_index_name)

        return subset[self.gem.sample_index_name].values.copy()

    @property
    def selection(self) -> dict:
        """Returns the currently selected gene and sample indexes as a dictionary.

        This is usefull for selecting data from `xarray` objects.

        :return: A dictionary {index_name: active_genes}.
        """
        return {self.gem.gene_index_name: self.get_gene_index(),
                self.gem.sample_index_name: self.get_sample_index()}

    @property
    def x_data(self) -> xr.Dataset:
        """Returns the currently selected 'x_data'. Usually this will be a subset of the
        active count array.

        :return: An Xarray.Dataset selection of the currently active 'x_data'.
        """
        # TODO: Consider adding a copy option.
        gene_index = self.get_gene_index()
        sample_index = self.get_sample_index()
        selection_dict = {self.gem.gene_index_name: gene_index,
                          self.gem.sample_index_name: sample_index}
        selection_dict = {k: v for k, v in selection_dict.items() if v is not None}
        selection = self.gem.data.sel(selection_dict)

        if self.x_variable is not None:
            count_variable = self.x_variable
        else:
            count_variable = self.gem.count_array_name

        if self.count_mask == "masked":
            return selection[count_variable].where(selection[self.gem.count_array_name] > 0)
        else:
            return selection[count_variable].copy(deep=True)

    @property
    def y_data(self):
        """Returns the currently selected 'y_data', or None, based on the `y_variables` parameter.

        :return: An Xarray.Dataset or Xarray.DataArray object of the currently selected y_data.
        """
        # TODO: Consider adding a copy option.
        if self.y_variables is None:
            return None

        sample_index = self.get_sample_index()
        subset = self.gem.data.sel({self.gem.sample_index_name: sample_index})
        return subset[self.y_variables].copy(deep=True)


# Python 3.8 will let us move this code into the class body, and add
# add register the @classmethods.
@functools.singledispatch
def _interface_dispatch(*args, **params):
    """Calls the appropriate classmethod."""
    print("dispatch called")
    raise TypeError(f"Source of type: {type(args[0])} not supported.")


_interface_dispatch.register(AnnotatedGEM, Interface.parse_annotated_gem)
_interface_dispatch.register(LineamentCollection, Interface.parse_lineament_collection)
_interface_dispatch.register(Lineament, Interface.parse_lineament)
