"""
The GSForge interface model forms a basis for interacting with data
stored in the `AnnotatedGEM` and `GeneSetCollection` objects.
"""
from typing import Union
import param
from textwrap import dedent
import numpy as np
import xarray as xr
import functools

from ._AnnotatedGEM import AnnotatedGEM, TypeAGEM
# from ._GeneSet import GeneSet
from ._GeneSetCollection import GeneSetCollection, TypeGSC


# TODO: Add GeneSet data access.
class Interface(param.Parameterized):
    """
    The Interface provides common API access for interacting with the `AnnotatedGEM` and
    `GeneSetCollection` objects. It also accepts an `AnnotatedGEM` and a single `GeneSet`
    for subset selection.

    For updating default parameters within subclasses, use the following, although it may
    cause 'watching' parameters to fire.

    ```
    self.set_param(key=value)
    ```

    """

    gem = param.ClassSelector(class_=AnnotatedGEM, doc=dedent("""\
    An AnnotatedGEM object."""), default=None, precedence=-1.0)

    gene_set_collection = param.ClassSelector(class_=GeneSetCollection, doc=dedent("""\
    A GeneSetCollection object."""), default=None, precedence=-1.0)

    selected_gene_sets = param.ListSelector(default=[None], doc=dedent("""\
    A list of keys from the provided GeneSetCollection (stored in gene_set_collection)
    that are to be used for selecting sets of genes from the count matrix."""))

    selected_genes = param.Parameter(default=None, doc=dedent("""\
    A list of genes to use in indexing from the count matrix. This parameter takes
    priority over all other gene selecting methods. That means that selected
    lineaments (or combinations thereof) will have no effect."""), precedence=-1.0)

    gene_set_mode = param.ObjectSelector(
        default="union",
        objects=["complete", "union", "intersection"],
        doc=dedent("""\
        Controls how any selected gene sets are returned by the interface.
        + complete
            Returns the entire gene set of the AnnotatedGEM.
        + union
            Returns the union of the selected gene sets support.
        + intersection
            Returns the intersection of the selected gene sets support.
        """)
    )

    sample_subset = param.Parameter(default=None, precedence=-1.0, doc=dedent("""\
    A list of samples to use in a given operation. These can be supplied
    directly as a list of genes, or can be drawn from a given GeneSet."""))

    count_variable = param.String(default=None, precedence=-1.0, doc=dedent("""\
    The name of the count matrix used."""))

    annotation_variables = param.List(doc=dedent("""\
    The name of the active annotation variable(s). These are the annotation columns that will
    be control the subset returned by `y_annotation_data`."""), precedence=-1.0, default=None)

    count_mask = param.ObjectSelector(doc=dedent("""\
    The type of mask to use for the count matrix.
        + 'complete' returns the entire count matrix as numbers.
        + 'masked' returns the entire count matrix with zero or missing as NaN values.
        + 'dropped' returns the count matrix without genes that have zero or missing values.
    """), default='complete', objects=["complete", "masked", "dropped"], precedence=-1.0)

    annotation_mask = param.ObjectSelector(doc=dedent("""\
    The type of mask to use for the target array.
        + 'complete' returns the entire target array.
        + 'masked' returns the entire target array with zero or missing as NaN values.
        + 'dropped' returns the target array without samples that have zero or missing values.
    """), default='complete', objects=["complete", "dropped"], precedence=1.0)

    count_transform = param.Callable(default=None, precedence=-1.0, doc=dedent("""\
    A transform that will be run on the x_data that is supplied by this Interface. 
    The transform runs on the subset of the matrix that has been selected."""))

    def __init__(self, *args, **params):
        if args:
            params = _interface_dispatch(*args, **params)
        super().__init__(**params)

        # TODO: Add and update an 'available count matrix' parameter.

        if self.gene_set_collection is not None:
            avail_mappings = list(self.gene_set_collection.gene_sets.keys())
            self.param["selected_gene_sets"].objects = avail_mappings + [None]

    @staticmethod
    def _parse_annotated_gem(annotated_gem: TypeAGEM, *args, **params) -> dict:
        """
        Parse arguments for creation of a new `Interface` instance from an `AnnotatedGEM`.

        :param annotated_gem:
            A `GSForge.AnnotatedGEM` object.

        :param args:
            Not used.

        :param params:
            Parameters to initialize this `Interface` with.

        :return:
            A parsed parameter dictionary.
        """
        params = {"gem": annotated_gem,
                  "count_variable": annotated_gem.count_array_name,
                  **params}
        # if args:
        #     params = _interface_dispatch(*args, **params)
        return params

    @staticmethod
    def _parse_gene_set_collection(gene_set_collection: TypeGSC, *args, **params) -> dict:
        """
        Parse arguments for creation of a new `Interface` instance from an `GeneSetCollection`.

        :param gene_set_collection:
            A `GSForge.GeneSetCollection` object.

        :param args:
            Not used.

        :param params:
            Parameters to initialize this `Interface` with.

        :return:
            A parsed parameter dictionary.
        """
        if gene_set_collection.gem is not None:
            params = {"gem": gene_set_collection.gem,
                      "count_variable": gene_set_collection.gem.count_array_name,
                      **params}
        params = {"gene_set_collection": gene_set_collection, **params}

        # if args:
        #     params = _interface_dispatch(*args, **params)

        return params

    def get_gene_index(self, count_variable=None) -> np.array:
        """
        Get the currently selected gene index as a numpy array.

        :param count_variable:
            Optional, for selecting alternate count arrays. The variable to be retrieved.

        :return:
            A numpy array of the currently selected genes.
        """
        gene_set_combinations = {
            "union": lambda sel_gs: self.gene_set_collection.union(sel_gs),
            "intersection": lambda sel_gs: self.gene_set_collection.intersection(sel_gs),
            "complete": lambda sel_gs: self.gem.gene_index,
        }

        mask_modes = {
            "complete": lambda counts: counts.fillna(0),
            "masked": lambda counts: counts.where(counts > 0.0),
            "dropped": lambda counts: counts.where(counts > 0.0).dropna(dim=self.gem.gene_index_name),
        }

        # Ensure the correct count array is selected.
        count_variable = self.count_variable if count_variable is None else count_variable

        # If an explicit list of genes is provided, use those and move on to masking.
        if self.selected_genes is not None:
            support = self.selected_genes

        # Check if a collection was provided, if not: return the entire gene index for masking.
        elif self.selected_gene_sets == [None] or not self.gene_set_collection:
            support = self.gem.gene_index

        # Otherwise, use  some combination of GeneSet supports should be used.
        else:
            support = gene_set_combinations[self.gene_set_mode](self.selected_gene_sets)

        # Now the counts can be selected based on the gene support.
        counts = self.gem.data.sel({self.gem.gene_index_name: support})[count_variable]

        # Mask the counts.
        masked_counts = mask_modes[self.count_mask](counts)

        # Get the gene index from this mask.
        # A .copy() of the array values should be used to prevent strange index behavior.
        genes = masked_counts[self.gem.gene_index_name].values.copy()
        return genes

    @property
    def active_count_variable(self) -> str:
        """Returns the name of the currently active count matrix."""
        if self.count_variable is not None:
            count_variable = self.count_variable
        else:
            count_variable = self.gem.count_array_name
        return count_variable

    @property
    def gene_index_name(self) -> str:
        """Returns the name of the gene index."""
        return self.gem.gene_index_name

    @property
    def sample_index_name(self) -> str:
        """Returns the name of the sample index."""
        return self.gem.sample_index_name

    def get_sample_index(self) -> np.array:
        """
        Get the currently selected sample index as a numpy array.

        :return:
            A numpy array of the currently selected samples.
        """

        if self.sample_subset is not None:
            subset = self.gem.data.sel({self.gem.sample_index_name: self.sample_subset})
        else:
            subset = self.gem.data

        if self.annotation_mask == "complete":
            pass
        elif self.annotation_mask == "dropped":
            if self.annotation_variables is not None:
                subset = subset[self.annotation_variables].dropna(dim=self.gem.sample_index_name)

            subset = subset.dropna(dim=self.gem.sample_index_name)

        return subset[self.gem.sample_index_name].values.copy()

    def get_selection_indexes(self) -> dict:
        """Returns the currently selected indexes as a dictionary."""
        return {self.gem.gene_index_name: self.get_gene_index(),
                self.gem.sample_index_name: self.get_sample_index()}

    @property
    def selection(self) -> xr.Dataset:
        """Returns the currently selected data as an ``xarray.Dataset`` object.."""
        selected_variables = [self.active_count_variable]

        if self.annotation_variables is not None:
            selected_variables += self.annotation_variables

        selection = self.gem.data[selected_variables].sel({
            self.gem.gene_index_name: self.get_gene_index(),
            self.gem.sample_index_name: self.get_sample_index()})

        # Optional transform.
        if self.count_transform is not None:
            selection[self.count_variable] = self.count_transform(selection[self.count_variable])

        return selection

    @property
    def x_count_data(self) -> xr.DataArray:
        """
        Returns the currently selected 'x_data'. Usually this will be a subset of the active count array.

        :return:
            An `xarray.Dataset` selection of the currently active 'x_data'.
        """
        # TODO: Consider adding a copy option.

        mask_modes = {
            "complete": lambda counts: counts.fillna(0),
            "masked": lambda counts: counts.where(counts > 0.0),
            "dropped": lambda counts: counts.where(counts > 0.0).dropna(dim=self.gem.gene_index_name),
        }
        # Ensure the correct count array is selected.
        count_variable = self.count_variable if self.count_variable is not None else self.gem.count_array_name

        # Get the gene and sample indexes.
        gene_index = self.get_gene_index()
        sample_index = self.get_sample_index()

        selection_dict = {self.gem.gene_index_name: gene_index,
                          self.gem.sample_index_name: sample_index}
        selection_dict = {k: v for k, v in selection_dict.items() if v is not None}
        selection = self.gem.data.sel(selection_dict)

        # Mask the counts.
        data = mask_modes[self.count_mask](selection[count_variable]).copy(deep=True)

        # Optional transform.
        if self.count_transform is not None:
            data = self.count_transform(data)

        return data

    @property
    def y_annotation_data(self) -> Union[xr.Dataset, None]:
        """
        Returns the currently selected 'y_data', or None, based on the `selected_annotation_variables` parameter.

        :return:
            An `xarray.Dataset` or `xarray.DataArray` object of the currently selected y_data.
        """
        # TODO: Consider enforcing list input for standardizing outputs to datasets.
        # TODO: Consider adding a copy option.
        if self.annotation_variables is None:
            return None

        sample_index = self.get_sample_index()
        subset = self.gem.data.sel({self.gem.sample_index_name: sample_index})
        # if len(self.annotation_variables) == 1:
        #     return subset[self.annotation_variables[0]].copy(deep=True)
        return subset[self.annotation_variables].copy(deep=True)


# Python 3.8 will let us move this code into the class body, and add
# add register the @classmethods.
@functools.singledispatch
def _interface_dispatch(*args, **params):
    """Calls the appropriate classmethod."""
    print("dispatch called")
    raise TypeError(f"Source of type: {type(args[0])} not supported.")


_interface_dispatch.register(AnnotatedGEM, Interface._parse_annotated_gem)
_interface_dispatch.register(GeneSetCollection, Interface._parse_gene_set_collection)
# _interface_dispatch.register(GeneSet, Interface._parse_gene_set)
