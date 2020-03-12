from __future__ import annotations

from textwrap import dedent
from typing import Union

import numpy as np
import param
import xarray as xr

from ._AnnotatedGEM import AnnotatedGEM
from ._GeneSetCollection import GeneSetCollection
from GSForge._singledispatchmethod import singledispatchmethod


# TODO: Add .keys() and other dict like functionality.
class Interface(param.Parameterized):
    """
    The Interface provides common API access for interacting with the ``AnnotatedGEM`` and
    ``GeneSetCollection`` objects.

    Parameters
    ----------
    gem : AnnotatedGEM
        An instance of an ``AnnotatedGEM`` object that stores gene expression matrix and
        any associated sample or gene annotations.

    gene_set_collection : GeneSetCollection
        An instance of a ``GeneSetCollection`` object.

    selected_gene_sets : List[str]
        A list of keys from the provided ``GeneSetCollection`` (stored in `gene_set_collection`)
        that are to be used for selecting sets of genes from the count matrix.

    selected_genes : Union[List, np.ndarray]
        A list of genes to use in indexing from the count matrix. This parameter takes
        priority over all other gene selecting methods. That means that selected
        GeneSets (or combinations thereof) will have no effect.

    gene_set_mode : str
        Controls how any selected gene sets are returned by the interface.
        **complete**
            Returns the entire gene set of the ``AnnotatedGEM``.
        **union**
            Returns the union of the selected gene sets support.
        **intersection**
            Returns the intersection of the selected gene sets support.

    sample_subset : Union[List, np.ndarray]
        A list of samples to use in a given operation. These can be supplied
        directly as a list of genes, or can be drawn from a given GeneSet.

    count_variable : str
        The name of the count matrix used.

    annotation_variables : List[str]
        The name of the active annotation variable(s). These are the annotation columns
        that will be control the subset returned by `y_annotation_data`.

    count_mask : str
        The type of mask to use for the count matrix.
        **complete**
            Returns the entire count matrix as numbers.
        **masked**
            Returns the entire count matrix with zero or missing as NaN values.
        **dropped**
            Returns the count matrix without genes that have zero or missing values.

    annotation_mask : str
        The type of mask to use for the target array.
        **complete**
            Returns the entire target array.
        **dropped**
            Returns the target array without samples that have zero or missing values.

    count_transform : Callable
        A transform that will be run on the `x_data` that is supplied by this Interface.
        The transform runs on the subset of the matrix that has been selected.

        Some common transforms:

        .. doctest::
            :options: +SKIP

            lambda counts: np.log2(counts.where(counts > 0)))


    For updating default parameters within subclasses, use the following:

    .. doctest::
        :options: +SKIP

        >>> self.set_param(key=value)

    Although it may cause 'watching' parameters to fire.
    """

    gem = param.ClassSelector(class_=AnnotatedGEM, doc=dedent("""\
    An ``AnnotatedGEM`` object."""), default=None, precedence=-1.0)

    gene_set_collection = param.ClassSelector(class_=GeneSetCollection, doc=dedent("""\
    A ``GeneSetCollection`` object."""), default=None, precedence=-1.0)

    selected_gene_sets = param.ListSelector(default=None, doc=dedent("""\
    A list of keys from the provided GeneSetCollection (stored in gene_set_collection)
    that are to be used for selecting sets of genes from the count matrix."""))

    selected_genes = param.Parameter(default=None, doc=dedent("""\
    A list of genes to use in indexing from the count matrix. This parameter takes
    priority over all other gene selecting methods. That means that selected
    GeneSets (or combinations thereof) will have no effect."""), precedence=-1.0)

    gene_set_mode = param.ObjectSelector(
        default="union",
        objects=["complete", "union", "intersection"],
        doc=dedent("""\
    Controls how any selected gene sets are returned by the interface.
    **complete**
        Returns the entire gene set of the ``AnnotatedGEM``.
    **union**
        Returns the union of the selected gene sets support.
    **intersection**
        Returns the intersection of the selected gene sets support.
    """))

    sample_subset = param.Parameter(default=None, precedence=-1.0, doc=dedent("""\
    A list of samples to use in a given operation. These can be supplied
    directly as a list of genes, or can be drawn from a given GeneSet."""))

    count_variable = param.ObjectSelector(default=None, precedence=1.0, doc="""\
    The name of the count matrix used.""", objects=[None], check_on_set=False)

    # TODO: Change to an object selector?
    annotation_variables = param.List(doc=dedent("""\
    The name of the active annotation variable(s). These are the annotation columns that will
    be control the subset returned by ``y_annotation_data``."""), precedence=-1.0, default=None)

    count_mask = param.ObjectSelector(doc=dedent("""\
    The type of mask to use for the count matrix.
    **complete**
        Returns the entire count matrix as numbers.
    **masked** 
        Returns the entire count matrix with zero or missing as NaN values.
    **dropped** 
        Returns the count matrix without genes that have zero or missing values.
    """), default='complete', objects=["complete", "masked", "dropped"], precedence=1.0)

    annotation_mask = param.ObjectSelector(doc=dedent("""\
    The type of mask to use for the target array.
    **complete** 
        Returns the entire target array.
    **dropped** 
        Returns the target array without samples that have zero or missing values.
    """), default='complete', objects=["complete", "dropped"], precedence=1.0)

    count_transform = param.Callable(default=None, precedence=-1.0, doc=dedent("""\
    A transform that will be run on the `x_data` that is supplied by this Interface. 
    The transform runs on the subset of the matrix that has been selected."""))

    @singledispatchmethod
    def __interface_dispatch(*args, **params):
        raise TypeError(f"Source of type: {type(args[0])} not supported.")

    def __init__(self, *args, **params):
        if args:
            params = self.__interface_dispatch(*args, **params)

        # If the user passes a string, place it as a single item within a list.
        if isinstance(params.get("annotation_variables"), str):
            params["annotation_variables"] = [params.get("annotation_variables")]

        super().__init__(**params)

        # Populate the count variable selector with valid count arrays.
        self.param["count_variable"].objects = [None] + sorted(self.gem.count_array_names)

        if self.count_variable is None:
            self.set_param(**{"count_variable": self.gem.count_array_name})

        if self.gene_set_collection is not None:
            avail_mappings = list(self.gene_set_collection.gene_sets.keys())
            self.param["selected_gene_sets"].objects = avail_mappings + [None]

    @__interface_dispatch.register(AnnotatedGEM)
    @staticmethod
    def _parse_annotated_gem(annotated_gem: AnnotatedGEM, *_args, **params) -> dict:
        """
        Parse arguments for creation of a new `Interface` instance from an `AnnotatedGEM`.

        Parameters
        ----------
        annotated_gem : AnnotatedGEM
            A `GSForge.AnnotatedGEM` object.

        _args :
            Not used.

        params :
            Parameters to initialize this `Interface` with.

        Returns
        -------
        params : dict
            A parsed parameter dictionary.
        """
        params = {"gem": annotated_gem,
                  "count_variable": annotated_gem.count_array_name,
                  **params}
        return params

    @__interface_dispatch.register(GeneSetCollection)
    @staticmethod
    def _parse_gene_set_collection(gene_set_collection: GeneSetCollection, *_args, **params) -> dict:
        """
        Parse arguments for creation of a new `Interface` instance from an `GeneSetCollection`.

        Parameters
        ----------
        annotated_gem : AnnotatedGEM
            A `GSForge.AnnotatedGEM` object.

        _args :
            Not used.

        params :
            Parameters to initialize this `Interface` with.

        Returns
        -------
        params : dict
            A parsed parameter dictionary.
        """
        if gene_set_collection.gem is not None:
            params = {"gem": gene_set_collection.gem,
                      "count_variable": gene_set_collection.gem.count_array_name,
                      **params}
        params = {"gene_set_collection": gene_set_collection, **params}

        return params

    def get_gene_index(self, count_variable=None) -> np.array:
        """
        Get the currently selected gene index as a numpy array.

        Parameters
        ----------
        count_variable : str
            Optional, for selecting alternate count arrays. The variable to be retrieved.

        Returns
        -------
        np.ndarray
            An array of the currently selected genes.
        """
        gene_set_combinations = {
            "union": lambda sel_gs: self.gene_set_collection.union(sel_gs),
            "intersection": lambda sel_gs: self.gene_set_collection.intersection(sel_gs),
            "complete": lambda sel_gs: self.gem.gene_index,
        }

        mask_modes = {
            "complete": lambda counts_: counts_.fillna(0),
            "masked": lambda counts_: counts_.where(counts_ > 0.0),
            "dropped": lambda counts_: counts_.where(counts_ > 0.0).dropna(dim=self.gem.gene_index_name),
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

    def get_sample_index(self) -> np.ndarray:
        """
        Get the currently selected sample index as a numpy array.

        Returns
        -------
        np.ndarray
            An array of the currently selected samples.
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

        Returns
        -------
        xarray.Dataset
            The selection of the currently active count data.
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
    def y_annotation_data(self) -> Union[xr.Dataset, xr.DataArray, None]:
        """
        Returns the currently selected 'y_data', or None, based on the `selected_annotation_variables` parameter.

        :return:
            An `xarray.Dataset` or `xarray.DataArray` object of the currently selected y_data.
        """
        if self.annotation_variables is None:
            return None

        sample_index = self.get_sample_index()
        subset = self.gem.data.sel({self.gem.sample_index_name: sample_index})
        if len(self.annotation_variables) == 1:
            return subset[self.annotation_variables[0]].copy()
        return subset[self.annotation_variables].copy()
