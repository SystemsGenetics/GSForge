from __future__ import annotations

import logging
from textwrap import dedent
from typing import Union

import numpy as np
import param
import xarray as xr

from GSForge._singledispatchmethod import singledispatchmethod
from ._AnnotatedGEM import AnnotatedGEM
from ._GeneSetCollection import GeneSetCollection
from ..utils import transient_log_handler


logger = logging.getLogger("GSForge")


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

        >>> self.param.set_param(key=value)

    """

    gem = param.ClassSelector(class_=AnnotatedGEM, doc="""\
    An ``AnnotatedGEM`` object.""", default=None, precedence=-1.0)

    gene_set_collection = param.ClassSelector(class_=GeneSetCollection, doc=dedent("""\
    A ``GeneSetCollection`` object."""), default=None, precedence=-1.0)

    selected_gene_sets = param.ListSelector(default=[None], doc=dedent("""\
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

    count_variable = param.ObjectSelector(default=None, precedence=1.0,
                                          doc="The name of the count matrix used.",
                                          objects=[None], check_on_set=False)

    # TODO: Change to an object selector?
    annotation_variables = param.List(doc=dedent("""\
    The name of the active annotation variable(s). These are the annotation columns that will
    be control the subset returned by ``y_annotation_data``."""), precedence=-1.0, default=[None])

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

    #     verbose = param.Boolean(default=False)

    @singledispatchmethod
    def _interface_dispatch(*args, **params):
        raise TypeError(f"Source of type: {type(args[0])} not supported.")

    def __init__(self, *args, **params):

        # If the user passes a string, place it as a single item within a list.
        if isinstance(params.get("annotation_variables"), str):
            params["annotation_variables"] = [params.get("annotation_variables")]

        if args:
            params = self._interface_dispatch(*args, **params)

        super().__init__(**params)

        # Populate the count variable selector with valid count arrays.
        self.param["count_variable"].objects = [None] + sorted(self.gem.count_array_names)

        if self.count_variable is None:
            self.param.set_param(**{"count_variable": self.gem.count_array_name})

        self.param["count_variable"].objects = self.gem.count_array_names# + [None]

        if self.gene_set_collection is not None:
            avail_mappings = list(self.gene_set_collection.gene_sets.keys())
            self.param["selected_gene_sets"].objects = avail_mappings + [None]

    @_interface_dispatch.register(AnnotatedGEM)
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

    @_interface_dispatch.register(GeneSetCollection)
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
        logger.info(f'Determining sample index.')

        if self.sample_subset is not None:
            subset = self.gem.data.sel({self.gem.sample_index_name: self.sample_subset})
            logger.info('No sample subset selected.')
        else:
            subset = self.gem.data

        if self.annotation_mask == "complete":
            pass
        elif self.annotation_mask == "dropped":
            if self.annotation_variables != [None]:
                subset = subset[self.annotation_variables].dropna(dim=self.gem.sample_index_name)

            subset = subset.dropna(dim=self.gem.sample_index_name)
            logger.info(
                f'Dropped samples with missing labels, {subset[self.gem.sample_index_name].shape[0]} samples remain.')

        selected_samples = subset[self.gem.sample_index_name].values.copy()
        return selected_samples

    @property
    def selection_indexes(self) -> dict:
        """Returns the currently selected indexes as a dictionary."""
        return {self.gem.gene_index_name: self.get_gene_index(),
                self.gem.sample_index_name: self.get_sample_index()}

    # TODO: Consider adding a copy option.
    @property
    def x_count_data(self) -> xr.DataArray:
        """
        Returns the currently selected 'x_data'. Usually this will be a subset of the active count array.
        
        Note: In constructing the a gene index, the count data is constructed first in order to infer
        coordinate selection based on masking.

        Returns
        -------
        xarray.Dataset
            The selection of the currently active count data.
        """

        gene_set_combinations = {
            "union": lambda sel_gs: self.gene_set_collection.union(sel_gs),
            "intersection": lambda sel_gs: self.gene_set_collection.intersection(sel_gs),
            "joint_difference": lambda sel_gs: self.gene_set_collection.joint_difference(sel_gs),
            "complete": lambda sel_gs: self.gem.gene_index,
        }

        mask_modes = {
            "complete": lambda counts: counts.fillna(0),
            "masked": lambda counts: counts.where(counts > 0.0),
            "dropped": lambda counts: counts.where(counts > 0.0).dropna(dim=self.gem.gene_index_name),
        }

        # Ensure the correct count array is selected.
        count_variable = self.count_variable if self.count_variable is not None else self.gem.count_array_name

        logger.info(f'Preparing count data from the variable {count_variable}.')

        if self.selected_genes is not None:
            support = self.selected_genes
            logger.info(f'Gene selection of: {support.shape[0]} genes provided.')

        # Check if a collection was provided, if not: return the entire gene index for masking.
        elif self.gene_set_collection is None:
            support = self.gem.gene_index
            logger.info(f'No collection or gene selection provided, using the entire gene index.')

        # If a collection has been provided, but no genesets have been selected,  use the entire index.
        elif self.selected_gene_sets == [None] or self.gene_set_collection is None:
            support = self.gem.gene_index
            logger.info(f'No genesets selected, using the entire gene index.')

        # Otherwise, use  some combination of GeneSet supports should be used.
        else:
            support = gene_set_combinations[self.gene_set_mode](self.selected_gene_sets)
            logger.info(
                f'Selected {len(self.selected_gene_sets)} GeneSets, using mode: {self.gene_set_mode} for a support of size: {support.shape[0]}.')

        # Now the counts can be selected based on the gene support, and the selected samples.
        sample_support = self.get_sample_index()

        counts = self.gem.data.sel({self.gem.gene_index_name: support,
                                    self.gem.sample_index_name: sample_support})[count_variable]
        logger.info(f'Selected count array of shape: {counts.shape}')

        logger.info(f'Preparing count data using mask mode {self.count_mask}.')
        counts = mask_modes[self.count_mask](counts)

        # Optional transform.
        if self.count_transform is not None:
            logger.info(f'Applying given transform to counts...')
            counts = self.count_transform(counts.copy(deep=True))

        return counts

    def get_gene_index(self) -> np.array:
        """
        Get the currently selected gene index as a numpy array.

        Returns
        -------
        np.ndarray
            An array of the currently selected genes.
        """
        logger.info(f'Preparing the gene index, this requires determining x_data.')
        return self.x_count_data[self.gene_index_name].values.copy()

    @property
    def y_annotation_data(self) -> Union[xr.Dataset, xr.DataArray, None]:
        """
        Returns the currently selected 'y_data', or None, based on the `selected_annotation_variables` parameter.

        Returns
        -------
            An ``xarray.Dataset`` of the currently selected y_data.
        """
        if (self.annotation_variables is None) or (self.annotation_variables == [None]):
            logger.info('No annotations selected.')
            return None

        logger.info(f'The following annotations where selected: {self.annotation_variables}.')
        sample_index = self.get_sample_index()
        subset = self.gem.data.sel({self.gem.sample_index_name: sample_index})
        # If only one label has been selected, return this as an xarray.DataArray.
        if len(self.annotation_variables) == 1:
            return subset[self.annotation_variables[0]].copy()
        return subset[self.annotation_variables].copy()

    # TODO: Should this be a private function?
    #       Users should call gsf.get_gem_data..., or should agem.get_gem_data...
    # For now make private.
    # Internal use should make use of x_ and y_ data.
    def get_gem_data(self, single_object=False, output_type='xarray', **params):
        """
        Returns count [and annotation] data based on the current parameters.

        Allows selection of the type and grouping of the output.

        Parameters
        ----------
        single_object
        output_type

        Returns
        -------

        """

        if params:
            logger.info(f'params {params}')
            self.param.set_param(**params)

        def xarray_single(self_):
            logger.info('Returning data as a single ``xarray.DataArray``.')
            if self_.y_annotation_data is not None:
                data = xr.merge([self_.x_count_data, self_.y_annotation_data])
            else:
                data = self_.x_count_data
            return data

        def xarray_tuple(self_):
            logger.info('Returning counts as an xarray.DataArray and annotations as an xarray.Dataset.')
            return self_.x_count_data, self_.y_annotation_data

        def pandas_single(self_):
            logger.info('Returning counts and annotations as a single pandas.DataFrame.')
            if self_.y_annotation_data is not None:
                data = xr.merge([self_.x_count_data, self_.y_annotation_data])
            else:
                data = self_.x_count_data
            return data.to_dataframe()

        def pandas_tuple(self_):
            logger.info('Returning a tuple of counts and annotations each as a pandas.DataFrame.')
            return self_.x_count_data.to_dataframe(), self_.y_annotation_data.to_dataframe()

        def numpy_single(self_):
            return np.dstack(self_.x_count_data.values, self_.y_annotation_data.values)

        def numpy_tuple(self_):
            return self_.x_count_data.values, self_.y_annotation_data.values

        modes = {
            ('xarray', True): xarray_single,
            ('xarray', False): xarray_tuple,
            ('pandas', True): pandas_single,
            ('pandas', False): pandas_tuple,
            ('numpy', True): numpy_single,
            ('numpy', False): numpy_tuple,
        }
        key = (output_type, single_object)
        # TODO: Clarify error message.
        if key not in modes.keys():
            raise ValueError(f'key given: {key} is not one of the available'
                             f'types: {list(modes.keys())}')

        return modes[key](self)


class CallableInterface(Interface, param.ParameterizedFunction):

    @transient_log_handler
    def __new__(cls, *args, **params):
        logger.debug('Creating a new GSForge.CallableInterface instance.')
        if args:
            params = Interface._interface_dispatch(*args, **params)
        inst = cls.instance(**params)  # See the param code for more on this `instance` function.
        return inst.__call__()

    def __call__(self):
        raise NotImplementedError
