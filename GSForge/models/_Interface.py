from __future__ import annotations

import logging
from textwrap import dedent
from typing import Union
import warnings

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
    """The Interface provides common API access for interacting with the ``AnnotatedGEM``
    and ``GeneSetCollection`` objects.
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
    """), default='complete', objects=["complete", "dropped"], precedence=-1.0)

    count_transform = param.Callable(default=None, precedence=-1.0, doc=dedent("""\
    A transform that will be run on the `x_data` that is supplied by this Interface. 
    The transform runs on the subset of the matrix that has been selected."""))

    @singledispatchmethod
    def _interface_dispatch(*args, **params):
        raise TypeError(f"Source of type: {type(args[0])} not supported.")

    def __init__(self, *args, **params):

        # If the user passes a string, place it as a single item within a list.
        if isinstance(params.get("annotation_variables"), str):
            params["annotation_variables"] = [params.get("annotation_variables")]

        if isinstance(params.get("selected_gene_sets"), str):
            params["selected_gene_sets"] = [params.get("selected_gene_sets")]

        if args:
            params = self._interface_dispatch(*args, **params)

        super().__init__(**params)

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
            # We need to load the data if a user-list is supplied to prevent some straneg issues
            # with nested numpy arrays being returned.
            self.gem.data.load()
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

        selected_samples = subset[self.gem.sample_index_name].copy(deep=True).values
        return selected_samples

    @property
    def get_selection_indices(self) -> dict:
        """Returns the currently selected indexes as a dictionary."""
        return {self.gem.gene_index_name: self.get_gene_index(),
                self.gem.sample_index_name: self.get_sample_index()}

    @property
    def x_count_data(self) -> Union[xr.DataArray, None]:
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
            "complete": lambda c: c.fillna(0),
            "masked": lambda c: c.where(c > 0.0),
            "dropped": lambda c: c.where(c > 0.0).dropna(dim=self.gem.gene_index_name),
        }

        # Ensure the correct count array is selected.
        count_variable = self.count_variable if self.count_variable is not None else self.gem.count_array_name

        logger.info(f'Preparing count data from the variable {count_variable}.')
        if self.selected_gene_sets == ['all']:
            self.param.set_param(selected_gene_sets=list(self.gene_set_collection.gene_sets.keys()))

        if self.selected_genes is not None:
            support = self.selected_genes
            logger.info(f'Gene selection of: {support.shape[0]} genes provided.')

        elif self.selected_gene_sets == [None]:
            support = self.gem.gene_index

        # Otherwise, use  some combination of GeneSet supports should be used.
        else:
            support = gene_set_combinations[self.gene_set_mode](self.selected_gene_sets)
            logger.info(
                f'Selected {len(self.selected_gene_sets)} GeneSets, '
                f'using mode: {self.gene_set_mode} for a support of size: {support.shape[0]}.')

        # Now the counts can be selected based on the gene support, and the selected samples.
        sample_support = self.get_sample_index()

        # Check the overlap of the selected genes with those in the GEM index.
        # Give the user a warning if any unavailable genes are requested, but still
        # return the available genes.
        avail_support = np.intersect1d(self.gem.gene_index, support)
        if len(avail_support) < len(support):
            diff = len(support) - len(avail_support)
            warnings.warn(
                f'{diff} Unavailable genes of the {len(support)} requested. '
                f'Using the {len(avail_support)} available.',
                UserWarning)

        counts = self.gem.data.sel({self.gem.gene_index_name: avail_support,
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
        # If only one label has been selected, return this as an xarray.DataArray.
        if len(self.annotation_variables) == 1:
            return self.gem.data[self.annotation_variables].sel(
                {self.gem.sample_index_name: sample_index})[self.annotation_variables[0]].copy(deep=True)

        return self.gem.data[self.annotation_variables].sel(
                {self.gem.sample_index_name: sample_index}).copy(deep=True)

    def get_gem_data(self, single_object=False, output_type='xarray', **params):
        """
        Returns count [and annotation] data based on the current parameters.

        Users should call gsf.get_gem_data
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
            # It may be faster to call .values and create a new dataframe instead of unstacking.
            logger.info('Returning a tuple of counts and annotations each as a pandas.DataFrame.')
            if self_.y_annotation_data is not None:
                return self_.x_count_data.to_dataframe().unstack().droplevel(0, axis=1), \
                       self_.y_annotation_data.to_dataframe()
            else:
                return self_.x_count_data.to_dataframe().unstack().droplevel(0, axis=1), self_.y_annotation_data

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
            # TODO: Convert to a helper function or consider the implementation in operations.core (commented out).
            if isinstance(params.get("annotation_variables"), str):
                params["annotation_variables"] = [params.get("annotation_variables")]

            if isinstance(params.get("selected_gene_sets"), str):
                params["selected_gene_sets"] = [params.get("selected_gene_sets")]
        inst = cls.instance(**params)  # See the param code for more on this `instance` function.
        return inst.__call__()

    def __call__(self):
        raise NotImplementedError
