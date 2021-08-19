from __future__ import annotations

import copy
import functools
import itertools
import logging
import os
from collections import defaultdict, UserDict
from functools import reduce
from pathlib import Path
from typing import Dict, Tuple, List, Union, Callable, IO, AnyStr, FrozenSet

import numpy as np
import pandas as pd
import param
import xarray as xr

from ._AnnotatedGEM import AnnotatedGEM
from ._GeneSet import GeneSet
from .._singledispatchmethod import singledispatchmethod

logger = logging.getLogger("GSForge")


# TODO: Add warning if .from_folder() or other creation function is empty.
class GeneSetDictionary(UserDict):
    """
    A dictionary with hooks to update support arrays.
    """
    # TODO: Add warnings if provided genes are not found within the given GEM index.

    @singledispatchmethod
    def __dispatch(*args, **params):
        raise TypeError(f"Source of type: {type(args[0])} not supported.")

    @__dispatch.register(list)
    @__dispatch.register(tuple)
    def __from_iterable(self, source):
        for gs in source:
            self.__setitem__(gs.name, gs)

    @__dispatch.register(dict)
    def __from_dict(self, source):
        for name, gs in source.items():
            self.__setitem__(name, gs)

    def __init__(self, parent_index=None, source=None):
        super().__init__()  # Gives us a dictionary as a .data attribute.
        if parent_index is not None:
            self.parent_index = parent_index
        if source:
            self.__dispatch(source)

    def __getitem__(self, item):
        return self.data[item]

    def __setitem__(self, key, gene_set):
        """Ensures that the support array references the correct index, given by parent_index."""
        # if gene_set.name != key:
        #     gene_set.param.set_param(name=key)
        # TODO: Add a check for a support array.
        #       If no support array is found and the shape is less than the complete index, use the
        #       provided index as an implicit support index.
        # Ensure all provided genes are within the parent index.
        # gene_set.data = gene_set.data.reindex({"Gene": self.parent_index})
        # Raise an error if this is not the case?
        # Update the geneset data index.
        # if gene_set.data.Gene.shape
        # TODO: Add a check so this only runs if the incoming geneset has a larger gene-index
        #   than the existing one. or perhaps through an error or warning instead.
        # updated_support_index = np.isin(self.parent_index, gene_set.get_support(), assume_unique=True)
        # gene_set.data[gene_set.support_index_name] = ((gene_set.gene_index_name,), updated_support_index)
        self.data[key] = gene_set


class GeneSetCollection(param.Parameterized):
    """
    An interface class which contains an AnnotatedGEM and a dictionary of GeneSet objects.
    """

    gem = param.ClassSelector(class_=AnnotatedGEM, allow_None=True, doc="A GSForge.AnnotatedGEM object.")

    def __init__(self, **params):
        logger.debug('Initializing a new gsforge.GeneSetCollection object...')
        gene_sets = None
        if 'gene_sets' in params:
            gene_sets = params.pop('gene_sets')
        super().__init__(**params)
        parent_index = self.gem.gene_index if self.gem else None
        self.gene_sets = GeneSetDictionary(parent_index=parent_index, source=gene_sets)
        logger.debug('GeneSetCollection initialization complete.')

    def __getitem__(self, item):
        return self.gene_sets[item]

    def __setitem__(self, key, value):
        self.gene_sets[key] = value

    def summarize_gene_sets(self) -> Dict[str, int]:
        """
        Summarize this GeneSetCollection, returns a dictionary of ``{gene_set_name: support_length}``.
        This is used to generate display used in the ``__repr__`` function.
        """
        counts = {key: len(gs.get_support()) for key, gs in self.gene_sets.items()}
        counts = {key: counts[key] for key in sorted(counts, key=counts.get, reverse=True)}
        return counts

    def __repr__(self) -> str:
        """
        Construct a user-friendly representation of this GeneSetCollection.
        """
        summary = [f"<GSForge.{type(self).__name__}>"]
        summary += [self.name]
        gene_set_info = self.summarize_gene_sets()
        summary += [f"GeneSets ({len(gene_set_info)} total): Support Count"]
        gene_summary = self.summarize_gene_sets()
        summary += [f"    {k}: {v}" for k, v in itertools.islice(self.summarize_gene_sets().items(), 5)]
        if len(gene_summary) > 5:
            summary += [f"    ... and {len(gene_summary) - 5} more."]
        return "\n".join(summary)

    def get_support(self, key: str) -> np.ndarray:
        """
        Get the support array for a given key.

        Parameters
        ----------
        key : str
            The GeneSet from which to get the gene support.

        Returns
        -------
        np.ndarray : An array of the genes that make up the support of this ``GeneSet``.
        """
        return self.gene_sets[key].get_support()

    def gene_sets_to_dataframes(self, keys: List[str] = None, only_supported: bool = True) -> Dict[str, pd.DataFrame]:
        """
        Returns a dictionary of {key: pd.DataFrame} of the ``GeneSet.data``. The DataFrame is limited
        to only those genes that are 'supported' within the GeneSet by default.

        Parameters
        ----------
        keys : List[str]
            An optional list of gene_set keys to return, by default all keys are selected.

        only_supported : bool
            Whether to return a subset defined by each GeneSet support, or the complete data frame.

        Returns
        -------
        dict : A dictionary of {key: pd.DataFrame} of the ``GeneSet.data`` attribute.
        """
        keys = self.gene_sets.keys() if keys is None else keys
        keys = [key for key in keys if self.gene_sets[key].support_exists]
        return {key: self.gene_sets[key].to_dataframe(only_supported) for key in keys}

    # TODO: Consider overwrite protection.
    def gene_sets_to_csv_files(self, target_dir: str = None, keys: List[str] = None,
                               only_supported: bool = True) -> None:
        """
        Writes GeneSet.data as .csv files.

        By default this creates creates a folder with the current working directory and saves the .csv
        files within. By default only genes that are "supported" by a GeneSet are included.

        Parameters
        ----------
        target_dir :
            The target directory to save the .csv files to. This defaults to the name of this
            GeneSetCollection, which creates a folder in the current working directory.

        keys : List[str]
            An optional list of gene_set keys to return, by default all keys are selected.

        only_supported : bool
            Whether to return a subset defined by each GeneSet support, or the complete data frame.

        Returns
        -------
        None
        """
        keys = self.gene_sets.keys() if keys is None else keys
        keys = [key for key in keys if self.gene_sets[key].support_exists]
        target_dir = self.name if target_dir is None else target_dir
        data_frames = self.gene_sets_to_dataframes(keys, only_supported)
        os.makedirs(target_dir, exist_ok=True)
        for name, df in data_frames.items():
            df.to_csv(f"{target_dir}/{name}.csv")

    # TODO: Consider overwrite protection.
    def gene_sets_to_excel_sheet(self, name: str = None, keys: List[str] = None, only_supported: bool = True) -> None:
        """
        Writes the GeneSet.data within this GeneSetCollection as a single Excel worksheet.

        By default this sheet is named using the ``.name`` of this GeneSetCollection. By default
        only genes that are "supported" by a GeneSet are included.

        Parameters
        ----------
        name : str
            The name of the Excel sheet. ``.xlsx`` will be appended to the given name.

        keys : List[str]
             An optional list of gene_set keys to return, by default all keys are selected.

        only_supported : bool
            Whether to return a subset defined by each GeneSet support, or the complete data frame.

        Returns
        -------
        None
        """
        keys = self.gene_sets.keys() if keys is None else keys
        keys = [key for key in keys if self.gene_sets[key].support_exists]
        name = f'{self.name}.xlsx' if name is None else f'{name}.xlsx'
        data_frames = self.gene_sets_to_dataframes(keys, only_supported)
        with pd.ExcelWriter(name) as writer:
            for set_name, df in data_frames.items():
                df.to_excel(writer, sheet_name=set_name)

    def _as_dict(self, keys: Tuple[AnyStr]) -> Dict[str, np.ndarray]:
        return copy.deepcopy({key: self.gene_sets[key].get_support() for key in keys})

    def _parse_keys(self, keys: List[str] = None, exclude: List[str] = None) -> Tuple:
        if isinstance(keys, np.ndarray):
            keys = keys.tolist()

        logger.debug(f'Parsing keys: {keys}, excluding: {exclude}')

        # Use all keys in the collection if none are provided.
        if keys == [None] or keys is None:
            # keys = list(self.gene_sets.keys()) if keys is None else list(keys)
            keys = list(self.gene_sets.keys())
        else:
            keys = list(keys)
        exclude = [] if exclude is None else exclude

        # Ensure all keys provided are actually within the collection.
        for key in keys + exclude:
            if key not in self.gene_sets.keys():
                raise ValueError(f"Key {key} not found in available keys:\n{list(self.gene_sets.keys())}")

        if exclude is not None:
            keys = [key for key in keys if key not in exclude]

        # Sort the remaining keys and pass them as a 'hash-able' tuple to the cached function.
        sorted_keys = tuple(sorted(keys))
        logger.debug(f'Parsed to:  {sorted_keys}')

        return sorted_keys

    def as_dict(self, keys: List[str] = None, exclude: List[str] = None,
                empty_supports: bool = False) -> Dict[str, np.ndarray]:
        """
        Returns a dictionary of {name: supported_genes} for each GeneSet, or those specified
        by the `keys` argument.

        Parameters
        ----------
        keys : List[str]
            An optional list of gene_set keys to return, by default all keys are selected.

        exclude : List[str]
            An optional list of `GeneSet` keys to exclude from the returned dictionary.

        empty_supports:
            Whether to include GeneSets that have no support array, or no genes supported within
            the support array.

        Returns
        -------
        dict : Dictionary of {name: supported_genes} for each GeneSet.
        """
        sorted_keys = self._parse_keys(keys, exclude)
        return self._as_dict(sorted_keys)

    def _intersection(self, keys: Tuple[AnyStr]) -> np.ndarray:
        gene_set_dict = self._as_dict(keys)
        return reduce(np.intersect1d, gene_set_dict.values())

    def intersection(self, keys: List[str] = None, exclude: List[str] = None) -> np.ndarray:
        """
        Return the intersection of supported genes in this GeneSet collection.

        Parameters
        ----------
        keys : List[str]
            An optional list of gene_set keys to return, by default all keys are selected.

        exclude : List[str]
            An optional list of `GeneSet` keys to exclude from the returned dictionary.

        Returns
        -------
        np.ndarray : Intersection of the supported genes within GeneSets.
        """
        sorted_keys = self._parse_keys(keys, exclude)
        return self._intersection(sorted_keys)

    def _union(self, keys: Tuple[AnyStr]) -> np.ndarray:
        gene_set_dict = self._as_dict(keys)
        return reduce(np.union1d, gene_set_dict.values())

    def union(self, keys: List[str] = None, exclude: List[str] = None) -> np.ndarray:
        """
        Get the union of supported genes in this GeneSet collection.

        Parameters
        ----------
        keys : List[str]
            An optional list of gene_set keys to return, by default all keys are selected.

        exclude : List[str]
            An optional list of `GeneSet` keys to exclude from the returned dictionary.

        Returns
        -------
        np.ndarray : Union of the supported genes within GeneSets.
        """
        sorted_keys = self._parse_keys(keys, exclude)
        return self._union(sorted_keys)

    def _difference(self, primary_key: str, other_keys: FrozenSet[str], mode: str) -> np.ndarray:
        modes = {'union': self.union, 'intersection': self.intersection}
        other_set = modes[mode](other_keys)
        primary_set = self.gene_sets[primary_key].get_support()
        return np.setdiff1d(primary_set, other_set)

    def difference(self, primary_key: str, other_keys: List[str] = None, mode: str = 'union') -> np.ndarray:
        """
        Finds the genes within `primary_key` that are not within the `mode` of the sets
        given in `other_keys`.

        If no `other_keys` are provided, all remaining keys are used.
        The default `mode` is `union`.

        Parameters
        ----------
        primary_key : List[str]
            The set

        other_keys : List[str]
            An optional list of `GeneSet` keys...

        mode : str
            Mode by which to join the GeneSets given by `other_keys`.

        Returns
        -------
        np.ndarray
            ...
        """

        if other_keys is None:
            other_keys = [key for key in self.gene_sets.keys() if key != primary_key]
        other_keys_set = frozenset(sorted(other_keys))

        valid_modes = {'union', 'intersection'}

        if mode not in valid_modes:
            raise ValueError(f"Given mode {mode} is not a valid mode. ({valid_modes})")

        return self._difference(primary_key, other_keys_set, mode)

    def joint_difference(self, primary_keys: List[str], other_keys: List[str] = None,
                         primary_join_mode: str = 'union',
                         others_join_mode: str = 'union'):
        """

        Parameters
        ----------
        primary_keys
        other_keys
        primary_join_mode
        others_join_mode

        Returns
        -------

        """
        if other_keys is None:
            other_keys = [key for key in self.gene_sets.keys() if key not in primary_keys]

        join_modes = {'union': self.union, 'intersection': self.intersection}
        for mode in [primary_join_mode, others_join_mode]:
            if mode not in join_modes.keys():
                raise ValueError(f"Given mode {mode} is not a valid mode. ({list(join_modes.keys())})")

        primary_join_function = join_modes[primary_join_mode]
        others_join_function = join_modes[others_join_mode]

        primary_join = primary_join_function(primary_keys)
        others_join = others_join_function(other_keys)

        joint_difference = np.setdiff1d(primary_join, others_join)

        return joint_difference

    def pairwise_unions(self, keys: List[str] = None, exclude: List[str] = None) -> Dict[Tuple[str, str], np.ndarray]:
        """
        Construct pairwise permutations of GeneSets within this collection, and return
        the union of each pair in a dictionary.

        Parameters
        ----------
        keys : List[str]
            An optional list of gene_set keys to return, by default all keys are selected.

        exclude : List[str]
            An optional list of `GeneSet` keys to exclude from the returned dictionary.

        Returns
        -------
        dict : A dictionary of ``{(GeneSet.name, GeneSet.name): gene support union}``.
        """
        gene_set_dict = self.as_dict(keys, exclude)
        return {(ak, bk): np.union1d(av, bv)
                for (ak, av), (bk, bv) in itertools.permutations(gene_set_dict.items(), 2)
                if ak != bk}

    def pairwise_intersection(self, keys: List[str] = None, exclude: List[str] = None
                              ) -> Dict[Tuple[str, str], np.ndarray]:
        """
        Construct pairwise combinations of GeneSets within this collection, and return
        the intersection of each pair in a dictionary.

        Parameters
        ----------
        keys : List[str]
            An optional list of gene_set keys to return, by default all keys are selected.

        exclude : List[str]
            An optional list of `GeneSet` keys to exclude from the returned dictionary.

        Returns
        -------
        dict : A dictionary of ``{GeneSet.Name, GeneSet.name): GeneSets.get_support() intersection}``.
        """
        gene_set_dict = self.as_dict(keys, exclude)

        return {(ak, bk): np.intersect1d(av, bv)
                for (ak, av), (bk, bv) in itertools.combinations(gene_set_dict.items(), 2)
                if ak != bk}

    def pairwise_percent_intersection(self, keys=None, exclude=None) -> List[Tuple[str, str, float]]:
        """
        Construct pairwise permutations of GeneSets within this collection, and return
        the intersection of each pair within a dictionary.

        Parameters
        ----------
        keys : List[str]
            An optional list of gene_set keys to return, by default all keys are selected.

        exclude : List[str]
            An optional list of `GeneSet` keys to exclude from the returned dictionary.

        Returns
        -------
        dict : A dictionary of ``{GeneSet.Name, GeneSet.name): percent gene intersection}``.
        """
        gene_set_dict = self.as_dict(keys, exclude)
        zero_filtered_dict = {k: v for k, v in gene_set_dict.items()
                              if len(v) > 0}
        return [(ak, bk, len(np.intersect1d(av, bv)) / len(av))
                for (ak, av), (bk, bv) in itertools.permutations(zero_filtered_dict.items(), 2)
                if ak != bk]

    def construct_standard_specification(self, include: List[str] = None, exclude=None) -> dict:
        """
        Construct a standard specification that can be used to view unions, intersections and
        differences (unique genes) of the sets within this collection.

        Parameters
        ----------
        include : List[str]
            An optional list of gene_set keys to return, by default all keys are selected.

        exclude : List[str]
            An optional list of `GeneSet` keys to exclude from the returned dictionary.

        Returns
        -------
        dict: A specification dictionary.
        """

        include = self._parse_keys(include, exclude)
        standard_spec = defaultdict(list)

        standard_spec['union'].append({'name': f'{self.name}__standard_union',
                                       'keys': include})

        standard_spec['intersection'].append({'name': f'{self.name}__standard_intersection',
                                              'keys': include})

        for primary_key in include:
            other_keys = [key for key in include if primary_key != key]
            standard_spec['difference'].append({'name': f'{self.name}__{primary_key}__unique',
                                                'primary_key': primary_key,
                                                'other_keys': other_keys})

        return standard_spec

    @staticmethod
    def merge_specifications(*specs):
        """
        Merges sets of defaultdict(list) objects with common keys.
        """
        # TODO: Document me.
        keys = functools.reduce(set.union, [set(d.keys()) for d in specs])
        output_spec = defaultdict(list)
        for key in keys:
            for spec in specs:
                if spec.get(key):
                    output_spec[key].extend(spec.get(key))
        return output_spec

    def process_set_operation_specification(self, specification: dict = None) -> dict:
        """
        Calls and stores the results from a specification. The specification must declare
        set operation functions and their arguments.

        Parameters
        ----------
        specification : Dict
        """
        # TODO: Add input validation for the specification.

        function_map = {
            'intersection': self.intersection,
            'union': self.union,
            'difference': self.difference,
            'joint_difference': self.joint_difference,
            'pairwise_unions': self.pairwise_unions,
            'pairwise_intersection': self.pairwise_intersection,
            'pairwise_percent_intersection': self.pairwise_percent_intersection,
        }

        processed_spec = dict()

        for key, function in function_map.copy().items():
            if specification.get(key):
                for entry in specification.get(key):
                    entry_kwargs = {k: v for k, v in entry.items() if k != 'name'}
                    name = entry.get('name')
                    processed_spec[name] = function(**entry_kwargs)

        return processed_spec

    @classmethod
    def from_specification(cls, source_collection, specification=None, name="processed_specification"):

        if specification is None:
            specification = source_collection.construct_standard_specification()

        processed_specification = source_collection.process_set_operation_specification(specification)

        collection = cls(gem=source_collection.gem, name=name)

        for key, genes in processed_specification.items():
            collection[key] = GeneSet.from_gene_array(
                name=key,
                selected_gene_array=genes,
                complete_gene_index=source_collection.gem.gene_index)

        collection.gene_sets = {**collection.gene_sets, **source_collection.gene_sets}

        return collection

    @classmethod
    def from_folder(cls, gem: AnnotatedGEM, target_dir: Union[str, Path, IO[AnyStr]], glob_filter: str = "*.nc",
                    filter_func: Callable = None, **params) -> GeneSetCollection:
        """
        Create a `GeneSetCollection` from a directory of saved GeneSet objects.

        The file name of each gene_set.nc file will be used as the key in the `gene_sets` dictionary.

        Parameters
        ----------
        gem : AnnotatedGEM
            A `GSForge.AnnotatedGEM` object.

        target_dir : Union[str, Path, IO[AnyStr]]
            The directory which contains the saved GeneSet .netcdf files.

        glob_filter : str
            A glob by which to restrict the files found within `target_dir`.

        filter_func : Callable
             A function by which to filter which `xarray.Dataset` objects are included.
             This function should take an `xarray.Dataset` and return a boolean.

        params :
            Parameters to configure the GeneSetCollection.

        Returns
        -------
        GeneSetCollection : A new GeneSetCollection.
        """
        # TODO: Add a warning if the glob returns nothing.
        gene_sets = dict()
        for file in Path(target_dir).expanduser().resolve().glob(glob_filter):
            data = xr.open_dataset(file)

            if filter_func is not None:
                if filter_func(data):
                    pass
                else:
                    continue

            data = xr.align(data, gem.gene_index, join="outer")[0]
            new_gene_set = GeneSet.from_xarray_dataset(data=data)

            if data.attrs.get("__GSForge.GeneSet.params.name") is None:
                name = os.path.basename(str(file)).rsplit(".nc")[0]
            else:
                name = new_gene_set.name

            gene_sets[name] = new_gene_set

        return cls(gem=gem, gene_sets=gene_sets, **params)

    def save(self, target_dir: str, keys: List[str] = None) -> None:
        """
        Save  this collection to ``target_dir``. Each GeneSet will be saved as a separate
        .netcdf file within this directory.

        Parameters
        ----------
        target_dir : str
            The path to which GeneSet ``xarray.Dataset`` .netcdf files will be written.

        keys : List[str]
            The list of GeneSet keys that should be saved. If this is not provided, all
            GeneSet objects are saved.

        Returns
        -------
        None
        """
        # Save all the gene sets in this collection in the target_dir path.
        if keys is None:
            keys = self.gene_sets.keys()

        # Create any needed intermediate directories.
        if not os.path.exists(target_dir):
            os.makedirs(target_dir)

        for key in keys:
            save_path = self.gene_sets[key].save_as_netcdf(target_dir)
            print(save_path)
