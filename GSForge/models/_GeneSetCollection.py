from __future__ import annotations

import itertools
import os
from pathlib import Path
from functools import reduce
from textwrap import dedent
from typing import Dict, Tuple, List, Union, Callable, IO, AnyStr, FrozenSet

import methodtools
import numpy as np
import pandas as pd
import param
import xarray as xr

from ._AnnotatedGEM import AnnotatedGEM
from ._GeneSet import GeneSet


class GeneSetCollection(param.Parameterized):
    """
    An interface class which contains an AnnotatedGEM and a dictionary of GeneSet objects.
    """

    ###############################################################################################
    # PARAMETERS
    # See the param documentation.
    ###############################################################################################

    gem = param.ClassSelector(class_=AnnotatedGEM, doc=dedent("""\
    A GSForge.AnnotatedGEM object."""))

    gene_sets = param.Dict(doc=dedent("""\
    A dictionary of `{key: GSForge.GeneSet}`."""))

    ###############################################################################################
    # PRIVATE FUNCTIONS
    ###############################################################################################

    def __init__(self, **params):
        super().__init__(**params)
        if self.gene_sets is None:
            self.set_param(gene_sets=dict())

    def __getitem__(self, item):
        return self.gene_sets[item]

    def __setitem__(self, key, value):
        self.gene_sets[key] = value

    def summarize_gene_sets(self) -> Dict[str, int]:
        """
        Summarize this GeneSetCollection, returns a dictionary of ``{gene_set_name: support_length}``.
        This is used to generate display used in the ``__repr__`` function.
        """
        counts = {key: len(gs.gene_support()) for key, gs in self.gene_sets.items()}
        counts = {key: counts[key] for key in sorted(counts, key=counts.get, reverse=True)}
        return counts

    def __repr__(self) -> str:
        summary = [f"<GSForge.{type(self).__name__}>"]
        summary += [self.name]
        # summary += [indent(self.gem.name, "    ")]
        gene_set_info = self.summarize_gene_sets()
        summary += [f"GeneSets ({len(gene_set_info)} total): Support Count"]
        summary += [f"    {k}: {v}" for k, v in itertools.islice(self.summarize_gene_sets().items(), 10)]
        return "\n".join(summary)

    ###############################################################################################
    # PUBLIC FUNCTIONS
    ###############################################################################################

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
        return self.gene_sets[key].gene_support()

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

    @methodtools.lru_cache()
    def _as_dict(self, keys: Tuple[AnyStr]) -> Dict[str, np.ndarray]:
        return {key: self.gene_sets[key].gene_support() for key in keys}

    def _parse_keys(self, keys: List[str] = None, exclude: List[str] = None,
                    empty_supports: bool = False) -> Tuple[AnyStr]:
        # Use all keys in the collection if none are provided.
        keys = self.gene_sets.keys() if keys is None else keys

        # Ensure all keys provided are actually within the collection.
        if not all(key in self.gene_sets.keys() for key in keys):
            raise ValueError(f"Not all keys given were found in the available keys: {list(self.gene_sets.keys())}")

        if empty_supports is False:
            keys = [key for key in keys if self.gene_sets[key].support_exists]

        if exclude is not None:
            keys = [key for key in keys if key not in exclude]

        # Sort the remaining keys and pass them as a 'hash-able' tuple to the cached function.
        sorted_keys = tuple(sorted(keys))
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
        sorted_keys = self._parse_keys(keys, exclude, empty_supports)
        return self._as_dict(sorted_keys)

    @methodtools.lru_cache()
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

    @methodtools.lru_cache()
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

    @methodtools.lru_cache()
    def _difference(self, primary_key: str, other_keys: FrozenSet[str], mode: str) -> np.ndarray:
        modes = {'union': self.union, 'intersection': self.intersection}
        other_set = modes[mode](other_keys)
        primary_set = self.gene_sets[primary_key].gene_support()
        return np.setdiff1d(primary_set, other_set)

    # def difference(self, keys: List[str] = None, exclude: List[str] = None) -> np.ndarray:
    def difference(self, primary_key: str, other_keys: List[str] = None, mode: str = 'union'):
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
        dict : A dictionary of ``{GeneSet.Name, GeneSet.name): GeneSets.gene_support() intersection}``.
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

    ###############################################################################################
    # CONSTRUCTOR FUNCTIONS
    ###############################################################################################
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
