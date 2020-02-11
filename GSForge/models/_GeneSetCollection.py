import itertools
import os
import pathlib
from functools import reduce
from textwrap import dedent
from typing import Dict, Tuple, List, TypeVar, Type, Union

import numpy as np
import pandas as pd
import param
import xarray as xr

from ._AnnotatedGEM import AnnotatedGEM
from ._GeneSet import GeneSet

# from os import PathLike

# Declare typing hints.
# Typing hints are not used by the Python runtime in any way.
# They are used by third party tools and linters.
# In python 3.8 we can use forward references.
TypeGSC = TypeVar("TypeGSC", bound="GeneSetCollection")
TypePath = Union[str, 'PathLike[Any]']


class GeneSetCollection(param.Parameterized):
    """
    An interface class which contains an AnnotatedGEM and a dictionary of GeneSet objects.
    """

    gem = param.ClassSelector(class_=AnnotatedGEM, doc=dedent("""\
    A GSForge.AnnotatedGEM object."""))

    gene_sets = param.Dict(doc=dedent("""\
    A dictionary of `{key: GSForge.GeneSet}`."""))

    def __init__(self, **params):
        super().__init__(**params)
        if self.gene_sets is None:
            self.set_param(gene_sets=dict())

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

    def get_support(self, key) -> np.ndarray:
        """
        Get the support array for a given key.

        :param key:
            The GeneSet from which to get the gene support.

        :return:
            A ``numpy`` array of the genes that make up the support of this ``GeneSet``.
        """
        return self.gene_sets[key].gene_support()

    def save(self, target_dir, keys=None) -> None:
        """
        Save  this collection to ``target_dir``. Each GeneSet will be saved as a separate
        .netcdf file within this directory.

        :param target_dir:
            The path to which GeneSet ``xarray.Dataset`` .netcdf files will be written.

        :param keys: The list of GeneSet keys that should be saved. If this is not provided, all
            GeneSet objects are saved.

        :return: None
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

    def gene_sets_to_dataframes(self, keys=None, only_supported: bool = True) -> Dict[str, pd.DataFrame]:
        """
        Returns a dictionary of {key: pd.DataFrame} of the ``GeneSet.data``. The DataFrame is limited
        to only those genes that are 'supported' within the GeneSet by default.

        :param keys:
            An optional list of gene_set keys to return, by default all keys are selected.

        :param only_supported:
            Whether to return a subset defined by each GeneSet support, or the complete data frame.

        :return:
            A dictionary of {key: pd.DataFrame} of the `GeneSet.data`.
        """
        keys = self.gene_sets.keys() if keys is None else keys
        keys = [key for key in keys if self.gene_sets[key].support_exists]
        return {key: self.gene_sets[key].to_dataframe(only_supported) for key in keys}

    # TODO: Consider overwrite protection.
    def gene_sets_to_csv_files(self, target_dir=None, keys=None, only_supported: bool = True) -> None:
        """
        Writes GeneSet.data as .csv files.

        By default this creates creates a folder with the current working directory and saves the .csv
        files within. By default only genes that are "supported" by a GeneSet are included.

        :param target_dir:
            The target directory to save the .csv files to. This defaults to the name of this
            GeneSetCollection, which creates a folder in the current working directory.

        :param keys:
            An optional list of gene_set keys to return, by default all keys are selected.

        :param only_supported:
            Whether to return a subset defined by each GeneSet support, or the complete data frame.

        :return: None
        """
        keys = self.gene_sets.keys() if keys is None else keys
        keys = [key for key in keys if self.gene_sets[key].support_exists]
        target_dir = self.name if target_dir is None else target_dir
        data_frames = self.gene_sets_to_dataframes(keys, only_supported)
        os.makedirs(target_dir, exist_ok=True)
        for name, df in data_frames.items():
            df.to_csv(f"{target_dir}/{name}.csv")

    # TODO: Consider overwrite protection.
    def gene_sets_to_excel_sheet(self, name: str = None, keys=None, only_supported: bool = True) -> None:
        """
        Writes the GeneSet.data within this GeneSetCollection as a single Excel worksheet.

        By default this sheet is named using the `.name` of this GeneSetCollection. By default
        only genes that are "supported" by a GeneSet are included.

        :param name:
            The name of the Excel sheet. `.xlsx` will be appended to the given name.

        :param keys:
             An optional list of gene_set keys to return, by default all keys are selected.

        :param only_supported:
            Whether to return a subset defined by each GeneSet support, or the complete data frame.

        :return: None
        """
        keys = self.gene_sets.keys() if keys is None else keys
        keys = [key for key in keys if self.gene_sets[key].support_exists]
        name = f'{self.name}.xlsx' if name is None else f'{name}.xlsx'
        data_frames = self.gene_sets_to_dataframes(keys, only_supported)
        with pd.ExcelWriter(name) as writer:
            for set_name, df in data_frames.items():
                df.to_excel(writer, sheet_name=set_name)

    def as_dict(self, keys: List[str] = None, exclude: List[str] = None,
                empty_supports: bool = False) -> Dict[str, np.ndarray]:
        """
        Returns a dictionary of {name: supported_genes} for each GeneSet, or those specified
        by the `keys` argument.

        :param keys:
            An optional list of gene_set keys to return, by default all keys are selected.

        :param exclude:
            An optional list of `GeneSet` keys to exclude from the returned dictionary.

        :param empty_supports:
            Whether to include GeneSets that have no support array, or no genes supported within
            the support array.

        :return:
            Dictionary of {name: supported_genes} for each GeneSet.
        """
        keys = self.gene_sets.keys() if keys is None else keys

        if empty_supports is False:
            keys = [key for key in keys if self.gene_sets[key].support_exists]

        if exclude is not None:
            keys = [key for key in keys if key not in exclude]

        return {key: self.gene_sets[key].gene_support() for key in keys}

    def intersection(self, keys: List[str] = None, exclude: List[str] = None) -> np.ndarray:
        """
        Return the intersection of supported genes in this GeneSet collection.

        :param keys:
            An optional list of gene_set keys to return, by default all keys are selected.

        :param exclude:
            An optional list of `GeneSet` keys to exclude from the returned dictionary.

        :return:
            Intersection of the supported genes within GeneSets.
        """
        gene_set_dict = self.as_dict(keys, exclude)
        return reduce(np.intersect1d, gene_set_dict.values())

    def union(self, keys: List[str] = None, exclude: List[str] = None) -> np.ndarray:
        """
        Get the union of supported genes in this GeneSet collection.

        :param keys:
            An optional list of gene_set keys to return, by default all keys are selected.

        :param exclude:
            An optional list of `GeneSet` keys to exclude from the returned dictionary.

        :return:
            Union of the supported genes within GeneSets.
        """
        gene_set_dict = self.as_dict(keys, exclude)
        return reduce(np.union1d, gene_set_dict.values())

    def difference(self, keys: List[str] = None, exclude: List[str] = None) -> np.ndarray:
        """
        Get the difference of supported genes in this GeneSet collection.

        :param keys:
            An optional list of gene_set keys to return, by default all keys are selected.

        :param exclude:
            An optional list of `GeneSet` keys to exclude from the returned dictionary.

        :return:
            Difference of the supported genes within GeneSets.
        """
        gene_set_dict = self.as_dict(keys, exclude)
        return reduce(np.setdiff1d, gene_set_dict.values())

    def pairwise_unions(self, keys: List[str] = None, exclude: List[str] = None) -> Dict[str, np.ndarray]:
        """
        Construct pairwise permutations of GeneSets within this collection, and return
        the union of each pair in a dictionary.

        :param keys:
            An optional list of gene_set keys to return, by default all keys are selected.

        :param exclude:
            An optional list of `GeneSet` keys to exclude from the returned dictionary.

        :return:
            A dictionary of `{(GeneSet.name, GeneSet.name): gene support union}`.
        """
        gene_set_dict = self.as_dict(keys, exclude)
        return {(ak, bk): np.union1d(av, bv)
                for (ak, av), (bk, bv) in itertools.permutations(gene_set_dict.items(), 2)
                if ak != bk}

    def pairwise_intersection(self, keys: List[str] = None, exclude: List[str] = None) -> Dict[str, np.ndarray]:
        """
        Construct pairwise combinations of GeneSets within this collection, and return
        the intersection of each pair in a dictionary.

        :param keys:
            An optional list of gene_set keys to return, by default all keys are selected.

        :param exclude:
            An optional list of `GeneSet` keys to exclude from the returned dictionary.

        :return:
            A dictionary of `{GeneSet.Name, GeneSet.name): GeneSets.gene_support() intersection}`.
        """
        gene_set_dict = self.as_dict(keys, exclude)

        return {(ak, bk): np.intersect1d(av, bv)
                for (ak, av), (bk, bv) in itertools.combinations(gene_set_dict.items(), 2)
                if ak != bk}

    def pairwise_percent_intersection(self, keys=None, exclude=None) -> List[Tuple[str, str, float]]:
        """
        Construct pairwise permutations of GeneSets within this collection, and return
        the intersection of each pair within a dictionary.

        :param keys:
            An optional list of gene_set keys to return, by default all keys are selected.

        :param exclude:
            An optional list of `GeneSet` keys to exclude from the returned dictionary.

        :return:
            A dictionary of `{GeneSet.Name, GeneSet.name): percent gene intersection}`.
        """
        gene_set_dict = self.as_dict(keys, exclude)
        zero_filtered_dict = {k: v for k, v in gene_set_dict.items()
                              if len(v) > 0}
        return [(ak, bk, len(np.intersect1d(av, bv)) / len(av))
                for (ak, av), (bk, bv) in itertools.permutations(zero_filtered_dict.items(), 2)
                if ak != bk]

    # def get_gene_sets_data(self, variables: list, keys=None) -> Dict[str, np.ndarray]:
    #     """
    #     Extracts variables from each Dataset, and returns them as a dictionary.
    #
    #     If you need these to be collated into a dataset, use `collate_gene_sets_data()`.
    #
    #     :param variables:
    #     :param keys:
    #     :return:
    #     """
    #     keys = self.gene_sets.keys() if keys is None else keys
    #     return {key: self.gene_sets[key].gene_support()[variables] for key in keys}

    # def collate_gene_sets_data(self, variables: list, keys=None) -> xr.Dataset:
    #     keys = self.gene_sets.keys() if keys is None else keys
    #     data_dict = self.get_gene_sets_data(variables, keys)
    #     # TODO: Change this call to ensure that it returns an xr.Dataset.
    #     key_ds = xr.concat(data_dict.values(), dim="gene_set")
    #     key_ds["gene_set"] = keys
    #     return key_ds

    @classmethod
    def from_folder(cls: Type[TypeGSC], gem: AnnotatedGEM, target_dir, glob_filter="*.nc", filter_func=None,
                    **params) -> TypeGSC:
        """
        Create a `GeneSetCollection` from a directory of saved GeneSet objects.

        The file name of each gene_set.nc file will be used as the key in the `gene_sets` dictionary.

        :param gem:
            A `GSForge.AnnotatedGEM` object.

        :param target_dir:
            The directory which contains the saved GeneSet .netcdf files.

        :param glob_filter:
            A glob by which to restrict the files found within `target_dir`.

        :param filter_func:
             A function by which to filter which `xarray.Dataset` objects are included.
             This function should take an `xarray.Dataset` and return a boolean.

        :param params:
            Parameters to configure the GeneSetCollection.

        :return:
            A new GeneSetCollection.
        """
        # TODO: Add a warning if the glob returns nothing.
        gene_sets = dict()
        for file in pathlib.Path(target_dir).expanduser().resolve().glob(glob_filter):
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
