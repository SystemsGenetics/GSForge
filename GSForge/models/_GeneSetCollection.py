import pathlib
import param
import os
import pandas as pd
import itertools
import xarray as xr
import numpy as np
from functools import reduce

from textwrap import dedent, indent

from ._AnnotatedGEM import AnnotatedGEM
from ._GeneSet import GeneSet


class GeneSetCollection(param.Parameterized):
    """
    An interface class which contains an AnnotatedGEM and a dictionary of GeneSet objects.
    """

    gem = param.ClassSelector(class_=AnnotatedGEM, doc=dedent("""\
    A Gene Expression Matrix (GEM) object."""))

    gene_sets = param.Dict(doc=dedent("""\
    A dictionary of {key: xarray.DataArray}, boolean arrays indicating
    support for a given gene."""))

    def __init__(self, **params):
        super().__init__(**params)
        if self.gene_sets is None:
            self.set_param(gene_sets=dict())

    def _summarize_gene_sets(self):
        """Summarize this GeneSetCollection."""
        counts = {key: len(gs.gene_support()) for key, gs in self.gene_sets.items()}
        counts = {key: counts[key] for key in sorted(counts, key=counts.get, reverse=True)}
        return counts

    def __repr__(self):
        summary = [f"<GSForge.{type(self).__name__}>"]
        summary += [self.name]
        # summary += [indent(self.gem.__repr__(), "    ")]
        gene_set_info = self._summarize_gene_sets()
        summary += [f"GeneSets ({len(gene_set_info)} total): Support Count"]
        summary += [f"    {k}: {v}" for k, v in itertools.islice(self._summarize_gene_sets().items(), 10)]
        return "\n".join(summary)

    def get_support(self, key) -> np.array:
        """Get the support array for a given key.

        :param key: The GeneSet from which to get the gene support.

        :return: A `numpy` array of the genes that make up the support of this `GeneSet`.
        """
        return self.gene_sets[key].gene_support()

    def save(self, target_dir, keys=None):
        """
        Save  this collection to `target_dir`. Each `GeneSet` will be saved as a separate
        .netcdf file within this directory.

        :param target_dir: The path to which the 'GeneSet' `xarray.Dataset` .netcdf files will be written.

        :param keys: The list of `GeneSet` keys that should be saved. If this is not provided, all
            `GeneSet` objects are saved.

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

    def gene_sets_to_dataframes(self, keys=None, only_supported: bool = True) -> dict:
        keys = self.gene_sets.keys() if keys is None else keys
        return {key: self.gene_sets[key].to_dataframe(only_supported) for key in keys}

    # TODO: Consider overwrite protection.
    def gene_sets_to_csv_files(self, target_dir=None, keys=None, only_supported: bool = True):
        keys = self.gene_sets.keys() if keys is None else keys
        target_dir = self.name if target_dir is None else target_dir
        data_frames = self.gene_sets_to_dataframes(keys, only_supported)
        os.makedirs(target_dir, exist_ok=True)
        for name, df in data_frames.items():
            df.to_csv(f"{target_dir}/{name}.csv")

    # TODO: Consider overwrite protection.
    def gene_sets_to_excel_sheet(self, name: str = None, keys=None, only_supported: bool = True):
        keys = self.gene_sets.keys() if keys is None else keys
        name = f'{self.name}.xlsx' if name is None else f'{name}.xlsx'
        data_frames = self.gene_sets_to_dataframes(keys, only_supported)
        with pd.ExcelWriter(name) as writer:
            for set_name, df in data_frames.items():
                df.to_excel(writer, sheet_name=set_name)

    def as_dict(self, keys=None, exclude=None):
        """
        Returns a dictionary of {name: supported_genes} for each gene set, or those specified
        by the `keys` argument.

        :param keys: The list of `GeneSet` keys to be included in the returned dictionary.

        :param exclude: A list of `GeneSet` keys to exclude from the returned dictionary.
        """
        if keys is None:
            keys = self.gene_sets.keys()

        keys = [key for key in keys if "support" in self.gene_sets[key].data]
        if exclude is not None:
            keys = [key for key in keys if key not in exclude]
        gene_sets = {k: self.gene_sets[k] for k in keys}
        return {k: v.gene_support() for k, v in gene_sets.items()}

    def intersection(self, keys=None, exclude=None):
        """
        Get the intersection of supported genes in this GeneSet collection.

        :param keys: The list of `GeneSet` keys to be included in the returned dictionary.

        :param exclude: A list of `GeneSet` keys to exclude from the returned dictionary.
        """
        gene_set_dict = self.as_dict(keys, exclude)
        return reduce(np.intersect1d, gene_set_dict.values())

    def union(self, keys=None, exclude=None):
        """
        Get the union of supported genes in this GeneSet collection.

        :param keys: The list of `GeneSet` keys to be included in the returned dictionary.

        :param exclude: A list of `GeneSet` keys to exclude from the returned dictionary.
        """
        gene_set_dict = self.as_dict(keys, exclude)
        return reduce(np.union1d, gene_set_dict.values())

    def difference(self, keys=None, exclude=None):
        """
        Get the difference of supported genes in this GeneSet collection.

        :param keys: The list of `GeneSet` keys to be included in the returned dictionary.

        :param exclude: A list of `GeneSet` keys to exclude from the returned dictionary.
        """
        gene_set_dict = self.as_dict(keys, exclude)
        return reduce(np.setdiff1d, gene_set_dict.values())

    def pairwise_unions(self, keys=None, exclude=None, size=2):
        gene_set_dict = self.as_dict(keys, exclude)
        return {(ak, bk): np.union1d(av, bv)
                for (ak, av), (bk, bv) in itertools.permutations(gene_set_dict.items(), size)
                if ak != bk}

    def pairwise_intersection(self, keys=None):
        gene_set_dict = self.as_dict(keys)
        return {(ak, bk): np.intersect1d(av, bv)
                for (ak, av), (bk, bv) in itertools.combinations(gene_set_dict.items(), 2)
                if ak != bk}

    def pairwise_percent_intersection(self, keys=None):
        """Get the normalized intersection length of each facet combination."""
        gene_set_dict = self.as_dict(keys)
        zero_filtered_dict = {k: v for k, v in gene_set_dict.items()
                              if len(v) > 0}
        return [(ak, bk, len(np.intersect1d(av, bv)) / len(av))
                for (ak, av), (bk, bv) in itertools.permutations(zero_filtered_dict.items(), 2)
                if ak != bk]

    def get_gene_sets_data(self, variables: list, keys=None):
        """
        Extracts variables from each Dataset, and returns them as a dictionary.

        If you need these to be collated into a dataset, use ...

        :param variables:
        :param keys:
        :return:
        """
        keys = self.gene_sets.keys() if keys is None else keys
        return {key: self.gene_sets[key].gene_support()[variables] for key in keys}

    def collate_gene_sets_data(self, variables: list, keys=None) -> xr.Dataset:
        keys = self.gene_sets.keys() if keys is None else keys
        data_dict = self.get_gene_sets_data(variables, keys)
        key_ds = xr.concat(data_dict.values(), dim="gene_set")
        key_ds["gene_set"] = keys
        return key_ds

    @classmethod
    def from_folder(cls, gem, target_dir, glob_filter="*.nc", filter_func=None, **params):
        """Create a `CompoundFacet` from a list of file paths. The base file names
        will be used as the key values.
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
