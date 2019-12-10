import pathlib
import param
import os
import itertools
import xarray as xr
import numpy as np

from textwrap import dedent, indent

from ._AnnotatedGEM import AnnotatedGEM
from ._GeneSet import GeneSet


# TODO:
#     + Create default GeneSets on initialization (These should
#       probably not be saved). "total" and "zero_dropped".
#       Consider a `basic_GeneSets` function.
#     + Add a `get_support(key)` function.
class GeneSetCollection(param.Parameterized):
    """
    Contains both an `AnnotatedGEM` and a dictionary of `GeneSet` objects, as well
    as functions for comparing and analyzing those objects.
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
        summary = [f"<GEMprospector.{type(self).__name__}>"]
        summary += [self.name]
        summary += [indent(self.gem.__repr__(), "    ")]
        summary += ["GeneSet Keys and # of Selected Genes"]
        summary += [f"    {k}: {v}" for k, v in self._summarize_gene_sets().items()]
        return "\n".join(summary)

    def get_support(self, key) -> np.array:
        """Get the support array for a given key.

        :param key: The GeneSet from which to get the gene support.

        :return: A `numpy` array of the genes that make up the support of this `GeneSet`.
        """
        return self.gene_sets[key].gene_support()

    def save(self, target_dir, keys=None):
        """Save  this collection to `target_dir`. Each `GeneSet` will be saved as a separate
        .netcdf file within this directory.

        :param target_dir: The path to which the 'GeneSet' `xarray.Dataset` .netcdf files will be written.

        :param keys: The list of `GeneSet` keys that should be saved. If this is not provided, all
            `GeneSet` objects are saved.

        :return:
        """
        # Save all the gene sets in this collection in the target_dir path.
        if keys is None:
            keys = self.gene_sets.keys()

        # Create any needed intermediate directories.
        if not os.path.exists(target_dir):
            os.makedirs(target_dir)

        for key in keys:
            save_path = self.gene_sets[key].save_as_netcdf(target_dir)
            # yield save_path

    def as_dict(self, keys=None, exclude=None):
        """Returns a dictionary of {name: supported_genes} for each gene set, or those specified
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
        """Get the intersection of supported genes in this GeneSet collection."""
        gene_set_dict = self.as_dict(keys, exclude)
        return set.intersection(*[set(x) for x in gene_set_dict.values()])

    def union(self, keys=None, exclude=None):
        """Get the union of supported genes in this GeneSet collection."""
        gene_set_dict = self.as_dict(keys, exclude)
        return set(itertools.chain.from_iterable(gene_set_dict.values()))

    def difference(self, keys=None, exclude=None):
        gene_set_dict = self.as_dict(keys, exclude)
        return set.difference(*[set(x) for x in gene_set_dict.values()])

    def pairwise_unions(self, keys=None, exclude=None, size=2):
        gene_set_dict = self.as_dict(keys, exclude)
        return {(ak, bk): set.union(set(av), set(bv))
                for (ak, av), (bk, bv) in itertools.permutations(gene_set_dict.items(), size)
                if ak != bk}

    def pairwise_intersection(self, keys=None):
        gene_set_dict = self.as_dict(keys)
        return {(ak, bk): set.intersection(set(av), set(bv))
                for (ak, av), (bk, bv) in itertools.combinations(gene_set_dict.items(), 2)
                if ak != bk}

    def pairwise_percent_intersection(self, keys=None):
        """Get the normalized intersection length of each facet combination."""
        gene_set_dict = self.as_dict(keys)
        zero_filtered_dict = {k: v for k, v in gene_set_dict.items()
                              if len(v) > 0}
        return [(ak, bk, len(set.intersection(set(av), set(bv))) / len(set(av)))
                for (ak, av), (bk, bv) in itertools.permutations(zero_filtered_dict.items(), 2)
                if ak != bk]

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
                name = key = os.path.basename(str(file)).rsplit(".nc")[0]
            else:
                name = new_gene_set.name

            gene_sets[name] = new_gene_set

        return cls(gem=gem, gene_sets=gene_sets, **params)

# # Python 3.8 will let us move this code into the class body, and add
# # add register the @classmethod functions.
# @functools.singledispatch
# def _lineament_collection_dispatch(source, *args, **params):
#     raise TypeError(f"Source of type: {type(source)} not supported.")
#
#
# _lineament_collection_dispatch.register(str, GeneSetCollection.from_folder)
# _lineament_collection_dispatch.register(pd.DataFrame, GeneSetCollection.from_pandas)
