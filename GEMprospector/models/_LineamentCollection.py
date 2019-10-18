import pathlib
import param
import os
import itertools
import xarray as xr
import numpy as np

from textwrap import dedent, indent

from ._AnnotatedGEM import AnnotatedGEM
from ._Lineament import Lineament


# TODO:
#     + Create default lineaments on initialization (These should
#       probably not be saved). "total" and "zero_dropped".
#       Consider a `basic_lineaments` function.
#     + Add a `get_support(key)` function.
class LineamentCollection(param.Parameterized):
    """
    Contains both an `AnnotatedGEM` and a dictionary of `Lineament` objects, as well
    as functions for comparing and analyzing those objects.
    """

    gem = param.ClassSelector(class_=AnnotatedGEM, doc=dedent("""\
    A Gene Expression Matrix (GEM) object."""))

    lineaments = param.Dict(doc=dedent("""\
    A dictionary of {key: xarray.DataArray}, boolean arrays indicating
    support for a given gene."""))

    def __init__(self, **params):
        super().__init__(**params)
        if self.lineaments is None:
            self.lineaments = dict()

    def _summarize_lineaments(self):
        """Summarize this LineamentCollection."""
        lin_counts = {key: len(lin.gene_support()) for key, lin in self.lineaments.items()}
        lin_counts = {key: lin_counts[key] for key in sorted(lin_counts, key=lin_counts.get, reverse=True)}
        return lin_counts

    def __repr__(self):
        summary = [f"<GeneSelector.{type(self).__name__}>"]
        summary += [indent(self.gem.__repr__(), "    ")]
        summary += ["Lineament Keys and # of Selected Genes"]
        summary += [f"    {k}: {v}" for k, v in self._summarize_lineaments().items()]
        return "\n".join(summary)

    def get_support(self, key) -> np.array:
        """Get the support array for a given key.

        :param key: The lineament from which to get the gene support.

        :return: A `numpy` array of the genes that make up the support of this `Lineament`.
        """
        return self.lineaments[key].gene_support()

    def save(self, target_dir, keys=None):
        """Save  this collection to `target_dir`. Each `Lineament` will be saved as a separate
        .netcdf file within this directory.

        :param target_dir: The path to which the 'Lineament' `xarray.Dataset` .netcdf files will be written.

        :param keys: The list of `Lineament` keys that should be saved. If this is not provided, all
            `Lineament` objects are saved.

        :return:
        """
        # Save all the lineaments in this collection in the target_dir path.
        if keys is None:
            keys = self.lineaments.keys()

        # Create any needed intermediate directories.
        if not os.path.exists(target_dir):
            os.makedirs(target_dir)

        for key in keys:
            save_path = self.lineaments[key].save_as_netcdf(target_dir)
            print(save_path)
            # xds = self.lineaments[key].data
            # xds.to_netcdf(os.path.join(target_dir, key + ".nc"))

    def as_dict(self, keys=None, exclude=None):
        if keys is None:
            keys = self.lineaments.keys()

        keys = [key for key in keys if "support" in self.lineaments[key].data]
        if exclude is not None:
            keys = [key for key in keys if key not in exclude]
        lineaments = {k: self.lineaments[k] for k in keys}
        return {k: v.gene_support() for k, v in lineaments.items()}

    def intersection(self, keys=None, exclude=None):
        """Get the intersection of supported genes in this lineament collection."""
        lineament_dict = self.as_dict(keys, exclude)
        return set.intersection(*[set(x) for x in lineament_dict.values()])

    def union(self, keys=None, exclude=None):
        """Get the union of supported genes in this lineament collection."""
        lineament_dict = self.as_dict(keys, exclude)
        return set(itertools.chain.from_iterable(lineament_dict.values()))

    def difference(self, keys=None, exclude=None):
        lineament_dict = self.as_dict(keys, exclude)
        return set.difference(*[set(x) for x in lineament_dict.values()])

    def union_combinations(self, keys=None, exclude=None, size=2):
        lineament_dict = self.as_dict(keys, exclude)
        return {(ak, bk): set.union(set(av), set(bv))
                for (ak, av), (bk, bv) in itertools.permutations(lineament_dict.items(), size)
                if ak != bk}

    def pairwise_intersection(self, keys=None):
        lineament_dict = self.as_dict(keys)
        return {(ak, bk): set.intersection(set(av), set(bv))
                for (ak, av), (bk, bv) in itertools.combinations(lineament_dict.items(), 2)
                if ak != bk}

    def pairwise_percent_intersection(self, keys=None):
        """Get the normalized intersection length of each facet combination."""
        lineament_dict = self.as_dict(keys)
        zero_filtered_dict = {k: v for k, v in lineament_dict.items()
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
        lineaments = dict()
        for file in pathlib.Path(target_dir).glob(glob_filter):
            # key = os.path.basename(file).rsplit(".nc")[0]
            data = xr.open_dataset(file)

            if filter_func is not None:
                if filter_func(data):
                    pass
                else:
                    continue

            data = xr.align(data, gem.gene_index, join="outer")[0]
            new_lineament = Lineament.from_xarray_dataset(data=data)
            lineaments[new_lineament.name] = new_lineament

        return cls(gem=gem, lineaments=lineaments, **params)

# # Python 3.8 will let us move this code into the class body, and add
# # add register the @classmethod functions.
# @functools.singledispatch
# def _lineament_collection_dispatch(source, *args, **params):
#     raise TypeError(f"Source of type: {type(source)} not supported.")
#
#
# _lineament_collection_dispatch.register(str, LineamentCollection.from_folder)
# _lineament_collection_dispatch.register(pd.DataFrame, LineamentCollection.from_pandas)
