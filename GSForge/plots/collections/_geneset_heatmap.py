import holoviews as hv
import itertools
import param

from ...models import GeneSetCollection, Interface
from ..utils import AbstractPlottingOperation


class WithinCollectionOverlapHeatMap(Interface, AbstractPlottingOperation):
    mode = param.ObjectSelector(default="overlap", objects=["overlap", "percent"])

    @staticmethod
    def bokeh_opts():
        return hv.opts.HeatMap(xrotation=45, width=450, height=450, labelled=[], colorbar=True)

    @staticmethod
    def matplotlib_opts():
        return hv.opts.HeatMap(xrotation=45, fig_size=250, labelled=[], colorbar=True)

    @staticmethod
    def within_collection_overlap(gene_dict, mode="overlap"):
        modes = {"overlap": lambda va, vb: len(set.intersection(set(va), set(vb))),
                 "percent": lambda va, vb: len(set.intersection(set(va), set(vb))) / len(set(va))}

        mode_formaters = {"overlap": lambda x: f"{x}",
                          "percent": lambda x: f"{x:.0%}"}

        if mode not in modes.keys():
            raise ValueError(f"{mode} is not a valid mode. Select from {list(modes.keys())}.)")

        mode_func = modes[mode]
        overlap_dim = hv.Dimension(f"Overlap {mode}", value_format=mode_formaters[mode])
        data = [(ak, bk, mode_func(av, bv))
                for (ak, av), (bk, bv) in itertools.permutations(gene_dict.items(), 2)
                if ak != bk]
        heatmap = hv.HeatMap(data, vdims=overlap_dim)
        return heatmap * hv.Labels(heatmap)

    def process(self):
        gene_dict = self.gene_set_collection.as_dict(self.selected_gene_sets)
        return self.within_collection_overlap(gene_dict, mode=self.mode)


# TODO: Re-implement below.
# class BetweenCollectionOverlap(AbstractPlottingOperation):
#
#     @staticmethod
#     def between_collection_overlap(alpha: dict, beta: dict):
#         """
#         View the overlap (as a heatmap) of two GeneSetCollection dictionaries,
#         as provided by gsc.as_dict().
#
#         :param alpha:
#         :param beta:
#         :return:
#         """
#         overlap_dim = hv.Dimension("Overlap Count", value_format=lambda x: f"{x}")
#         data = [(f"{ak}: {len(av)}", f"{bk}: {len(bv)}", len(set.intersection(set(av), set(bv))))
#                 for (ak, av), (bk, bv) in itertools.product(alpha.items(), beta.items())]
#         heatmap = hv.HeatMap(data, vdims=overlap_dim)
#         return heatmap * hv.Labels(heatmap)
#
#     def process(self):
#         gene_dict = self.gene_set_collection.as_dict(self.selected_gene_sets)
#         return self.within_collection_overlap(gene_dict, mode=self.mode)

# def dictionary_pair_overlap(alpha, beta):
#     overlap_dim = hv.Dimension("Overlap Count", value_format=lambda x: f"{x}")
#
#     data = [(f"{ak}: {len(av)}", f"{bk}: {len(bv)}", len(set.intersection(set(av), set(bv))))
#             for (ak, av), (bk, bv) in itertools.product(alpha.items(), beta.items())]
#
#     heatmap = hv.HeatMap(data, vdims=overlap_dim).opts(labelled=[], colorbar=True)
#     return heatmap * hv.Labels(heatmap)
