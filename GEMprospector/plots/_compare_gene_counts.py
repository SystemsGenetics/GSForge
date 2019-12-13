import holoviews as hv
import numpy as np
import pandas as pd
import functools
import param

from ..models import OperationInterface


class GenesVsCounts(OperationInterface):
    hue = param.Parameter(default=None)
    mode = param.ObjectSelector(default="scatter", objects=["scatter"])

    @staticmethod
    def genewise_scatter(data, hue=None,
                         gene_dim="Gene", sample_dim="Sample", count_dim="counts"):

        hvds = hv.Dataset(data.to_dataframe().reset_index())

        kdims = [gene_dim, count_dim]
        vdims = [sample_dim]

        if hue:
            vdims += [hue]

        scatter = hv.Scatter(hvds, kdims=kdims, vdims=vdims)

        if hue is not None:
            overlays = {key: scatter.select(**{hue: key})
                        for key in scatter.data[hue].unique()}
            return hv.NdOverlay(overlays)

        return scatter

    # @staticmethod
    # def genewise_violin(data, hue=None,
    #                     gene_dim="Gene", sample_dim="Sample", count_dim="counts"):
    #     hvds = hv.Dataset(data.to_dataframe().reset_index())
    #
    #     kdims = [gene_dim]
    #     vdims = [count_dim]
    #
    #     if hue:
    #         kdims += [hue]
    #
    #     violin = hv.Violin(hvds, kdims=kdims, vdims=vdims)
    #
    #     if hue is not None:
    #         overlays = {key: violin.select(**{hue: key})
    #                     for key in violin.data[hue].unique()}
    #         return hv.NdOverlay(overlays)
    #
    #     return violin

    def process(self):

        self.set_param(gene_set_mode="intersection")

        if self.gene_set_collection is not None:
            self.set_param(selected_gene_sets=list(self.gene_set_collection.gene_sets.keys()))

        selected_subset = self.selection

        modes = {
            "scatter": self.genewise_scatter,
            # "violin": self.genewise_violin,
        }

        return modes[self.mode](selected_subset, self.hue)
