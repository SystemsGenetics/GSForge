import holoviews as hv
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
        vdims = [sample_dim] if hue is None else [sample_dim, hue]

        scatter = hv.Scatter(hvds, kdims=kdims, vdims=vdims)

        if hue is not None:
            overlays = {key: scatter.select(**{hue: key}) for key in scatter.data[hue].unique()}
            return hv.NdOverlay(overlays)

        return scatter

    def process(self):
        modes = {"scatter": self.genewise_scatter,
                 # "violin": self.genewise_violin,
                 }
        self.set_param(annotation_variables=[self.hue])
        return modes[self.mode](self.selection, self.hue)
