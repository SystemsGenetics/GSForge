import holoviews as hv
import param
import pandas as pd
from textwrap import dedent

from ...models import Interface
from ..utils import AbstractPlottingOperation


class GeneVsCountsScatter(Interface, AbstractPlottingOperation):
    """
    Display the counts of a small selection of genes on a scatter plot (genes vs counts).

    This function is a `GSForge.OperationInterface` method, and shares those common parameters.

    A selection of genes must be supplied to this function, either specifically via the
    `selected_genes` parameter, or implicitly by the selection or combination of some
    GeneSet support arrays via a GeneSetCollection.

    A warning is generated if the number of genes surpasses a 'soft_max' parameter.
    This is only there to prevent useless, accidental plots of many thousands of genes.
    You would need to override it to show more than 50 genes by default.

    Parameters specific to this function:

    :param hue:
        Color by which to shade the observations.

    :param soft_max:
        The number of genes above which this function will return a ValueError rather than
        attempt to plot an unreasonable number of genes.

    :param backend:
        The selected plotting backend to use for display. Options are ["bokeh", "matplotlib"].

    :param apply_default_opts:
        Whether to apply the default styling.

    :returns:
        A `holoviews.Layout` object. The display of this object depends on the currently
        selected backend.

    """
    hue = param.String(default=None, doc="Color by which to shade the observations.")

    soft_max = param.Integer(default=250, precedence=-1.0, doc=dedent("""\
    The number of genes above which this function will return a ValueError rather than attempt to 
    plot an unreasonable number of genes."""))

    @staticmethod
    def genewise_scatter(data: pd.DataFrame, hue: str = None, gene_dim: str = "Gene", sample_dim: str = "Sample",
                         count_dim: str = "counts"):
        hvds = hv.Dataset(data)
        kdims = [gene_dim, count_dim]
        vdims = [sample_dim] if hue is None else [sample_dim, hue]

        scatter = hv.Scatter(hvds, kdims=kdims, vdims=vdims)

        if hue is not None:
            overlays = {key: scatter.select(**{hue: key}) for key in scatter.data[hue].unique()}
            return hv.NdOverlay(overlays)

        return scatter

    @staticmethod
    def bokeh_opts():
        return hv.opts.Scatter(jitter=0.2, width=800, height=500, legend_position="right", xrotation=90, padding=0.1,
                               backend="bokeh")

    @staticmethod
    def matplotlib_opts():
        return hv.opts.Scatter(fig_size=250, aspect=1.8, xrotation=90, padding=0.1, backend="matplotlib")

    def process(self):
        self.set_param(annotation_variables=[self.hue])
        if self.get_gene_index().shape[0] > self.soft_max:
            genes_selected = self.get_gene_index().shape[0]
            raise ValueError(f"Selected number of genes: {genes_selected} is likely too much."
                             f"Provide an array of genes to `selected_genes` less than (or override) "
                             f"the `soft_max` parameter of {self.soft_max}")

        layout = self.genewise_scatter(data=self.selection.to_dataframe().reset_index(),
                                       hue=self.hue,
                                       gene_dim=self.gene_index_name,
                                       sample_dim=self.sample_index_name,
                                       count_dim=self.count_variable)
        return layout
