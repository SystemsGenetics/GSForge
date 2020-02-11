import param
import xarray as xr
import holoviews as hv
import colorcet as cc
from holoviews.operation.stats import univariate_kde
import itertools

from ...models import Interface
from ..utils import AbstractPlottingOperation


class sample_wise_distributions(Interface, AbstractPlottingOperation):
    # TODO: Document me.
    # TODO: Add matplotlib options.

    @staticmethod
    def bokeh_opts():
        return hv.opts.Area(fill_alpha=0.0, width=600, height=500)

    # @staticmethod
    # def matplotlib_opts():
    #     return

    @staticmethod
    def sample_wise_count_distributions(counts: xr.DataArray, sample_dim: str = "Sample",
                                        cmap: list = cc.glasbey) -> hv.Overlay:
        distributions = [
            univariate_kde(hv.Distribution(counts.sel(**{sample_dim: sample}).values)).opts(line_color=color)
            for sample, color in zip(counts[sample_dim].values, itertools.cycle(cmap))]
        return hv.Overlay(distributions).opts(show_legend=False)

    def process(self):
        layout = self.sample_wise_count_distributions(counts=self.x_count_data,
                                                      sample_dim=self.sample_index_name)
        return layout
