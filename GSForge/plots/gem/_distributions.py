import param
import xarray as xr
import numpy as np
import holoviews as hv
import colorcet as cc
from holoviews.operation.stats import univariate_kde
import itertools
from textwrap import dedent

from GSForge.models import OperationInterface


class sample_wise_distributions(OperationInterface):

    backend = param.ObjectSelector(default="bokeh", objects=["bokeh", "matplotlib"], doc=dedent("""\
    The selected plotting backend to use for display. Options are ["bokeh", "matplotlib"]."""))

    apply_default_opts = param.Boolean(default=True, precedence=-1.0, doc=dedent("""\
    Whether to apply the default styling based on the current backend."""))

    @staticmethod
    def bokeh_options():
        return hv.opts.Area(fill_alpha=0.0, width=600, height=500)

    @staticmethod
    def matplotlib_options():
        return

    @staticmethod
    def sample_wise_count_distributions(counts: xr.DataArray, sample_dim: str = "Sample",
                                        cmap: list = cc.glasbey) -> hv.Overlay:
        distributions = [
            univariate_kde(hv.Distribution(counts.sel(**{sample_dim: sample}).values)).opts(line_color=color)
            for sample, color in zip(counts[sample_dim].values, itertools.cycle(cmap))]
        return hv.Overlay(distributions).opts(show_legend=False)

    def process(self):

        if self.backend == "matplotlib":
            raise NotImplemented("Matplotlib backend not currently working for this.")

        layout = self.sample_wise_count_distributions(
            counts=self.x_count_data, sample_dim=self.sample_index_name)
        options = {"bokeh": self.bokeh_options, "matplotlib": self.matplotlib_options}
        if self.apply_default_opts:
            default_options = options[self.backend]()
            return layout.opts(default_options)

        return layout
