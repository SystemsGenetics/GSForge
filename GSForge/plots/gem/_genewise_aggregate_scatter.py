import param
import xarray as xr
import numpy as np
import holoviews as hv

from holoviews.operation.stats import univariate_kde
from holoviews.operation import datashader
import warnings

from ...models import Interface
from ..utils import AbstractPlottingOperation


class GenewiseAggregateScatter(Interface, AbstractPlottingOperation):
    """
    Displays the output of selected aggregations upon the count array on a scatter plot with optional
    adjoined kernel density estimates. e.g. mean counts vs mean variance etc.

    This function is a `GSForge.OperationInterface` method, and shares those common parameters.

    Axis aggregation functions:
        + frequency
        + mean
        + variance
        + standard_dev
        + fano
        + mean_rank
        + cv_squared

    Parameters specific to this function:

    :param x_axis_selector:
        Select from the available axis aggregation functions.

    :param y_axis_selector:
        Select from the available axis aggregation functions.

    :param datashade:
        Whether to apply the datashader.datashade operation.

    :param dynspread:
        Whether to apply the datashader.dynspread operation.

    :param backend:
        The selected plotting backend to use for display. Options are ["bokeh", "matplotlib"].

    :param apply_default_opts:
        Whether to apply the default styling.

    :returns:
        A `holoviews.Layout` object. The display of this object depends on the currently
        selected backend.
    """
    axis_functions = {
        "frequency": lambda counts, dim: (counts > 0).sum(dim=dim) / counts[dim].shape,
        "mean": lambda counts, dim: counts.mean(dim=dim),
        "variance": lambda counts, dim: counts.var(dim=dim),
        "standard_dev": lambda counts, dim: counts.std(dim=dim),
        "fano": lambda counts, dim: counts.var(dim=dim) / counts.mean(dim=dim),
        "mean_rank": lambda counts, dim: counts.mean(dim=dim).rank(dim="Gene"),
        "cv_squared": lambda counts, dim: counts.var(dim=dim) / counts.mean(dim=dim) ** 2
    }

    datashade = param.Boolean(default=False)
    dynspread = param.Boolean(default=False)

    x_axis_selector = param.ObjectSelector(default="mean", doc="Select from the available axis aggregation functions.",
                                           objects=axis_functions.keys())

    y_axis_selector = param.ObjectSelector(default="variance",
                                           doc="Select from the available axis aggregation functions.",
                                           objects=axis_functions.keys())

    axis_transform = param.Parameter(default=lambda ds: np.log2(ds.where(ds > 0)),
                                     doc="A transform (usually log) for getting a viewable spread of the results.")

    @staticmethod
    def bokeh_opts():
        return [
            hv.opts.Points(width=500, height=500, bgcolor="lightgrey", size=1.2, muted_alpha=0.05,
                           show_grid=True),
            hv.opts.Area(bgcolor="lightgrey", show_grid=True, show_legend=False, alpha=0.25),
            hv.opts.Area("dist_x", width=150),
            hv.opts.Area("dist_y", height=150),
            hv.opts.RGB(width=500, height=500, bgcolor="lightgrey", show_grid=True),
        ]

    @staticmethod
    def matplotlib_opts():
        return [
            hv.opts.Points(fig_size=250, bgcolor="lightgrey", s=1.2, muted_alpha=0.05,
                           show_grid=True),
            hv.opts.Area(bgcolor="lightgrey", show_grid=True, show_legend=False, alpha=0.25),
            hv.opts.Area("dist_x", width=150),
            hv.opts.Area("dist_y", height=150),
            hv.opts.RGB(width=500, height=500, bgcolor="lightgrey", show_grid=True),
        ]

    @staticmethod
    def scatter_dist(dataset, x_kdims, y_kdims,
                     datashade=False,
                     dynspread=False):

        if dynspread and not datashade:
            warnings.warn("Dynspread can only be used with datashade, setting both to true.")
            datashade = True

        df = dataset[[x_kdims, y_kdims]].to_dataframe()
        points = hv.Points(df, kdims=[x_kdims, y_kdims])

        dist_x = univariate_kde(hv.Distribution(points, kdims=[y_kdims], group="dist_x"), n_samples=1000)
        dist_y = univariate_kde(hv.Distribution(points, kdims=[x_kdims], group="dist_y"), n_samples=1000)

        if datashade:
            points = datashader.datashade(points)
            if dynspread:
                points = datashader.dynspread(points)

        return points << dist_x << dist_y

    @staticmethod
    def scatter_dist_by_mappings(dataset, x_kdims, y_kdims,
                                 mappings,
                                 selection_dim="Gene",
                                 datashade=False,
                                 dynspread=False,
                                 ):

        data_groups = {name: dataset.sel({selection_dim: genes}) for name, genes in mappings.items()}
        data_group_dfs = {k: v[[x_kdims, y_kdims]].to_dataframe() for k, v in data_groups.items()}

        points = {k: hv.Points(val, kdims=[x_kdims, y_kdims]) for k, val in data_group_dfs.items()}

        dist_x = {k: univariate_kde(hv.Distribution(p, kdims=[y_kdims], group="dist_x"), n_samples=1000)
                  for k, p in points.items()}
        dist_y = {k: univariate_kde(hv.Distribution(p, kdims=[x_kdims], group="dist_y"), n_samples=1000)
                  for k, p in points.items()}

        if datashade:
            points_overlay = datashader.datashade(hv.NdOverlay(points))
            if dynspread:
                points_overlay = datashader.dynspread(points_overlay)
        else:
            points_overlay = hv.NdOverlay(points)

        return points_overlay << hv.NdOverlay(dist_x) << hv.NdOverlay(dist_y)

    def process(self):
        counts = self.x_count_data
        dataset = xr.Dataset({
            self.x_axis_selector: self.axis_functions[self.x_axis_selector](counts, self.sample_index_name),
            self.y_axis_selector: self.axis_functions[self.y_axis_selector](counts, self.sample_index_name),
        })

        # TODO:
        if self.axis_transform:
            dataset = self.axis_transform(dataset)

        if self.selected_gene_sets == [None] or not self.gene_set_collection:

            layout = self.scatter_dist(
                dataset=dataset,
                x_kdims=self.x_axis_selector,
                y_kdims=self.y_axis_selector,
                datashade=self.datashade,
                dynspread=self.dynspread,
            )

        else:
            mappings = self.gene_set_collection.as_dict(self.selected_gene_sets)

            layout = self.scatter_dist_by_mappings(
                dataset=dataset,
                x_kdims=self.x_axis_selector,
                y_kdims=self.y_axis_selector,
                mappings=mappings,
                datashade=self.datashade,
                dynspread=self.dynspread,
            )

        return layout
