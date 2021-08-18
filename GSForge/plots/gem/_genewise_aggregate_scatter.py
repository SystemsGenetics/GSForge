import param
import xarray as xr
import numpy as np
import holoviews as hv

from holoviews.operation.stats import univariate_kde
from holoviews.operation import datashader
import warnings

from ..abstract_plot_models import InterfacePlottingBase


# TODO: Display for the user information on how the plot was created.
# TODO: Add n_samples for univariate kde function.
class GenewiseAggregateScatter(InterfacePlottingBase):
    """
    Displays the output of selected aggregations upon the count array on a scatter plot with optional
    adjoined kernel density estimates. e.g. mean counts vs mean variance etc. By default such outputs are
    log2 transformed as well.

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
        "frequency": lambda counts: (counts > 0).sum(axis=0) / counts.shape[1],
        "mean": lambda counts: counts.mean(axis=0),
        "variance": lambda counts: counts.var(axis=0),
        "standard_dev": lambda counts: counts.std(axis=0),
        "fano": lambda counts: counts.var(axis=0) / counts.mean(axis=0),
        "mean_rank": lambda counts: np.rank(counts.mean(axis=0)),
        "cv_squared": lambda counts: (counts.std(axis=0) / counts.mean(axis=0)) ** 2
    }

    datashade = param.Boolean(default=False)
    dynspread = param.Boolean(default=False)
    
    adjoint_distributions = param.Boolean(default=True)

    x_axis_selector = param.ObjectSelector(default="mean", doc="Select from the available axis aggregation functions.",
                                           objects=axis_functions.keys())

    y_axis_selector = param.ObjectSelector(default="variance",
                                           doc="Select from the available axis aggregation functions.",
                                           objects=axis_functions.keys())

    # axis_transform = param.Parameter(default=('log 2', lambda ds: np.log2(ds.where(ds > 0))),
    #                                  doc="A transform (usually log2) for getting a viewable spread of the results.")
    # axis_transform_name = param.String(default='log 2')

    @staticmethod
    def bokeh_opts():
        return [
            hv.opts.Points(width=500, height=500, bgcolor="lightgrey", size=1.2, muted_alpha=0.05, show_grid=True),
            hv.opts.Area(bgcolor="lightgrey", show_grid=True, show_legend=False, alpha=0.25),
            hv.opts.Area("dist_x", width=100),
            hv.opts.Area("dist_y", height=100),
            hv.opts.RGB(width=500, height=500, bgcolor="lightgrey", show_grid=True),
        ]

    @staticmethod
    def matplotlib_opts():
        return [
            hv.opts.Points(aspect=1.0, bgcolor="lightgrey", s=1.2, show_grid=True),
            hv.opts.Area(bgcolor="lightgrey", show_grid=True, show_legend=False, alpha=0.25),
            hv.opts.RGB(bgcolor="lightgrey", show_grid=True, aspect=1.0),
        ]
            
    @staticmethod
    def scatter_dist(dataset, x_kdims, y_kdims,
                     datashade=False,
                     dynspread=False,
                     adjoint_distributions=True,
                     options=None):

        if dynspread and not datashade:
            warnings.warn("Dynspread can only be used with datashade, setting both to true.")
            datashade = True

        df = dataset[[x_kdims, y_kdims]].to_dataframe()
        points = hv.Points(df, kdims=[x_kdims, y_kdims])

        # The distributions must be completed before the datashader operation occurs.
        if adjoint_distributions:
            dist_x = univariate_kde(hv.Distribution(points, kdims=[y_kdims], group="dist_x"), n_samples=1000)
            dist_y = univariate_kde(hv.Distribution(points, kdims=[x_kdims], group="dist_y"), n_samples=1000)

        if datashade:
            points = datashader.datashade(points)
            if dynspread:
                points = datashader.dynspread(points)

        layout = points
        # This must be done here to prevent the options getting lost upon creation of the
        # holoviews.AdjointLayout object.
        if options:
            points = points.opts(options)

            if adjoint_distributions:
                dist_x = dist_x.opts(options)
                dist_y = dist_y.opts(options)
                layout = points << dist_x << dist_y

        return layout

    @staticmethod
    def scatter_dist_by_mappings(dataset, x_kdims, y_kdims,
                                 mappings,
                                 selection_dim="Gene",
                                 datashade=False,
                                 dynspread=False,
                                 adjoint_distributions=True,
                                 options=None,
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

        if options:
            points_overlay = points_overlay.opts(options)

            if adjoint_distributions:
                dist_x = hv.NdOverlay(dist_x).opts(options)
                dist_y = hv.NdOverlay(dist_y).opts(options)
                points_overlay = points_overlay << dist_x << dist_y

        return points_overlay

        # return points_overlay << hv.NdOverlay(dist_x) << hv.NdOverlay(dist_y)

    def __call__(self):
        counts = self.x_count_data

        x_axis_name = f'{self.x_axis_selector}'
        y_axis_name = f'{self.y_axis_selector}'

        data = xr.Dataset({
            x_axis_name: self.axis_functions[self.x_axis_selector](counts),
            y_axis_name: self.axis_functions[self.y_axis_selector](counts),
        })
        
        options = None
        if self.apply_default_opts is True:
            options = self.get_default_options()

        if self.selected_gene_sets == [None] or not self.gene_set_collection:

            layout = self.scatter_dist(
                dataset=data,
                x_kdims=x_axis_name,
                y_kdims=y_axis_name,
                datashade=self.datashade,
                dynspread=self.dynspread,
                adjoint_distributions=self.adjoint_distributions,
                options=options,
            )

        else:
            mappings = self.gene_set_collection.as_dict(self.selected_gene_sets)

            layout = self.scatter_dist_by_mappings(
                dataset=data,
                x_kdims=x_axis_name,
                y_kdims=y_axis_name,
                mappings=mappings,
                datashade=self.datashade,
                dynspread=self.dynspread,
                adjoint_distributions=self.adjoint_distributions,
                options=options,
            )
        if self.apply_default_opts is True:
            layout = layout.opts(self.get_default_options())
        return layout

