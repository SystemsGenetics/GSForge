import param
import xarray as xr
import numpy as np
import holoviews as hv
from ..models import OperationInterface
from holoviews.operation.stats import univariate_kde
from holoviews.operation.datashader import datashade, dynspread
import warnings

# hv.extension("bokeh", "matplotlib")


__all__ = [
    "ScatterDistributionBase",
]

# Define default options here so they can be used in several places.
default_hv_opts = {
    "points": dict(width=500, height=500, bgcolor="lightgrey", size=1.2, muted_alpha=0.05, show_grid=True),
    "dists": dict(bgcolor="lightgrey", show_grid=True, show_legend=False, alpha=0.25),
    "x_dist": dict(width=150),
    "y_dist": dict(height=150),
    "shaded": dict(width=500, height=500, bgcolor="lightgrey", show_grid=True),
}


class ScatterDistributionBase(OperationInterface):

    # transform = param.Callable(default=None)

    hv_point_options = param.Dict(default=default_hv_opts["points"])
    hv_dist_opts = param.Dict(default=default_hv_opts["dists"])
    hv_x_dist_opts = param.Dict(default=default_hv_opts["x_dist"])
    hv_y_dist_opts = param.Dict(default=default_hv_opts["y_dist"])

    datashade_ = param.Boolean(default=False)
    dynspread_ = param.Boolean(default=False)

    x_axis_selector = param.ObjectSelector(default="mean")
    y_axis_selector = param.ObjectSelector(default="variance")

    mappings = param.Parameter()

    backend = param.ObjectSelector(default="bokeh", objects=["bokeh", "matplotlib"])

    # scaler = param.ObjectSelector(default=np.log2)

    axis_functions = {
        "frequency": lambda counts, dim: (counts > 0).sum(dim=dim) / counts[dim].shape,
        "mean": lambda counts, dim: counts.mean(dim=dim),
        "variance": lambda counts, dim: counts.var(dim=dim),
    }

    @staticmethod
    def hv_scatter_dist(dataset, x_kdims, y_kdims,
                        points_opts=None,
                        dist_x_opts=None,
                        dist_y_opts=None,
                        datashade_=False,
                        dynspread_=False,
                        shaded_opts=None):

        if dynspread_ and not datashade_:
            warnings.warn("Dynspread can only be used with datashade, setting both to true.")
            datashade_= True

        df = dataset[[x_kdims, y_kdims]].to_dataframe()
        points = hv.Points(df, kdims=[x_kdims, y_kdims]).opts(**points_opts)

        dist_x = univariate_kde(hv.Distribution(points, kdims=[y_kdims]), n_samples=1000)
        dist_y = univariate_kde(hv.Distribution(points, kdims=[x_kdims]), n_samples=1000)

        if datashade_:
            points = datashade(points).opts(**shaded_opts)
            if dynspread_:
                points = dynspread(points)

        return points << dist_x.opts(**dist_x_opts) << dist_y.opts(**dist_y_opts)

    @staticmethod
    def scatter_dist_by_mappings():
        pass

    def process(self):
        counts = self.x_data
        dataset = xr.Dataset({
            self.x_axis_selector: self.axis_functions[self.x_axis_selector](counts, self.sample_index_name),
            self.y_axis_selector: self.axis_functions[self.y_axis_selector](counts, self.sample_index_name),
        })

        plot = self.hv_scatter_dist(
            dataset=dataset,
            x_kdims=self.x_axis_selector,
            y_kdims=self.y_axis_selector,
            points_opts=self.hv_point_options,
            dist_x_opts={**self.hv_dist_opts, **self.hv_x_dist_opts},
            dist_y_opts={**self.hv_dist_opts, **self.hv_y_dist_opts},
            datashade_=self.datashade_,
            dynspread_=self.dynspread_,
        )

        return plot
