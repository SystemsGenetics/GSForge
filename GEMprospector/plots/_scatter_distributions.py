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
# bokeh_options = hv.opts(
#     hv.opts.
# )
default_hv_opts = {
    "points": dict(width=500, height=500, bgcolor="lightgrey", size=1.2, muted_alpha=0.05, show_grid=True),
    "dists": dict(bgcolor="lightgrey", show_grid=True, show_legend=False, alpha=0.25),
    "x_dist": dict(width=150),
    "y_dist": dict(height=150),
    "shaded": dict(width=500, height=500, bgcolor="lightgrey", show_grid=True),
}


class ScatterDistributionBase(OperationInterface):
    """


    TODO: Remake

    + create dataset / aggregate function.
    + scatter function
    + dist function
    + scatter + dist function
    + option application function (layout(?), mode)

    process:
    + creates dataset
    + creates scatter
    + creates dists
    + applies options
    + returns layout


    """
    # transform = param.Callable(default=None)

    hv_point_options = param.Dict(default=default_hv_opts["points"])
    hv_dist_opts = param.Dict(default=default_hv_opts["dists"])
    hv_x_dist_opts = param.Dict(default=default_hv_opts["x_dist"])
    hv_y_dist_opts = param.Dict(default=default_hv_opts["y_dist"])
    shaded_opts = param.Dict(default=default_hv_opts["shaded"])

    datashade_ = param.Boolean(default=False)
    dynspread_ = param.Boolean(default=False)
    
    log_scale = param.Boolean(default=True)

    x_axis_selector = param.ObjectSelector(default="mean")
    y_axis_selector = param.ObjectSelector(default="variance")

    mappings = param.Parameter()

    backend = param.ObjectSelector(default="bokeh", objects=["bokeh", "matplotlib"])

    axis_functions = {
        "frequency": lambda counts, dim: (counts > 0).sum(dim=dim) / counts[dim].shape,
        "mean": lambda counts, dim: counts.mean(dim=dim),
        "variance": lambda counts, dim: counts.var(dim=dim),
        "standard_dev": lambda counts, dim: counts.std(dim=dim),
        "fano": lambda counts, dim: counts.var(dim=dim) / counts.mean(dim=dim),
        "mean_rank": lambda counts, dim: counts.mean(dim=dim).rank(dim="Gene"),
        "cv_squared": lambda counts, dim: counts.var(dim=dim) / counts.mean(dim=dim)**2
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
            datashade_ = True

        df = dataset[[x_kdims, y_kdims]].to_dataframe()
        points = hv.Points(df, kdims=[x_kdims, y_kdims])#.opts(**points_opts)

        dist_x = univariate_kde(hv.Distribution(points, kdims=[y_kdims], group="dist_x"), n_samples=1000)#.opts(**dist_x_opts)
        dist_y = univariate_kde(hv.Distribution(points, kdims=[x_kdims], group="dist_y"), n_samples=1000)#.opts(**dist_y_opts)

        if datashade_:
            points = datashade(points)#.opts(**shaded_opts)
            if dynspread_:
                points = dynspread(points)

        return points << dist_x << dist_y

    @staticmethod
    def scatter_dist_by_mappings(dataset, x_kdims, y_kdims,
                                 mappings,
                                 selection_dim="Gene",
                                 dist_x_opts=None,
                                 dist_y_opts=None,
                                 points_opts=None,
                                 datashade_=False,
                                 dynspread_=False,
                                 shaded_opts=None,
                                 ):

        data_groups = {name: dataset.sel({selection_dim: genes}) for name, genes in mappings.items()}
        data_group_dfs = {k: v[[x_kdims, y_kdims]].to_dataframe() for k, v in data_groups.items()}

        points = {k: hv.Points(val, kdims=[x_kdims, y_kdims]) for k, val in data_group_dfs.items()}

        dist_x = {k: univariate_kde(hv.Distribution(p, kdims=[y_kdims], group="dist_x"), n_samples=1000)#.opts(**dist_x_opts)
                  for k, p in points.items()}
        dist_y = {k: univariate_kde(hv.Distribution(p, kdims=[x_kdims], group="dist_y"), n_samples=1000)#.opts(**dist_y_opts)
                  for k, p in points.items()}

        if datashade_:
            points_overlay = datashade(hv.NdOverlay(points))#.opts(shaded_opts)
            if dynspread_:
                points_overlay = dynspread(points_overlay)
        else:
            points_overlay = hv.NdOverlay(points)#.opts({"Points": points_opts})

        return points_overlay << hv.NdOverlay(dist_x) << hv.NdOverlay(dist_y)

    def process(self):
        counts = self.x_count_data
        dataset = xr.Dataset({
            self.x_axis_selector: self.axis_functions[self.x_axis_selector](counts, self.sample_index_name),
            self.y_axis_selector: self.axis_functions[self.y_axis_selector](counts, self.sample_index_name),
        })
        
        if self.log_scale:
            dataset = np.log2(dataset.where(dataset > 0))

        if self.lineament_collection is not None and self.lineament_keys is not None:
            print("by mappings")
            mappings = self.lineament_collection.as_dict(self.lineament_keys)
            
#             for name, mapping in mappings.items():
#                 selected_gene_set = set(dataset[self.gene_index_name].values)
#                 print(len(selected_gene_set))
#                 active_genes = set.intersection(selected_gene_set, set(mappings))
#                 print(active_genes)
                
#                 mappings[name] = np.array(list(active_genes))
                
            plot = self.scatter_dist_by_mappings(
                dataset=dataset,
                x_kdims=self.x_axis_selector,
                y_kdims=self.y_axis_selector,
                mappings=mappings,
                points_opts=self.hv_point_options,
                dist_x_opts={**self.hv_dist_opts, **self.hv_x_dist_opts},
                dist_y_opts={**self.hv_dist_opts, **self.hv_y_dist_opts},
                datashade_=self.datashade_,
                dynspread_=self.dynspread_,
                shaded_opts=self.shaded_opts,
            )

        else:
            plot = self.hv_scatter_dist(
                dataset=dataset,
                x_kdims=self.x_axis_selector,
                y_kdims=self.y_axis_selector,
                # points_opts=self.hv_point_options,
                # dist_x_opts={**self.hv_dist_opts, **self.hv_x_dist_opts},
                # dist_y_opts={**self.hv_dist_opts, **self.hv_y_dist_opts},
                datashade_=self.datashade_,
                dynspread_=self.dynspread_,
                # shaded_opts=self.shaded_opts,
            )

        return plot
