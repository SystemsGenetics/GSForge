import param
import xarray as xr
import numpy as np
import holoviews as hv
from GSForge.models import OperationInterface

from holoviews.operation.stats import univariate_kde
from holoviews.operation.datashader import datashade, dynspread
import warnings

__all__ = [
    "ScatterDistributionBase",
]


class ScatterDistributionBase(OperationInterface):

    datashade_ = param.Boolean(default=False)
    dynspread_ = param.Boolean(default=False)

    log_scale = param.Boolean(default=True)

    x_axis_selector = param.ObjectSelector(default="mean")
    y_axis_selector = param.ObjectSelector(default="variance")

    axis_transform = param.Parameter(default=lambda ds: np.log2(ds.where(ds > 0)))

    backend = param.ObjectSelector(default="bokeh", objects=["bokeh", "matplotlib"])

    axis_functions = {
        "frequency": lambda counts, dim: (counts > 0).sum(dim=dim) / counts[dim].shape,
        "mean": lambda counts, dim: counts.mean(dim=dim),
        "variance": lambda counts, dim: counts.var(dim=dim),
        "standard_dev": lambda counts, dim: counts.std(dim=dim),
        "fano": lambda counts, dim: counts.var(dim=dim) / counts.mean(dim=dim),
        "mean_rank": lambda counts, dim: counts.mean(dim=dim).rank(dim="Gene"),
        "cv_squared": lambda counts, dim: counts.var(dim=dim) / counts.mean(dim=dim) ** 2
    }

    @staticmethod
    def hv_scatter_dist(dataset, x_kdims, y_kdims,
                        datashade_=False,
                        dynspread_=False):

        if dynspread_ and not datashade_:
            warnings.warn("Dynspread can only be used with datashade, setting both to true.")
            datashade_ = True

        df = dataset[[x_kdims, y_kdims]].to_dataframe()
        points = hv.Points(df, kdims=[x_kdims, y_kdims])

        dist_x = univariate_kde(hv.Distribution(points, kdims=[y_kdims], group="dist_x"), n_samples=1000)
        dist_y = univariate_kde(hv.Distribution(points, kdims=[x_kdims], group="dist_y"), n_samples=1000)

        if datashade_:
            points = datashade(points)
            if dynspread_:
                points = dynspread(points)

        return points << dist_x << dist_y

    @staticmethod
    def scatter_dist_by_mappings(dataset, x_kdims, y_kdims,
                                 mappings,
                                 selection_dim="Gene",
                                 datashade_=False,
                                 dynspread_=False,
                                 ):

        data_groups = {name: dataset.sel({selection_dim: genes}) for name, genes in mappings.items()}
        data_group_dfs = {k: v[[x_kdims, y_kdims]].to_dataframe() for k, v in data_groups.items()}

        points = {k: hv.Points(val, kdims=[x_kdims, y_kdims]) for k, val in data_group_dfs.items()}

        dist_x = {k: univariate_kde(hv.Distribution(p, kdims=[y_kdims], group="dist_x"), n_samples=1000)
                  for k, p in points.items()}
        dist_y = {k: univariate_kde(hv.Distribution(p, kdims=[x_kdims], group="dist_y"), n_samples=1000)
                  for k, p in points.items()}

        if datashade_:
            points_overlay = datashade(hv.NdOverlay(points))
            if dynspread_:
                points_overlay = dynspread(points_overlay)
        else:
            points_overlay = hv.NdOverlay(points)

        return points_overlay << hv.NdOverlay(dist_x) << hv.NdOverlay(dist_y)

    def process(self):
        counts = self.x_count_data
        dataset = xr.Dataset({
            self.x_axis_selector: self.axis_functions[self.x_axis_selector](counts, self.sample_index_name),
            self.y_axis_selector: self.axis_functions[self.y_axis_selector](counts, self.sample_index_name),
        })

        if self.axis_transform:
            dataset = self.axis_transform(dataset)

        if self.selected_gene_sets == [None] or not self.gene_set_collection:

            plot = self.hv_scatter_dist(
                dataset=dataset,
                x_kdims=self.x_axis_selector,
                y_kdims=self.y_axis_selector,
                datashade_=self.datashade_,
                dynspread_=self.dynspread_,
            )

        else:
            mappings = self.gene_set_collection.as_dict(self.selected_gene_sets)

            plot = self.scatter_dist_by_mappings(
                dataset=dataset,
                x_kdims=self.x_axis_selector,
                y_kdims=self.y_axis_selector,
                mappings=mappings,
                datashade_=self.datashade_,
                dynspread_=self.dynspread_,
            )

        return plot
