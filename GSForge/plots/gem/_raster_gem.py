import holoviews as hv
# from holoviews.operation.datashader import datashade
import xarray as xr
from holoviews.operation import datashader
# from holoviews.operation.datashader import datashade, dynspread
# from datashader.reductions import count_cat
import itertools
import param
import warnings

from ..abstract_plot_models import InterfacePlottingBase


class RasterGEM(InterfacePlottingBase):
    # TODO: Document me.
    # TODO: Fix colorized_raster. It's not displaying colors.

    hue = param.String(default=None, doc="Color by which to shade the observations.")

    @staticmethod
    def bokeh_opts():
        return hv.opts.RGB(width=900, height=500, xaxis=None, yaxis=None, labelled=[])

    @staticmethod
    def matplotlib_opts():
        return hv.opts.RGB(fig_size=250)

    @staticmethod
    def gem_raster(counts, use_datashader=True):
        if use_datashader:
            return hv.operation.datashader.datashade(hv.Image(counts.values))
        else:
            return hv.Image(counts.values)

    @staticmethod
    def colorized_raster(counts: xr.DataArray, labels:xr.DataArray, colors=None):
        warnings.warn("Group color maps are not yet working correctly.")
        if colors is None:
            colors = ["Blues", "Greens", "Greys", "Oranges", "Purples", "Reds"]

        label_df = labels.to_dataframe()
        groups = label_df.groupby(labels.name).groups

        colors = {name: color for name, color in zip(groups.keys(), itertools.cycle(colors))}

        data_groups = {name: counts.where(labels == name)
                       for name in groups.keys()}

        images = [hv.Image(values.values).opts(cmap=colors[name])
                  for name, values in data_groups.items()]

        return hv.Overlay(images)

    def __call__(self):
        if self.hue is None:
            raster = self.gem_raster(self.x_count_data)
        else:
            raster = self.colorized_raster(self.x_count_data, self.y_annotation_data[self.hue])

        return raster.opts(self.get_default_options())
