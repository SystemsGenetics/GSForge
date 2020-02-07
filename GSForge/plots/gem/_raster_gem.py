import holoviews as hv
from holoviews.operation.datashader import datashade
import pandas as pd
import xarray as xr

import itertools
import param

from GSForge.models import OperationInterface


class RasterGEM(OperationInterface):

    # use_datashader = param.Boolean(default=True)
    hue = param.String(default=None, doc="Color by which to shade the observations.")

    @staticmethod
    def gem_raster(counts, use_datashader=True):
        if use_datashader:
            return datashade(hv.Image(counts.values))
        else:
            return hv.Image(counts.values)

    @staticmethod
    def colorized_raster(counts: xr.DataArray, labels:xr.DataArray, colors=None):
        if colors is None:
            colors = ["Blues", "Greens", "Greys", "Oranges", "Purples", "Reds"]

        label_df = labels.to_dataframe()
        # label_series = pd.Series(labels)
        groups = label_df.groupby(labels.name).groups

        colors = {name: color for name, color in zip(groups.keys(), itertools.cycle(colors))}

        data_groups = {name: counts.where(labels == name)
                       for name in groups.keys()}

        images = [hv.Image(values.values).opts(cmap=colors[name])
                  for name, values in data_groups.items()]
        # images = {name: hv.Image(values.values).opts(cmap=colors[name])
        #           for name, values in data_groups.items()}

        return hv.Overlay(images)
        # return hv.NdOverlay(images)

    def process(self):
        if self.hue is None:
            raster = self.gem_raster(self.x_count_data)
        else:
            raster = self.colorized_raster(self.x_count_data, self.y_annotation_data[self.hue])

        return raster.opts()
