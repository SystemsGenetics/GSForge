import holoviews as hv
from holoviews.operation.datashader import datashade
import pandas as pd

import itertools
import param

from ..models import OperationInterface


class RasterGEM(OperationInterface):

    # use_datashader = param.Boolean(default=True)

    @staticmethod
    def gem_raster(counts, use_datashader=True):
        if use_datashader:
            return datashade(hv.Image(counts.values))
        else:
            return hv.Image(counts.values)

    @staticmethod
    def colorized_raster(counts, labels, colors=None):
        if colors is None:
            colors = ["Blues", "Greens", "Greys", "Oranges", "Purples", "Reds"]

        label_series = pd.Series(labels)
        groups = pd.DataFrame({"label": label_series}).groupby("label").groups

        colors = {name: color for name, color in zip(groups.keys(), itertools.cycle(colors))}

        data_groups = {name: counts.where(labels == name)
                       for name in groups.keys()}

        images = {name: hv.Image(values.values).opts(cmap=colors[name])
                  for name, values in data_groups.items()}

        return hv.NdOverlay(images)  #.options(width=600, height=400, xaxis=None, yaxis=None)

    def process(self):
        if self.selected_annotation_variables is None:
            raster = self.gem_raster(self.x_count_data)
        else:
            raster = self.colorized_raster(self.x_count_data, self.y_annotation_data)

        return raster.opts()
