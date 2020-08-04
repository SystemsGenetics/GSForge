import holoviews as hv
# from holoviews.operation.datashader import datashade
import xarray as xr
# from holoviews.operation import datashader
from holoviews.operation.datashader import datashade, dynspread
# from datashader.reductions import count_cat
import itertools
import param
import warnings
import numpy as np
import pandas as pd
import datashader as ds
from datashader import transfer_functions

from ..abstract_plot_models import InterfacePlottingBase


class RasterGEM(InterfacePlottingBase):
    # TODO: Document me.

    color_key = param.String(default=None, doc="Color by which to shade the observations.")
    cmap = param.Parameter(default=None)
    plot_options = param.Parameter(default=dict(plot_width=800, plot_height=400))


    @staticmethod
    def gem_raster(counts, canvas_opts=None):


        sample_index = np.arange(counts.shape[0])
        gene_index = np.arange(counts.shape[1])
        reset_counts = xr.DataArray(data=counts,
                                    name='counts',
                                    coords=[('sample_index', sample_index), ('gene_index', gene_index)])
        if canvas_opts is None:
            canvas_opts = dict(plot_width=800, plot_height=400)
        cvs = ds.Canvas(**canvas_opts)

        return ds.transfer_functions.shade(cvs.raster(reset_counts))

    @staticmethod
    def colorized_gem_raster(counts, labels, cmap, canvas_opts=None):
        if canvas_opts is None:
            canvas_opts = dict(plot_width=800, plot_height=400)

        sample_index = np.arange(counts.shape[0])
        gene_index = np.arange(counts.shape[1])
        reset_counts = xr.DataArray(data=counts,
                                    name='counts',
                                    coords=[('sample_index', sample_index), ('gene_index', gene_index)])
        reset_counts = reset_counts.to_dataset()
        label_broadcast = xr.DataArray(data=np.broadcast_to(np.expand_dims(labels, axis=1), counts.shape),
                                       coords=[('sample_index', sample_index), ('gene_index', gene_index)])
        reset_counts['labels'] = label_broadcast

        unique_keys = pd.Series(labels).unique()
        cvs = ds.Canvas(**canvas_opts)
        images = []

        for group in unique_keys:
            count_subset = reset_counts.where(reset_counts['labels'] == group)
            group_image = ds.transfer_functions.shade(cvs.raster(count_subset['counts']), cmap=cmap[group])
            images.append(group_image)

        image = ds.transfer_functions.stack(*images)
        return image

    def __call__(self):

        if self.color_key is None:
            return self.gem_raster(self.x_count_data.values, self.get_default_options())

        else:
            # Check for cmap.
            self.param.set_param(annotation_variables=[self.color_key])
            return self.colorized_gem_raster(
                counts=self.x_count_data.values,
                labels=self.y_annotation_data.values,
                cmap=self.cmap,
                canvas_opts=self.cmap)

# cmap = {
#     'CONTROL': hv.plotting.util.process_cmap('gray'),
#     'HEAT': hv.plotting.util.process_cmap('blues'),
#     'RECOV_HEAT': hv.plotting.util.process_cmap('fire'),
#     'DROUGHT': hv.plotting.util.process_cmap('blues'),
#     'RECOV_DROUGHT': hv.plotting.util.process_cmap('kg'),
# }