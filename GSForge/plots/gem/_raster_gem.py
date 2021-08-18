import datashader as datashader
import holoviews as hv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import param
import xarray as xr
from datashader import transfer_functions
from matplotlib.patches import Patch

from ..abstract_plot_models import InterfacePlottingBase


class RasterGEM(InterfacePlottingBase):
    # TODO: Document me.

    hue = param.String(default=None, doc="Color by which to shade the observations.")
    cmap = param.Parameter(default=None)
    plot_options = param.Parameter(default=dict(plot_width=800, plot_height=400))
    canvas_opts = param.Dict(default=None)

    @staticmethod
    def gem_raster(counts, canvas_opts=None):

        sample_index = np.arange(counts.shape[0])
        gene_index = np.arange(counts.shape[1])
        reset_counts = xr.DataArray(data=counts,
                                    name='counts',
                                    coords=[('sample_index', sample_index), ('gene_index', gene_index)])
        if canvas_opts is None:
            canvas_opts = dict(plot_width=800, plot_height=400)
        cvs = datashader.Canvas(**canvas_opts)

        return datashader.transfer_functions.shade(cvs.raster(reset_counts))

    @staticmethod
    def colorized_gem_raster(counts, labels, cmap=None, canvas_opts=None):

        if cmap is None:
            labels = pd.Series(labels)
            n_unique_labels = len(labels.unique())
            colors = hv.plotting.util.process_cmap(cmap='glasbey', ncolors=n_unique_labels, categorical=True)
            cmap = {label: color for label, color in zip(labels.unique(), colors)}

        # Create a legend.
        legend_items = [Patch(facecolor=color, edgecolor=color, label=label)
                        for label, color in cmap.items()]

        if canvas_opts is None:
            canvas_opts = dict(plot_width=800, plot_height=400)

        sample_index = np.arange(counts.shape[0])
        gene_index = np.arange(counts.shape[1])
        reset_counts = xr.DataArray(data=np.asarray(counts),
                                    name='counts',
                                    coords=[('sample_index', sample_index), ('gene_index', gene_index)])
        reset_counts = reset_counts.to_dataset()
        label_broadcast = xr.DataArray(data=np.broadcast_to(np.expand_dims(labels, axis=1), counts.shape),
                                       coords=[('sample_index', sample_index), ('gene_index', gene_index)])
        reset_counts['labels'] = label_broadcast

        unique_keys = pd.Series(labels).unique()
        cvs = datashader.Canvas(**canvas_opts)
        images = []

        for group in unique_keys:
            count_subset = reset_counts.where(reset_counts['labels'] == group)
            group_image = datashader.transfer_functions.shade(cvs.raster(count_subset['counts']), cmap=cmap[group])
            images.append(group_image)

        image = datashader.transfer_functions.stack(*images)

        # Create a legend.
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7), gridspec_kw={'width_ratios': [7, 1]})
        ax1.imshow(image.data.view(np.uint8).reshape(image.shape + (4,)))
        ax2.legend(handles=legend_items, loc='center')
        ax2.axis('off')
        plt.close()  # Prevents double-plotting.
        return fig

    def __call__(self):

        if self.hue is None:
            return self.gem_raster(self.x_count_data.values, self.get_default_options())

        else:
            self.param.set_param(annotation_variables=[self.hue])
            return self.colorized_gem_raster(
                counts=self.x_count_data.values,
                labels=self.y_annotation_data.values,
                cmap=self.cmap,
                canvas_opts=self.canvas_opts)
