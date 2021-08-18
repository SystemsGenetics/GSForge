import logging

import datashader
import holoviews as hv
import holoviews.operation.datashader
import holoviews.plotting
import numpy as np
import pandas as pd
import param
import scipy.stats
import xarray as xr

from ..abstract_plot_models import InterfacePlottingBase

logger = logging.getLogger("GSForge")


def generate_color_key_array(label_series: pd.Series, cmap='glasbey'):
    unique_keys = label_series.unique()
    colors = hv.plotting.util.process_cmap(cmap=cmap, ncolors=len(unique_keys))
    hue_colors = {label: color for label, color in zip(unique_keys, colors)}
    color_key = list(label_series.map(hue_colors).values)
    return color_key


class SamplewiseDistributions(InterfacePlottingBase):
    """Provides a base that iterates through the selected samples."""
    hue_key = param.Parameter(default=None)
    hue_colors = param.Parameter(default=None)

    cmap = param.Parameter('glasbey')

    sample_colormap = param.Parameter(doc=""" """)
    # x_axis_transform = param.Parameter(default=('log 2', lambda ds: np.log2(ds + 0.25)),
    #                                    doc="A transform (usually log2) for getting a viewable spread of the results.")

    datashade = param.Boolean(default=False)
    dynspread = param.Boolean(default=False)
    kde_kwargs = param.Dict(default={'filled': False})
    kde_sample_size = param.Integer(default=200)

    @staticmethod
    def bokeh_opts():
        return [hv.opts.NdOverlay(width=600, height=600, show_legend=False),
                hv.opts.Overlay(width=600, height=600, show_legend=False),
                hv.opts.Curve(width=600, height=600, line_width=1, bgcolor="lightgrey", show_grid=True, padding=0.02),
                hv.opts.RGB(width=600, height=600, show_grid=True, bgcolor="lightgrey")]

    @staticmethod
    def matplotlib_opts():
        return hv.opts.Curve(fig_size=200, linewidth=1, bgcolor="lightgrey", show_grid=True, padding=0.02)

    @staticmethod
    def kde_linespace(count_array, bin_range, sample_size, cut=3):
        """
        Apply the math behind the univariate_kde as implemented by Holoviews so that
        it can be applied with `numpy.apply_along_axis`.
        """
        kde = scipy.stats.gaussian_kde(dataset=count_array, bw_method=None, weights=None)
        bw = kde.scotts_factor() * count_array.std(ddof=1)
        kmin, kmax = bin_range[0] - bw * cut, bin_range[1] + bw * cut
        xs = np.linspace(kmin, kmax, sample_size)
        return xs

    @staticmethod
    def evaluate_kde(values, x_space):
        kde_model = scipy.stats.gaussian_kde(dataset=values, bw_method=None, weights=None)
        return kde_model.evaluate(x_space)

    @classmethod
    def plot_sample_wise_kde(cls, count_array, sample_dim='Sample', color_key=None, x_axis_name='counts', sample_size=100):
        bin_range = (count_array.min().values, count_array.max().values)

        x_space = cls.kde_linespace(count_array.fillna(0).values,
                                    bin_range=bin_range,
                                    sample_size=sample_size)
        kde_dist_y = np.apply_along_axis(func1d=cls.evaluate_kde,
                                         axis=1,
                                         arr=count_array.fillna(0).values,
                                         x_space=x_space)

        df = pd.DataFrame(kde_dist_y, columns=x_space, index=count_array[sample_dim].values)

        overlay = hv.NdOverlay(kdims=sample_dim)
        for index, (sample_name, row_series) in enumerate(df.iterrows()):
            plot = hv.Curve((df.columns.values, row_series.values), kdims=x_axis_name, vdims='Distribution')
            if color_key is not None:
                plot = plot.opts(color=color_key[index])
            overlay[sample_name] = plot

        return overlay

    def __call__(self):
        # Check options.
        # if self.dynspread is True and (self.datashader is False): ...
        # self.param.set_param(count_mask='masked')
        # Prepare and optionally transform the data and axis labels.

        x_axis_name = self.count_variable
        counts = self.x_count_data

        # if self.x_axis_transform:
        #     x_axis_name = f'{self.x_axis_transform[0]} {x_axis_name}'
        #     counts = self.x_axis_transform[1](counts)

        if self.hue_key is not None:
            self.param.set_param(annotation_variables=[self.hue_key])
            color_key = generate_color_key_array(self.y_annotation_data.to_series())
        else:
            color_key = hv.plotting.util.process_cmap(cmap='glasbey', ncolors=counts[self.sample_index_name].shape[0])
            # color_key = None

        # Construct the plot.
        plot = self.plot_sample_wise_kde(counts, sample_dim=self.sample_index_name,
                                         color_key=color_key, x_axis_name=x_axis_name, sample_size=self.kde_sample_size)

        if self.datashade is True:
            plot = hv.operation.datashader.datashade(plot,
                                                     color_key=color_key,
                                                     aggregator=datashader.count_cat(self.sample_index_name))

            if self.dynspread is True:
                plot = hv.operation.datashader.dynspread(plot)

        if self.apply_default_opts is True:
            plot = plot.opts(self.get_default_options())

        return plot


class EmpiricalCumulativeDistribution(InterfacePlottingBase):
    hue_key = param.Parameter()
    hue_colors = param.Parameter()

    cmap = param.Parameter('glasbey')

    sample_colormap = param.Parameter(doc=""" """)
    x_axis_transform = param.Parameter(default=('log 2', lambda ds: np.log2(ds + 0.25)),
                                       doc="A transform (usually log2) for getting a viewable spread of the results.")
    datashade = param.Boolean(default=False)
    dynspread = param.Boolean(default=False)

    @staticmethod
    def plot_sample_wise_ecdf(counts: xr.DataArray,
                              sample_dim: str = "Sample",
                              color_key: np.ndarray = None,
                              x_axis_name='counts'):
        """Draw sample-wise empirical cumulative distribution functions."""
        sorted_counts = np.apply_along_axis(np.sort, 1, counts.values)
        y_axis_values = np.arange(1, counts.shape[1] + 1 / counts.shape[1])
        # print((y_axis_values.shape))
        overlay = hv.NdOverlay(kdims=sample_dim)
        for index, (sample, sample_counts) in enumerate(zip(counts[sample_dim].values, sorted_counts)):
            plot = hv.Curve((sample_counts, y_axis_values), kdims=x_axis_name, vdims='ECDF')
            if color_key is not None:
                plot = plot.opts(color=color_key[index])
            overlay[sample] = plot
        return overlay

    @staticmethod
    def bokeh_opts():
        return [hv.opts.NdOverlay(width=600, height=400, show_legend=False),
                hv.opts.Curve(width=600, height=400, line_width=1, bgcolor="lightgrey", show_grid=True, padding=0.02),
                hv.opts.RGB(width=600, height=400, show_grid=True, bgcolor="lightgrey")]

    @staticmethod
    def matplotlib_opts():
        return [hv.opts.Curve(aspect=3.0, fig_size=300, show_legend=False, logx=True)]

    def __call__(self):

        self.param.set_param(count_mask='masked')

        x_axis_name = self.count_variable
        counts = self.x_count_data

        if self.x_axis_transform is not None:
            x_axis_name = f'{self.x_axis_transform[0]} {x_axis_name}'
            counts = self.x_axis_transform[1](counts)

        if self.hue_key is not None:
            self.param.set_param(annotation_variables=[self.hue_key])
            color_key = generate_color_key_array(self.y_annotation_data.to_series())
        else:
            color_key = hv.plotting.util.process_cmap(cmap='glasbey', ncolors=counts[self.sample_index_name].shape[0])

        plot = self.plot_sample_wise_ecdf(counts=counts,
                                          sample_dim=self.sample_index_name,
                                          color_key=color_key,
                                          x_axis_name=x_axis_name)

        if self.datashade is True:
            plot = hv.operation.datashader.datashade(plot,
                                                     color_key=color_key,
                                                     aggregator=datashader.count_cat(self.sample_index_name))

            if self.dynspread is True:
                plot = hv.operation.datashader.dynspread(plot)

        if self.apply_default_opts is True:
            plot = plot.opts(self.get_default_options())

        return plot
