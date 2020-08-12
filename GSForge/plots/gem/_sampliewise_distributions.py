import holoviews as hv
import holoviews.plotting
import numpy as np
import pandas as pd
import param
import scipy.stats
import xarray as xr
import datashader
import holoviews.operation.datashader
from ..abstract_plot_models import InterfacePlottingBase


def generate_color_key_array(label_series: pd.Series):
    unique_keys = label_series.unique()
    colors = hv.plotting.util.process_cmap(cmap='glasbey', ncolors=len(unique_keys))
    hue_colors = {label: color for label, color in zip(unique_keys, colors)}
    color_key = label_series.map(hue_colors)
    return color_key


class SamplewiseDistributions(InterfacePlottingBase):
    """Provides a base that iterates through the selected samples."""
    hue_key = param.Parameter(default=None)
    hue_colors = param.Parameter(default=None)

    cmap = param.Parameter('glasbey')

    sample_colormap = param.Parameter(doc=""" """)
    x_axis_transform = param.Parameter(default=('log 2', lambda ds: np.log2(ds + 0.25)),
                                       doc="A transform (usually log2) for getting a viewable spread of the results.")
    y_axis_transform = param.Parameter(default=None,  # ('log 2', lambda ds: np.log2(ds.where(ds > 0))),
                                       doc="A transform (usually log2) for getting a viewable spread of the results.")

    datashade = param.Boolean(default=False)
    dynspread = param.Boolean(default=False)
    kde_kwargs = param.Dict(default={'filled': False})

    @staticmethod
    def bokeh_opts():
        return [hv.opts.NdOverlay(width=600, height=400, show_legend=False),
                hv.opts.Overlay(width=600, height=400, show_legend=False),
                hv.opts.Curve(width=600, height=400, line_width=1, bgcolor="lightgrey", show_grid=True, padding=0.02),
                hv.opts.RGB(width=600, height=400, show_grid=True, bgcolor="lightgrey")]

    # @staticmethod
    # def matplotlib_opts():
    #     return hv.opts.Curve()

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
    def plot_sample_wise_kde(cls, count_array, sample_dim='Sample', color_key=None):
        bin_range = (count_array.min().values, count_array.max().values)

        x_space = cls.kde_linespace(count_array.values, bin_range=bin_range, sample_size=50)
        kde_dist_y = np.apply_along_axis(func1d=cls.evaluate_kde, axis=1, arr=count_array.values, x_space=x_space)

        df = pd.DataFrame(kde_dist_y, columns=x_space, index=count_array[sample_dim].values)

        overlay = hv.NdOverlay(kdims=sample_dim)
        for index, (sample_name, row_series) in enumerate(df.iterrows()):
            plot = hv.Curve((df.columns.values, row_series.values))
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
        if self.x_axis_transform:
            x_axis_name = f'{self.x_axis_transform[0]} {x_axis_name}'
            counts = self.x_axis_transform[1](counts)
        if self.hue_key is not None:
            self.param.set_param(annotation_variables=[self.hue_key])
            color_key = generate_color_key_array(self.y_annotation_data.to_series())
        else:
            color_key = hv.plotting.util.process_cmap(cmap='glasbey', ncolors=counts[self.sample_index_name].shape[0])

        # Construct the plot.
        plot = self.plot_sample_wise_kde(counts, sample_dim=self.sample_index_name, color_key=color_key)

        if self.datashade is True:
            if self.hue_key is not None:
                plot = hv.operation.datashader.datashade(plot,
                                                         aggregator=datashader.count_cat(self.sample_index_name))

            else:
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
    x_axis_transform = param.Parameter(default=None,  # ('log 2', lambda ds: np.log2(ds.where(ds > 0))),
                                       doc="A transform (usually log2) for getting a viewable spread of the results.")
    datashade = param.Boolean(default=False)
    dynspread = param.Boolean(default=False)

    @staticmethod
    def plot_sample_wise_ecdf(counts: xr.DataArray,
                              sample_dim: str = "Sample",
                              color_key: np.ndarray = None,
                              x_axis_name='Expression'):
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

    # @staticmethod
    # def plot_ecdf(values: np.ndarray, xaxis_name='expression', yaxis_name='ECDF'):
    #     df = pd.DataFrame({
    #         xaxis_name: np.sort(values),
    #         yaxis_name: np.arange(1, values.shape[0] + 1) / values.shape[0]
    #     })
    #     return hv.Curve(df)

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

        if self.hue_key:
            self.param.set_param(annotation_variables=[self.hue_key])
            color_key = generate_color_key_array(self.y_annotation_data.to_series())

        x_axis_name = self.count_variable

        counts = self.x_count_data
        if self.x_axis_transform is not None:
            x_axis_name = f'{self.x_axis_transform[0]} {x_axis_name}'
            counts = self.x_axis_transform[1](counts)

        plot = self.plot_sample_wise_ecdf(counts=counts, sample_dim=self.sample_index_name,
                                          color_key=color_key if self.hue_key else None,
                                          x_axis_name=x_axis_name)

        if self.datashade is True:
            plot = hv.operation.datashader.datashade(plot,
                                                     color_key=self.hue_key if self.hue_key is not None else None,
                                                     aggregator=datashader.count_cat(self.sample_index_name))
            if self.dynspread is True:
                plot = hv.operation.datashader.dynspread(plot)

        if self.apply_default_opts is True:
            plot = plot.opts(self.get_default_options())

        return plot
