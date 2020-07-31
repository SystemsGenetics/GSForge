import itertools
import param
import colorcet as cc
import holoviews as hv
import numpy as np
import pandas as pd
import xarray as xr
from holoviews.operation.stats import univariate_kde

from ..abstract_plot_models import InterfacePlottingBase


class SamplewiseDistributions(InterfacePlottingBase):
    # TODO: Document me.
    # TODO: Add matplotlib options.

    @staticmethod
    def bokeh_opts():
        return hv.opts.Curve()

    @staticmethod
    def matplotlib_opts():
        return hv.opts.Curve()

    @staticmethod
    def sample_wise_count_distributions(counts: xr.DataArray, sample_dim: str = "Sample",
                                        cmap: list = cc.glasbey) -> hv.Overlay:
        distributions = [
            univariate_kde(hv.Distribution(counts.sel(**{sample_dim: sample}).values,
                                           ), filled=False)#.opts(line_color=color)
            for sample, color in zip(counts[sample_dim].values, itertools.cycle(cmap))]
        return hv.Overlay(distributions).opts(show_legend=False)

    def process(self):
        layout = self.sample_wise_count_distributions(counts=self.x_count_data,
                                                      sample_dim=self.sample_index_name)
        return layout


class EmpiricalCumulativeDistribution(InterfacePlottingBase):
    hue_key = param.Parameter()
    hue_colors = param.Parameter()

    cmap = param.Parameter('glasbey')

    sample_colormap = param.Parameter(doc=""" """)
    axis_transform = param.Parameter(default=None,  # ('log 2', lambda ds: np.log2(ds.where(ds > 0))),
                                     doc="A transform (usually log2) for getting a viewable spread of the results.")

    @staticmethod
    def build_colormap(samples, labels=None,
                       hue_key=None, hue_colors=None,
                       cmap='glasbey'):
        """
        Construct a colormap for the samples based on ``hue_key``, ``hue_colormap``
        and ``sample_colormap``.

        ``sample_colormap`` has the highest priority, and the other arguments are ignored
        if it is provided.

        """

        if (labels is not None) and hue_key and not hue_colors:
            unique_keys = pd.Series(labels).unique()
            colors = hv.plotting.util.process_cmap(cmap=cmap, ncolors=len(unique_keys))
            hue_colors = {label: color for label, color in zip(unique_keys, colors)}
            sample_colors = {sample: hue_colors[label] for sample, label in zip(samples, labels)}

        elif (labels is not None) and hue_key and hue_colors:
            sample_colors = {sample: hue_colors[label] for sample, label in zip(samples, labels)}

        else:
            sample_colors = {k: v for k, v in zip(samples, itertools.cycle(hv.plotting.util.process_cmap(cmap=cmap)))}

        return sample_colors

    @staticmethod
    def plot_ecdf(values: np.ndarray, xaxis_name='expression', yaxis_name='ECDF'):
        df = pd.DataFrame({
            xaxis_name: np.sort(values),
            yaxis_name: np.arange(1, values.shape[0] + 1) / values.shape[0]
        })
        return hv.Curve(df)

    @staticmethod
    def bokeh_opts():
        return [
            hv.opts.Curve(
                logx=True,
                width=600,
                height=500,
            )
        ]

    @staticmethod
    def matplotlib_opts():
        return [
            hv.opts.Curve(aspect=3.0,
                          fig_size=300,
                          show_legend=False,
                          logx=True),
        ]

    def __call__(self):

        self.param.set_param(count_mask='masked')

        samples = self.get_sample_index()

        if self.hue_key:
            self.param.set_param(annotation_variables=[self.hue_key])

        labels = self.y_annotation_data if (self.y_annotation_data is not None) else None

        if self.hue_key:
            unique_keys = pd.Series(labels).unique()
            colors = hv.plotting.util.process_cmap(cmap=self.cmap, ncolors=len(unique_keys))
            hue_colors = {label: color for label, color in zip(unique_keys, colors)}

        x_axis_name = 'expression'
        y_axis_name = 'ECDF'

        if self.axis_transform:
            x_axis_name = f'{self.axis_transform[0]} {x_axis_name}'
            y_axis_name = f'{self.axis_transform[0]} {y_axis_name}'
            counts = self.axis_transform[1](self.x_count_data)
        else:
            counts = self.x_count_data

        kdims = [x_axis_name]
        vdims = [y_axis_name, self.hue_key] if (self.hue_key is not None) else [y_axis_name]

        # layout = []
        if self.hue_key:
            layout = hv.NdOverlay(kdims=[self.hue_key])
        else:
            layout = hv.NdOverlay()

        for sample_counts in counts:

            sample = str(counts['Sample'].values)

            plot_data = pd.DataFrame({
                x_axis_name: np.sort(sample_counts.values),
                y_axis_name: np.arange(1, sample_counts.values.shape[0] + 1) / sample_counts.values.shape[0],
            })

            if self.hue_key:
                plot_data[self.hue_key] = labels.sel(Sample=sample).values

            hvds = hv.Dataset(plot_data, kdims=kdims)

            ecdf_plot = hvds.to(hv.Curve, kdims, vdims)  # , groupby=self.hue_key)

            if self.hue_key:
                ecdf_plot = ecdf_plot.groupby(self.hue_key)
                ecdf_plot = ecdf_plot.overlay(self.hue_key).opts(
                    hv.opts.Curve(color=hv.dim(self.hue_key).categorize(hue_colors)))

            layout[sample] = ecdf_plot

        return layout.opts(self.get_default_options())
