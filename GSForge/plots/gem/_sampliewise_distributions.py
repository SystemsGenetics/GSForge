import holoviews as hv
import holoviews.plotting
import numpy as np
import pandas as pd
import param
import xarray as xr
from holoviews.operation.stats import univariate_kde

from ..abstract_plot_models import InterfacePlottingBase


class SamplewiseDistributions(InterfacePlottingBase):
    """Provides a base that iterates through the selected samples."""
    hue_key = param.Parameter(default=None)
    hue_colors = param.Parameter(default=None)

    cmap = param.Parameter('glasbey')

    sample_colormap = param.Parameter(doc=""" """)
    x_axis_transform = param.Parameter(default=('log 2', lambda ds: np.log2(ds.where(ds > 0))),
                                       doc="A transform (usually log2) for getting a viewable spread of the results.")
    y_axis_transform = param.Parameter(default=None,  # ('log 2', lambda ds: np.log2(ds.where(ds > 0))),
                                       doc="A transform (usually log2) for getting a viewable spread of the results.")

    datashade = param.Boolean(default=False)
    dynspread = param.Boolean(default=False)
    kde_kwargs = param.Dict(default={'filled': False})

    @staticmethod
    def bokeh_opts():
        return [hv.opts.Overlay(width=600, height=400)]

    @staticmethod
    def matplotlib_opts():
        return hv.opts.Curve()

    @staticmethod
    def sample_wise_count_distributions(counts: xr.DataArray,
                                        sample_dim: str = "Sample",
                                        x_axis_name='Expression',
                                        **kde_kwargs,
                                        ) -> hv.Overlay:
        """Draw sample-wise count distributions using a univariate kernel density estimate."""
        # kde_kwargs = {**{'filled': filled}, **kde_kwargs}

        kde_plots = list()
        for sample in counts[sample_dim].values:
            sample_count_values = counts.sel(**{sample_dim: sample}).values
            sample_dist = hv.Distribution(sample_count_values, kdims=[x_axis_name])
            sample_kde = univariate_kde(sample_dist, **kde_kwargs)
            kde_plots.append(sample_kde)

        layout = hv.Overlay(kde_plots)
        return layout

    @classmethod
    def color_keyed_sample_wise_count_distribution(
            cls,
            counts: xr.DataArray,
            labels: xr.DataArray,
            hue_colors: dict = None,
            sample_dim: str = "Sample",
            x_axis_name='Expression',
            cmap='glasbey',
            **kde_kwargs,
    ):
        if hue_colors is None:
            unique_keys = pd.Series(labels).unique()
            colors = hv.plotting.util.process_cmap(cmap=cmap, ncolors=len(unique_keys))
            hue_colors = {label: color for label, color in zip(unique_keys, colors)}

        group_overlays = []
        for group, group_df in labels.to_dataframe().groupby(labels.name):
            group_counts = counts.sel(**{sample_dim: group_df.index.values})
            group_overlay = cls.sample_wise_count_distributions(
                counts=group_counts,
                sample_dim=sample_dim,
                x_axis_name=x_axis_name,
                **kde_kwargs
            )
            group_overlay.opts(hv.opts.Curve(color=hue_colors[group]))
            group_overlay.param.set_param(group=group)
            group_overlays.append(group_overlay)

        layout = hv.Overlay(group_overlays)
        return layout

    def __call__(self):
        # Check options.
        # if self.dynspread is True and (self.datashader is False): ...

        # Prepare and optionally transform the data and axis labels.
        x_axis_name = 'counts'
        counts = self.x_count_data
        if self.x_axis_transform:
            x_axis_name = f'{self.x_axis_transform[0]} {x_axis_name}'
            counts = self.x_axis_transform[1](counts)

        # Check if the default configuration should be used.
        if (self.hue_key is None) and (self.datashade is False):
            layout = self.sample_wise_count_distributions(
                counts=counts,
                sample_dim=self.sample_index_name,
                x_axis_name=x_axis_name,
                **self.kde_kwargs
            )

        # Check if the datashaded option should be used.
        elif self.datashade is True:
            pass
            # Note that hue_key has no effect with datashader yet.

        # Check if a hue_key is provided.
        elif self.hue_key is not None:

            self.param.set_param(annotation_variables=[self.hue_key])

            layout = self.color_keyed_sample_wise_count_distribution(
                counts=counts,
                labels=self.y_annotation_data,
                sample_dim=self.sample_index_name,
                x_axis_name=x_axis_name,
                **self.kde_kwargs
            )

        if self.apply_default_opts:
            layout = layout.opts(self.get_default_options())
        return layout



class EmpiricalCumulativeDistribution(InterfacePlottingBase):
    hue_key = param.Parameter()
    hue_colors = param.Parameter()

    cmap = param.Parameter('glasbey')

    sample_colormap = param.Parameter(doc=""" """)
    axis_transform = param.Parameter(default=None,  # ('log 2', lambda ds: np.log2(ds.where(ds > 0))),
                                     doc="A transform (usually log2) for getting a viewable spread of the results.")
    datashade = param.Boolean(default=False)
    dynspread = param.Boolean(default=False)

    # @staticmethod
    # def build_colormap(samples, labels=None,
    #                    hue_key=None, hue_colors=None,
    #                    cmap='glasbey'):
    #     """
    #     Construct a colormap for the samples based on ``hue_key``, ``hue_colormap``
    #     and ``sample_colormap``.
    #
    #     ``sample_colormap`` has the highest priority, and the other arguments are ignored
    #     if it is provided.
    #
    #     """
    #
    #     if (labels is not None) and hue_key and not hue_colors:
    #         unique_keys = pd.Series(labels).unique()
    #         colors = hv.plotting.util.process_cmap(cmap=cmap, ncolors=len(unique_keys))
    #         hue_colors = {label: color for label, color in zip(unique_keys, colors)}
    #         sample_colors = {sample: hue_colors[label] for sample, label in zip(samples, labels)}
    #
    #     elif (labels is not None) and hue_key and hue_colors:
    #         sample_colors = {sample: hue_colors[label] for sample, label in zip(samples, labels)}
    #
    #     else:
    #         sample_colors = {k: v for k, v in zip(samples, itertools.cycle(hv.plotting.util.process_cmap(cmap=cmap)))}
    #
    #     return sample_colors

    @staticmethod
    def plot_ecdf(values: np.ndarray, xaxis_name='expression', yaxis_name='ECDF'):
        df = pd.DataFrame({
            xaxis_name: np.sort(values),
            yaxis_name: np.arange(1, values.shape[0] + 1) / values.shape[0]
        })
        return hv.Curve(df)

    @staticmethod
    def bokeh_opts():
        return [hv.opts.Curve(logx=True, width=600, height=500, )]

    @staticmethod
    def matplotlib_opts():
        return [hv.opts.Curve(aspect=3.0, fig_size=300, show_legend=False, logx=True)]

    def __call__(self):

        self.param.set_param(count_mask='masked')

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

            ecdf_plot = hvds.to(hv.Curve, kdims, vdims)

            if self.hue_key:
                ecdf_plot = ecdf_plot.groupby(self.hue_key)
                ecdf_plot = ecdf_plot.overlay(self.hue_key).opts(
                    hv.opts.Curve(color=hv.dim(self.hue_key).categorize(hue_colors)))

            layout[sample] = ecdf_plot

        return layout.opts(self.get_default_options())
