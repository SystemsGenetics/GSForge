"""
Plotting functions for GEMprospector.

"""

import holoviews as hv

from ._gene_mass import plot_count_sum_mass, plot_sample_mass, plot_sample
from ._scatter_distributions import ScatterDistributionBase
from ._distributions import SampleWiseDistribution
from ._lineament_heatmap import lineament_overlap_heatmap

# hv.extension("bokeh", "matplotlib")


__all__ = [
    "plot_count_sum_mass",
    "plot_sample_mass",
    "ScatterDistributionBase",
    "plot_sample",
    # "datashade_gem",
    # "colorized_raster",
    "lineament_overlap_heatmap",
]


def plot_label_bars(label_df, max_x=10):
    bar_list = list()

    for col in label_df.columns:

        if label_df[col].nunique() > max_x:
            continue

        label_names = label_df.groupby([col]).nunique().index.values
        label_counts = label_df.groupby([col]).count().max(axis=1).values

        new_bar = hv.Bars(zip(label_names, label_counts),
                          kdims=[col], vdims=["count"],
                          group='Label Counts', label=col)

        bar_list.append(new_bar)

    return hv.Layout(bar_list)
