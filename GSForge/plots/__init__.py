"""
Plotting functions for GSForge.

"""

import holoviews as hv

from ._gene_mass import plot_count_sum_mass, plot_sample_mass, plot_sample
from ._scatter_distributions import ScatterDistributionBase
from ._distributions import SampleWiseDistribution
from ._geneset_heatmap import geneset_overlap_heatmap
from ._raster_gem import RasterGEM
from ._compare_gene_counts import GenesVsCounts


__all__ = [
    "plot_count_sum_mass",
    "plot_sample_mass",
    "ScatterDistributionBase",
    "plot_sample",
    "SampleWiseDistribution",
    "RasterGEM",
    "geneset_overlap_heatmap",
    "GenesVsCounts"
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
