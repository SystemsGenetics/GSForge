"""
Plotting functions for GEMprospector.

"""

# import pandas as pd
import numpy as np
import xarray as xr
import itertools
import networkx as nx
import holoviews as hv
from holoviews.operation.stats import univariate_kde
from holoviews.operation.datashader import datashade
import joypy
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm.autonotebook import tqdm
import matplotlib.patches as mpatches


from ._gene_mass import plot_count_sum_mass, plot_sample_mass

__all__ = [
    "lineament_connection_network",
    "scatter_dist_by_mappings",
    "ridge_plot",
    "mean_vs_variance",
    "np_sample_distributions",
    "plot_count_sum_mass",
    "plot_sample_mass",
]


def np_sample_distributions(counts: np.ndarray, labels: np.ndarray = None):
    """
    Calculate and overlay kernel density estimates of the given count matrix on a per-sample basis.

    If this gives a ValueError you may need to install statsmodels for a more robust kernel estimate.

    :param numpy.ndarray counts:  The count matrix to be displayed.
    :param numpy.ndarray labels: If provided the output density estimates will be colored by their
        label membership.

    :returns: matplotlib.axes with overlayed kernel density estimates.
    """

    fig, ax = plt.subplots(1, figsize=(15, 8))

    if labels is not None:
        label_set = list(np.unique(labels))
        colors = {label: color for label, color in zip(label_set, sns.color_palette(n_colors=len(label_set)))}

    for index, sample_row in tqdm(enumerate(counts), total=counts.shape[0]):

        if labels is not None:
            label = labels[index]
            color = colors[label]
            sns.kdeplot(sample_row, ax=ax, legend=True, shade=False, gridsize=250, color=color)

        else:
            sns.kdeplot(sample_row, ax=ax, legend=False, shade=False, gridsize=250)
    if labels is not None:
        patches = [mpatches.Patch(color=color, label=label)
                   for label, color in colors.items()]
        plt.legend(handles=patches)

    return ax


def lineament_connection_network(mappings):
    G = nx.Graph()

    for key, genes in mappings.items():
        G.add_node(key, type_="Feature")
        G.add_nodes_from(genes, type_="Gene")
        new_edges = itertools.product([key, ], genes)
        G.add_edges_from(new_edges)

    #         subset = lcoll.lineaments[key].data.sel(Gene=genes)
    #         weights = subset["feature_importances"].values

    #         new_edges = [(key, gene, {"weight": weight})
    #                      for gene, weight in zip(genes, weights)]

    graph = hv.Graph.from_networkx(G, nx.layout.kamada_kawai_layout)

    return graph.opts(tools=['hover'], padding=0.2, height=600, width=700,
                      color_index="node_type", node_size=5, edge_line_width=.5,
                      xaxis=None, yaxis=None)


def _freq_vs_means(counts: xr.DataArray, dim="Sample"):
    return xr.Dataset({'Mean': counts.mean(dim=dim),
                       'Frequency': ((counts > 0).sum(dim=dim) / counts[dim].shape)})


def _mean_vs_var(values, dim="Sample"):
    """Returns a dataframe of means and variances from the given values."""
    return xr.Dataset({'Mean': values.mean(dim=dim), 'Variance': values.var(dim=dim)})


def scatter_dist_by_mappings(data, xkdim, ykdim, mappings=None,
                             selection_dim="Gene",
                             use_datashade=False,
                             backend="bokeh"):
    # Setup backend-specific plotting options.
    # These options are specific to the backend -- i.e. matplotlib or bokeh.
    # We also need to create options for the optional datashader output.
    # Note that 'Area' is used instead of 'Distribution' as that is what the
    # 'univariate_kde' function returns.
    if backend == "bokeh":
        hv.output(backend='bokeh')
        points_opts = hv.opts.Points(width=500, height=500, bgcolor="lightgrey",
                                     size=1.2, muted_alpha=0.05, show_grid=True)
        dist_x_opts = hv.opts.Area(width=150, bgcolor="lightgrey", show_grid=True, show_legend=False, alpha=0.25)
        dist_y_opts = hv.opts.Area(height=150, bgcolor="lightgrey", show_grid=True, show_legend=False, alpha=0.25)
        shaded_opts = hv.opts.RGB(width=500, height=500, bgcolor="lightgrey", show_grid=True)
    elif backend == "matplotlib":
        hv.output(backend='matplotlib')
        points_opts = hv.opts.Points(s=1, aspect=1.75, fig_size=300, show_grid=True,
                                     bgcolor="#eeeeee")
        dist_x_opts = hv.opts.Area(bgcolor="lightgrey", show_grid=True, show_legend=False, alpha=0.25)
        dist_y_opts = hv.opts.Area(bgcolor="lightgrey", show_grid=True, show_legend=False, alpha=0.25)
        shaded_opts = hv.opts.RGB(width=500, height=500, bgcolor="lightgrey", show_grid=True)
    else:
        raise ValueError(f"{backend} is not a valid backend. Select from 'bokeh', or 'matplotlib'")

    # If no mappings are provided then the task is simplified.
    if mappings is None:
        df = data[[xkdim, ykdim]].to_dataframe()
        points = hv.Points(df, kdims=[xkdim, ykdim]).opts(points_opts)
        dist_x = univariate_kde(hv.Distribution(points, kdims=[ykdim]), n_samples=1000)
        dist_y = univariate_kde(hv.Distribution(points, kdims=[xkdim]), n_samples=1000)
        if use_datashade:
            points = datashade(points).opts(shaded_opts)
        return points << dist_x.opts(dist_x_opts) << dist_y.opts(dist_y_opts)

    # TODO: Repair options usage as done in the above if-block.
    # Otherwise we must create dictionaries of our mappings and overlay them.
    data_groups = {name: data.sel({selection_dim: genes}) for name, genes in mappings.items()}
    data_group_dfs = {k: v[[xkdim, ykdim]].to_dataframe() for k, v in data_groups.items()}

    points = {k: hv.Points(val, kdims=[xkdim, ykdim]) for k, val in data_group_dfs.items()}

    dist_x = {k: univariate_kde(hv.Distribution(p, kdims=[ykdim]), n_samples=1000).opts(dist_x_opts)
              for k, p in points.items()}
    dist_y = {k: univariate_kde(hv.Distribution(p, kdims=[xkdim]), n_samples=1000).opts(dist_y_opts)
              for k, p in points.items()}

    if use_datashade:
        points_overlay = datashade(hv.NdOverlay(points)).opts(shaded_opts)
    else:
        points_overlay = hv.NdOverlay(points).opts(points_opts)

    return points_overlay << hv.NdOverlay(dist_x) << hv.NdOverlay(dist_y)


# TODO: Add a numpy interface to this function.
def mean_vs_variance(count_array: xr.DataArray,
                     mappings=None,
                     use_datashade=False,
                     calc_dim="Sample",
                     mapping_dim="Gene",
                     backend="bokeh") -> hv.NdOverlay:
    """Plots the mean vs the variance of the given count array."""
    ds = xr.Dataset({'Mean': count_array.mean(dim=calc_dim),
                     'Variance': count_array.var(dim=calc_dim)})
    return scatter_dist_by_mappings(data=ds,
                                    xkdim="Mean",
                                    ykdim="Variance",
                                    mappings=mappings,
                                    selection_dim=mapping_dim,
                                    use_datashade=use_datashade,
                                    backend=backend)


def ridge_plot(counts, **kwargs):
    df = counts.to_dataframe()
    default_kwargs = dict(
        by="Sample",
        column="counts",
        bins=100,
        kind="counts",
        linecolor="black",
        linewidth=1,
        overlap=0.5,
        figsize=(10, 10),
        ylabels=False,
        fade=True)
    args = {**default_kwargs, **kwargs}
    return joypy.joyplot(data=df, **args)


# def volcano(deseq_result, kdims=None, vdims=None):
#     if kdims is None:
#         kdims = ["baseMean", "log2FoldChange"]
#
#     if vdims is None:
#         vdims = ["padj"]
#
#     levels = [0.0, 0.1, 1]
#     colors = ['blue', 'red']
#     opts = hv.opts.Points(logx=True,
#                           color='padj',
#                           color_levels=levels,
#                           cmap=colors,
#                           width=500,
#                           height=500)
#
#     points = hv.Points(data=deseq_result,
#                        kdims=kdims,
#                        vdims=vdims)
#
#     return points.options(opts)


def lineament_overlap_heatmap(lineament_collection, keys=None, mode="overlap"):
    lineament_dict = lineament_collection.as_dict(keys)

    if mode == "overlap":
        overlap_dim = hv.Dimension("Overlap Count", value_format=lambda x: f"{x}")
        data = [(ak, bk, len(set.intersection(set(av), set(bv))))
                for (ak, av), (bk, bv) in itertools.permutations(lineament_dict.items(), 2)
                if ak != bk]

    elif mode == "percent":
        overlap_dim = hv.Dimension("Overlap Percent", value_format=lambda x: f"{x:.0%}")

        zero_filtered_dict = {k: v for k, v in lineament_dict.items()
                              if len(v) > 0}
        data = [(ak, bk, len(set.intersection(set(av), set(bv))) / len(set(av)))
                for (ak, av), (bk, bv) in itertools.permutations(zero_filtered_dict.items(), 2)
                if ak != bk]

    heatmap = hv.HeatMap(data, vdims=overlap_dim).options(xrotation=45, width=450, height=450, labelled=[], colorbar=True)
    return heatmap * hv.Labels(heatmap)


def heat_map(counts):
    return hv.HeatMap(counts).options(xrotation=90, tools=['hover'])


def lineament_overlap_heat_map(overlaps):
    # TODO: add hover tooltip
    return hv.HeatMap(overlaps).options(xrotation=90)


def gene_vs_gene_scatter(counts, dim="Gene"):
    hv.output(backend='bokeh')

    def gvg_points(alpha, beta, dim="Gene"):
        alpha_vals = counts.sel({dim: alpha})
        beta_vals = counts.sel({dim: beta})
        return hv.Points((alpha_vals, beta_vals)).options(padding=0.1)

    genes = counts[dim].values
    dmap = hv.DynamicMap(gvg_points, kdims=['alpha', 'beta'])
    return dmap.redim.values(alpha=genes, beta=genes)


def gem_raster(counts):
    return datashade(hv.Image(counts.values)).options(width=600, height=400)


def plot_label_bars(label_df, max_x=10):
    bar_list = list()

    for col in label_df.columns:

        if label_df[col].nunique() > max_x:
            continue

        label_names = label_df.groupby([col]).nunique().index.values
        label_counts = label_df.groupby([col]).count().max().values

        new_bar = hv.Bars(zip(label_names, label_counts),
                          kdims=[col], vdims=["count"],
                          group='Label Counts', label=col).options(xrotation=45)

        bar_list.append(new_bar)

    return hv.Layout(bar_list)
