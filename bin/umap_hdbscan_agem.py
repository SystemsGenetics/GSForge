#!/usr/bin/env python3

import GSForge as gsf
import holoviews as hv
import datashader
import hdbscan
from matplotlib.patches import Patch
import matplotlib.pyplot as plt

hv.extension("matplotlib")


def export_cmap_legend(cmap, filename="legend.png"):
    fig, ax = plt.subplots()
    ax.axis('off')
    custom_lines = [Patch(facecolor=cmap[k], label=k) for k in cmap.keys()]
    ncols = len(cmap.keys())
    if ncols > 50:
        ncols = 50
    legend = ax.legend(custom_lines, list(cmap.keys()), frameon=False, loc='lower center', ncol=ncols)

    fig = ax.figure
    fig.canvas.draw()

    bbox = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(filename, dpi="figure", bbox_inches=bbox)


agem = gsf.AnnotatedGEM("${params.gem_netcdf}")
agem.data.load()

datashade_kwargs = dict(
    height=${params.datashade_kwargs.height},
            width =${params.datashade_kwargs.width})

umap_panel = gsf.panels.UMAP_Interface(
    agem,
    datashade=True,
    random_state=${random_state},
                  hue = "${hue}",
                        n_neighbors =${n_neighbors},
                                      datashade_kwargs = datashade_kwargs)

df = umap_panel.build_embedding_data_frame()
df["${hue}"] = df["${hue}"].astype('category')

unique_hues = df["${hue}"].unique()
hue_colors = hv.plotting.util.process_cmap(
    cmap='glasbey_hv',
    ncolors=len(unique_hues),
    categorical=True
)

hue_cmap = {label: color for label, color in
            zip(unique_hues, hue_colors)}

df['${hue}_color'] = df["${hue}"].map(hue_cmap).astype('category')

clusterer = hdbscan.HDBSCAN(
    algorithm='best', alpha=1.0, approx_min_span_tree=True,
    gen_min_span_tree=False, leaf_size=40,
    metric='euclidean', min_cluster_size=${min_cluster_size},
                                          min_samples = None, p = None)

clusterer.fit(df[['x', 'y']].values)

df['cluster_label'] = clusterer.labels_
df['cluster_label'] = df['cluster_label'].astype('category')

cluster_colors = hv.plotting.util.process_cmap(
    cmap='glasbey',
    ncolors=len(df['cluster_label'].unique()),
    categorical=True)

cluster_cmap = {label: color for label, color in
                zip(df['cluster_label'].unique(), cluster_colors)}

df['cluster_color'] = df.cluster_label.map(cluster_cmap)
df['cluster_color'] = df['cluster_color'].astype('category')

for (label, color_key) in [('cluster_label', cluster_cmap), ('${hue}', hue_cmap)]:
    export_cmap_legend(color_key,
                       f'umap_nn_${n_neighbors}_color_{label}_cluster_${hue}_rs_${random_state}_mcs_${min_cluster_size}.legend.png')

    canvas = datashader.Canvas(plot_height=$params.datashade_kwargs.height,
                                            plot_width =$params.datashade_kwargs.width)
    agg = canvas.points(df, 'x', 'y', agg=datashader.count_cat(label))

    image = datashader.transfer_functions.shade(agg, alpha=255, min_alpha=60, color_key=color_key)
    image = datashader.transfer_functions.set_background(image, 'black')

    image.to_pil().save(
        f'umap_nn_${n_neighbors}_color_{label}_cluster_${hue}_rs_${random_state}_mcs_${min_cluster_size}.png')

    image = datashader.transfer_functions.dynspread(image, threshold=0.75, max_px=1)
    image.to_pil().save(
        f'dynn_umap_nn_${n_neighbors}_color_{label}_cluster_${hue}_rs_${random_state}_mcs_${min_cluster_size}.png')
