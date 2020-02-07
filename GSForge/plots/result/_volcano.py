import numpy as np
import pandas as pd
import xarray as xr

import holoviews as hv

# hv.extension('bokeh', 'matplotlib')

default_options = [
    # Point options.
    hv.opts.Points(backend='bokeh', show_grid=True, invert_yaxis=True, padding=(0, 0.05), width=500, height=500),
    hv.opts.Points(backend='matplotlib', invert_yaxis=True, padding=(0.01, 0.05), fig_size=200, show_grid=True,
                   alpha=0.5),
    # Horizontal Line options.
    hv.opts.HLine(backend='bokeh', color="black", line_width=0.75, line_dash='dashed'),
    hv.opts.HLine(backend='matplotlib', linewidth=0.5, linestyle="dashdot", color="black"),
    # Vertical Line options.
    hv.opts.VLine(backend='bokeh', color="black", line_width=0.75, line_dash='dashed'),
    hv.opts.VLine(backend='matplotlib', color="black", linewidth=0.5, linestyle="dashdot"),
    # Label / annotation options.
    hv.opts.Labels(backend='bokeh', xoffset=0.6, yoffset=3, text_font_size='7pt'),
    hv.opts.Labels(backend='matplotlib', xoffset=0.6, yoffset=3, size=8),
    # General group options.
    hv.opts.Points("No Signal", color="black"),
    hv.opts.Points("LFC within pval", color="red"),
    hv.opts.Points("LFC outside pval", color="green"),
    hv.opts.Points("No LFC within pval", color="blue"),
]


def volcano(data: xr.Dataset, lfc_var: str, pval_var: str, lfc_cutoff: float = 2.0, pval_cutoff: float = 1e-6,
            label_selected_genes: bool = False, apply_default_opts: bool = True):
    """
    A volcano plot for examining differential gene expression results.

    :param data:
        Data containing the log-fold change and p-value variables.

    :param lfc_var:
        The name of the log-fold change column.

    :param pval_var:
        The name of the p-value column.

    :param lfc_cutoff:
        Cutoff to use in grouping and coloring genes. Defaults to 2.0.

    :param pval_cutoff:
        Cutoff to use in grouping and coloring genes. Defaults to 1e-6.

    :param label_selected_genes:
        Whether to apply annotations of genes that pass both the log-fold-change and p-value cutoff values.

    :param apply_default_opts:
        Whether to apply the default styling.

    :return:
        A scatter plot of log-fold-change versus -log10(p-values).
    """

    data = data.where(data[pval_var] > 0).dropna(dim="Gene")

    # Determine the group that each gene belongs to based on the cutoff values given.
    above_abs_lfd = np.abs(data[lfc_var]) >= lfc_cutoff
    below_pval_cutoff = data[pval_var] <= pval_cutoff

    gene_groups = {
        "No Signal": ~above_abs_lfd & ~below_pval_cutoff,
        "LFC within pval": above_abs_lfd & below_pval_cutoff,
        "LFC outside pval": above_abs_lfd & ~below_pval_cutoff,
        "No LFC within pval": ~above_abs_lfd & below_pval_cutoff,
    }

    df = pd.DataFrame({
        "Gene": data["Gene"].values,
        "lfc": data[lfc_var].values,
        "p-values": np.log10(data[pval_var].values),
    }).reset_index(drop=True)

    # Convert the gene groups into a single, label column.
    for key, selection in gene_groups.items():
        df.loc[selection.values, "Gene_group"] = key

    kdims = [("lfc", "$log_2$ fold change"), ("p-values", "$-log_{10}$ p-values")]
    vdims = ["Gene_group", "Gene"]

    groups = df.groupby("Gene_group").groups

    scatter_dict = {group: hv.Points(df.iloc[genes.values], kdims=kdims, vdims=vdims, group=group)
                    # .opts(color=color_map[group])
                    for group, genes in groups.items()}

    scatter_overlay = hv.NdOverlay(scatter_dict)

    if label_selected_genes:
        labels = hv.Labels(scatter_dict["LFC within pval"], kdims=kdims, vdims="Gene")
        scatter_overlay = scatter_overlay * labels

    layout = scatter_overlay * hv.HLine(np.log10(pval_cutoff)) * hv.VLine(lfc_cutoff) * hv.VLine(-lfc_cutoff)
    if apply_default_opts:
        return layout.opts(default_options)

    return layout