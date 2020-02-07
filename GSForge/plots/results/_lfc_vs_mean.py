import xarray as xr
import numpy as np
import pandas as pd

import holoviews as hv

from ..utils import infer_default_kwargs


def bokeh_options():
    return [
        hv.opts.Points(backend='bokeh', show_grid=True, invert_yaxis=True, padding=(0, 0.05), width=500, height=500),
        hv.opts.HLine(backend='bokeh', color="black", line_width=0.75, line_dash='dashed'),
        hv.opts.VLine(backend='bokeh', color="black", line_width=0.75, line_dash='dashed'),
        hv.opts.Labels(backend='bokeh', xoffset=0.6, yoffset=3, text_font_size='7pt'),
        hv.opts.Points("No Signal", color="black"),
        hv.opts.Points("LFC within pval", color="red"),
        hv.opts.Points("LFC outside pval", color="green"),
        hv.opts.Points("No LFC within pval", color="blue"),
    ]


def matplotlib_options():
    return [
        hv.opts.Points(backend='matplotlib', invert_yaxis=True, padding=(0.01, 0.05), fig_size=200, show_grid=True,
                       alpha=0.5),
        hv.opts.HLine(backend='matplotlib', linewidth=0.5, linestyle="dashdot", color="black"),
        hv.opts.VLine(backend='matplotlib', color="black", linewidth=0.5, linestyle="dashdot"),
        hv.opts.Labels(backend='matplotlib', xoffset=0.6, yoffset=3, size=8),
        hv.opts.Points("No Signal", color="black"),
        hv.opts.Points("LFC within pval", color="red"),
        hv.opts.Points("LFC outside pval", color="green"),
        hv.opts.Points("No LFC within pval", color="blue"),
    ]


__options = {"bokeh": bokeh_options, "matplotlib": matplotlib_options}


@infer_default_kwargs
def mean_vs_lfc(data: xr.Dataset, lfc_var: str, pval_var: str, mean_var, lfc_cutoff: float = 2.0,
                pval_cutoff: float = 1e-6, label_genes: bool = False, backend: str = "bokeh",
                apply_default_opts: bool = True) -> hv.NdOverlay:
    # Omit missing data.
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
        "mean": np.log10(data[mean_var].values),
        "p-values": np.log10(data[pval_var].values),
    }).reset_index(drop=True)

    # Convert the gene groups into a single, label column.
    for key, selection in gene_groups.items():
        df.loc[selection.values, "Gene_group"] = key

    kdims = [("mean", "mean"), ("lfc", "log2 fold change")]
    vdims = ["Gene_group", "Gene"]

    groups = df.groupby("Gene_group").groups

    scatter_dict = {group: hv.Points(df.iloc[genes.values], kdims=kdims, vdims=vdims, group=group)
                    for group, genes in groups.items()}

    scatter = hv.NdOverlay(scatter_dict)

    layout = scatter * hv.HLine(lfc_cutoff) * hv.HLine(-lfc_cutoff)

    if label_genes:
        labels = hv.Labels(scatter_dict["LFC within pval"], kdims=kdims, vdims="Gene")
        layout = layout * labels

    if apply_default_opts:
        return layout.opts(__options[backend]())

    return layout
