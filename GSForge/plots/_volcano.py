import numpy as np
import pandas as pd
import xarray as xr

import param
import panel as pn
import holoviews as hv


def volcano(data: xr.Dataset, lfc_var: str, pval_var: str, lfc_cutoff: float = 2.0, pval_cutoff: float = 1e-6):
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

    :return:
        A scatter plot of log-fold-change versus -log10(p-values).
    """

    data = data.where(data[pval_var] > 0).dropna(dim="Gene")

    # Determine the group that each gene belongs to based on the cutoff values given.
    above_abs_lfd = np.abs(data[lfc_var]) >= lfc_cutoff
    below_pval_cutoff = data[pval_var] <= pval_cutoff

    gene_groups = {
        "No Signal": above_abs_lfd & below_pval_cutoff,
        "LFC within pval": above_abs_lfd & ~below_pval_cutoff,
        "LFC outside pval": ~above_abs_lfd & below_pval_cutoff,
        "No LFC within pval": ~above_abs_lfd & ~below_pval_cutoff,
    }

    df = pd.DataFrame({
        "Gene": data["Gene"].values,
        "lfc": data[lfc_var].values,
        "p-values": np.log10(data[pval_var].values),
    }).reset_index(drop=True)

    # Convert the gene groups into a single, label column.
    for key, selection in gene_groups.items():
        df.loc[selection.values, "Gene_group"] = key

    color_map = {
        "No Signal": "black",
        "LFC within pval": "red",
        "LFC outside pval": "green",
        "No LFC within pval": "blue",
    }

    kdims = [("lfc", "$log_2$ fold change"), ("pvalues", "$-log_{10}$ p-values")]
    vdims = ["Gene_group", "Gene"]

    groups = df.groupby("Gene_group").groups

    scatter_dict = {group: hv.Points(df.iloc[genes], kdims=kdims, vdims=vdims, label=group)  #.opts(color=color_map[group])
                    for group, genes in groups.items()}

    labels = hv.Labels(scatter_dict["LFC within pval"], kdims=kdims, vdims="Gene")

    scatter = hv.NdOverlay(scatter_dict) * labels

    layout = scatter * hv.HLine(np.log10(pval_cutoff)) * hv.VLine(lfc_cutoff) * hv.VLine(-lfc_cutoff)

    return layout
