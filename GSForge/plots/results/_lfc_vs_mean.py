import xarray as xr
import numpy as np
import pandas as pd
import param

import holoviews as hv

from ..abstract_plot_models import ResultPlottingOperation


class MeanVsLFC(ResultPlottingOperation):
    """
    Mean vs log-fold-change scatter plot used for visualizing differential gene expression results.

    Parameters
    ----------
    source : Union[GSForge.GeneSet, xarray.Dataset, pandas.DataFrame]
        Data source containing the log-fold change, and p-value variables.

    log_fold_change_var : str
        The name of the log-fold change column. Must be a variable within `source`.

    p_value_var : str
        The name of the p-value column. Must be a variable within `source`.

    mean_value_var : str
        The name of the base mean column. Must be a variable within `source`.

    log_fold_change_cutoff : float
        Cutoff to use in grouping and coloring genes. Defaults to 2.0.

    p_value_cutoff : float
        Cutoff to use in grouping and coloring genes. Defaults to 1e-6.

    label_selected_genes : bool
        Apply (if True) annotations of genes that pass both the log-fold-change and p-value cutoff values.

    apply_default_opts : bool
        Whether to apply the default styling.

    Returns
    -------
    holoviews.Overlay
        A holoviews scatter plot of log-fold-change versus mean values.
    """

    log_fold_change_var = param.String(default=None, doc="""
    The name of the log-fold change column. Must be a variable within `source`.""")
    p_value_var = param.String(default=None, doc="""
    The name of the p-value column. Must be a variable within `source`.""")
    log_fold_change_cutoff = param.Number(default=2.0, doc="""
    Cutoff to use in grouping and coloring genes. Defaults to 2.0.""")
    p_value_cutoff = param.Number(default=1e-6, doc="""
    Cutoff to use in grouping and coloring genes. Defaults to 1e-6.""")
    mean_value_var = param.String(default=None, doc="""
    The name of the base mean column. Must be a variable within `source`.""")
    label_selected_genes = param.Boolean(default=False, doc="""
    Apply (if True) annotations of genes that pass both the log-fold-change and p-value cutoff values.""")

    @staticmethod
    def bokeh_opts():
        return [
            hv.opts.Points(backend='bokeh', show_grid=True, invert_yaxis=True, padding=(0, 0.05), width=500,
                           height=500),
            hv.opts.HLine(backend='bokeh', color="black", line_width=0.75, line_dash='dashed'),
            hv.opts.VLine(backend='bokeh', color="black", line_width=0.75, line_dash='dashed'),
            hv.opts.Labels(backend='bokeh', xoffset=0.6, yoffset=3, text_font_size='7pt'),
            hv.opts.Points("No Signal", color="black"),
            hv.opts.Points("LFC within pval", color="red"),
            hv.opts.Points("LFC outside pval", color="green"),
            hv.opts.Points("No LFC within pval", color="blue"),
        ]

    @staticmethod
    def matplotlib_opts():
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

    @staticmethod
    def mean_vs_lfc(source: xr.Dataset,
                    log_fold_change_var: str,
                    p_value_var: str,
                    mean_value_var: str,
                    log_fold_change_cutoff: float = 2.0,
                    p_value_cutoff: float = 1e-6,
                    label_selected_genes: bool = False,
                    gene_dim="Gene") -> hv.NdOverlay:
        # Omit missing data.
        data = source.where(source[p_value_var] > 0).dropna(dim=gene_dim)

        # Determine the group that each gene belongs to based on the cutoff values given.
        above_abs_lfd = np.abs(data[log_fold_change_var]) >= log_fold_change_cutoff
        below_pval_cutoff = data[p_value_var] <= p_value_cutoff

        gene_groups = {
            "No Signal": ~above_abs_lfd & ~below_pval_cutoff,
            "LFC within pval": above_abs_lfd & below_pval_cutoff,
            "LFC outside pval": above_abs_lfd & ~below_pval_cutoff,
            "No LFC within pval": ~above_abs_lfd & below_pval_cutoff,
        }

        df = pd.DataFrame({
            "Gene": data["Gene"].values,
            "lfc": data[log_fold_change_var].values,
            "mean": np.log10(data[mean_value_var].values),
            "p-values": np.log10(data[p_value_var].values),
        }).reset_index(drop=True)

        # Convert the gene groups into a single, label column.
        for key, selection in gene_groups.items():
            df.loc[selection.values, "Gene_group"] = key

        kdims = [("mean", "mean"), ("lfc", "log2 fold change")]
        vdims = ["Gene_group", "Gene"]

        groups = df.groupby("Gene_group").groups

        scatter_dict = {group: hv.Points(df.iloc[genes.values], kdims=kdims, vdims=vdims, group=group)
                        for group, genes in groups.items()}

        layout = hv.NdOverlay(scatter_dict) \
                 * hv.HLine(log_fold_change_cutoff) \
                 * hv.HLine(-log_fold_change_cutoff)

        if label_selected_genes:
            layout = layout * hv.Labels(scatter_dict["LFC within pval"], kdims=kdims, vdims="Gene")

        return layout

    def __call__(self, *args, **params):
        kwargs = {**self.infer_kwarg_defaults_from_data(self.source, self.mean_vs_lfc),
                  **self.get_param_process_overlap_kwargs(self.mean_vs_lfc)}
        plot = self.mean_vs_lfc(**kwargs)
        if self.apply_default_opts is True:
            plot = plot.opts(self.get_default_options())
        return plot
