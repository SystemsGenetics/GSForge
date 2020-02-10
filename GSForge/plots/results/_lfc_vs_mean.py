import xarray as xr
import numpy as np
import pandas as pd
import param

import holoviews as hv

from ..utils import ResultPlottingOperation, infer_kwarg_defaults_from_data, get_param_process_overlap_kwargs


class Mean_vs_LFC(ResultPlottingOperation):
    log_fold_change_var = param.String(default=None)
    p_value_var = param.String(default=None)

    log_fold_change_cutoff = param.Number(default=2.0)
    p_value_cutoff = param.Number(default=1e-6)
    mean_value_var = param.String(default=None)
    label_selected_genes = param.Boolean(default=False)

    @staticmethod
    def bokeh_options():
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

    def process(self):
        kwargs = {**infer_kwarg_defaults_from_data(self.source, self.mean_vs_lfc),
                  **get_param_process_overlap_kwargs(self, self.mean_vs_lfc)}
        return self.mean_vs_lfc(**kwargs)
