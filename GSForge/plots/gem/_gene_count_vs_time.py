


class GeneCountOverTime(Interface):
    """
    For each treatment group:
        + Curve: mean / median / trend
        + Points: given by count_variable.
        + Area: log fold change
        + Area: variance, error.

    """

    selected_gene = param.Parameter()
    time_variable = param.Parameter()
    treatment_variable = param.Parameter()
    dispersion_variable = param.Parameter()
    log_fold_change = param.Parameter(default=1.0)
    treatment_colormap = param.Parameter()
    log2_output = param.Boolean(default=True)

    @staticmethod
    def create_plot_dataframe(
            dataset,
            selected_gene,
            time_var,
            treatment_var,
            dispersion_var,
            count_var,
            gene_var,
    ):
        active_variables = list(filter(None, [time_var, treatment_var, dispersion_var, count_var]))
        gene_subset = dataset.sel({gene_var: selected_gene})[active_variables]
        df = gene_subset.to_dataframe().reset_index()
        return df.copy(deep=True)

    def get_plot_dataframe(self):
        return self.create_plot_dataframe(
            dataset=self.gem.data,
            selected_gene=self.selected_gene,
            time_var=self.time_variable,
            treatment_var=self.treatment_variable,
            dispersion_var=self.dispersion_variable,
            count_var=self.count_variable,
            gene_var=self.gem.gene_index_name,
        )

    @staticmethod
    def create_spread_dataframe(points_dataframe, time_var,
                                log_fold_spread=1.0,
                                count_var='counts', treatment_var=None, dispersion_var=None):
        group_vars = list(filter(None, [treatment_var, time_var]))
        spread_df = points_dataframe.groupby(group_vars).mean().reset_index()

        spread_df['lfc'] = log_fold_spread

        spread_df['variance_spread'] = spread_df[dispersion_var] * spread_df[count_var]

        if dispersion_var is not None:
            spread_df[f'log2_{dispersion_var}'] = spread_df[dispersion_var] / np.log(2)

        return spread_df

    def get_spread_dataframe(self):
        return self.create_spread_dataframe(self.get_plot_dataframe(),
                                            time_var=self.time_variable,
                                            treatment_var=self.treatment_variable,
                                            dispersion_var=self.dispersion_variable
                                            )

    #     def plot_as_scatter(self):
    #         pdf = self.get_plot_dataframe()
    #         kdims = [self.time_variable, self.count_variable]
    #         vdims = self.treatment_variable if self.treatment_variable else None
    #         scatter = hv.Scatter(pdf, kdims, vdims)

    #         return scatter.opts(self.bokeh_opts)

    def plot_as_scatter(self):
        pdf = self.get_plot_dataframe()
        kdims = [self.time_variable, self.count_variable]
        vdims = self.treatment_variable if self.treatment_variable else None

        if self.treatment_variable:
            overlay = hv.NdOverlay(kdims=[self.treatment_variable])

            for treatment, group in pdf.groupby(self.treatment_variable):
                scatter = hv.Scatter(group, kdims, vdims, group='Counts', label=treatment)

                if self.treatment_colormap:
                    scatter.opts(color=self.treatment_colormap[treatment])

                overlay[treatment] = scatter

            return overlay.opts(self.bokeh_opts)

        else:
            return hv.Scatter(pdf, kdims, vdims).opts(self.bokeh_opts)

    def plot_mean_curve(self):
        sdf = self.get_spread_dataframe()
        kdims = [self.time_variable]
        vdims = list(filter(None, [self.count_variable, self.treatment_variable]))

        if self.treatment_variable:
            curve_overlay = hv.NdOverlay(kdims=[self.treatment_variable])

            for treatment, group in sdf.groupby(self.treatment_variable):
                curve = hv.Curve(group, kdims, vdims, group='Mean', label=treatment)

                if self.treatment_colormap:
                    curve.opts(color=self.treatment_colormap[treatment])

                curve_overlay[treatment] = curve

            return curve_overlay.opts(self.bokeh_opts)

        else:
            return hv.Curve(sdf, kdims, vdims).opts(self.bokeh_opts)

    def plot_dispersion_spread(self):

        sdf = self.get_spread_dataframe()
        kdims = [self.time_variable, self.count_variable]
        vdims = ['variance_spread']

        if self.treatment_variable:
            overlay = hv.NdOverlay(kdims=[self.treatment_variable])

            for treatment in sdf[self.treatment_variable].unique():

                group = sdf[sdf[self.treatment_variable] == treatment]
                spread = hv.Spread(group, kdims, vdims, group='Dispersion', label=treatment)

                if self.treatment_colormap:
                    spread.opts(color=self.treatment_colormap[treatment])

                overlay[treatment] = spread

            return overlay.opts(self.bokeh_opts)

        else:
            return hv.Spread(sdf, kdims, vdims).opts(self.bokeh_opts)

    @property
    def bokeh_opts(self):
        options = [
            hv.opts.Scatter(width=600, height=400, legend_position='right'),
            hv.opts.NdOverlay(width=600, height=400, legend_position='right')
        ]

        #         if self.treatment_variable:
        #             options.append(hv.opts.Scatter(color=self.treatment_variable, cmap=self.treatment_colormap))

        return options

    def process(self):
        return

# tdf = GeneCountOverTime(
#     processed_collection,
#     selected_gene='TraesCS5D02G318700.1',
#     time_variable='Exposed_Hours',
#     treatment_variable='Bio_Sample_ID',
#     dispersion_variable='model_1_edgeR_trended_dispersion',
# ).get_spread_dataframe()