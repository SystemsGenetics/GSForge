import param
import holoviews as hv

from ...models import Interface
from ..utils import AbstractPlottingOperation


class GroupedGeneCovariance(AbstractPlottingOperation, Interface, param.ParameterizedFunction):
    # TODO: Document me.
    group_variable = param.String()
    x_group_label = param.String()
    y_group_label = param.String()
    grouped_data = param.Parameter(precedence=-1)

    @staticmethod
    def grouped_mean_scatter(count_xarray, labels, group_variable, x_group_label, y_group_label):
        groups = labels.groupby(group_variable).groups
        x_group_mean = count_xarray.sel({"Sample": groups[x_group_label]}).mean(dim="Sample")
        y_group_mean = count_xarray.sel({"Sample": groups[y_group_label]}).mean(dim="Sample")
        return hv.Points((x_group_mean, y_group_mean), [x_group_label, y_group_label], "Gene")

    @staticmethod
    def bokeh_opts():
        return hv.opts.Points(backend="bokeh", padding=0.05, width=500, height=500, show_grid=True)

    @staticmethod
    def matplotlib_opts():
        return hv.opts.Points(backend="matplotlib", padding=0.05, fig_size=200, s=5, show_grid=True)

    def process(self):
        self.set_param(annotation_variables=[self.group_variable])
        plot = self.grouped_mean_scatter(self.x_count_data, self.y_annotation_data.to_dataframe(),
                                         self.group_variable, self.x_group_label, self.y_group_label)

        return plot
