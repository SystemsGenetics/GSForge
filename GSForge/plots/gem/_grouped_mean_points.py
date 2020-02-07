from ...models import OperationInterface
# from GSForge.models import OperationInterface
import param
import holoviews as hv
from textwrap import dedent


class grouped_mean_scatter_operation(OperationInterface):
    group_variable = param.String()
    x_group_label = param.String()
    y_group_label = param.String()

    apply_default_opts = param.Boolean(default=True, precedence=-1.0, doc=dedent("""\
    Whether to apply the default styling based on the current backend."""))

    grouped_data = param.Parameter(precedence=-1)

    backend = param.ObjectSelector(default="bokeh", objects=["bokeh", "matplotlib"], doc=dedent("""\
    The selected plotting backend to use for display. Options are ["bokeh", "matplotlib"]."""))

    @staticmethod
    def grouped_mean_scatter(count_xarray, labels, group_variable, x_group_label, y_group_label):
        groups = labels.groupby(group_variable).groups
        x_group_mean = count_xarray.isel({"Sample": groups[x_group_label]}).mean(dim="Sample")
        y_group_mean = count_xarray.isel({"Sample": groups[y_group_label]}).mean(dim="Sample")
        return hv.Points((x_group_mean, y_group_mean), [x_group_label, y_group_label], "Gene")

    @staticmethod
    def bokeh_options():
        return hv.opts.Points(backend="bokeh", padding=0.05, width=500, height=500, show_grid=True)


    @staticmethod
    def matplotlib_options():
        return hv.opts.Points(backend="matplotlib",padding=0.05, fig_size=200, s=5, show_grid=True)

    def process(self):
        self.set_param(annotation_variables=[self.group_variable])
        plot = self.grouped_mean_scatter(self.x_count_data, self.y_annotation_data, self.group_variable,
                                         self.x_group_label, self.y_group_label)
        if self.apply_default_opts:
            options = {"bokeh": self.bokeh_options, "matplotlib": self.matplotlib_options}
            default_options = options[self.backend]()
            return plot.opts(default_options)

        return plot
