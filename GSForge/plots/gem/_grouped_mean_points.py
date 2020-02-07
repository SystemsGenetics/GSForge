from ...models import OperationInterface
# from GSForge.models import OperationInterface
import param
import holoviews as hv


class grouped_mean_scatter_operation(OperationInterface):
    group_variable = param.String()
    x_group_label = param.String()
    y_group_label = param.String()

    apply_options = param.Boolean(default=True)  # Eventually a dict or a hook list.

    grouped_data = param.Parameter(precedence=-1)

    @staticmethod
    def grouped_mean_scatter(count_xarray, labels, group_variable, x_group_label, y_group_label):
        groups = labels.groupby(group_variable).groups
        x_group_mean = count_xarray.isel({"Sample": groups[x_group_label]}).mean(dim="Sample")
        y_group_mean = count_xarray.isel({"Sample": groups[y_group_label]}).mean(dim="Sample")
        return hv.Points((x_group_mean, y_group_mean), [x_group_label, y_group_label], "Gene")

    @staticmethod
    def bokeh_options():
        return [
            hv.opts.Points(padding=0.05),
            hv.opts.Points(backend="matplotlib", fig_size=200, padding=0.05, s=5, show_grid=True),
            hv.opts.Points(backend="bokeh", width=500, height=500, show_grid=True),
        ]

    @staticmethod
    def matplotlib_options():
        return [
            hv.opts.Points(padding=0.05),
            hv.opts.Points(backend="matplotlib", fig_size=200, padding=0.05, s=5, show_grid=True),
            hv.opts.Points(backend="bokeh", width=500, height=500, show_grid=True),
        ]

    def process(self):
        self.set_param(annotation_variables=[self.group_variable])
        plot = self.grouped_mean_scatter(self.x_count_data, self.y_annotation_data, self.group_variable,
                                         self.x_group_label, self.y_group_label)
        if self.apply_options:
            return self.apply_options(plot)
        return plot