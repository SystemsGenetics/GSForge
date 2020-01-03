import param
from textwrap import dedent
import inspect
import panel as pn
import holoviews as hv
from bokeh.models import HoverTool
from methodtools import lru_cache
from sklearn.decomposition import PCA

from ..models import Interface


class PCA_Panel(Interface):
    """
    For lineamenet collections:
        By default all keys will be made available in the selector.
        By default the mode should be complete(?) i.e. no lineaments.

    """

    mapping_selector = param.ListSelector(default=["Complete"], doc=dedent("""\
    The active mapping to be used to subset the data."""))

    hue = param.ObjectSelector(doc=dedent("""\
    Select the column by which the points drawn should be colored by."""), default=None)

    variable_categories = param.Dict(default=None)

    update_view = param.Action(lambda self: self.param.trigger('update_view'),
                               doc="Update panel view.")

    def __init__(self, *args, **params):
        super().__init__(*args, **params)

        # Mapping_selector options must include any lineaments.
        if self.lineament_collection is not None:
            lineament_names = list(self.lineament_collection.lineaments.keys())
            self.param["mapping_selector"].objects = lineament_names + ["Complete"]

        # Try to infer the variables of the data.
        if self.variable_categories is None:
            inferred_variables = gsf.utils.infer_xarray_variables(self.gem.data)
            self.set_param(variable_categories=inferred_variables)

        # If any lineament_keys are provided, they should be selected by default.
        if self.lineament_keys is not None:
            self.set_param(mapping_selector=list(self.lineament_keys))

        avail_hues = []
        if self.variable_categories.get("discrete") is not None:
            avail_hues += self.variable_categories.get("discrete")
        if self.variable_categories.get("quantile") is not None:
            avail_hues += self.variable_categories.get("quantile")

        if avail_hues:
            self.param["hue"].objects = [None] + avail_hues

    def _get_transform_kwargs(self):
        key_set = set(inspect.signature(PCA).parameters.keys()
                      ).intersection(set(self.param.objects().keys()))
        return {key: getattr(self, key) for key in key_set}

    @lru_cache()
    def transform(self, frozen_mapping_keys, frozen_transform_kwargs):
        transform_kwargs = dict(frozen_transform_kwargs)
        mapping_keys = list(frozen_mapping_keys)

        if self.mapping_selector == ["Complete"]:
            subset = self.gem.data
        else:
            gene_selection = list(self.lineament_collection.union(mapping_keys))
            subset = self.gem.data.sel({self.mapping_index_name: gene_selection})

        zero_filled_data = subset[self.gem.count_array_name].fillna(0).values

        transform = PCA(**transform_kwargs).fit_transform(zero_filled_data)

        subset["x"] = ((self.gem.sample_index_name,), transform[:, 0])
        subset["y"] = ((self.gem.sample_index_name,), transform[:, 1])

        return subset

    @param.depends('update_view', 'hue')
    def view(self):
        frozen_umap_kwargs = frozenset(self._get_transform_kwargs().items())
        frozen_map_selector = frozenset(self.mapping_selector)

        transform_ds = self.transform(frozen_map_selector, frozen_umap_kwargs)

        plotting_dims = ['x', 'y'] + [self.sample_index_name] + [self.gene_index_name]
        vdims = [self.sample_index_name, self.gene_index_name]
        tooltips = None

        if self.variable_categories.get("all_labels") is not None:
            plotting_dims += self.variable_categories.get("all_labels")
            vdims += [item for item in self.variable_categories.get("all_labels")]
            tooltips = [(name, "@" + f"{name}") for name in self.variable_categories.get("all_labels")]
            hover = HoverTool(tooltips=tooltips)

        df = transform_ds[plotting_dims].to_dataframe().reset_index()

        plot = hv.Points(df, kdims=["x", "y"], vdims=vdims)

        if hover:
            plot = plot.opts(tools=[hover])

        color = self.hue if self.hue else "#3288bd"

        return plot.options(color=color, cmap="Set1", legend_position='bottom', xaxis=None, yaxis=None,
                            padding=0.1, show_grid=True, bgcolor="lightgrey", width=500, height=500)

    def panel(self):
        controls = pn.Param(self.param, widgets={
            'update_view': {'type': pn.widgets.Button, 'button_type': 'primary'}})
        layout = pn.Column("### PCA Exploration Panel",
                           pn.Row(self.view, controls))
        return layout