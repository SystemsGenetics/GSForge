import param
from textwrap import dedent
import inspect
import panel as pn
import holoviews as hv
from bokeh.models import HoverTool
from methodtools import lru_cache

import umap

from ..models import Interface
from ..utils import infer_xarray_variables
from .helpers import generate_help_pane


class UMAP_Panel(Interface):
    """
    UMAP Panel short description.
    """

    # UMAP model parameters.
    n_neighbors = param.Integer(default=15, doc=dedent("""\
    The size of local neighborhood (in terms of number of neighboring sample points) used 
    for manifold approximation. Larger values result in more global views of the manifold, 
    while smaller values result in more local data being preserved. In general values 
    should be in the range 2 to 100.\n
    Sourced from the UMAP documentation.\n
    This parameter is limited by the number of samples in a given dataset."""))

    n_components = param.Integer(default=2, precedence=-2, doc=dedent("""\
    The number of components (dimensions) to reduce to. Maybe one day we will go 3D.
    For now this should not be changed."""))

    min_dist = param.Number(default=0.1, bounds=[0, 1.0], doc=dedent("""\
    The effective minimum distance between embedded points. Smaller values will 
    result in a more clustered/clumped embedding where nearby points on the manifold 
    are drawn closer together, while larger values will result on a more even dispersal
    of points. The value should be set relative to the spread value, which determines 
    the scale at which embedded points will be spread out.\n UMAP Documentation."""))

    metric = param.ObjectSelector(
        default="manhattan", objects=["euclidean", "manhattan", "chebyshev", "minkowski"],
        doc=dedent("""The metric to use to compute distances in high dimensional space."""))
    random_state = param.Integer(default=42, doc="Random state seed.")

    # Display Parameters.
    hue = param.ObjectSelector(doc=dedent("""\
    Select the column by which the points drawn should be colored by."""), default=None)

    variable_categories = param.Dict(default=None, precedence=-2)

    # Use a secret, undocumented, method to update the view with a button.
    # https://github.com/pyviz/panel/issues/389
    update_umap = param.Action(lambda self: self.param.trigger('update_umap'),
                               doc="Update panel view.")

    def __init__(self, *args, **params):
        super().__init__(*args, **params)

        # UMAP n_neighbors is limited by the number of samples.
        self.param["n_neighbors"].bounds = [1, len(self.gem.data[self.gem.sample_index_name]) - 1]

        # Try to infer the variables of the data.
        if self.variable_categories is None:
            inferred_variables = infer_xarray_variables(self.gem.data,
                                                        skip=[self.gem.gene_index_name] + self.gem.count_array_names)
            self.set_param(variable_categories=inferred_variables)

        avail_hues = []
        if self.variable_categories.get("discrete") is not None:
            avail_hues += self.variable_categories.get("discrete")
        if self.variable_categories.get("quantile") is not None:
            avail_hues += self.variable_categories.get("quantile")

        avail_hues = list(set(avail_hues))
        if avail_hues:
            self.param["hue"].objects = [None] + sorted(avail_hues)

    def _get_umap_kwargs(self):
        """Returns the arguments needed for the UMAP reduction."""
        key_set = set(inspect.signature(umap.umap_.UMAP).parameters.keys()
                      ).intersection(set(self.param.objects().keys()))
        return {key: getattr(self, key) for key in key_set}

    @lru_cache()
    def transform(self, frozen_mapping_key, frozen_umap_kwargs):
        """A cached UMAP transform function. This lets the user to update the plot quickly
        for those parameters that do not require a new UMAP calculation."""
        # The arguments come in as frozen sets so that the lru_cache works correctly.
        # Convert them to usable forms.
        umap_kwargs = dict(frozen_umap_kwargs)
        _ = list(frozen_mapping_key)  # Just used to hash the result of the data key.
        subset = self.selection

        zero_filled_data = subset[self.active_count_variable].fillna(0).values
        transform = umap.UMAP(**umap_kwargs).fit_transform(zero_filled_data)

        subset["x"] = ((self.gem.sample_index_name,), transform[:, 0])
        subset["y"] = ((self.gem.sample_index_name,), transform[:, 1])

        return subset

    @param.depends('update_umap', 'hue')
    def view(self):
        """Draw a 2d umap embedding."""
        frozen_umap_kwargs = frozenset(self._get_umap_kwargs().items())
        frozen_map_selector = frozenset((self.selected_gene_sets + [self.gene_set_mode]))
        umap_ds = self.transform(frozen_map_selector, frozen_umap_kwargs)
        plotting_dims = ['x', 'y']

        if self.variable_categories.get("all_labels") is not None:
            plotting_dims += self.variable_categories.get("all_labels")

        df = umap_ds[plotting_dims].to_dataframe().reindex()

        # Set quantileable or group categories to the 'string' datatype.
        # for var_type in ["discrete", "quantile"]:
        #     if self.variable_categories.get(var_type) is not None:
        #         df[self.variable_categories.get(var_type)] = df[self.variable_categories.get(var_type)].astype("str")

        vdims = [item for item in self.variable_categories.get("all_labels")
                 if item in list(df.columns)]

        plot = hv.Points(df, kdims=["x", "y"], vdims=vdims)

        tooltips = [(name, "@" + f"{name}") for name in vdims
                    if name in list(df.columns)]

        hover = HoverTool(tooltips=tooltips)

        # plot = plot.opts(hv.opts.Points(tools=[hover]))
        plot = plot.opts(tools=[hover])
        color = self.hue if self.hue else "#3288bd"

        return plot.options(color=color, cmap="Set1", legend_position='bottom', xaxis=None, yaxis=None,
                            padding=0.1, show_grid=True, bgcolor="lightgrey", width=500, height=500)

    def panel(self):
        controls = pn.Param(self.param, widgets={
            'update_umap': {'type': pn.widgets.Button, 'button_type': 'primary'}, }
                            )

        umap_layout = pn.Column("### UMAP Exploration Panel",
                                pn.Row(self.view, controls))
        help_text = generate_help_pane(self)

        return pn.Tabs(("UMAP", umap_layout), ("Documentation", help_text))
