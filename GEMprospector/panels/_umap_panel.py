import param
from textwrap import dedent
import inspect
import panel as pn
import holoviews as hv
from bokeh.models import HoverTool
from methodtools import lru_cache

import umap

from ..models import Interface
from .helpers import generate_help_pane


class UMAP_Panel(Interface):
    """
    UMAP Panel short description.
    """

    # mappings = param.Dict(precedence=-1.0, doc=dedent("""\
    # A dictionary of Genes to be used in subseting the data."""))

    mapping_index_name = param.String(default="Gene", precedence=-1.0, doc=dedent("""\
    The index name to which the Genes in the `mappings` parameter refer."""))

    mapping_selector = param.ListSelector(default=["Complete"], doc=dedent("""\
    The active mapping to be used to subset the data."""))

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

    # Use a secret, undocumented, method to update the view with a button.
    # https://github.com/pyviz/panel/issues/389
    update_umap = param.Action(lambda self: self.param.trigger('update_umap'),
                               doc="Update panel view.")

    def __init__(self, *args, **params):
        super().__init__(*args, **params)
        if self.lineament_collection is not None:
            avail_mappings = list(self.lineament_collection.lineaments.keys())
            self.param["mapping_selector"].objects = avail_mappings + ["Complete"]

        # UMAP n_neighbors is limited by the number of samples.
        self.param["n_neighbors"].bounds = [1, len(self.gem.data[self.gem.sample_index_name]) - 1]

        avail_hues = sorted(list(set(self.gem.data.attrs.get('quantile')
                                     + self.gem.data.attrs.get('discrete'))))
        self.param["hue"].objects = [None] + avail_hues
        if self.hue is None:
            self.set_param(hue=avail_hues[0])

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
        mapping_key = list(frozen_mapping_key)
        if mapping_key == ["Complete"]:
            subset = self.gem.data
        else:
            gene_selection = list(self.lineament_collection.union(mapping_key))
            subset = self.gem.data.sel({self.mapping_index_name: gene_selection})

        zero_filled_data = subset[self.gem.count_array_name].fillna(0).values

        transform = umap.UMAP(**umap_kwargs).fit_transform(zero_filled_data)

        subset["x"] = ((self.gem.sample_index_name,), transform[:, 0])
        subset["y"] = ((self.gem.sample_index_name,), transform[:, 1])

        return subset

    @param.depends('update_umap', 'hue')
    def view(self):
        """Draw a 2d umap embedding."""
        frozen_umap_kwargs = frozenset(self._get_umap_kwargs().items())
        frozen_map_selector = frozenset(self.mapping_selector)
        umap_ds = self.transform(frozen_map_selector, frozen_umap_kwargs)

        plotting_dims = ['x', 'y'] + self.gem.data.attrs.get("all_labels")
        df = umap_ds[plotting_dims].to_dataframe()

        vdims = [item for item in self.gem.data.attrs.get("all_labels")]

        plot = hv.Points(df, kdims=["x", "y"], vdims=vdims)

        tooltips = [(name, "@" + f"{name}") for name in self.gem.data.attrs.get("all_labels")]
        hover = HoverTool(tooltips=tooltips)

        # plot = plot.opts(hv.opts.Points(tools=[hover]))
        plot = plot.opts(tools=[hover])
        color = self.hue if self.hue else "#3288bd"

        return plot.options(color=color, cmap="Set1", legend_position='bottom', xaxis=None, yaxis=None,
                            padding=0.2, show_grid=True, bgcolor="lightgrey", width=500, height=500)

    def panel(self):
        controls = pn.Param(self.param, widgets={
            'update_umap': {'type': pn.widgets.Button, 'button_type': 'primary'}, }
                            )

        umap_layout = pn.Column("### UMAP Exploration Panel",
                                pn.Row(self.view, controls))
        help_text = generate_help_pane(self)

        return pn.Tabs(("UMAP", umap_layout), ("Documentation", help_text))
