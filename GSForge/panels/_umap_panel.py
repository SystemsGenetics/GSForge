import param
from textwrap import dedent
import inspect
import panel as pn
import holoviews as hv
from bokeh.models import HoverTool
from methodtools import lru_cache

import umap

from GSForge import Interface
from GSForge.utils._panel_utils import generate_help_pane


class UMAP_Panel(param.Parameterized):
    """A UMAP Panel Exploration Tool.

    """

    interface = param.Parameter(
        precedence=-1.0,
        doc="An instance of a GSForge.Interface object.",
    )

    data_var_cats = param.Dict(
        precedence=-1.0,
        doc="Categories of variables. Controls which values are selectable as a 'hue'.",
    )
    ###############################################################################################
    # Transform Parameters.
    ###############################################################################################

    n_neighbors = param.Integer(
        precedence=1.0,
        default=15,
        doc=dedent("""\
        The size of local neighborhood (in terms of number of neighboring sample points) used 
        for manifold approximation. Larger values result in more global views of the manifold, 
        while smaller values result in more local data being preserved. In general values 
        should be in the range 2 to 100.\n  
        *This parameter is limited by the number of samples in a given dataset.*""")
    )

    n_components = param.Integer(
        default=2,
        precedence=-2,
        doc=dedent("""\
        The number of components (dimensions) to reduce to. Maybe one day we will go 3D.
        For now this should not be changed.""")
    )

    min_dist = param.Number(
        precedence=1.0,
        default=0.1,
        bounds=[0, 1.0],
        doc=dedent("""\
        The effective minimum distance between embedded points. Smaller values will 
        result in a more clustered/clumped embedding where nearby points on the manifold 
        are drawn closer together, while larger values will result on a more even dispersal
        of points. The value should be set relative to the spread value, which determines 
        the scale at which embedded points will be spread out.""")
    )

    metric = param.ObjectSelector(
        precedence=1.0,
        default="manhattan",
        objects=["euclidean", "manhattan", "chebyshev", "minkowski"],
        doc=dedent("""The metric to use to compute distances in high dimensional space."""))

    random_state = param.Integer(
        precedence=1.0,
        default=None,
        doc="Random state seed. Set for (exactly) reproducible plots."
    )

    ###############################################################################################
    # Plot display parameters.
    ###############################################################################################

    hue = param.ObjectSelector(
        precedence=1.0,
        doc="Select the column by which the points drawn should be colored by.",
        default=None
    )

    # Create a button so that the transform only occurs when clicked.
    update = param.Action(lambda self: self.param.trigger('update'))

    def __init__(self, source, interface_opts=None, **params):
        # Set up the Interface object.
        default_interface_opts = {"gene_set_mode": "complete",
                                  "count_mask": "complete"}
        if interface_opts is not None:
            interface_opts = {**default_interface_opts, **interface_opts}
        else:
            interface_opts = default_interface_opts

        interface = Interface(source, **interface_opts)
        super().__init__(interface=interface, **params)

        # Infer the variable categories of the supplied dataframe.
        if self.data_var_cats is None:
            self.set_param(data_var_cats=self.interface.gem.infer_variables())

        # Limit the number of neighbors to the number of samples.
        self.param["n_neighbors"].bounds = [1, len(self.interface.gem.data[self.interface.gem.sample_index_name]) - 1]

        avail_hues = []
        if self.data_var_cats.get("discrete") is not None:
            avail_hues += self.data_var_cats.get("discrete")
        if self.data_var_cats.get("quantile") is not None:
            avail_hues += self.data_var_cats.get("quantile")

        avail_hues = list(set(avail_hues))
        if avail_hues:
            self.param["hue"].objects = [None] + sorted(avail_hues)

    def get_transform_kwargs(self, transform=umap.umap_.UMAP):
        """Gets the overlapping arguments of the transform and parameters of this panel class,
        and returns them as a dictionary."""
        key_set = set(inspect.signature(transform).parameters.keys()
                      ).intersection(set(self.param.objects().keys()))
        return {key: getattr(self, key) for key in key_set}

    @staticmethod
    def static_transform(array, **kwargs):
        """Runs the transform with the selected panel parameters on the given array."""
        return umap.UMAP(**kwargs).fit_transform(array)

    def transform(self):
        """Runs the transform with the selected panel parameters and data."""
        array = self.interface.x_count_data
        kwargs = self.get_transform_kwargs()
        return self.static_transform(array, **kwargs)

    @lru_cache()
    def cached_transform(self, gene_set, transform_state):
        """A chached transform based on the genes and transform arguments selected."""
        return self.transform()

    @param.depends('update', 'hue')
    def hvplot_view(self):
        """A holoviews.Points plot of the selected transform."""
        df = self.interface.gem.data[self.data_var_cats["all_labels"]].to_dataframe().reset_index()
        gene_set = frozenset(self.interface.get_gene_index())
        transform_state = frozenset(self.get_transform_kwargs().items())
        hover = HoverTool(tooltips=[(name, "@" + f"{name}") for name in list(df.columns)])

        transform = self.cached_transform(gene_set, transform_state)
        df["x"] = transform[:, 0]
        df["y"] = transform[:, 1]

        plot = hv.Points(df, kdims=["x", "y"])

        hue = self.hue if self.hue else "#3288bd"

        return plot.opts(tools=[hover], color=hue, cmap="Set1", legend_position='bottom', xaxis=None, yaxis=None,
                         padding=0.05, show_grid=True, bgcolor="lightgrey", width=500, height=500)

    def panel(self):
        """Interactive panel application for transform exploration."""
        # Update styles or types of buttons for the interface.
        transform_controls = pn.Param(self.param, widgets={
            'update': {'type': pn.widgets.Button, 'button_type': 'primary'},
        })

        transform_layout = pn.Row(self.interface, self.hvplot_view, transform_controls)
        docs = generate_help_pane(self)

        tab_layout = pn.Tabs(("UMAP", transform_layout), ("Documentation", docs))

        return tab_layout
