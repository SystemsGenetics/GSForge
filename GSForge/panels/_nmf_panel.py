# import param
# import panel as pn
# import inspect
# import holoviews as hv
# from bokeh.models import HoverTool
# from methodtools import lru_cache
#
# from sklearn.decomposition import NMF
#
# from textwrap import dedent
#
#
# from GSForge.utils._panel_utils import generate_help_pane
#
#
# class NMF_Panel(param.Parameterized):
#     """A Negative Matrix Factorization Panel Exploration Tool.
#
#     """
#
#     interface = param.Parameter(
#         precedence=-1.0,
#         doc="An instance of a GSForge.Interface object.",
#     )
#
#     data_var_cats = param.Dict(
#         precedence=-1.0,
#         doc="Categories of variables. Controls which values are selectable as a 'hue'.",
#     )
#     ###############################################################################################
#     # Transform Parameters.
#     ###############################################################################################
#
#     n_components = param.Integer(
#         default=2,
#         precedence=-2,
#         doc=dedent("""\
#         The number of components (dimensions) to reduce to. Maybe one day we will go 3D.
#         For now this should not be changed.""")
#     )
#
#     random_state = param.Integer(
#         precedence=1.0,
#         default=None,
#         doc="Random state seed. Set for (exactly) reproducible plots."
#     )
#
#     ###############################################################################################
#     # Plot display parameters.
#     ###############################################################################################
#
#     hue = param.ObjectSelector(
#         precedence=1.0,
#         doc="Select the column by which the points drawn should be colored by.",
#         default=None
#     )
#
#     # Create a button so that the transform only occurs when clicked.
#     update = param.Action(lambda self: self.param.trigger('update'))
#
#     def __init__(self, *args, **params):
#         super().__init__(*args, **params)
#
#         # Infer the variable categoires of the supplied dataframe.
#         if self.data_var_cats is None:
#             self.set_param(data_var_cats=self.interface.gem.infer_variables())
#
#         avail_hues = []
#         if self.data_var_cats.get("discrete") is not None:
#             avail_hues += self.data_var_cats.get("discrete")
#         if self.data_var_cats.get("quantile") is not None:
#             avail_hues += self.data_var_cats.get("quantile")
#
#         avail_hues = list(set(avail_hues))
#         if avail_hues:
#             self.param["hue"].objects = [None] + sorted(avail_hues)
#
#     def get_transform_kwargs(self, transform=NMF):
#         """Gets the overlapping arguments of the transform and parameters of this panel class,
#         and returns them as a dictionary."""
#         key_set = set(inspect.signature(transform).parameters.keys()
#                       ).intersection(set(self.param.objects().keys()))
#         return {key: getattr(self, key) for key in key_set}
#
#     @staticmethod
#     def static_transform(array, **kwargs):
#         """Runs the transform on the given array."""
#         return NMF(**kwargs).fit_transform(array)
#
#     def transform(self):
#         """Runs the transform with the selected panel parameters and data."""
#         array = self.interface.x_count_data
#         kwargs = self.get_transform_kwargs()
#         return self.static_transform(array, **kwargs)
#
#     @lru_cache()
#     def cached_transform(self, gene_set, transform_state):
#         """A chached transform based on the genes and transform arguments selected."""
#         return self.transform()
#
#     @param.depends('update', 'hue')
#     def view(self):
#         """A holoviews.Points plot of the selected transform."""
#         df = self.interface.gem.data[self.data_var_cats["all_labels"]].to_dataframe().reset_index()
#         gene_set = frozenset(self.interface.get_gene_index())
#         transform_state = frozenset(self.get_transform_kwargs().items())
#         hover = HoverTool(tooltips=[(name, "@" + f"{name}") for name in list(df.columns)])
#
#         transform = self.cached_transform(gene_set, transform_state)
#         df["x"] = transform[:, 0]
#         df["y"] = transform[:, 1]
#
#         plot = hv.Points(df, kdims=["x", "y"])
#
#         hue = self.hue if self.hue else "#3288bd"
#
#         return plot.opts(tools=[hover], color=hue, cmap="Set1", legend_position='bottom', xaxis=None, yaxis=None,
#                          padding=0.05, show_grid=True, bgcolor="lightgrey", width=500, height=500)
#
#     def panel(self):
#         """Interactive panel application for transform exploration."""
#         # Update styles or types of buttons for the interface.
#         transform_controls = pn.Param(self.param, widgets={
#             'update': {'type': pn.widgets.Button, 'button_type': 'primary'},
#         })
#
#         transform_layout = pn.Row(self.interface, self.view, transform_controls)
#         docs = generate_help_pane(self)
#
#         tab_layout = pn.Tabs(("UMAP", transform_layout), ("Documentation", docs))
#
#         return tab_layout
