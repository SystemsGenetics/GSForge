import inspect
import logging

import holoviews as hv
import holoviews.plotting
import pandas as pd
import panel as pn
import param
import umap
from bokeh.models import HoverTool
from datashader.reductions import count_cat
from holoviews.operation import datashader
from methodtools import lru_cache

from .utils import generate_help_pane
from ..models._Interface import Interface

logger = logging.getLogger("GSForge")


# TODO: Allow size selection of the points drawn.
# TODO: Disable initial transform.
# TODO: Disable selected gene sets when only an AnnotatedGEM is supplied.
class UMAP_Interface(Interface):
    """A UMAP Panel Exploration Tool."""

    annotation_categories = param.Dict(
        precedence=-1.0,
        doc="Categories of variables. Controls which values are selectable as a 'hue'.",
    )

    ###############################################################################################
    # Transform Parameters.
    ###############################################################################################

    n_neighbors = param.Integer(
        precedence=1.0,
        default=15,
        doc="""\
        The size of local neighborhood (in terms of number of neighboring sample points) used 
        for manifold approximation. Larger values result in more global views of the manifold, 
        while smaller values result in more local data being preserved. In general values 
        should be in the range 2 to 100.\n  
        *This parameter is limited by the number of samples in a given dataset.*"""
    )

    n_components = param.Integer(
        default=2,
        precedence=-2,
        doc="""\
        The number of components (dimensions) to reduce to. Maybe one day we will go 3D.
        For now this should not be changed."""
    )

    min_dist = param.Number(
        precedence=1.0,
        default=0.1,
        bounds=[0, 1.0],
        doc="""\
        The effective minimum distance between embedded points. Smaller values will 
        result in a more clustered/clumped embedding where nearby points on the manifold 
        are drawn closer together, while larger values will result on a more even dispersal
        of points. The value should be set relative to the spread value, which determines 
        the scale at which embedded points will be spread out."""
    )

    metric = param.ObjectSelector(
        precedence=1.0,
        default="manhattan",
        objects=["euclidean", "manhattan", "chebyshev", "minkowski"],
        doc="""The metric to use to compute distances in high dimensional space.""")

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

    apply_default_opts = param.Boolean(default=True, precedence=-1.0, doc="""
        Whether to apply the default styling based on the current backend.""")

    plot_options = param.Parameter(default=None, doc="""
    User supplied options to the plotting functions. If provided (and is not None), these will take
    precedence over a functions built-in defaults.""")

    datashade = param.Boolean(default=None)
    dynspread = param.Boolean(default=None)
    datashade_kwargs = param.Dict(default=None, precedence=-1.0)

    # Create a button so that the transform only occurs when clicked.
    update = param.Action(lambda self: self.param.trigger('update'))

    @staticmethod
    def bokeh_opts():
        return [
            hv.opts.Points(cmap="glasbey_cool", legend_position='top_right', axiswise=True, xaxis=None, yaxis=None,
                           padding=0.05, show_grid=True, bgcolor="lightgrey", width=600, height=600),
            hv.opts.RGB(bgcolor="black", xaxis=None, yaxis=None, width=600, height=600, backend='bokeh'),
        ]

    # @staticmethod
    # def datashader_opts():
    #     return [
    #         hv.opts.RGB(bgcolor="black", xaxis=None, yaxis=None, width=900, height=600, backend='bokeh'),
    #         hv.opts.Points(bgcolor="black", xaxis=None, yaxis=None, backend='matplotlib'),
    #         hv.opts.Image(bgcolor="black", xaxis=None, yaxis=None, backend='matplotlib'),
    #     ]

    @staticmethod
    def matplotlib_opts():
        def hook(plot, element):
            for i, _ in enumerate(plot.handles['legend'].legendHandles):
                plot.handles['legend'].legendHandles[i]._sizes = [30]
        return [
            # show_grid may not be functional, see this issue:
            # https://github.com/holoviz/holoviews/issues/3729#issue-447808528
            hv.opts.Points(s=6, padding=0.05, aspect=1, bgcolor="lightgrey", cmap="glasbey_cool"),
            hv.opts.RGB(bgcolor="black", padding=0.05, xaxis=None, yaxis=None, aspect=1),
            hv.opts.NdOverlay(xaxis=None, yaxis=None, legend_position='right', fig_inches=8, hooks=[hook]),
        ]

    def __init__(self, *args, **params):
        logger.info('Creating UMAP Interface...')
        super().__init__(*args, **params)

        if self.hue is not None:
            self.param.set_param(annotation_variables=[self.hue])

        if self.count_variable is None:
            self.param.set_param(**{"count_variable": self.gem.count_array_name})
        self.param["count_variable"].objects = sorted(self.gem.count_array_names)

        avail_mappings = [None]
        if self.gene_set_collection is not None:
            avail_mappings = list(sorted(self.gene_set_collection.gene_sets.keys()))

        self.param["selected_gene_sets"].objects = avail_mappings

        # # Infer the variable categories of the supplied annotations.
        # if self.annotation_categories is None:
        #     logger.info('No annotations selected, this causes all available annotations to be available by default.')
        #     self.param.set_param(annotation_categories=self.gem.infer_variables())
        #     self.param.set_param(annotation_variables=self.annotation_categories["all_labels"])

        # Limit the number of neighbors to the number of samples.
        self.param["n_neighbors"].bounds = [1, len(self.gem.data[self.gem.sample_index_name]) - 1]

        # logger.info('Setting hue options...')
        # avail_hues = []
        # if self.annotation_categories.get("discrete") is not None:
        #     avail_hues += self.annotation_categories.get("discrete")
        # if self.annotation_categories.get("quantile") is not None:
        #     avail_hues += self.annotation_categories.get("quantile")

        # avail_hues = list(set(avail_hues))
        # if avail_hues:
        #     self.param["hue"].objects = [None] + sorted(avail_hues)

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
        logger.info('Running UMAP transform...')
        array = self.x_count_data
        kwargs = self.get_transform_kwargs()
        return self.static_transform(array, **kwargs)

    @lru_cache()
    def cached_transform(self, _count_array_name, _gene_set, _transform_state):
        """A cached transform based on the genes and transform arguments selected."""
        logger.info('Cached UMAP transform...')
        return self.transform()

    # TODO: Prevent extra coordinates?
    def build_embedding_data_frame(self):
        gene_set = frozenset(self.get_gene_index())
        transform_state = frozenset(self.get_transform_kwargs().items())
        count_array_state = frozenset(self.count_variable)

        logger.info('Construction annotations...')
        if self.hue:
            labels = list(set([self.hue] + self.annotation_variables)) if self.annotation_variables else [self.hue]
            self.param.set_param(annotation_variables=labels)
            df = self.y_annotation_data.to_dataframe()
            # colors = hv.plotting.util.process_cmap(cmap='glasbey', ncolors=len(df[self.hue].unique()))
        else:
            df = pd.DataFrame()

        transform = self.cached_transform(count_array_state, gene_set, transform_state)
        df["x"] = transform[:, 0]
        df["y"] = transform[:, 1]
        df["Sample"] = self.get_sample_index()

        return df

    def build_plot(self, df, opts):
        vdims = [col for col in df.columns if not any(col == kdim for kdim in ['x', 'y'])]

        points = hv.Points(df, kdims=["x", "y"], vdims=vdims)#.opts(opts)

        if hv.Store.current_backend == 'bokeh':
            hover = HoverTool(tooltips=[(name, "@" + f"{name}") for name in vdims])
            points = points.opts(tools=[hover])

        if self.hue:
            points = points.groupby([self.hue], container_type=hv.NdOverlay)
            # points = points.opts(color=hv.Cycle())

        return points.opts(opts)

    def datashaded_view(self, df, opts):
        points = hv.Points(df, kdims=["x", "y"]).opts(opts)

        kwargs = {}

        if self.datashade_kwargs:
            kwargs = {**self.datashade_kwargs}

        if self.hue:
            points = points.groupby([self.hue], container_type=hv.NdOverlay)
            colors = hv.plotting.util.process_cmap(cmap='glasbey', ncolors=len(df[self.hue].unique()))
            kwargs = {**kwargs, 'aggregator': count_cat(self.hue), 'color_key': colors}

        return datashader.datashade(points, **kwargs).opts(opts)

    @param.depends('update')
    def view(self, **params):
        """A `holoviews.Points` plot of the selected transform."""
        if params:
            self.param.set_param(**params)

        if self.hue is not None:
            self.param.set_param(annotation_variables=[self.hue])

        if self.selected_gene_sets is [None]:
            return pn.pane.HTML('Select Genes and press "Update".')

        if self.dynspread is True and self.datashade is False:
            logger.warning('Setting dynspread=True has no effect without also setting datashader=True.')

        df = self.build_embedding_data_frame()

        if self.plot_options is not None:
            opts = self.plot_options
        else:
            backend_options = {"bokeh": self.bokeh_opts, "matplotlib": self.matplotlib_opts}
            backend = hv.Store.current_backend
            # if backend not in backend_options.keys():
            #     raise ValueError(f"{backend} is not a valid backend selection. Select from 'bokeh' or 'matplotlib'.")
            opts = backend_options[backend]()

        if self.datashade is True:
            plot = self.datashaded_view(df, opts)
            if self.dynspread is True:
                plot = datashader.dynspread(plot)
        else:
            plot = self.build_plot(df, opts)

        return plot

    def panel(self):
        """Interactive panel application for transform exploration."""
        # Set the selected sets to None so that the user does not have to wait for the UMAP transform to complete.
        if self.selected_gene_sets is None:
            self.param.set_param(selected_gene_sets=[None])
        # Update styles or types of buttons for the interface.
        transform_controls = pn.Param(self.param, widgets={
            'update': {'type': pn.widgets.Button, 'button_type': 'primary'},
        })

        # tab_layout = pn.Tabs(
        #     ("UMAP", pn.Row(self.view, transform_controls)),
        #     ("Documentation", generate_help_pane(self))
        # )
        return pn.Row(transform_controls, self.view)
        # return tab_layout
