import inspect
import logging

import holoviews as hv
import panel as pn
import param
import umap
from bokeh.models import HoverTool
from methodtools import lru_cache
import pandas as pd
from holoviews.operation import datashader
# from holoviews.operation.datashader import datashade, dynspread
from datashader.reductions import count_cat

from .utils import generate_help_pane
from ..models._Interface import Interface

logger = logging.getLogger("GSForge")


# TODO: Allow size selection of the points drawn.
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

    color_map = param.Parameter()

    # marker_variable = param.Parameter(default=None)

    datashade = param.Boolean(default=None)
    dynspread = param.Boolean(default=None)

    # Create a button so that the transform only occurs when clicked.
    update = param.Action(lambda self: self.param.trigger('update'))

    @staticmethod
    def bokeh_opts():
        return [
            hv.opts.Points(
                cmap="Set1",
                legend_position='right',
                axiswise=True,
                xaxis=None,
                yaxis=None,
                padding=0.05,
                show_grid=True,
                bgcolor="lightgrey",
                width=900,
                height=600,
                backend="bokeh"),
            hv.opts.RGB(
                show_grid=True,
                bgcolor="lightgrey",
                xaxis=None,
                yaxis=None,
                width=900,
                height=600,
            )
        ]

    @staticmethod
    def matplotlib_opts():
        raise NotImplementedError("'matplotlib' options are not supported for this plotting function.")

    def __init__(self, *args, **params):
        logger.info('Creating UMAP Interface...')
        super().__init__(*args, **params)  # This may be fragile.

        # Infer the variable categories of the supplied dataframe.
        if self.annotation_categories is None:
            logger.info('No annotations selected, this causes all available annotations to be available by default.')
            self.param.set_param(annotation_categories=self.gem.infer_variables())
            self.param.set_param(annotation_variables=self.annotation_categories["all_labels"])

        # Limit the number of neighbors to the number of samples.
        self.param["n_neighbors"].bounds = [1, len(self.gem.data[self.gem.sample_index_name]) - 1]

        logger.info('Setting hue options...')
        avail_hues = []
        if self.annotation_categories.get("discrete") is not None:
            avail_hues += self.annotation_categories.get("discrete")
        if self.annotation_categories.get("quantile") is not None:
            avail_hues += self.annotation_categories.get("quantile")

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
        logger.info('Running UMAP transform...')
        array = self.x_count_data
        kwargs = self.get_transform_kwargs()
        return self.static_transform(array, **kwargs)

    @lru_cache()
    def cached_transform(self, _count_array_name, _gene_set, _transform_state):
        """A cached transform based on the genes and transform arguments selected."""
        logger.info('Cached UMAP transform...')
        return self.transform()

    def build_embedding_data_frame(self):
        gene_set = frozenset(self.get_gene_index())
        transform_state = frozenset(self.get_transform_kwargs().items())
        count_array_state = frozenset(self.count_variable)

        if len(gene_set) == 0:
            return pn.pane.Markdown('No genes in selected set.')

        logger.info('Construction annotations...')
        if self.hue:
            labels = list(set([self.hue] + self.annotation_variables)) if self.annotation_variables else [self.hue]
            self.param.set_param(annotation_variables=labels)
            df = self.y_annotation_data.to_dataframe()
        else:
            df = pd.DataFrame()

        transform = self.cached_transform(count_array_state, gene_set, transform_state)
        df["x"] = transform[:, 0]
        df["y"] = transform[:, 1]

        return df

    def bokeh_view(self):
        df = self.build_embedding_data_frame()
        vdims = [col for col in df.columns if not any(col == kdim for kdim in ['x', 'y'])]
        hover = HoverTool(tooltips=[(name, "@" + f"{name}") for name in vdims])
        points = hv.Points(df, kdims=["x", "y"], vdims=vdims).opts(tools=[hover])

        # if self.marker_variable:
        #     markers = hv.Cycle(['o', 's', '^', 'v', '*', 'D', 'h', 'x', '+', '8', 'p', '<', '>', 'd', 'H'])
        # marker_mapping = {k: v for k, v in zip(df[self.marker_variable].unique(), markers)}
        # points = points.opts(marker=hv.dim(self.marker_variable))

        if self.hue:
            # points = points.groupby([self.hue], container_type=hv.NdOverlay)
            points = points.opts(color=self.hue)

        return points.opts(self.bokeh_opts())

    def datashaded_view(self):
        df = self.build_embedding_data_frame()
        vdims = [col for col in df.columns if not any(col == kdim for kdim in ['x', 'y'])]
        points = hv.Points(df, kdims=["x", "y"], vdims=vdims).opts(self.bokeh_opts())
        if self.hue:
            points = points.groupby([self.hue], container_type=hv.NdOverlay)
            plot = datashader.datashade(points, aggregator=count_cat(self.hue))
        else:
            plot = datashader.datashade(points, cmap='darkblue')
        return plot.opts(self.bokeh_opts())

    @param.depends('update', 'hue')
    def view(self):
        """A `holoviews.Points` plot of the selected transform."""
        # Call to bokeh or datashade?
        if self.dynspread is True and self.datashade is False:
            logger.warning('Setting dynspread=True has no effect without also setting datashader=True.')

        if self.datashade is True:
            plot = self.datashaded_view()
            if self.dynspread is True:
                plot = datashader.dynspread(plot)
        else:
            plot = self.bokeh_view()

        return plot

    def panel(self):
        """Interactive panel application for transform exploration."""
        # Update styles or types of buttons for the interface.
        transform_controls = pn.Param(self.param, widgets={
            'update': {'type': pn.widgets.Button, 'button_type': 'primary'},
        })

        tab_layout = pn.Tabs(
            ("UMAP", pn.Row(transform_controls, self.view)),
            ("Documentation", generate_help_pane(self))
        )
        return tab_layout
