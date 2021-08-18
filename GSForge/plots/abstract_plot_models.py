import inspect

import holoviews as hv
import pandas as pd
import param
import xarray as xr

from .._singledispatchmethod import singledispatchmethod
from ..models._GeneSet import GeneSet
from ..models._Interface import CallableInterface


class AbstractPlottingOperation(param.ParameterizedFunction):
    # backend = param.ObjectSelector(default=None, objects=["bokeh", "matplotlib"], doc="""
    #     The selected plotting backend to use for display. Options are ["bokeh", "matplotlib"].""")

    apply_default_opts = param.Boolean(default=True, precedence=-1.0, doc="""
        Whether to apply the default styling based on the current backend.""")

    plot_options = param.Parameter(default=None, doc="""
    User supplied options to the plotting functions. If provided (and is not None), these will take
    precedence over a functions built-in defaults.""")

    @staticmethod
    def bokeh_opts():
        raise NotImplementedError("'bokeh' options are not supported for this plotting function.")

    @staticmethod
    def matplotlib_opts():
        raise NotImplementedError("'matplotlib' options are not supported for this plotting function.")

    def get_default_options(self):
        """Apply default styling options by default."""
        if self.plot_options is not None:
            return self.plot_options
        else:
            backend_options = {"bokeh": self.bokeh_opts, "matplotlib": self.matplotlib_opts}
            backend = hv.Store.current_backend if hv.Store.current_backend else "bokeh"
            if backend not in backend_options.keys():
                raise ValueError(f"{backend} is not a valid backend selection. Select from 'bokeh' or 'matplotlib'.")
            return backend_options[backend]()

    def get_param_process_overlap_kwargs(self, process) -> dict:
        """Gets overlapping kwargs of the given process and parameters of a Parameterized class,
        and returns them as a dictionary."""
        key_set = set.intersection(set(inspect.signature(process).parameters.keys()),
                                   set(self.param.objects().keys()))
        return {key: getattr(self, key) for key in key_set if getattr(self, key) is not None}

    def __call__(self, *args, **params):
        raise NotImplementedError("Sub-classes of PlottingOperation must define a `__call__` function.")


class InterfacePlottingBase(CallableInterface, AbstractPlottingOperation):
    """
    This abstract base class should be used for plotting operations that act upon a GSForge interface object,
    and which directly return a plot.
    """
    pass


class ResultPlottingOperation(AbstractPlottingOperation):
    source = param.Parameter()
    DGE_DEFAULT_KWARGS = dict(
        log_fold_change_var=[
            "logFC",  # from EdgeR.
            "log2FoldChange",  # from DESeq2.
        ],
        mean_value_var=[
            "baseMean",
            "logCPM",
        ],
        p_value_var=[
            "pvalue",
            "PValue",
        ]
    )

    def __new__(cls, *args, **params):
        inst = cls.instance()
        if args:
            params = inst.__dispatch_input_args(*args, **params)
        inst.__init__(**params)
        return inst.__call__(*args, **params)

    def __call__(self, *args, **params):
        raise NotImplementedError("Sub-classes of PlottingOperation must define a `__call__` function.")

    @classmethod
    def infer_kwarg_defaults_from_data(cls, source: xr.Dataset, function) -> dict:
        """Try to infer variable names for plotting functions.
        e.g. Try to find what the 'p-value' column.
        """
        kwargs = dict()
        key_overlap = set.intersection(set(inspect.signature(function).parameters.keys()),
                                       set(cls.DGE_DEFAULT_KWARGS.keys()))
        for key in key_overlap:
            default_argument = set.intersection(set(source.variables.keys()), cls.DGE_DEFAULT_KWARGS[key])

            if len(default_argument) == 0:
                continue

            if len(default_argument) > 1:
                raise ValueError("More than one potential default found. Explicitly set the arguments"
                                 "to this function.")

            default_argument = list(default_argument)[0]
            if default_argument is not None:
                kwargs[key] = default_argument

        return kwargs

    @singledispatchmethod
    def __dispatch_input_args(self, source, *args, **params):
        """Dispatches *args based on the first time to be parsed and joined with an updated params dict."""
        raise NotImplementedError(f"Cannot parse {source} of type {type(source)}.")

    @__dispatch_input_args.register
    def __parse_xarray_datasets(self, source: xr.Dataset, *args, **params) -> dict:
        return {"source": source, **params}

    @__dispatch_input_args.register
    def __parse_gene_set(self, source: GeneSet, **params) -> dict:
        return {"source": source.data, **params}

    @__dispatch_input_args.register
    def __parse_pandas_dataframe(self, source: pd.DataFrame, **params) -> dict:
        print("Converting `pandas.DataFrame` to an `xarray` object.")
        return {"source": source.to_xarray(), **params}

# class CollectionOperation(AbstractPlottingOperation):
#     """
#     Handles plotting options and input for functions that act upon GeneSetCollection objects or some equivalent.
#
#     Collection functions are assumed to require only GeneSet membership information, e.g. a dictionary that
#     contains arrays of gene names is assumed to be sufficient.
#
#     Supported methods of providing `GeneSetCollection` or equivalent data:
#     + *GeneSetCollection, **params
#     + *Tuple[GeneSetCollection, List[str]], **params
#     + *Dict[str: np.ndarray]
#
#     """
#
#     mappings = param.Dict()
#
#     @functools.singledispatchmethod
#     def __dispatch_input_args(self, source, *args, **params):
#         """Dispatches *args based on the first time to be parsed and joined with an updated params dict."""
#         raise NotImplementedError(f"Cannot parse {source} of type {type(source)}.")
#
#     @__dispatch_input_args.register
#     def __parse_tuples(self, source: tuple, **params) -> dict:
#         """Tuple input is assumed to be one or more pairs of GeneSetCollection objects
#         and keys corresponding to GeneSets.."""
#         return {"mappings": {gsc.name: gsc.as_dict(keys)
#                              for (gsc, keys) in source},
#                 **params}
#
#     @__dispatch_input_args.register
#     def __parse_gene_set_collection(self, source: GeneSetCollection, *args, **params) -> dict:
#         collections = [source, *args]
#         # TODO: Ensure that all *args are GeneSetCollection.
#         mappings = dict()
#         for gsc in collections:
#             coll_mappings = gsc.as_dict()
#
#         return {"source": source, **params}
