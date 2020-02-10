import xarray as xr
import param
import inspect
import holoviews as hv
import functools
import pandas as pd

from ..models import Interface, GeneSet

DGE_DEFAULT_KWARGS = dict(
    log_fold_change_var=[
        "logFC",  # from EdgeR.
        "log2FoldChange",  # from DESeq2.
    ],
    mean_value_var=["baseMean"],
    p_value_var=["pvalue"]
)


def infer_kwarg_defaults_from_data(source: xr.Dataset, function) -> dict:
    kwargs = dict()
    key_overlap = set.intersection(set(inspect.signature(function).parameters.keys()),
                                   set(DGE_DEFAULT_KWARGS.keys()))

    for key in key_overlap:
        default_argument = set.intersection(set(source.variables.keys()), DGE_DEFAULT_KWARGS[key])
        if len(default_argument) > 1:
            raise ValueError("More than one potential default found. Explicitly set the arguments"
                             "to this function.")
        default_argument = list(default_argument)[0]
        if default_argument is not None:
            kwargs[key] = default_argument

    return kwargs


def get_param_process_overlap_kwargs(paramed, process) -> dict:
    """Gets overlapping kwargs of the given process and parameters of a Parameterized class,
    and returns them as a dictionary."""
    key_set = set.intersection(set(inspect.signature(process).parameters.keys()),
                               set(paramed.param.objects().keys()))
    return {key: getattr(paramed, key) for key in key_set if getattr(paramed, key) is not None}


class AbstractPlottingOperation(param.ParameterizedFunction):
    backend = param.ObjectSelector(default=None, objects=["bokeh", "matplotlib"], doc="""
        The selected plotting backend to use for display. Options are ["bokeh", "matplotlib"].""")

    apply_default_opts = param.Boolean(default=True, precedence=-1.0, doc="""
        Whether to apply the default styling based on the current backend.""")

    @staticmethod
    def bokeh_opts():
        raise NotImplementedError("'bokeh' options are not supported for this plotting function.")

    @staticmethod
    def matplotlib_opts():
        raise NotImplementedError("'matplotlib' options are not supported for this plotting function.")

    def process(self):
        """The core plotting function that must be defined for subclasses."""
        raise NotImplementedError("Sub-classes of PlottingOperation must define a `process` function.")

    def get_default_options(self):
        """Apply default styling options by default."""
        backend_options = {"bokeh": self.bokeh_opts, "matplotlib": self.matplotlib_opts}
        backend = hv.Store.current_backend if self.backend is None else "bokeh"
        if backend not in backend_options.keys():
            raise ValueError(f"{backend} is not a valid backend selection. Select from 'bokeh' or 'matplotlib'.")
        return backend_options[backend]()

    def __call__(self, *args, **params):
        super().__init__(**params)
        layout = self.process()
        if self.apply_default_opts is False:
            return layout
        return layout.opts(self.get_default_options())


class ResultPlottingOperation(AbstractPlottingOperation):
    source = param.Parameter()

    @functools.singledispatchmethod
    def __dispatch_input_args(self, source, **params):
        """Dispatches *args based on the first time to be parsed and joined with an updated params dict."""
        raise NotImplementedError(f"Cannot parse {source} of type {type(source)}.")

    @__dispatch_input_args.register
    def __parse_xarray_datasets(self, source: xr.Dataset, **params) -> dict:
        return {"source": source, **params}

    @__dispatch_input_args.register
    def __parse_xarray_datasets(self, source: GeneSet, **params) -> dict:
        return {"source": source.data, **params}

    @__dispatch_input_args.register
    def __parse_xarray_datasets(self, source: pd.DataFrame, **params) -> dict:
        print("Converting `pandas.DataFrame` to an `xarray` object.")
        return {"source": source.to_xarray(), **params}

    def __call__(self, *args, **params):
        if args:
            params = {**self.__dispatch_input_args(*args), **params}
        super().__init__(**params)
        layout = self.process()
        if self.apply_default_opts is False:
            return layout
        return layout.opts(self.get_default_options())
