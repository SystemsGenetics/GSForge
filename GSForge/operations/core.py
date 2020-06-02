import param

from ..models import Interface

__all__ = [
    "get_gem_data"
]


class get_gem_data(param.ParameterizedFunction, Interface):
    """
    Gets the GEM matrix and an optional annotation column.
    """

    tuple_output = param.Boolean(default=True)
    output_type = param.ObjectSelector(default="xarray", objects=["xarray", "pandas", "numpy"])

    def __new__(cls, *args, **params):
        inst = cls.instance(**params)
        inst.__init__(*args, **params)
        # return inst
        return inst.__call__()

    def __call__(self):
        # super().__init__(*args, **params)
        # return self
        return self.get_gem_data(self.tuple_output, self.output_type)
