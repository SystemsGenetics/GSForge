import param

from ..models import Interface

__all__ = [
    "get_gem_data"
]


class get_gem_data(Interface, param.ParameterizedFunction):
    """
    Gets the GEM matrix and an optional annotation column.
    """

    tuple_output = param.Boolean(default=True)
    output_type = param.ObjectSelector(default="xarray", objects=["xarray", "pandas", "numpy"])

    def __call__(self, *args, **params):
        super().__init__(*args, **params)
        return self.get_gem_data(self.tuple_output, self.output_type)
