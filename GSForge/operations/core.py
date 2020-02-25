from typing import Union

import param
from ..models import Interface

__all__ = [
    "get_gem_data"
]


class get_gem_data(Interface, param.ParameterizedFunction):
    """
    Gets the GEM matrix and an optional annotation column.
    """

    output_mode = param.ObjectSelector(default="tuple", objects=["tuple", "single"])
    output_type = param.ObjectSelector(default="xarray", objects=["xarray", "pandas", "numpy"])

    def __call__(self, *args, **params):

        super().__init__(*args, **params)

        mode_funcs = {"tuple": lambda self_: (self_.x_count_data, self_.y_annotation_data),
                      "single": lambda self_: self_.selection}

        def __pandas_output(data_):
            if data_ is not None:
                return data_.to_dataframe()
            else:
                return None

        def __numpy_output(data_):
            if data_ is not None:
                return data_.values
            else:
                return None

        type_funcs = {"xarray": lambda data_: data_,
                      "pandas": __pandas_output,
                      "numpy": __numpy_output}

        data = mode_funcs[self.output_mode](self)
        data = tuple(map(type_funcs[self.output_type], data)) if isinstance(data, tuple) \
            else type_funcs[self.output_type](data)
        return data
