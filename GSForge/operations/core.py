from typing import Union

import param
from ..models import Interface


class get_gem_data(Interface, param.ParameterizedFunction):
    """
    Gets the GEM matrix and an optional annotation column.
    """

    output_mode = param.ObjectSelector(default="tuple", objects=["tuple", "single"])
    output_type = param.ObjectSelector(default="xarray", objects=["xarray", "pandas", "numpy"])

    def __call__(self, *args, **params):
        super().__init__(*args, **params)
        return self.process()

    def process(self):
        mode_funcs = {"tuple": lambda self_: (self_.x_count_data, self_.y_annotation_data),
                      "single": lambda self_: self_.selection}

        type_funcs = {"xarray": lambda data: data,
                      "pandas": lambda data: data.to_dataframe(),
                      "numpy": lambda data: data.values}

        data = mode_funcs[self.output_mode]()
        data = map(type_funcs[self.output_type], data) if isinstance(data, tuple) \
            else type_funcs[self.output_type]
        return data
