import logging

import param

from ..models import Interface

__all__ = [
    "get_gem_data"
]

logger = logging.getLogger(__name__)


class get_gem_data(param.ParameterizedFunction, Interface):
    """
    Gets the GEM matrix and an optional annotation column.
    """

    single_object = param.Boolean(default=False)
    output_type = param.ObjectSelector(default="xarray", objects=["xarray", "pandas", "numpy"])

    def __new__(cls, *args, **params):
        inst = cls.instance(**params)
        inst.__init__(*args, **params)
        return inst.__call__()

    def __call__(self):
        return self.get_gem_data(self.single_object, self.output_type)
