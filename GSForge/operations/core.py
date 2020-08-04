import logging

import param

from ..models import CallableInterface
# from ..models import Interface

__all__ = [
    "get_gem_data"
]

logger = logging.getLogger("GSForge")


class get_gem_data(CallableInterface):
    """
    Gets the GEM matrix and an optional annotation column from either an AnnotatedGEM or GeneSetCollection.

    For both cases the user may:
        + Have gene N/A values masked, dropped gene-wise, or zero-filled.
        + Have a transform (e.g. log) run on the resultant count array.
        + Return the count [and annotation] data be returned as a ``numpy.ndarray``, ``pandas.DataFrame``, or
          a ``xarray.DataArray`` [and ``xarray.DataSet``].

    ```
    agem = gsf.AnnotatedGEM('my_gem.nc')
    count_xr, empty_labels = gsf.get_gem_data(agem, count_mask='complete')
    count_xr, labels = gsf.get_gem_data(agem, count_mask='complete', annotation_variables=['class', 'treatment'])
    ```

    For the ``GeneSetCollection`` case users may select genes by their membership in one or more ``GeneSets``, and
    some set operation:
        + Union of the genes within the selected sets.
        + Intersection of the genes within the selected sets.
        + Unique to selected sets.


    """

    single_object = param.Boolean(default=False)
    output_type = param.ObjectSelector(default="xarray", objects=["xarray", "pandas", "numpy"])

    # def __new__(cls, *args, **params):
    #     """
    #     When this class is used created by a user calling ``gsf.get_gem_data``, we create an new instance
    #     of the Interface class, and return the result of ``__call__`` immediately.
    #     """
    #     logger.debug('Creating a new GSForge.Interface instance.')
    #     if args:
    #         params = CallableInterface._interface_dispatch(*args, **params)
    #     inst = cls.instance(**params)  # See the param code for more on this `instance` function.
    #     return inst.__call__()

    def __call__(self):
        return self.get_gem_data(self.single_object, self.output_type)
