"""
Some simple wrappers for maintaining ``xarray`` coordinates through ``sklearn`` functions.

This module may eventually be supplanted by use of `sklearn xarray <https://phausamann.github.io/sklearn-xarray/index.html>`_.
"""

from sklearn.model_selection import cross_val_score, train_test_split
from ..models import OperationInterface

__all__ = [
    "train_test_split_wrapper",
]


class train_test_split_wrapper(OperationInterface):
    """
    Performs an ``sklearn.preprocessing.train_test_split()`` call on the subset of data
    specified by the interface options (the same options passed to ``get_data()``.

    :returns: x_train, x_test, y_train, y_test
    """
    # TODO: Add links and reference to the sklearn function and docs.
    train_test_split_options = param.Parameter(default=dict())

    def process(self):
        # Get the subset of data selected by this operation.
        y_index = self.get_sample_index()

        # Get the sample index and make the train and test indexes.
        train_idx, test_idx = train_test_split(y_index, **self.train_test_split_options)

        x_train = self.x_count_data.sel({self.gem.sample_index_name: train_idx})
        x_test = self.x_count_data.sel({self.gem.sample_index_name: test_idx})

        y_train = self.y_annotation_data.sel({self.gem.sample_index_name: train_idx})
        y_test = self.y_annotation_data.sel({self.gem.sample_index_name: test_idx})

        return x_train, x_test, y_train, y_test
