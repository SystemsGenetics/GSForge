"""
GSForge operations can be broken down into three categories:

Analytics
  For discrete operations, *i.e.* chi-squared tests, differential gene expression, etc.

Normalizations
  For those operations that are meant to create an entire transform of the GEM.

Prospectors
  For non-deterministic operations, used in ranking and comparing gene selections.

"""
from ..models import OperationInterface

from .analytics import *
from .normalizations import *
from .prospectors import *


class get_data(OperationInterface):
    """
    Gets the GEM matrix and an optional annotation column.
    """

    # TODO: Expand comment, describe how the sample and gene indexes are built.

    def process(self):
        return self.x_count_data, self.y_annotation_data
