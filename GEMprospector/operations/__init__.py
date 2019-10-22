"""
GEMprospector operations can be broken down into three categories:

`Analytics`
    For discrete operations, *i.e.* chi-squared tests, differential gene expression, etc.
`Normalizations`
    For those operations that are meant to create an entire transform of the GEM.
`Prospectors`
    For non-deterministic operations, used in ranking and comparing gene selections.

"""

from .analytics import *
# from .normalize import *
from .prospectors import *
