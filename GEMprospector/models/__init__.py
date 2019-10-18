"""
There are two 'core' data models in *GEMprospector*, both of which
store their associated data in `xarray.Dataset` object under a `data`
attribute. You are encouraged to consult the xarray documentation:
http://xarray.pydata.org/en/stable/ for how to perform any transform
or selection not provided by *GEMprospector*. The two 'core' data
classes are:

**AnnotatedGEM**
    Contains the gene expression matrix, which is indexed by a 'Gene'
    and 'Sample' coordinates. This `xarray.Dataset` object also contains
    (but is not limited to) phenotype information as well.

**Lineament**
    A lineament is a set of genes and any associated values. A lineament
    can a set of 'supported' genes, *i.e.* genes that are 'within' a
    given lineament.

The `interface` classes provide patterns of data access and common
transformations that researchers may need from the `core` classes. They are:

**LineamentCollection**
    The work-horse of the *GEMprospector* package. This object contains an
    *AnnotatedGEM* and a python dictionary of {name: Lineament} objects.
    This class contains functions for comparing and analyzing *Lineaments*,
    as well as tools to pass of *Lineament*-derived subsets to other functions.

**Interface**
    The *Interface* object provides a common API to interacting with *AnnotatedGEM*
    or *LineamentCollection*. It provides functions that facilitate pulling gene
    or sample subsets and access to any transforms of the count matrix.

**OperationInterface**
    Aside from being abstract, this is the same as the above *Interface*, except
    this calls a single function as defined by `process` in a subclass.
"""

from ._AnnotatedGEM import AnnotatedGEM
from ._Lineament import Lineament
from ._LineamentCollection import LineamentCollection
from ._Interface import Interface
from ._OperationInterface import OperationInterface


__all__ = [
    "AnnotatedGEM",
    "Lineament",
    "LineamentCollection",
    "Interface",
    "OperationInterface"
]
