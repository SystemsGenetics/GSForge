"""
There are two 'core' data models in *GSForge*, both of which store their associated data
in ``xarray.Dataset`` object under a ``data`` attribute. You are encouraged to consult the
`xarray documentation <http://xarray.pydata.org/en/stable/>`_
for how to perform any transform or selection not provided by *GSForge*.
The two 'core' data classes are:

**AnnotatedGEM**
    Contains the gene expression matrix, which is indexed by a 'Gene'
    and 'Sample' coordinates. This ``xarray.Dataset`` object also contains
    (but is not limited to) phenotype information as well.

**GeneSet**
    A GeneSet is a set of genes and any associated values. A GeneSet
    can a set of 'supported' genes, *i.e.* genes that are 'within' a
    given GeneSet.

The interface classes provide patterns of data access and common
transformations that researchers may need from the core data classes. They are:

**GeneSetCollection**
    The work-horse of the *GSForge* package. This object contains an
    *AnnotatedGEM* and a python dictionary of {name: GeneSet} objects.
    This class contains functions for comparing and analyzing *GeneSet*,
    as well as tools to pass of *GeneSet*-derived subsets to other functions.

**Interface**
    The *Interface* object provides a common API to interacting with *AnnotatedGEM*
    or *GeneSetCollection*. It provides functions that facilitate pulling gene
    or sample subsets and access to any transforms of the count matrix.

**OperationInterface**
    Aside from being abstract, this is the same as the above *Interface*, except
    this calls a single function as defined by ``process`` function in a subclass.


.. image:: ../../doc/_static/GSForge_index_selection.svg
  :width: 400
  :alt: GSForge Index Selection

"""

from ._AnnotatedGEM import AnnotatedGEM
from ._GeneSet import GeneSet
from ._GeneSetCollection import GeneSetCollection
from ._Interface import Interface
from ._OperationInterface import OperationInterface


__all__ = [
    "AnnotatedGEM",
    "GeneSet",
    "GeneSetCollection",
    "Interface",
    "OperationInterface"
]
