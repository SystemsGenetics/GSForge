"""
There are two core data models in **GSForge**, both of which store their associated data
in ``xarray.Dataset`` object under a ``data`` attribute. You are encouraged to consult the
`xarray documentation <http://xarray.pydata.org/en/stable/>`_
for how to perform any transform or selection not provided by *GSForge*.

Core Data Classes
-----------------

**AnnotatedGEM**
    Contains the gene expression matrix, which is indexed by a 'Gene'
    and 'Sample' coordinates. This ``xarray.Dataset`` object also contains
    (but is not limited to) phenotype information as well.

**GeneSet**
    A GeneSet is a set of genes and any associated values. A GeneSet
    can a set of 'supported' genes, *i.e.* genes that are 'within' a
    given GeneSet.

These core data classes are constructed with a limited set of packages:

- ``numpy``
- ``pandas``
- ``xarray``
- ``param``

This allows the creation of container images without interactive visualization libraries.
"""
from ._AnnotatedGEM import AnnotatedGEM
from ._GeneSet import GeneSet
from ._GeneSetCollection import GeneSetCollection
from ._Interface import Interface, CallableInterface


__all__ = [
    "AnnotatedGEM",
    "GeneSet",
    "GeneSetCollection",
    "Interface",
    "CallableInterface",
]
