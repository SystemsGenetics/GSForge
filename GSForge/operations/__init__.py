"""
GSForge operations have been separated into categories:

Core
  Core operations provide a data-access interface to AnnotatedGEM and GeneSetCollection objects.

Analytics
  For discrete operations, *i.e.* chi-squared tests, differential gene expression, etc. As well
  as tools for analytics based on machine-learning models.

Normalizations
  For those operations that are meant to create an entire transform of the GEM.

Prospectors
  For non-deterministic operations, used in ranking and comparing gene selections.

"""

from .analytics import *
from .core import *
from .normalizations import *
from .prospectors import *
