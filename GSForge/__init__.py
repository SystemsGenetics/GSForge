"""
****************
Reference Manual
****************

Welcome to the API reference for GSForge.

`GSForge.models`_
 Primary data models for GSForge.
`GSForge.operations`_
 Transforms and feature selection tools.
`GSForge.panels`_
 Interactive visualizations.
`GSForge.plots`_
 Plotting tools.
`GSForge.utils`_
 Utilities.

"""


from .models import *
from .operations.core import *

import os

if "GSFORGE_MINIMAL" not in os.environ:
    from . import plots
    from . import panels
