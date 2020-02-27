"""
Welcome to the GSForge API reference documentation.

:ref:`GSForge.models`
 Primary data models for GSForge.

:ref:`GSForge.operations`
 Transforms and feature selection tools.

:ref:`GSForge.plots`
 Plotting tools.

:ref:`GSForge.panels`
 Interactive visualizations.

:ref:`GSForge.utils`
 Utilities.
"""

from .models import *
from .operations.core import *

import os

if "GSFORGE_MINIMAL" not in os.environ:
    from . import plots
    from . import panels
