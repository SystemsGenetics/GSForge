"""
Welcome to the GSForge API reference documentation.

:doc:`GSForge.models`
 Primary data models for GSForge.

:doc:`GSForge.operations`
 Transforms and feature selection tools.

:doc:`GSForge.plots`
 Plotting tools.

:doc:`GSForge.panels`
 Interactive visualizations.

:doc:`GSForge.utils`
 Utilities.
"""
from .models import *
from .operations.core import *

import os
# import logging
# import sys


# logging.basicConfig(stream=sys.stdout, level=logging.WARNING)


if "GSFORGE_INSTALL_MODE" not in os.environ:
    from . import plots
    from . import panels
