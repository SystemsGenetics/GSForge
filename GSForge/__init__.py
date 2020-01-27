from .models import *
from .operations import get_data

import os

if "GSFORGE_MINIMAL" not in os.environ:
    from . import plots
    from . import panels
