"""
Utility functions for GSForge.
"""
import logging
import sys

from . import R_interface
from . import converters
from ._models import *

# TODO: implement __all__?

logger = logging.getLogger('GSForge')


# TODO: Check for an existing handler that outputs to stdout.
def transient_log_handler(func):

    def call(*args, **kwargs):

        # if not len(logger.handlers):

        if kwargs.get('verbose'):
            log_level = getattr(logging, kwargs.pop('verbose'))

            handler = logging.StreamHandler(sys.stdout)
            handler.setFormatter(logging.Formatter('%(name)s - %(levelname)s - %(message)s'))

            logger.setLevel(log_level)
            logger.addHandler(handler)
            logger.info(f'Log level set to {log_level} for call to function name.')
            result = func(*args, **kwargs)
            logger.removeHandler(handler)

        else:
            result = func(*args, **kwargs)
        return result

    return call
