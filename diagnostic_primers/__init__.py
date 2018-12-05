# -*- coding: utf-8 -*-
__version__ = "0.2.0.dev"

import sys
import traceback


# Define base PDPException
class PDPException(Exception):
    """General exception for PDP/find_differential_primers"""

    def __init__(self, msg="Error in pdp/diagnostic_primers"):
        Exception.__init__(self, msg)


# Report last exception as string
def last_exception():
    """Return last exception as a string, or use in logging."""
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return "".join(traceback.format_exception(exc_type, exc_value, exc_traceback))
