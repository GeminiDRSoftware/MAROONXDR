"""
A basic init file so python knows to treat maroonx_fit as a sub-package
"""

from gempy.utils import logutils

# maroonx_fit level logger variable
_logger = None

def set_logger(logger):
    """Set the logger for this module"""
    global _logger
    _logger = logger

def get_logger():
    """Get the configured logger"""
    global _logger
    if _logger is None:
        logutils.config(file_name="maroonx_fit.log", mode="debug")
        _logger = logutils.get_logger(__name__)
    return _logger