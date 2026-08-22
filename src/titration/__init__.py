"""Anaerobic Digester 5-Point Titration Analysis

Converts legacy Pascal algorithms (1992) to modern Python for
water chemistry analysis of anaerobic digesters.
"""

__version__ = "0.1.0"
__author__ = "Charles Mellish"

from .core import calculate_alkalinity, find_endpoint, parse_pascal_equivalent

__all__ = [
    "calculate_alkalinity",
    "find_endpoint", 
    "parse_pascal_equivalent",
    "__version__",
]
