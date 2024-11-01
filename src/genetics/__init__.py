"""
Computational biology tools for genetic analysis.
"""

__version__ = "0.1.0"

from .core import calc_pi, calc_watterson, tajimas_d
from .frontend import main

__all__ = [
    "calc_pi",
    "calc_watterson",
    "tajimas_d",
    "main"
]