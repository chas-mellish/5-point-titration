"""5-Point Titration Analysis Package"""

__version__ = "0.1.0"
__author__ = "Charles Mellish"

from .core import run_titration
from .constants import THESIS_DEFAULTS
from .models import TitrationInput, TitrationResult

__all__ = [
    "run_titration",
    "TitrationInput",
    "TitrationResult",
    "THESIS_DEFAULTS",
    "__version__",
]
