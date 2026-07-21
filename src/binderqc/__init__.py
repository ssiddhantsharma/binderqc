"""Score binder termini for tag/conjugation suitability from predicted complexes."""

from .core import score_structure
from .paths import gather_paths

__version__ = "0.1.0"
__all__ = ["score_structure", "gather_paths", "__version__"]
