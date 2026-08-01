"""Score binder termini for tag/conjugation suitability from predicted complexes."""

from .core import score_structure, grippability_consensus
from .paths import gather_paths

__version__ = "0.2.0"
__all__ = ["score_structure", "grippability_consensus", "gather_paths", "__version__"]
