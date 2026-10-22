"""Score binder termini for tag/conjugation suitability from predicted complexes."""

from .core import grippability_consensus, score_structure, self_association
from .paths import gather_paths

__version__ = "0.6.0"
__all__ = ["score_structure", "grippability_consensus", "self_association",
           "gather_paths", "__version__"]
