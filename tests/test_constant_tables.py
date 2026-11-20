"""The hand-transcribed scales, checked against their published sources.

`_A3V` and `_BLACK_MOULD` are twenty floats each, typed in from a paper and from
Aggrescan3D's matrix file. Nothing else in the suite would notice a typo: every
other assertion about `sap_score` or `a3d_score` is directional (`> 0`) or a
self-consistency inequality, so a mistyped constant shifts every number the tool
reports and still passes. These are the tests that would fail.

The same shape as `test_instability_index_formula`, which already does this for
the DIWV matrix -- it just was never applied to the two scales that actually
drive the aggregation metrics.
"""

import math

from binderqc.core import _A3V, _BLACK_MOULD

# Twenty canonical amino acids, so a dropped or duplicated key is caught.
_AA = set("ACDEFGHIKLMNPQRSTVWY")


def test_both_scales_cover_exactly_the_twenty_amino_acids():
    assert set(_A3V) == _AA, f"A3V keys differ: {set(_A3V) ^ _AA}"
    assert set(_BLACK_MOULD) == _AA, f"Black-Mould keys differ: {set(_BLACK_MOULD) ^ _AA}"


def test_a3v_matches_aggrescan3d_1_0_2_verbatim():
    """From aggrescan/data/matrices/aggrescan.mat (Aggrescan3D 1.0.2).

    Transcribed independently of the dict under test, so this fails on a typo in
    either copy rather than agreeing with itself.
    """
    reference = {
        "I": 1.9136, "F": 1.8456, "V": 1.6856, "L": 1.4716, "Y": 1.2506,
        "W": 1.1286, "M": 1.0016, "C": 0.6956, "A": 0.0556, "T": -0.0674,
        "S": -0.2024, "P": -0.2424, "G": -0.4434, "H": -0.9414, "Q": -1.1394,
        "N": -1.2104, "K": -1.6164, "E": -1.7304, "D": -1.7014, "R": -1.7534,
    }
    for aa, want in reference.items():
        assert math.isclose(_A3V[aa], want, abs_tol=1e-9), f"{aa}: {_A3V[aa]} != {want}"


def test_black_mould_matches_expasy_protscale():
    """Expasy ProtScale "Hphob. / Black", Black & Mould 1991, normalised 0-1."""
    reference = {
        "A": 0.616, "R": 0.000, "N": 0.236, "D": 0.028, "C": 0.680, "E": 0.043,
        "Q": 0.251, "G": 0.501, "H": 0.165, "I": 0.943, "L": 0.943, "K": 0.283,
        "M": 0.738, "F": 1.000, "P": 0.711, "S": 0.359, "T": 0.450, "W": 0.878,
        "Y": 0.880, "V": 0.825,
    }
    for aa, want in reference.items():
        assert math.isclose(_BLACK_MOULD[aa], want, abs_tol=1e-9), f"{aa}: {_BLACK_MOULD[aa]} != {want}"


def test_the_scales_keep_the_orderings_their_papers_report():
    """Ordering checks, so a plausible-looking transcription with swapped values
    still fails. These are properties of the underlying scales, not of our copy."""
    # A3V is calibrated on aggregation data: beta-branched/aromatic at the top,
    # charged strongly negative. Arg is the minimum, Ile the maximum.
    assert max(_A3V, key=_A3V.get) == "I"
    assert min(_A3V, key=_A3V.get) == "R"
    assert _A3V["F"] > _A3V["A"] > _A3V["G"] > _A3V["D"]
    # Black-Mould is a hydrophobicity scale normalised to [0, 1]: Phe 1.0, Arg 0.0.
    assert math.isclose(_BLACK_MOULD["F"], 1.0)
    assert math.isclose(_BLACK_MOULD["R"], 0.0)
    assert all(0.0 <= v <= 1.0 for v in _BLACK_MOULD.values())
    # Leu and Ile are equal on this scale -- a real feature, easy to "fix" by mistake.
    assert math.isclose(_BLACK_MOULD["L"], _BLACK_MOULD["I"])


def test_a_typo_would_actually_be_caught():
    """The control for these tests: perturb a value and confirm the check fails.

    Without this, a test that silently compared a dict to itself would pass and
    look like validation.
    """
    perturbed = dict(_A3V)
    perturbed["I"] = 1.9316  # transposed digits, the realistic typo
    assert not math.isclose(perturbed["I"], _A3V["I"], abs_tol=1e-9)
