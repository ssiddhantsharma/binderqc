"""Tests for the terminal-accessibility scorer.

Example data: PDB 7JZU -- the de novo designed minibinder LCB1 (chain A, 55 aa)
bound to the SARS-CoV-2 RBD (chain B), from Cao et al., Science 2020. LCB1 is a
three-helix bundle whose termini both sit near the binding face, so it is a good
worked example of the "ambiguous" case. RCSB coordinate files are public domain.
"""

import math
from pathlib import Path

import pytest

from terminal_accessibility import score_structure

FIXTURE = Path(__file__).parent / "data" / "7JZU_LCB1_RBD.pdb"

EXPECTED_COLUMNS = {
    "pdb", "binder_chain", "target_chains", "binder_len", "n_interface_res", "binder_bsa",
    "nterm_resnum", "nterm_resname", "nterm_relsasa", "nterm_dist_to_interface",
    "nterm_orientation", "nterm_sg_sasa",
    "cterm_resnum", "cterm_resname", "cterm_relsasa", "cterm_dist_to_interface",
    "cterm_orientation", "cterm_sg_sasa",
    "recommended_tag", "warnings",
}


@pytest.fixture(scope="module")
def row():
    """The single scored row for LCB1 (binder A) against the RBD (target B)."""
    rows = score_structure(str(FIXTURE), binder_chains=["A"], target_chains=["B"])
    assert len(rows) == 1
    return rows[0]


def test_schema(row):
    assert set(row) == EXPECTED_COLUMNS
    assert row["binder_chain"] == "A"
    assert row["target_chains"] == "B"
    assert row["binder_len"] == 55


def test_finds_an_interface(row):
    # Pure CA/heavy-atom geometry -> deterministic across biotite versions.
    assert row["n_interface_res"] > 10


def test_interface_bsa_is_substantial(row):
    # LCB1 is a real picomolar binder -> a large buried interface, well clear of
    # the small-interface (300 A^2) warning threshold.
    assert row["binder_bsa"] > 500
    assert "small interface" not in row["warnings"]


def test_recommends_terminus_farther_from_interface(row):
    # LCB1: N-term ~5.3 A, C-term ~9.9 A from the paratope -> C is recommended.
    assert row["recommended_tag"] == "C"
    assert row["cterm_dist_to_interface"] >= row["nterm_dist_to_interface"]


def test_ambiguous_warning_raised(row):
    # The two termini are within 5 A of each other in interface distance.
    assert "ambiguous" in row["warnings"]


def test_sg_sasa_is_nan_for_noncysteine_termini(row):
    # Termini are Asp/Glu, not Cys: SG-SASA must be NaN, never silently 0.
    assert row["nterm_resname"] != "CYS" and row["cterm_resname"] != "CYS"
    assert math.isnan(row["nterm_sg_sasa"])
    assert math.isnan(row["cterm_sg_sasa"])


def test_metrics_are_in_range(row):
    for key in ("nterm_relsasa", "cterm_relsasa"):
        assert 0.0 <= row[key] <= 1.5           # relSASA, occasionally slightly >1
    for key in ("nterm_orientation", "cterm_orientation"):
        assert -1.0 <= row[key] <= 1.0          # a cosine


def test_auto_guess_picks_the_small_chain():
    rows = score_structure(str(FIXTURE), verbose=False)  # no binder chains given
    assert len(rows) == 1
    assert rows[0]["binder_chain"] == "A"       # 55 aa vs 193 aa RBD


def test_missing_binder_chain_returns_error_row_not_crash():
    rows = score_structure(str(FIXTURE), binder_chains=["Z"])
    assert len(rows) == 1
    assert "error" in rows[0]
    assert "recommended_tag" not in rows[0]


def test_nonstandard_residue_relsasa_is_nan():
    # A residue absent from the Tien reference table must yield NaN relSASA,
    # never 0 -- the "unknown, don't guess" contract.
    import numpy as np
    import biotite.structure as struc
    from terminal_accessibility.core import _residue_relsasa

    atom = struc.Atom([0.0, 0.0, 0.0], chain_id="A", res_id=1,
                      res_name="PCA", atom_name="CA", element="C")
    arr = struc.array([atom])
    rel = _residue_relsasa(arr, np.array([12.3]))
    assert math.isnan(rel[("A", 1)])
