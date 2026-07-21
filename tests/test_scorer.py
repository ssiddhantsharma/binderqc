"""Tests for the binderqc scorer.

Example data: PDB 7JZU -- the de novo designed minibinder LCB1 (chain A, 55 aa)
bound to the SARS-CoV-2 RBD (chain B), from Cao et al., Science 2020. LCB1 is a
three-helix bundle whose termini both sit near the binding face, so it is a good
worked example of the "ambiguous" case. RCSB coordinate files are public domain.
"""

import math
from pathlib import Path

import pytest

from binderqc import score_structure

FIXTURE = Path(__file__).parent / "data" / "7JZU_LCB1_RBD.pdb"

EXPECTED_COLUMNS = {
    "pdb", "binder_chain", "target_chains", "binder_len", "n_interface_res", "binder_bsa",
    "approach_angle", "epitope_planarity", "epitope_hydrophobic_frac", "epitope_aromatic_n",
    "nterm_resnum", "nterm_resname", "nterm_relsasa", "nterm_dist_to_interface",
    "nterm_orientation", "nterm_sg_sasa",
    "cterm_resnum", "cterm_resname", "cterm_relsasa", "cterm_dist_to_interface",
    "cterm_orientation", "cterm_sg_sasa",
    "recommended_tag", "mw", "gravy", "net_charge_ph74", "pi", "ext_coeff_280",
    "sequence_liabilities", "warnings", "qc_pass", "binder_sequence",
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
    assert len(row["binder_sequence"]) == 55


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


def test_pose_and_grippability(row):
    assert 0.0 <= row["approach_angle"] <= 90.0     # undirected axis angle
    assert row["epitope_planarity"] >= 0.0          # RMSD to best-fit plane


def test_sequence_liabilities_is_a_string(row):
    # May be empty; must never be missing.
    assert isinstance(row["sequence_liabilities"], str)


def test_expression_signals_in_range(row):
    assert -4.5 <= row["gravy"] <= 4.5              # Kyte-Doolittle bounds
    assert 0.0 <= row["pi"] <= 14.0


def test_charge_and_pi_logic():
    # Polybasic -> high pI + positive charge; polyacidic -> low pI + negative.
    from binderqc.core import _net_charge, _isoelectric_point
    assert _isoelectric_point("K" * 10) > 9.0
    assert _isoelectric_point("E" * 10) < 5.0
    assert _net_charge("K" * 10) > 0
    assert _net_charge("E" * 10) < 0


def test_protparam_formulas():
    # Check the formulas self-consistently + a physical MW range (no memorized
    # external numbers, which are the classic fabrication trap).
    from binderqc.core import _protparam
    seq = "FVNQHLCGSHLVEALYLVCGERGFFYTPKT"   # 30-mer
    pp = _protparam(seq)
    assert 30 * 100 < pp["mw"] < 30 * 140                       # ~110 Da/residue
    assert pp["ext_coeff_280"] == 5500 * seq.count("W") + 1490 * seq.count("Y")


def test_epitope_composition_reported(row):
    assert 0.0 <= row["epitope_hydrophobic_frac"] <= 1.0
    assert row["epitope_aromatic_n"] >= 0


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
    from binderqc.core import _residue_relsasa

    atom = struc.Atom([0.0, 0.0, 0.0], chain_id="A", res_id=1,
                      res_name="PCA", atom_name="CA", element="C")
    arr = struc.array([atom])
    rel = _residue_relsasa(arr, np.array([12.3]))
    assert math.isnan(rel[("A", 1)])


def test_qc_pass_ignores_tag_site_warnings(row):
    # LCB1's only warning is the ambiguous tag site (an advisory about where to
    # tag, not a quality problem), so it must still pass QC.
    assert "ambiguous" in row["warnings"]
    assert row["qc_pass"] is True


def test_cif_input_matches_pdb(tmp_path, row):
    # The bundled fixture is a PDB; confirm the CIF path gives the same numbers.
    import biotite.structure.io as strucio
    arr = strucio.load_structure(str(FIXTURE))
    cif = tmp_path / "7jzu.cif"
    strucio.save_structure(str(cif), arr)
    cif_row = score_structure(str(cif), binder_chains=["A"], target_chains=["B"])[0]
    assert cif_row["binder_bsa"] == row["binder_bsa"]
    assert cif_row["recommended_tag"] == row["recommended_tag"]
    assert cif_row["binder_sequence"] == row["binder_sequence"]


def test_multiple_binder_chains_give_multiple_rows():
    rows = score_structure(str(FIXTURE), binder_chains=["A", "B"])
    assert len(rows) == 2
    assert {r["binder_chain"] for r in rows} == {"A", "B"}


def test_cli_writes_csv_and_fasta(tmp_path):
    from binderqc.cli import main
    out, fa = tmp_path / "out.csv", tmp_path / "clean.fasta"
    main(["--binder-chains", "A", "--target-chains", "B",
          "--out", str(out), "--fasta", str(fa), str(FIXTURE)])
    assert out.exists() and out.read_text().count("\n") >= 2      # header + >=1 row
    assert fa.read_text().startswith(">")                          # LCB1 passes QC
