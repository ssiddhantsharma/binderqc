"""self_association: does the binder's self-interface sit on its paratope?

Builds two homodimers from the bundled LCB1 binder by translating a second copy onto
the paratope face and onto the opposite face, and checks the two are told apart.
"""
import biotite.structure.io.pdb as pdb
import numpy as np
import pytest

from binderqc import score_structure, self_association
from binderqc.core import _interface_residue_ids, _load_protein, _paratope_centroid

COMPLEX = "tests/data/7JZU_LCB1_RBD.pdb"


def _dimer(binder, direction, out, min_res=8):
    """Second copy slid in along `direction` until it first makes a small contact patch."""
    for d in np.arange(60, 6, -0.25):
        copy = binder.copy()
        copy.coord = binder.coord + direction * d
        copy.chain_id = np.full(copy.array_length(), "Z")
        both = binder + copy
        if len(_interface_residue_ids(both, binder.chain_id[0], ["Z"], 5.0)) >= min_res:
            f = pdb.PDBFile()
            f.set_structure(both)
            f.write(out)
            return out
    raise RuntimeError("no contact found")


@pytest.fixture(scope="module")
def row():
    return score_structure(COMPLEX, verbose=False)[0]


@pytest.fixture(scope="module")
def faces(row, tmp_path_factory):
    array = _load_protein(COMPLEX)
    chain = row["binder_chain"]
    binder = array[array.chain_id == chain]
    paratope = {int(r) for r in row["paratope_res"].split(",")}
    axis = _paratope_centroid(array, chain, paratope) - binder.coord.mean(axis=0)
    axis /= np.linalg.norm(axis)
    d = tmp_path_factory.mktemp("dimers")
    return (_dimer(binder, axis, str(d / "para.pdb")),
            _dimer(binder, -axis, str(d / "far.pdb")))


def test_row_exposes_paratope_and_normalised_sap(row):
    assert row["paratope_res"] and len(row["paratope_res"].split(",")) == row["n_interface_res"]
    assert row["sap_per_res"] == pytest.approx(row["sap_total"] / len(row["binder_sequence"]), abs=1e-3)


def test_paratope_face_is_enriched(row, faces):
    r = self_association(faces[0], row)
    assert r["residues_matched"]
    assert r["self_paratope_overlap_n"] > 0
    assert r["self_paratope_enrichment"] > 1.0
    assert r["self_bsa"] > 0


def test_opposite_face_is_clear(row, faces):
    r = self_association(faces[1], row)
    assert r["self_paratope_overlap_frac"] < self_association(faces[0], row)["self_paratope_overlap_frac"]
    assert r["self_verdict"] == "paratope-clear"


def test_needs_two_copies(row):
    assert "error" in self_association(COMPLEX, row)   # binder + target, not binder + binder


def test_unmatched_numbering_is_flagged(faces):
    r = self_association(faces[0], {"paratope_res": "9001,9002"})
    assert not r["residues_matched"] and r["self_verdict"] == "unknown"
