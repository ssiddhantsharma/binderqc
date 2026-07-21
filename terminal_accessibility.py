#!/usr/bin/env python
"""Terminal (N-/C-terminus) accessibility scorer for protein binder complexes.

Given predicted binder--target complex structures (PDB or CIF), this decides,
for each binder chain in each structure, which terminus is the better place to
attach a purification / immobilization / conjugation tag. Pure geometry off the
structures you already have -- no folding, no GPU, no network.

A tag site should be (1) solvent-exposed so it does not disrupt the fold, and
(2) FAR from the binding interface so that tethering through it does not occlude
the paratope when the binder is immobilized on a sensor surface / displayed.
So three numbers per terminus matter:

  * relative SASA of the terminal residue (exposure; Tien et al. 2013 max-ASA ref)
  * distance from the terminal CA to the nearest binder interface (paratope) CA
  * orientation: does the chain extend TOWARD the paratope (a tag would head at
    the target, >0) or AWAY from it (<0)? -- cosine of the terminus's outward
    chain-extension direction against the direction to the paratope centroid.

`recommended_tag` picks the terminus farther from the interface, with a warning
when that terminus is buried, when the two termini are ~equally placed, or when
the recommended terminus points back toward the interface.

It also reports the terminal residue's exposure and the SG-SASA of any terminal
cysteine, so you can judge site-directed Cys conjugatability at a glance.

Chain convention
----------------
Binder/target chains are given EXPLICITLY as comma-separated sets. This is
deliberate -- no length heuristic is reliable across target types (in a pMHC
complex the peptide is the shortest chain; in other complexes the binder can be
the larger chain). So:

  * `--binder-chains` / `--target-chains` are authoritative when given.
  * `--target-chains` defaults to "all chains that are not binder chains".
  * If `--binder-chains` is omitted the script AUTO-GUESSES and prints the
    guess for every file (never silent), and additionally warns whenever a
    declared binder chain is the largest chain in the complex (a strong sign
    the binder/target assignment is flipped).

Output is one row per (structure, binder_chain).

Usage
-----
    uv run --python 3.12 --with biotite --with pandas --with numpy \
        python terminal_accessibility.py \
        --binder-chains B --target-chains A \
        --out tag_metrics.csv \
        path/to/preds/*.cif a_directory/

  --binder-chains     comma-separated binder chain ids (default: auto-guess)
  --target-chains     comma-separated target chains (default: all non-binder chains)
  --interface-cutoff  heavy-atom dist (A) to call a binder residue "interface" (default 5.0)
  --exposure-cutoff   relSASA below which a terminus is "buried" (default 0.25)
  --out               output CSV path (default: terminal_accessibility.csv)

Inputs may be files, globs, or directories (recursively scanned for *.pdb/*.cif).
"""

import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd
import biotite.structure as struc
import biotite.structure.io as strucio

# Tien et al. 2013, "theoretical" maximum accessible surface area (A^2) per
# residue. Used to normalize total-residue SASA into a 0-1 relative exposure.
_REF_MAX_ASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLU": 223.0, "GLN": 225.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}

# Auto-guess: a binder chain is a single chain in this residue-length window.
# Excludes short peptides below and large targets above. Heuristic only --
# the printed guess + "largest chain" warning are the real safety net.
_AUTO_MIN_LEN = 20
_AUTO_MAX_LEN = 250


def _load_protein(path):
    """Load a structure file, keep amino acids of the first model only."""
    array = strucio.load_structure(path)
    if isinstance(array, struc.AtomArrayStack):
        array = array[0]
    return array[struc.filter_amino_acids(array)]


def _chain_order(array):
    """Chain ids in order of first appearance (N->C along the file)."""
    seen = []
    for c in array.chain_id:
        if c not in seen:
            seen.append(c)
    return seen


def _chain_lengths(array):
    """{chain_id: residue count}, in chain order."""
    starts = struc.get_residue_starts(array)
    lengths = {}
    for s in starts:
        c = str(array.chain_id[s])
        lengths[c] = lengths.get(c, 0) + 1
    return lengths


def _guess_binder_chains(chain_lens):
    """Best single-chain binder guess: shortest chain within the length window."""
    candidates = {c: n for c, n in chain_lens.items() if _AUTO_MIN_LEN <= n <= _AUTO_MAX_LEN}
    if not candidates:
        return []
    return [min(candidates, key=candidates.get)]


def _chain_res_ids(array, chain_id):
    """Ordered res_ids for a chain, following chain order (N->C), not numbering."""
    starts = struc.get_residue_starts(array)
    return [int(array.res_id[s]) for s in starts if array.chain_id[s] == chain_id]


def _termini_res_ids(array, chain_id):
    """(N-term res_id, C-term res_id) following chain order, not numbering."""
    ids = _chain_res_ids(array, chain_id)
    return (ids[0], ids[-1]) if ids else (None, None)


def _residue_relsasa(array, atom_sasa):
    """Per-atom -> per-residue total SASA, normalized by Tien max-ASA.

    Returns {(chain_id, res_id): relative SASA}. SASA is computed on the whole
    complex so interface burial is reflected.
    """
    res_sasa = struc.apply_residue_wise(array, atom_sasa, np.sum)
    res_starts = struc.get_residue_starts(array)
    rel = {}
    for res_total, start in zip(res_sasa, res_starts):
        ref = _REF_MAX_ASA.get(array.res_name[start], 0.0)
        rel[(array.chain_id[start], int(array.res_id[start]))] = (
            (res_total / ref) if ref > 0 else float("nan")
        )
    return rel


def _sg_sasa(array, atom_sasa, chain_id, res_id):
    """Absolute SASA (A^2) of a cysteine SG atom, for conjugatability. NaN if none."""
    mask = (array.chain_id == chain_id) & (array.res_id == res_id) & (array.atom_name == "SG")
    return float(atom_sasa[mask][0]) if mask.any() else float("nan")


def _ca_coord(array, chain_id, res_id):
    mask = (array.chain_id == chain_id) & (array.res_id == res_id) & (array.atom_name == "CA")
    return array.coord[mask][0] if mask.any() else None


def _interface_residue_ids(array, binder_chain, target_chains, cutoff):
    """Binder res_ids with any heavy atom within `cutoff` A of any target atom."""
    binder_mask = array.chain_id == binder_chain
    target_mask = np.isin(array.chain_id, list(target_chains))
    if not binder_mask.any() or not target_mask.any():
        return set()
    b_coord, t_coord = array.coord[binder_mask], array.coord[target_mask]
    b_resid = array.res_id[binder_mask]
    interface = set()
    # Chunk over binder atoms to bound memory on large targets.
    for i in range(0, len(b_coord), 2048):
        chunk = b_coord[i : i + 2048]
        d = np.linalg.norm(chunk[:, None, :] - t_coord[None, :, :], axis=-1)
        hit = d.min(axis=1) <= cutoff
        interface.update(int(r) for r in b_resid[i : i + 2048][hit])
    return interface


def _paratope_centroid(array, binder_chain, interface_ids):
    """Mean CA of the binder interface (paratope) residues. None if unavailable."""
    cas = [c for c in (_ca_coord(array, binder_chain, r) for r in interface_ids) if c is not None]
    return np.mean(cas, axis=0) if cas else None


def _terminus_orientation(term_ca, adj_ca, paratope_centroid):
    """Cosine of the angle between the terminus's outward chain-extension direction
    and the direction from the terminus toward the paratope centroid.

    A tag continues the chain outward, along (term_ca - adjacent_ca). If that
    direction points toward the paratope the appendage heads at the target.

      > 0  terminus extends TOWARD the interface (bad, even if far)
      < 0  terminus extends AWAY from the interface (good)
      NaN  geometry unavailable (single-residue chain, missing CA, no interface)
    """
    if term_ca is None or adj_ca is None or paratope_centroid is None:
        return float("nan")
    ext = term_ca - adj_ca
    to_iface = paratope_centroid - term_ca
    n_ext, n_if = np.linalg.norm(ext), np.linalg.norm(to_iface)
    if n_ext == 0 or n_if == 0:
        return float("nan")
    return float(np.dot(ext, to_iface) / (n_ext * n_if))


def _min_dist_to_interface(array, binder_chain, term_res_id, interface_ids):
    """Min CA-CA distance from terminal residue to any interface residue."""
    if not interface_ids:
        return float("nan")
    term_ca = _ca_coord(array, binder_chain, term_res_id)
    if term_ca is None:
        return float("nan")
    iface_cas = np.array([c for c in (_ca_coord(array, binder_chain, r) for r in interface_ids) if c is not None])
    if len(iface_cas) == 0:
        return float("nan")
    return float(np.linalg.norm(iface_cas - term_ca, axis=1).min())


def _score_binder_chain(array, atom_sasa, name, binder_chain, target_chains, chain_lens,
                        relsasa, interface_cutoff, exposure_cutoff):
    ordered_ids = _chain_res_ids(array, binder_chain)
    nterm_id, cterm_id = (ordered_ids[0], ordered_ids[-1]) if ordered_ids else (None, None)
    # Residue one in from each terminus, for the chain-extension direction.
    n_adj_id = ordered_ids[1] if len(ordered_ids) >= 2 else None
    c_adj_id = ordered_ids[-2] if len(ordered_ids) >= 2 else None

    def _resname(rid):
        m = (array.chain_id == binder_chain) & (array.res_id == rid)
        return str(array.res_name[m][0]) if m.any() else ""

    interface_ids = _interface_residue_ids(array, binder_chain, target_chains, interface_cutoff)
    n_rel = relsasa.get((binder_chain, nterm_id), float("nan"))
    c_rel = relsasa.get((binder_chain, cterm_id), float("nan"))
    n_dist = _min_dist_to_interface(array, binder_chain, nterm_id, interface_ids)
    c_dist = _min_dist_to_interface(array, binder_chain, cterm_id, interface_ids)

    # Orientation: does the terminus extend toward (>0) or away from (<0) the paratope?
    paratope = _paratope_centroid(array, binder_chain, interface_ids)
    n_orient = _terminus_orientation(_ca_coord(array, binder_chain, nterm_id),
                                     _ca_coord(array, binder_chain, n_adj_id), paratope)
    c_orient = _terminus_orientation(_ca_coord(array, binder_chain, cterm_id),
                                     _ca_coord(array, binder_chain, c_adj_id), paratope)

    # SG-SASA only meaningful when the terminal residue is a cysteine.
    n_sg = _sg_sasa(array, atom_sasa, binder_chain, nterm_id) if _resname(nterm_id) == "CYS" else float("nan")
    c_sg = _sg_sasa(array, atom_sasa, binder_chain, cterm_id) if _resname(cterm_id) == "CYS" else float("nan")

    warnings = []
    if chain_lens.get(binder_chain) == max(chain_lens.values()):
        warnings.append("binder is the LARGEST chain -- binder/target may be flipped")

    if not interface_ids:
        recommended = "N/A"
        warnings.append("no interface residues found (check target chains / cutoff)")
    else:
        recommended = "C" if (np.nan_to_num(c_dist) >= np.nan_to_num(n_dist)) else "N"
        rec_dist, other_dist = (c_dist, n_dist) if recommended == "C" else (n_dist, c_dist)
        rec_rel = c_rel if recommended == "C" else n_rel
        rec_orient = c_orient if recommended == "C" else n_orient
        if np.isfinite(rec_rel) and rec_rel < exposure_cutoff:
            warnings.append(f"recommended {recommended}-term is buried (relSASA={rec_rel:.2f})")
        if np.isfinite(rec_dist) and np.isfinite(other_dist) and abs(rec_dist - other_dist) < 5.0:
            warnings.append("both termini ~equidistant from interface (ambiguous)")
        if np.isfinite(rec_dist) and rec_dist < 8.0:
            warnings.append(f"recommended terminus is close to interface ({rec_dist:.1f} A)")
        if np.isfinite(rec_orient) and rec_orient > 0.5:
            warnings.append(f"recommended {recommended}-term points toward interface (orientation={rec_orient:.2f})")

    return {
        "pdb": name,
        "binder_chain": binder_chain,
        "target_chains": ",".join(target_chains),
        "binder_len": chain_lens.get(binder_chain, 0),
        "n_interface_res": len(interface_ids),
        "nterm_resnum": nterm_id,
        "nterm_resname": _resname(nterm_id),
        "nterm_relsasa": round(n_rel, 3) if np.isfinite(n_rel) else float("nan"),
        "nterm_dist_to_interface": round(n_dist, 2) if np.isfinite(n_dist) else float("nan"),
        "nterm_orientation": round(n_orient, 2) if np.isfinite(n_orient) else float("nan"),
        "nterm_sg_sasa": round(n_sg, 2) if np.isfinite(n_sg) else float("nan"),
        "cterm_resnum": cterm_id,
        "cterm_resname": _resname(cterm_id),
        "cterm_relsasa": round(c_rel, 3) if np.isfinite(c_rel) else float("nan"),
        "cterm_dist_to_interface": round(c_dist, 2) if np.isfinite(c_dist) else float("nan"),
        "cterm_orientation": round(c_orient, 2) if np.isfinite(c_orient) else float("nan"),
        "cterm_sg_sasa": round(c_sg, 2) if np.isfinite(c_sg) else float("nan"),
        "recommended_tag": recommended,
        "warnings": "; ".join(warnings),
    }


def score_structure(path, binder_chains, target_chains, interface_cutoff, exposure_cutoff):
    array = _load_protein(path)
    name = os.path.basename(path)
    chain_lens = _chain_lengths(array)
    present = set(chain_lens)

    if binder_chains:
        binders = [c for c in binder_chains if c in present]
        if len(binders) < len(binder_chains):
            missing = sorted(set(binder_chains) - present)
            return [{"pdb": name, "error": f"binder chains {missing} not in {sorted(present)}"}]
    else:
        binders = _guess_binder_chains(chain_lens)
        if not binders:
            return [{"pdb": name, "error": f"could not auto-guess a binder chain among {chain_lens}"}]
        print(f"[auto] {name}: binder guess = {binders[0]} (chain lengths {chain_lens})")

    # Compute complex SASA once; reuse for both per-residue relSASA and SG-SASA.
    atom_sasa = np.nan_to_num(struc.sasa(array), nan=0.0)
    relsasa = _residue_relsasa(array, atom_sasa)
    rows = []
    for bc in binders:
        if target_chains:
            targets = [c for c in target_chains if c in present]
        else:
            targets = sorted(present - set(binders))
        rows.append(_score_binder_chain(
            array, atom_sasa, name, bc, targets, chain_lens, relsasa, interface_cutoff, exposure_cutoff
        ))
    return rows


def gather_paths(inputs):
    """Expand files, globs, and directories into a sorted unique list of .pdb/.cif paths."""
    paths = []
    for item in inputs:
        if os.path.isdir(item):
            for root, _, files in os.walk(item):
                paths.extend(os.path.join(root, f) for f in files if f.lower().endswith((".pdb", ".cif")))
        else:
            paths.extend(glob.glob(item))
    return sorted(set(paths))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="PDB/CIF files, globs, or directories")
    ap.add_argument("--binder-chains", default="", help="comma-separated; default = auto-guess")
    ap.add_argument("--target-chains", default="", help="comma-separated; default = all non-binder chains")
    ap.add_argument("--interface-cutoff", type=float, default=5.0)
    ap.add_argument("--exposure-cutoff", type=float, default=0.25)
    ap.add_argument("--out", default="terminal_accessibility.csv")
    args = ap.parse_args()

    paths = gather_paths(args.inputs)
    if not paths:
        sys.exit("No .pdb/.cif files found in inputs.")
    binder_chains = [c for c in args.binder_chains.split(",") if c]
    target_chains = [c for c in args.target_chains.split(",") if c]

    rows = []
    for p in paths:
        try:
            rows.extend(score_structure(p, binder_chains, target_chains, args.interface_cutoff, args.exposure_cutoff))
        except Exception as e:  # noqa: BLE001 - keep batch going, record failure
            rows.append({"pdb": os.path.basename(p), "error": str(e)})

    df = pd.DataFrame(rows)
    df.to_csv(args.out, index=False)
    with pd.option_context("display.max_columns", None, "display.width", 220):
        print(df.to_string(index=False))
    print(f"\nWrote {len(df)} rows -> {args.out}")


if __name__ == "__main__":
    main()
