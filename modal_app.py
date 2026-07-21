"""Standalone Modal app for the terminal-accessibility scorer.

Single file, no companion module needed. Structure files stay local; each file's
bytes are shipped to a lightweight CPU container that runs the scorer and returns
rows. Useful when you have a big directory of predictions and want the
SASA/interface geometry fanned out in parallel instead of chewing through them
one at a time on a laptop. No GPU, no secrets, no volumes.

    modal run modal_app.py --inputs "path/to/preds" --binder-chains B --target-chains A
    modal run modal_app.py --inputs "a.pdb,b.cif,globs/*.cif" --out tags.csv

The scoring logic between the BEGIN/END SCORER markers below is copied verbatim
from `terminal_accessibility.py` (the standalone local CLI). This file does NOT
import that one -- it is fully self-contained. If you change the science, change
it in both files so the local and Modal paths stay identical.
"""

import csv
import glob
import os

import modal

app = modal.App("terminal-accessibility")

image = (
    modal.Image.debian_slim(python_version="3.12")
    # Pinned to the versions the scorer was tested against, so a fresh build on
    # anyone's machine reproduces the same numbers.
    .pip_install("biotite==1.0.1", "numpy==2.1.0")
)

# Heavy imports load ONLY in the container (skipped locally), so `modal run`
# needs nothing but `modal` installed on your machine.
with image.imports():
    import numpy as np
    import biotite.structure as struc
    import biotite.structure.io as strucio


# ============================ BEGIN SCORER ============================
# Verbatim from terminal_accessibility.py -- keep in sync.

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

# ============================= END SCORER =============================


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


@app.function(image=image)
def score_one(name, data, binder_chains, target_chains, interface_cutoff, exposure_cutoff):
    """Score a single structure from its raw bytes. Returns a list of row dicts."""
    import tempfile

    suffix = os.path.splitext(name)[1] or ".pdb"
    with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as fh:
        fh.write(data)
        path = fh.name
    try:
        rows = score_structure(path, binder_chains, target_chains, interface_cutoff, exposure_cutoff)
    except Exception as e:  # noqa: BLE001 - keep the batch going, record the failure
        rows = [{"pdb": name, "error": str(e)}]
    finally:
        os.unlink(path)
    # Preserve the original filename (the temp file has a random name).
    for r in rows:
        r["pdb"] = name
    return rows


@app.local_entrypoint()
def main(
    inputs: str,
    binder_chains: str = "",
    target_chains: str = "",
    interface_cutoff: float = 5.0,
    exposure_cutoff: float = 0.25,
    out: str = "terminal_accessibility.csv",
):
    paths = gather_paths([s for s in inputs.split(",") if s])
    if not paths:
        raise SystemExit("No .pdb/.cif files found in inputs.")

    bc = [c for c in binder_chains.split(",") if c]
    tc = [c for c in target_chains.split(",") if c]

    args = []
    for p in paths:
        with open(p, "rb") as fh:
            args.append((os.path.basename(p), fh.read(), bc, tc, interface_cutoff, exposure_cutoff))

    print(f"Scoring {len(args)} structures on Modal...")
    rows = [r for sub in score_one.starmap(args) for r in sub]

    # Union of keys, preserving first-seen order, so error rows don't drop columns.
    fieldnames = []
    for r in rows:
        for k in r:
            if k not in fieldnames:
                fieldnames.append(k)

    with open(out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, restval="")
        writer.writeheader()
        writer.writerows(rows)

    for r in rows:
        tag = r.get("recommended_tag", "")
        note = r.get("error") or r.get("warnings") or ""
        print(f"  {r.get('pdb',''):<40} tag={tag:<4} {note}")
    print(f"\nWrote {len(rows)} rows -> {out}")
