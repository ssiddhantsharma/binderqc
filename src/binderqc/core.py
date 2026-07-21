"""Where to put a tag on a designed binder, plus a few QC numbers, read straight
off a predicted binder-target complex.

For each binder chain we work out which terminus (N or C) is the safer place to
hang a purification/immobilization/conjugation tag. A good tag site is exposed
(so the tag doesn't wreck the fold) and sits well away from the binding interface
(so tethering there doesn't block the paratope once the binder is on a chip or
displayed). So three numbers per terminus:

  - relative SASA of the terminal residue (how exposed it is; Tien 2013 max-ASA)
  - distance from the terminal CA to the nearest paratope CA
  - orientation: does the chain point back toward the paratope (>0) or away (<0)

recommended_tag is just the terminus farther from the interface, and we raise a
warning if that terminus turns out to be buried, roughly tied with the other, or
pointing back at the target. We also report the SG SASA of a terminal cysteine so
you can eyeball whether it's a usable conjugation handle.

Everything here is geometry/sequence off a single structure - no folding, no GPU,
no network.
"""

import os
import re

import numpy as np
import biotite.structure as struc
import biotite.structure.io as strucio

_THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# Tien et al. 2013 theoretical max ASA per residue, used to turn a residue's
# total SASA into a 0-1 exposure.
_REF_MAX_ASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLU": 223.0, "GLN": 225.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}

# When the binder chain isn't given we guess it as the shortest chain in this
# residue window. It's only a fallback - the printed guess and the "largest
# chain" warning are what actually keep you honest.
_AUTO_MIN_LEN = 20
_AUTO_MAX_LEN = 250


def _load_protein(path):
    """Load a structure, first model only, amino acids only."""
    array = strucio.load_structure(path)
    if isinstance(array, struc.AtomArrayStack):
        array = array[0]
    return array[struc.filter_amino_acids(array)]


def _chain_lengths(array):
    starts = struc.get_residue_starts(array)
    lengths = {}
    for s in starts:
        c = str(array.chain_id[s])
        lengths[c] = lengths.get(c, 0) + 1
    return lengths


def _guess_binder_chains(chain_lens):
    candidates = {c: n for c, n in chain_lens.items() if _AUTO_MIN_LEN <= n <= _AUTO_MAX_LEN}
    if not candidates:
        return []
    return [min(candidates, key=candidates.get)]


def _chain_res_ids(array, chain_id):
    """res_ids in chain order (N to C), not by numbering."""
    starts = struc.get_residue_starts(array)
    return [int(array.res_id[s]) for s in starts if array.chain_id[s] == chain_id]


def _residue_relsasa(array, atom_sasa):
    """{(chain, res_id): relative SASA}. SASA is over the whole complex, so a
    residue buried at the interface reads as buried."""
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
    """SASA of a Cys SG atom (conjugatability). NaN if there's no SG here."""
    mask = (array.chain_id == chain_id) & (array.res_id == res_id) & (array.atom_name == "SG")
    return float(atom_sasa[mask][0]) if mask.any() else float("nan")


def _binder_bsa(array, atom_sasa, binder_chain):
    """Area the binder buries on binding = its SASA alone minus its SASA in the
    complex. A tiny interface is the clearest tell of a junk binder. No PLIP,
    just biotite SASA."""
    mask = array.chain_id == binder_chain
    if not mask.any():
        return float("nan")
    iso = np.nan_to_num(struc.sasa(array[mask]), nan=0.0)
    buried = np.clip(iso - atom_sasa[mask], 0.0, None).sum()
    return float(buried)


def _ca_coord(array, chain_id, res_id):
    mask = (array.chain_id == chain_id) & (array.res_id == res_id) & (array.atom_name == "CA")
    return array.coord[mask][0] if mask.any() else None


def _interface_residue_ids(array, binder_chain, target_chains, cutoff):
    """Binder residues with any heavy atom within cutoff of the target."""
    binder_mask = array.chain_id == binder_chain
    target_mask = np.isin(array.chain_id, list(target_chains))
    if not binder_mask.any() or not target_mask.any():
        return set()
    b_coord, t_coord = array.coord[binder_mask], array.coord[target_mask]
    b_resid = array.res_id[binder_mask]
    interface = set()
    # chunk the binder atoms so the pairwise distance matrix stays small on big targets
    for i in range(0, len(b_coord), 2048):
        chunk = b_coord[i : i + 2048]
        d = np.linalg.norm(chunk[:, None, :] - t_coord[None, :, :], axis=-1)
        hit = d.min(axis=1) <= cutoff
        interface.update(int(r) for r in b_resid[i : i + 2048][hit])
    return interface


def _paratope_centroid(array, binder_chain, interface_ids):
    cas = [c for c in (_ca_coord(array, binder_chain, r) for r in interface_ids) if c is not None]
    return np.mean(cas, axis=0) if cas else None


def _terminus_orientation(term_ca, adj_ca, paratope_centroid):
    """Cosine between the way the chain extends past the terminus and the
    direction to the paratope. A tag grows outward along (term_ca - adj_ca); if
    that points at the paratope (>0) the tag heads straight for the target.
    Positive is bad even when the terminus is far away; NaN if we can't tell."""
    if term_ca is None or adj_ca is None or paratope_centroid is None:
        return float("nan")
    ext = term_ca - adj_ca
    to_iface = paratope_centroid - term_ca
    n_ext, n_if = np.linalg.norm(ext), np.linalg.norm(to_iface)
    if n_ext == 0 or n_if == 0:
        return float("nan")
    return float(np.dot(ext, to_iface) / (n_ext * n_if))


def _min_dist_to_interface(array, binder_chain, term_res_id, interface_ids):
    if not interface_ids:
        return float("nan")
    term_ca = _ca_coord(array, binder_chain, term_res_id)
    if term_ca is None:
        return float("nan")
    iface_cas = np.array([c for c in (_ca_coord(array, binder_chain, r) for r in interface_ids) if c is not None])
    if len(iface_cas) == 0:
        return float("nan")
    return float(np.linalg.norm(iface_cas - term_ca, axis=1).min())


def _binder_sequence(array, chain_id):
    """One-letter sequence in chain order; anything non-standard becomes X."""
    starts = struc.get_residue_starts(array)
    return "".join(
        _THREE_TO_ONE.get(str(array.res_name[s]), "X")
        for s in starts if array.chain_id[s] == chain_id
    )


def _sequence_liabilities(seq):
    """Common developability red flags, straight from the sequence. Motifs and
    thresholds are the ones in Adaptyv's protein-qc skill (MIT). These are things
    to look at, not automatic rejections."""
    flags = []
    n_cys = seq.count("C")
    if n_cys % 2 == 1:
        flags.append(f"odd Cys count ({n_cys}): possible unpaired thiol")
    if re.search(r"N[^P][ST]", seq):
        flags.append("N-glycosylation sequon (NxS/T)")
    if re.search(r"N[GST]", seq):
        flags.append("deamidation/isomerization motif (NG/NS/NT)")
    if re.search(r"[KR]{3,}", seq):
        flags.append("polybasic run (>=3 K/R): proteolysis")
    if re.search(r"[AILMFWV]{6,}", seq):
        flags.append("hydrophobic run (>=6): aggregation")
    return flags


# Kyte & Doolittle 1982 hydropathy.
_KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# pKa for side chains and termini (a common simple set) - used for charge and pI.
_PKA_POS = {"K": 10.5, "R": 12.5, "H": 6.0}
_PKA_NEG = {"D": 3.9, "E": 4.1, "C": 8.5, "Y": 10.1}
_PKA_NTERM, _PKA_CTERM = 9.0, 3.1


def _gravy(seq):
    """Mean Kyte-Doolittle hydropathy; higher means more hydrophobic."""
    vals = [_KD[a] for a in seq if a in _KD]
    return float(np.mean(vals)) if vals else float("nan")


# Average residue masses (Da), the usual Expasy values.
_RES_MASS = {
    "A": 71.0788, "R": 156.1875, "N": 114.1038, "D": 115.0886, "C": 103.1388,
    "E": 129.1155, "Q": 128.1307, "G": 57.0519, "H": 137.1411, "I": 113.1594,
    "L": 113.1594, "K": 128.1741, "M": 131.1926, "F": 147.1766, "P": 97.1167,
    "S": 87.0782, "T": 101.1051, "W": 186.2132, "Y": 163.1760, "V": 99.1326,
}
_WATER = 18.01528
_HYDROPHOBIC = set("AILMFWV")
_AROMATIC = set("FWY")


def _protparam(seq):
    """The two ProtParam numbers you actually reach for at the bench: molecular
    weight and the reduced 280 nm extinction coefficient (Pace 1995), i.e. what
    you need to express a construct and read its concentration."""
    if not seq:
        return {"mw": float("nan"), "ext_coeff_280": float("nan")}
    mw = sum(_RES_MASS.get(a, 0.0) for a in seq) + _WATER
    ext = 5500 * seq.count("W") + 1490 * seq.count("Y")   # reduced (no cystines)
    return {"mw": mw, "ext_coeff_280": float(ext)}


def _epitope_composition(array, target_iface):
    """Hydrophobic fraction and aromatic count of the epitope residues - the
    anchor chemistry a binder has to grip. This is about residue identity, so it
    complements epitope_planarity, which is about shape."""
    if not target_iface:
        return float("nan"), 0
    types = []
    for ch, rid in target_iface:
        m = (array.chain_id == ch) & (array.res_id == rid)
        if m.any():
            types.append(_THREE_TO_ONE.get(str(array.res_name[m][0]), "X"))
    if not types:
        return float("nan"), 0
    hyd = sum(1 for t in types if t in _HYDROPHOBIC) / len(types)
    arom = sum(1 for t in types if t in _AROMATIC)
    return float(hyd), int(arom)


def _charge_at_ph(seq, ph):
    """Net charge at a given pH (Henderson-Hasselbalch)."""
    pos = 1.0 / (1.0 + 10 ** (ph - _PKA_NTERM))
    for a, pk in _PKA_POS.items():
        pos += seq.count(a) / (1.0 + 10 ** (ph - pk))
    neg = 1.0 / (1.0 + 10 ** (_PKA_CTERM - ph))
    for a, pk in _PKA_NEG.items():
        neg += seq.count(a) / (1.0 + 10 ** (pk - ph))
    return pos - neg


def _net_charge(seq, ph=7.4):
    return _charge_at_ph(seq, ph) if seq else float("nan")


def _isoelectric_point(seq):
    """pI by bisection. It leans on the pKa set above, so treat it as a guide."""
    if not seq:
        return float("nan")
    lo, hi = 0.0, 14.0
    for _ in range(60):
        mid = (lo + hi) / 2.0
        if _charge_at_ph(seq, mid) > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def _sequence_metrics(seq):
    """Sequence-derived row fields (developability + expression), pre-rounded."""
    pp = _protparam(seq)
    gravy = _gravy(seq)
    return {
        "mw": round(pp["mw"], 1) if np.isfinite(pp["mw"]) else float("nan"),
        "gravy": round(gravy, 3) if np.isfinite(gravy) else float("nan"),
        "net_charge_ph74": round(_net_charge(seq), 2) if seq else float("nan"),
        "pi": round(_isoelectric_point(seq), 2) if seq else float("nan"),
        "ext_coeff_280": pp["ext_coeff_280"],
        "sequence_liabilities": "; ".join(_sequence_liabilities(seq)),
        "binder_sequence": seq,
    }


def _principal_axis(coords):
    """Unit vector along the long axis (direction of most spread)."""
    centered = coords - coords.mean(axis=0)
    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    return vh[0]


def _target_interface_res(array, binder_chain, target_chains, cutoff):
    """(chain, res_id) of target residues in contact with the binder."""
    binder_mask = array.chain_id == binder_chain
    target_mask = np.isin(array.chain_id, list(target_chains))
    if not binder_mask.any() or not target_mask.any():
        return []
    b_coord = array.coord[binder_mask]
    t_coord, t_chain, t_resid = array.coord[target_mask], array.chain_id[target_mask], array.res_id[target_mask]
    hits = set()
    for i in range(0, len(t_coord), 2048):
        chunk = t_coord[i : i + 2048]
        d = np.linalg.norm(chunk[:, None, :] - b_coord[None, :, :], axis=-1)
        hit = d.min(axis=1) <= cutoff
        for ch, rr, h in zip(t_chain[i : i + 2048], t_resid[i : i + 2048], hit):
            if h:
                hits.add((str(ch), int(rr)))
    return sorted(hits)


def _approach_angle(binder_ca, paratope_centroid, epitope_centroid):
    """Angle (0-90 deg) between the binder's long axis and the paratope->epitope
    axis. Near 0 the binder points end-on into the target; near 90 it lies across
    the surface. STCRpy measures a docking angle for TCRs against a fixed MHC
    frame; this is the same idea without needing a reference frame, so it works
    for any binder. NaN if we can't compute it."""
    if len(binder_ca) < 2 or paratope_centroid is None or epitope_centroid is None:
        return float("nan")
    axis = _principal_axis(binder_ca)
    binding = epitope_centroid - paratope_centroid
    n_bind = np.linalg.norm(binding)
    if n_bind == 0:
        return float("nan")
    cos = abs(float(np.dot(axis, binding / n_bind)))   # axis has no direction
    return float(np.degrees(np.arccos(min(cos, 1.0))))


def _planarity_rmsd(coords):
    """RMSD of points to their best-fit plane. Small = a flat epitope patch
    (little to grip); bigger = concave/knobby. NaN below 3 points."""
    if len(coords) < 3:
        return float("nan")
    centered = coords - coords.mean(axis=0)
    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    normal = vh[-1]
    return float(np.sqrt(np.mean((centered @ normal) ** 2)))


def _score_binder_chain(array, atom_sasa, name, binder_chain, target_chains, chain_lens,
                        relsasa, interface_cutoff, exposure_cutoff):
    ordered_ids = _chain_res_ids(array, binder_chain)
    nterm_id, cterm_id = (ordered_ids[0], ordered_ids[-1]) if ordered_ids else (None, None)
    # one residue in from each end, for the chain-extension direction
    n_adj_id = ordered_ids[1] if len(ordered_ids) >= 2 else None
    c_adj_id = ordered_ids[-2] if len(ordered_ids) >= 2 else None

    def _resname(rid):
        m = (array.chain_id == binder_chain) & (array.res_id == rid)
        return str(array.res_name[m][0]) if m.any() else ""

    interface_ids = _interface_residue_ids(array, binder_chain, target_chains, interface_cutoff)
    binder_bsa = _binder_bsa(array, atom_sasa, binder_chain)
    n_rel = relsasa.get((binder_chain, nterm_id), float("nan"))
    c_rel = relsasa.get((binder_chain, cterm_id), float("nan"))
    n_dist = _min_dist_to_interface(array, binder_chain, nterm_id, interface_ids)
    c_dist = _min_dist_to_interface(array, binder_chain, cterm_id, interface_ids)

    paratope = _paratope_centroid(array, binder_chain, interface_ids)
    n_orient = _terminus_orientation(_ca_coord(array, binder_chain, nterm_id),
                                     _ca_coord(array, binder_chain, n_adj_id), paratope)
    c_orient = _terminus_orientation(_ca_coord(array, binder_chain, cterm_id),
                                     _ca_coord(array, binder_chain, c_adj_id), paratope)

    # only meaningful if the terminal residue is actually a cysteine
    n_sg = _sg_sasa(array, atom_sasa, binder_chain, nterm_id) if _resname(nterm_id) == "CYS" else float("nan")
    c_sg = _sg_sasa(array, atom_sasa, binder_chain, cterm_id) if _resname(cterm_id) == "CYS" else float("nan")

    # pose + grippability, both computed straight from coordinates
    binder_ca = array.coord[(array.chain_id == binder_chain) & (array.atom_name == "CA")]
    target_iface = _target_interface_res(array, binder_chain, target_chains, interface_cutoff)
    epitope_ca = np.array([c for c in (_ca_coord(array, ch, r) for ch, r in target_iface) if c is not None])
    epitope_centroid = epitope_ca.mean(axis=0) if len(epitope_ca) else None
    approach = _approach_angle(binder_ca, paratope, epitope_centroid)
    planarity = _planarity_rmsd(epitope_ca)
    epi_hyd_frac, epi_aromatic_n = _epitope_composition(array, target_iface)

    # sequence-side developability + expression
    binder_seq = _binder_sequence(array, binder_chain)
    seqm = _sequence_metrics(binder_seq)

    # quality warnings say the binder/interface itself looks bad; tag-site
    # warnings are only about where to put a tag. qc_pass ignores the latter.
    quality, tagsite = [], []
    if np.isfinite(planarity) and planarity < 1.0:
        quality.append(f"flat epitope (planarity RMSD={planarity:.2f} A): low grippability")
    if np.isfinite(epi_hyd_frac) and epi_hyd_frac < 0.2 and epi_aromatic_n == 0:
        quality.append("polar epitope (few hydrophobic/aromatic anchors): hard to grip")
    if np.isfinite(seqm["gravy"]) and seqm["gravy"] > 0.4:
        quality.append(f"hydrophobic (GRAVY={seqm['gravy']:.2f}): solubility/aggregation risk")
    if chain_lens.get(binder_chain) == max(chain_lens.values()):
        quality.append("binder is the largest chain, so binder/target may be flipped")
    if np.isfinite(binder_bsa) and binder_bsa < 300.0:
        quality.append(f"small interface (binder BSA={binder_bsa:.0f} A^2): possibly weak/spurious")

    if not interface_ids:
        recommended = "N/A"
        quality.append("no interface residues found (check target chains / cutoff)")
    else:
        recommended = "C" if (np.nan_to_num(c_dist) >= np.nan_to_num(n_dist)) else "N"
        rec_dist, other_dist = (c_dist, n_dist) if recommended == "C" else (n_dist, c_dist)
        rec_rel = c_rel if recommended == "C" else n_rel
        rec_orient = c_orient if recommended == "C" else n_orient
        if np.isfinite(rec_rel) and rec_rel < exposure_cutoff:
            tagsite.append(f"recommended {recommended}-term is buried (relSASA={rec_rel:.2f})")
        if np.isfinite(rec_dist) and np.isfinite(other_dist) and abs(rec_dist - other_dist) < 5.0:
            tagsite.append("both termini ~equidistant from interface (ambiguous)")
        if np.isfinite(rec_dist) and rec_dist < 8.0:
            tagsite.append(f"recommended terminus is close to interface ({rec_dist:.1f} A)")
        if np.isfinite(rec_orient) and rec_orient > 0.5:
            tagsite.append(f"recommended {recommended}-term points toward interface (orientation={rec_orient:.2f})")

    warnings = quality + tagsite
    qc_pass = not quality

    return {
        "pdb": name,
        "binder_chain": binder_chain,
        "target_chains": ",".join(target_chains),
        "binder_len": chain_lens.get(binder_chain, 0),
        "n_interface_res": len(interface_ids),
        "binder_bsa": round(binder_bsa, 1) if np.isfinite(binder_bsa) else float("nan"),
        "approach_angle": round(approach, 1) if np.isfinite(approach) else float("nan"),
        "epitope_planarity": round(planarity, 2) if np.isfinite(planarity) else float("nan"),
        "epitope_hydrophobic_frac": round(epi_hyd_frac, 2) if np.isfinite(epi_hyd_frac) else float("nan"),
        "epitope_aromatic_n": epi_aromatic_n,
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
        "mw": seqm["mw"],
        "gravy": seqm["gravy"],
        "net_charge_ph74": seqm["net_charge_ph74"],
        "pi": seqm["pi"],
        "ext_coeff_280": seqm["ext_coeff_280"],
        "sequence_liabilities": seqm["sequence_liabilities"],
        "warnings": "; ".join(warnings),
        "qc_pass": qc_pass,
        "binder_sequence": seqm["binder_sequence"],
    }


def score_structure(path: str, binder_chains: "list[str] | None" = None,
                    target_chains: "list[str] | None" = None,
                    interface_cutoff: float = 5.0, exposure_cutoff: float = 0.25,
                    verbose: bool = True) -> "list[dict]":
    """Score every binder chain in one structure file.

    binder_chains / target_chains are lists of chain ids. Leave binder_chains
    empty to auto-guess it (shortest chain in the length window, printed when
    verbose); leave target_chains empty to use every other chain as target.
    interface_cutoff is the heavy-atom contact distance in A; exposure_cutoff is
    the relSASA below which a terminus counts as buried.

    Returns one dict per binder chain, or a single {"pdb", "error"} dict if the
    requested chains aren't there.
    """
    binder_chains = list(binder_chains) if binder_chains else []
    target_chains = list(target_chains) if target_chains else []

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
        if verbose:
            print(f"[auto] {name}: binder guess = {binders[0]} (chain lengths {chain_lens})")

    # SASA of the whole complex once, reused for relSASA and the SG lookups
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
