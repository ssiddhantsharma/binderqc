<h1 align="center">binderqc</h1>

<p align="center">
  Quality control and tag-site scoring for designed protein binders — geometry and
  sequence only, one CSV row per binder, straight from a predicted complex.
</p>

<p align="center">
  <a href="https://github.com/ssiddhantsharma/binderqc/actions/workflows/ci.yml"><img src="https://github.com/ssiddhantsharma/binderqc/actions/workflows/ci.yml/badge.svg" alt="CI"></a>
  <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT"></a>
  <a href="https://www.python.org/"><img src="https://img.shields.io/badge/python-3.10%2B-blue.svg" alt="Python 3.10+"></a>
</p>

<p align="center">
  <img src="docs/schematic.png" width="860"
       alt="From a predicted binder–target complex, binderqc reports interface, pose, grippability, tag site, and developability">
</p>

## Install

```bash
pip install -e .          # or: uv pip install -e .
```

Python 3.10+. Pulls in `biotite`, `numpy`, and `pandas`.

## Usage

```bash
binderqc --binder-chains A --target-chains B --out out.csv complex.cif some_dir/
```

```python
from binderqc import score_structure
rows = score_structure("complex.pdb", binder_chains=["A"], target_chains=["B"])
```

Inputs are PDB/CIF files, globs, or directories. Leave `--binder-chains` off to
guess the binder as the shortest chain (20–250 aa, printed for each file);
`--target-chains` defaults to the remaining chains.

| flag | default | meaning |
|---|---|---|
| `--binder-chains` | auto-guess | comma-separated binder chain ids |
| `--target-chains` | all non-binder | comma-separated target chain ids |
| `--interface-cutoff` | `5.0` | heavy-atom contact distance (Å) |
| `--exposure-cutoff` | `0.25` | relSASA below which a terminus is buried |
| `--out` | `binderqc.csv` | output CSV path |
| `--fasta` | off | also write the QC-passing binders to this FASTA |

Example output for the bundled LCB1 minibinder (a few of the columns):

| recommended_tag | binder_bsa | epitope_planarity | epitope_aromatic_n | pi | qc_pass |
|---|---|---|---|---|---|
| C | 1021.4 | 3.21 | 11 | 4.17 | True |

Its `warnings` field reads *both termini ~equidistant from interface (ambiguous)* —
a tag-site advisory, so `qc_pass` stays `True`.

## What it reports

Per binder chain:

- **Interface** — buried surface area and interface residue count.
- **Pose** — approach angle (end-on vs. lying across the surface).
- **Grippability** — epitope planarity, hydrophobic fraction, aromatic anchors.
- **Tag site** — recommended terminus (N/C) and the numbers behind it: relative
  SASA, CA–CA distance to the paratope, orientation, and a terminal cysteine's SG SASA.
- **Developability** — sequence liabilities, GRAVY, net charge, pI, MW, ε₂₈₀.

A `warnings` column flags problems (small, flat, or anchorless interfaces;
buried, ambiguous, or interface-facing tag sites; hydrophobic sequences).
`qc_pass` is true when there are no *quality* warnings — tag-site advisories like
an ambiguous terminus don't count — and `--fasta` dumps exactly those binders.

<details>
<summary>Full column list</summary>

`pdb, binder_chain, target_chains, binder_len, n_interface_res, binder_bsa,
approach_angle, epitope_planarity, epitope_hydrophobic_frac, epitope_aromatic_n,
nterm_resnum, nterm_resname, nterm_relsasa, nterm_dist_to_interface,
nterm_orientation, nterm_sg_sasa, cterm_resnum, cterm_resname, cterm_relsasa,
cterm_dist_to_interface, cterm_orientation, cterm_sg_sasa, recommended_tag, mw,
gravy, net_charge_ph74, pi, ext_coeff_280, sequence_liabilities, warnings,
qc_pass, binder_sequence`
</details>

## Tests

```bash
pip install -e ".[test]"
pytest
```

Runs against a bundled example, PDB 7JZU (the LCB1 minibinder on the SARS-CoV-2
RBD). `tests/pisa_correctness.py` is a separate script (not part of the unit
tests) that downloads 18 public complexes from RCSB and PDBePISA and checks the
interface area against PISA (r ≈ 1.0, ~1% median error):

```bash
pip install -e ".[validation]"
python tests/pisa_correctness.py
```

## License

MIT
