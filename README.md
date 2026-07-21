# binderqc

[![CI](https://github.com/ssiddhantsharma/binderqc/actions/workflows/ci.yml/badge.svg)](https://github.com/ssiddhantsharma/binderqc/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)

Quality control and tag-site scoring for designed protein binders, computed
straight from a predicted binder–target complex. Geometry and sequence only, one
CSV row per binder chain.

## Install

```bash
pip install -e .          # or: uv pip install -e .
```

Python 3.10+. Pulls in `biotite`, `numpy`, and `pandas`.

## Usage

Command line:

```bash
binderqc --binder-chains A --target-chains B --out out.csv complex.cif some_dir/
```

Python:

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
| `--fasta` | off | also write binders with no QC warnings to this FASTA |

## What it reports

Per binder chain:

- **Tag site** — recommended terminus (N/C) plus the numbers behind it: terminal
  relative SASA, CA–CA distance to the paratope, orientation, and the SG SASA of
  a terminal cysteine.
- **Interface** — buried surface area and interface residue count.
- **Pose** — approach angle (end-on vs. lying across the surface).
- **Grippability** — epitope planarity, hydrophobic fraction, aromatic anchor count.
- **Developability** — sequence liabilities, GRAVY, net charge, pI, MW, ε₂₈₀.

A `warnings` column flags buried/ambiguous/interface-facing tag sites, flat or
anchorless epitopes, small interfaces, hydrophobic sequences, and swapped chains.

Columns: `pdb, binder_chain, target_chains, binder_len, n_interface_res,
binder_bsa, approach_angle, epitope_planarity, epitope_hydrophobic_frac,
epitope_aromatic_n, nterm_*, cterm_*, recommended_tag, mw, gravy, net_charge_ph74,
pi, ext_coeff_280, sequence_liabilities, warnings, qc_pass, binder_sequence`.
`qc_pass` is true when there are no *quality* warnings (tag-site advisories like an
ambiguous terminus don't count). Pass `--fasta out.fasta` to also dump the
`qc_pass` binders.

## Tests

```bash
pip install -e ".[test]"
pytest
```

Runs against a bundled example, PDB 7JZU (the LCB1 minibinder on the SARS-CoV-2
RBD). `tests/pisa_correctness.py` is a separate script (not part of the unit
tests) that downloads 18 public complexes from RCSB and PDBePISA and checks the
interface area against PISA (r ≈ 1.0, ~1% median error); run it with
`pip install -e ".[validation]" && python tests/pisa_correctness.py`.

## License

MIT
