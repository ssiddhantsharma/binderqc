# terminal-accessibility

Work out where to put a tag on a designed protein binder, and get a few
developability QC numbers while you're at it, straight off a predicted
binder-target complex. No folding, no GPU, no network: just geometry and the
sequence.

Give it a predicted **binder-target complex** (PDB or CIF) and for each binder
chain it reports, per terminus (N and C):

- **relative SASA** of the terminal residue (how exposed it is; normalized by the
  Tien et al. 2013 max-ASA reference). A buried terminus is a poor tag site.
- **CA-CA distance to the nearest interface residue.** Tagging next to the
  paratope can block binding once the binder is immobilized or displayed.
- **orientation**: the cosine between the way the chain extends past the terminus
  and the direction to the paratope. Positive means the chain (and any tag) heads
  back toward the target even if the terminus is far away; negative means it
  points away, which is what you want.
- **SG-SASA** of a terminal cysteine, so you can see whether it's a usable handle
  for site-directed conjugation.

It recommends the terminus farther from the interface, and warns when that
terminus is buried, when the two are roughly tied, when it points back at the
target, or when the binder/target chains look swapped.

Per binder chain it also reports a handful of whole-complex QC signals:

- **`binder_bsa`** (Å²): buried surface area. A tiny interface is the clearest
  sign of a junk binder; anything under 300 Å² gets flagged.
- **`approach_angle`** (0-90°): binder long axis vs. the paratope-to-epitope axis.
  Near 0 it points end-on into the target; near 90 it lies across the surface.
- **`epitope_planarity`** (Å): how flat the epitope Cα patch is (RMSD to its
  best-fit plane). Flat means little to grip (flagged under 1 Å); larger is more
  concave.
- **`epitope_hydrophobic_frac`, `epitope_aromatic_n`**: the anchor chemistry on
  the epitope (hydrophobic fraction and aromatic count). A flat, polar, anchorless
  patch is hard to bind, and gets flagged.
- **`sequence_liabilities`**: motifs worth a look (odd cysteine count, N-glyc
  sequon, deamidation, polybasic run, hydrophobic run).
- **expression/solubility**, from the sequence alone (no BioPython): `gravy`
  (Kyte-Doolittle, flagged above 0.4), `net_charge_ph74` and an approximate `pi`
  (Henderson-Hasselbalch), plus the two ProtParam numbers you actually reach for,
  `mw` and `ext_coeff_280` (Pace 1995).

## Install

```bash
pip install -e .          # or: uv pip install -e .
```

Python 3.10+. It pulls in `biotite`, `numpy`, and `pandas`.

## Usage

```bash
terminal-accessibility \
    --binder-chains B --target-chains A \
    --out tag_metrics.csv \
    path/to/preds/*.cif some_directory/
```

Inputs can be files, globs, or directories (walked for `*.pdb` / `*.cif`). From
Python:

```python
from terminal_accessibility import score_structure
rows = score_structure("complex.pdb", binder_chains=["A"], target_chains=["B"])
```

**Chains.** Give the binder and target chains explicitly; no length rule is
reliable across target types. If you leave `--binder-chains` off, it guesses the
binder as the shortest chain in a 20-250 residue window and prints that guess for
every file. `--target-chains` defaults to every other chain.

| flag | default | meaning |
|---|---|---|
| `--binder-chains` | auto-guess | comma-separated binder chain ids |
| `--target-chains` | all non-binder | comma-separated target chain ids |
| `--interface-cutoff` | `5.0` | heavy-atom distance (Å) that counts as a contact |
| `--exposure-cutoff` | `0.25` | relSASA below which a terminus is "buried" |
| `--out` | `terminal_accessibility.csv` | output CSV path |

## Modal (optional, and you probably don't need it)

Each structure is a sub-second CPU calculation, so a few hundred binders finish
in a minute or two on a laptop. Modal is only there for very large batches, where
it runs the same CPU scorer across parallel containers. `modal_app.py` imports the
installed package, so put it in the same environment as `modal`:

```bash
pip install -e ".[modal]"
modal run modal_app.py --inputs "path/to/preds" --binder-chains B --target-chains A --out tags.csv
```

`--inputs` takes a directory, a glob, or a comma-separated list. The structures
stay on your machine; only their bytes go to the container.

## Tests and validation

```bash
pip install -e ".[test]"
pytest
```

The unit tests run against a bundled example, PDB **7JZU** (the de novo minibinder
LCB1 on the SARS-CoV-2 RBD, Cao et al., Science 2020).

`tests/pisa_correctness.py` is a separate, network-fetching benchmark (not part of
the pytest run) that checks the interface area against PISA across 18 public
complexes. It agrees closely (Pearson r ≈ 1.0, slope ≈ 1.0, ~1% median error); see
`tests/pisa_correctness.png`.

```bash
pip install -e ".[validation]"
python tests/pisa_correctness.py
```

## Layout

```
src/terminal_accessibility/   core.py (the scorer) · paths.py · cli.py
modal_app.py                  Modal wrapper (imports the package)
tests/                        test_scorer.py · pisa_correctness.py · data/7JZU_LCB1_RBD.pdb
```

## Output

One row per `(structure, binder_chain)`:

`pdb, binder_chain, target_chains, binder_len, n_interface_res, binder_bsa,`
`approach_angle, epitope_planarity, epitope_hydrophobic_frac, epitope_aromatic_n,`
`nterm_resnum, nterm_resname, nterm_relsasa, nterm_dist_to_interface, nterm_orientation, nterm_sg_sasa,`
`cterm_resnum, cterm_resname, cterm_relsasa, cterm_dist_to_interface, cterm_orientation, cterm_sg_sasa,`
`recommended_tag ("N" | "C" | "N/A"), mw, gravy, net_charge_ph74, pi, ext_coeff_280,`
`sequence_liabilities, warnings`

## License

MIT
