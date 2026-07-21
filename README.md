# terminal-accessibility

Pick the right terminus to tag on a designed protein binder.

Given predicted **binder–target complex** structures (PDB/CIF), for each binder
chain it reports, per terminus (N and C):

- **relative SASA** of the terminal residue (solvent exposure; normalized by the
  Tien et al. 2013 max-ASA reference) — a buried terminus is a bad tag site;
- **CA–CA distance to the nearest interface residue** — tagging near the paratope
  can occlude binding once the binder is immobilized/displayed;
- **orientation** — cosine of the terminus's outward chain-extension direction
  against the direction to the paratope centroid: `>0` the chain (and any tag)
  extends *toward* the target even if far away, `<0` it extends *away* (good);
- **SG-SASA** if the terminal residue is a cysteine — for site-directed Cys conjugation.

It then recommends the terminus **farther from the interface**, and warns when
that terminus is buried, when the two termini are roughly equidistant, when the
recommended terminus points back toward the interface, or when the binder/target
chain assignment looks flipped. Pure geometry — no folding, no GPU, no network.

It also reports the binder's **buried surface area** (`binder_bsa`, Å²) as an
interface-size sanity check — a tiny interface is the strongest single sign of a
spurious binder — and flags interfaces under 300 Å².

## Install

```bash
pip install -e .          # or: uv pip install -e .
```

Python ≥ 3.10. Dependencies (`biotite`, `numpy`, `pandas`) are pulled in
automatically.

## Local usage

```bash
terminal-accessibility \
    --binder-chains B --target-chains A \
    --out tag_metrics.csv \
    path/to/preds/*.cif some_directory/
```

Inputs can be files, globs, or directories (recursively scanned for `*.pdb` / `*.cif`).
Or from Python:

```python
from terminal_accessibility import score_structure
rows = score_structure("complex.pdb", binder_chains=["A"], target_chains=["B"])
```

**Chain convention.** Binder/target chains are given explicitly — no length
heuristic is reliable across target types. If `--binder-chains` is omitted the
binder is auto-guessed (shortest chain in a 20–250-residue window) and the guess
is printed for every file. `--target-chains` defaults to every chain that isn't a
binder chain.

| flag | default | meaning |
|---|---|---|
| `--binder-chains` | auto-guess | comma-separated binder chain ids |
| `--target-chains` | all non-binder | comma-separated target chain ids |
| `--interface-cutoff` | `5.0` | heavy-atom distance (Å) defining an interface residue |
| `--exposure-cutoff` | `0.25` | relSASA below which a terminus is "buried" |
| `--out` | `terminal_accessibility.csv` | output CSV path |

## Run on Modal (optional)

For a large directory, fan the work out over parallel CPU containers with
[Modal](https://modal.com). `modal_app.py` imports the installed package, so
install it into the same environment as `modal`:

```bash
pip install -e ".[modal]"
modal run modal_app.py --inputs "path/to/preds" --binder-chains B --target-chains A --out tags.csv
```

`--inputs` accepts a directory, a glob, or a comma-separated list. Structure
files stay on your machine; only their bytes are sent to the container. No GPU,
no secrets, no volumes.

## Tests

```bash
pip install -e ".[test]"
pytest
```

Tests run against a bundled example: PDB **7JZU** — the de novo designed
minibinder **LCB1** bound to the SARS-CoV-2 RBD (Cao et al., *Science* 2020).

## Layout

```
src/terminal_accessibility/   core.py (scorer) · paths.py · cli.py
modal_app.py                  Modal wrapper (imports the package)
tests/                        test_scorer.py · data/7JZU_LCB1_RBD.pdb
```

## Output

One row per `(structure, binder_chain)`:

`pdb, binder_chain, target_chains, binder_len, n_interface_res, binder_bsa,`
`nterm_resnum, nterm_resname, nterm_relsasa, nterm_dist_to_interface, nterm_orientation, nterm_sg_sasa,`
`cterm_resnum, cterm_resname, cterm_relsasa, cterm_dist_to_interface, cterm_orientation, cterm_sg_sasa,`
`recommended_tag ("N" | "C" | "N/A"), warnings`

## Related work & design notes

The interface-size metric (buried surface area) is a deliberately lightweight,
**PLIP-free** reimplementation of the interface-profiling idea from
[STCRpy](https://doi.org/10.1093/bioinformatics/btaf566) (Bioinformatics, 2025),
which computes richer per-contact interactions but is TCR-specific and depends on
[PLIP](https://github.com/pharmai/plip) (GPL-2.0). This tool keeps to buried-area
geometry computed with biotite so it stays MIT-licensed and dependency-light, and
works for any binder rather than only TCRs. The terminus / tag-site scoring is
distinct from design-time terminus losses (e.g.
[BindCraft](https://github.com/martinpacesa/BindCraft)) in that it is a post-hoc
QC score of an existing predicted complex, not an optimization objective.

## License

MIT

