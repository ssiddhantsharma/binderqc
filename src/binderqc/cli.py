"""Command-line entry point: `binderqc`.

    binderqc --binder-chains B --target-chains A \
        --out tag_metrics.csv path/to/preds/*.cif a_directory/

Inputs may be files, globs, or directories (recursively scanned for *.pdb/*.cif).
"""

import argparse
import os
import sys

from .core import score_structure
from .paths import gather_paths

_DESCRIPTION = """\
Pick the terminus to tag on a designed protein binder.

For each binder chain in a predicted binder-target complex it reports, per
terminus, the relative SASA (exposure), the CA-CA distance to the paratope, the
orientation (does the chain point back at the interface?) and the Cys-SG SASA,
then recommends the terminus farther from the interface and flags buried,
ambiguous or interface-facing cases. Just geometry: no folding, no GPU, no network.
"""


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="binderqc",
        description=_DESCRIPTION,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("inputs", nargs="+", help="PDB/CIF files, globs, or directories")
    ap.add_argument("--binder-chains", default="", help="comma-separated; default = auto-guess")
    ap.add_argument("--target-chains", default="", help="comma-separated; default = all non-binder chains")
    ap.add_argument("--interface-cutoff", type=float, default=5.0,
                    help="heavy-atom dist (A) to call a binder residue interface (default 5.0)")
    ap.add_argument("--exposure-cutoff", type=float, default=0.25,
                    help="relSASA below which a terminus is buried (default 0.25)")
    ap.add_argument("--out", default="binderqc.csv", help="output CSV path")
    ap.add_argument("--fasta", default="",
                    help="also write binders with no QC warnings to this FASTA path")
    args = ap.parse_args(argv)

    import pandas as pd  # imported here so `--help` works without pandas

    paths = gather_paths(args.inputs)
    if not paths:
        sys.exit("No .pdb/.cif files found in inputs.")
    binder_chains = [c for c in args.binder_chains.split(",") if c]
    target_chains = [c for c in args.target_chains.split(",") if c]

    rows = []
    for p in paths:
        try:
            rows.extend(score_structure(p, binder_chains, target_chains,
                                        args.interface_cutoff, args.exposure_cutoff))
        except Exception as e:  # noqa: BLE001 - keep the batch going, record the failure
            rows.append({"pdb": p, "error": str(e)})

    df = pd.DataFrame(rows)
    df.to_csv(args.out, index=False)
    shown = df.drop(columns=["binder_sequence"], errors="ignore")  # too wide to print
    with pd.option_context("display.max_columns", None, "display.width", 220):
        print(shown.to_string(index=False))
    print(f"\nWrote {len(df)} rows -> {args.out}")

    if args.fasta:
        clean = [r for r in rows
                 if not r.get("error") and not r.get("warnings") and r.get("binder_sequence")]
        with open(args.fasta, "w") as fh:
            for r in clean:
                fh.write(f">{os.path.splitext(r['pdb'])[0]}|{r['binder_chain']}\n{r['binder_sequence']}\n")
        scored = sum(1 for r in rows if not r.get("error"))
        print(f"Wrote {len(clean)}/{scored} binders with no warnings -> {args.fasta}")


if __name__ == "__main__":
    main()
