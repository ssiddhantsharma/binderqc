"""Command-line entry point: `binderqc`.

    binderqc --binder-chains B --target-chains A \
        --out tag_metrics.csv path/to/preds/*.cif a_directory/

Inputs may be files, globs, or directories (recursively scanned for *.pdb/*.cif).
"""

import argparse
import os
import sys

from .core import score_structure, grippability_consensus
from .paths import gather_paths


def _score_one(task):
    """Worker: score one file, returning rows (or a one-row error). Top-level so it
    is picklable for the process pool."""
    path, binder_chains, target_chains, interface_cutoff, exposure_cutoff, verbose = task
    try:
        return score_structure(path, binder_chains, target_chains,
                               interface_cutoff, exposure_cutoff, verbose=verbose)
    except Exception as e:  # noqa: BLE001 - keep the batch going, record the failure
        return [{"pdb": path, "error": str(e)}]

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
    ap.add_argument("-j", "--jobs", type=int, default=1,
                    help="parallel worker processes over a batch of files (default 1)")
    ap.add_argument("--fasta", default="",
                    help="also write QC-passing binders (no quality warnings) to this FASTA path")
    ap.add_argument("--iara-score", type=float, default=None,
                    help="learned target-side grippability of the epitope (0-100 mean hotspot "
                         "prob, e.g. from IARA); when given, adds a physical-vs-learned "
                         "grippability_consensus column (grippable/flat/disagree)")
    args = ap.parse_args(argv)

    import pandas as pd  # imported here so `--help` works without pandas

    paths = gather_paths(args.inputs)
    if not paths:
        sys.exit("No .pdb/.cif files found in inputs.")
    binder_chains = [c for c in args.binder_chains.split(",") if c]
    target_chains = [c for c in args.target_chains.split(",") if c]

    verbose = args.jobs <= 1  # auto-guess printout only makes sense when sequential
    tasks = [(p, binder_chains, target_chains, args.interface_cutoff, args.exposure_cutoff, verbose)
             for p in paths]
    rows = []
    if args.jobs > 1 and len(tasks) > 1:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=args.jobs) as ex:
            for res in ex.map(_score_one, tasks):   # map preserves input order
                rows.extend(res)
    else:
        for t in tasks:
            rows.extend(_score_one(t))

    # Optional: fold in a learned target-side grippability (e.g. IARA epitope mean,
    # computed by the caller so binderqc stays dependency- and license-clean) and
    # report where physical and learned agree. Only touches the CLI output, never
    # the core score_structure schema.
    if args.iara_score is not None:
        for r in rows:
            if "error" in r:
                continue
            r["iara_grippability"] = round(args.iara_score, 1)
            r["grippability_consensus"] = grippability_consensus(r, args.iara_score)["consensus"]

    df = pd.DataFrame(rows)
    df.to_csv(args.out, index=False)
    shown = df.drop(columns=["binder_sequence"], errors="ignore")  # too wide to print
    with pd.option_context("display.max_columns", None, "display.width", 220):
        print(shown.to_string(index=False))
    print(f"\nWrote {len(df)} rows -> {args.out}")

    if args.fasta:
        clean = [r for r in rows
                 if not r.get("error") and r.get("qc_pass") and r.get("binder_sequence")]
        with open(args.fasta, "w") as fh:
            for r in clean:
                fh.write(f">{os.path.splitext(r['pdb'])[0]}|{r['binder_chain']}\n{r['binder_sequence']}\n")
        scored = sum(1 for r in rows if not r.get("error"))
        print(f"Wrote {len(clean)}/{scored} QC-passing binders -> {args.fasta}")


if __name__ == "__main__":
    main()
