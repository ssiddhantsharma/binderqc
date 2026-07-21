"""Run the terminal-accessibility scorer on Modal.

Imports the installed `terminal_accessibility` package (no inlined copy -- the
science lives in one place). Structure files stay local; each file's bytes are
shipped to a lightweight CPU container that runs the scorer and returns rows.
Useful for fanning a big directory of predictions out in parallel. No GPU, no
secrets, no volumes.

Setup (once):   pip install -e .
Run:            modal run modal_app.py --inputs "path/to/preds" --binder-chains B --target-chains A
                modal run modal_app.py --inputs "a.pdb,b.cif,globs/*.cif" --out tags.csv
"""

import csv
import os

import modal

from terminal_accessibility.paths import gather_paths

app = modal.App("terminal-accessibility")

image = (
    modal.Image.debian_slim(python_version="3.12")
    # Pinned to the versions the scorer was tested against, so a fresh build on
    # anyone's machine reproduces the same numbers.
    .pip_install("biotite==1.0.1", "numpy==2.1.0")
    # Ship the local package source into the container so it can be imported.
    .add_local_python_source("terminal_accessibility")
)


@app.function(image=image)
def score_one(name, data, binder_chains, target_chains, interface_cutoff, exposure_cutoff):
    """Score a single structure from its raw bytes. Returns a list of row dicts."""
    import tempfile

    from terminal_accessibility.core import score_structure

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
        print(f"  {r.get('pdb', ''):<40} tag={tag:<4} {note}")
    print(f"\nWrote {len(rows)} rows -> {out}")
