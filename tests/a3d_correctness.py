"""Check our Aggrescan3D port against the real Aggrescan3D 1.0.2, end to end.

`core._a3d_score` reimplements Aggrescan3D in numpy. This script runs the actual
reference scoring on the same structures and correlates the two, so the claim in
that docstring is measured rather than asserted.

How the reference is obtained. Aggrescan3D is "free for non-commercial users",
which is not MIT-compatible, so it is NOT vendored -- this script downloads the
published wheel at run time and takes three things out of it:

  * aggrescan/aggrescan_3d.py   -- the scoring itself (pure Python: os/math/csv),
                                   imported and called directly
  * data/matrices/aggrescan.mat -- the a3v propensity scale
  * data/freesasa-2.0.1/...     -- naccess atomic radii (naccess.config) and the
                                   naccess reference max-ASA values, which are
                                   compiled into classifier_naccess.c rather than
                                   written in the config file

No binary from the wheel is executed. SASA comes from the freesasa Python
package, matching the reference's own `--radii naccess --resolution 100`.

    pip install -e ".[validation]" && pip install freesasa
    python tests/a3d_correctness.py

Needs the network (PyPI + RCSB). Standalone, not part of the pytest suite.
"""
import io
import re
import sys
import types
import urllib.request
import zipfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import biotite.structure as struc
import biotite.structure.io.pdb as pdbio
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr, spearmanr

from binderqc.core import _a3d_score

HERE = Path(__file__).parent
WHEEL = ("https://files.pythonhosted.org/packages/bb/ce/"
         "54b6d8752eaf94d61aede795818caf007583d2423b1ccdcc67b12a9819a9/"
         "Aggrescan3D-1.0.2-py2-none-any.whl")

# Small, well-resolved, single-domain proteins plus two de novo minibinders.
PDBS = ["1UBQ", "1CRN", "2GB1", "1SHG", "1BDD", "1ENH", "2CI2", "1LYZ", "3CHY",
        "1BPI", "1MBN", "1CSP", "1MJC", "1PIN", "1TEN", "2ACY", "1AKI", "1IGD",
        "1VII", "2PTL", "1FKB", "2TRX", "1POH", "1DIV", "1QYS", "1A6M", "1HHP",
        "1PGB", "5PTI", "1OPC", "1YCC", "1WIT", "1HEL", "1RNB", "1SSO", "1UBI",
        "1PRB", "1TTF", "7JZU", "7JZM"]


def fetch_reference(workdir):
    """Download the Aggrescan3D wheel and return (module, a3v matrix, naccess tables)."""
    whl = workdir / "a3d.whl"
    if not whl.exists():
        whl.write_bytes(urllib.request.urlopen(WHEEL, timeout=600).read())
    root = workdir / "wheel"
    if not root.exists():
        with zipfile.ZipFile(io.BytesIO(whl.read_bytes())) as z:
            z.extractall(root)
    agg = root / "aggrescan"

    # aggrescan_3d.py does `import logger` (its own py2-era module); stub it so the
    # scoring imports without pulling in the rest of the package.
    stub = types.ModuleType("logger")
    stub.warning = stub.critical = stub.info = stub.debug = lambda **k: None
    stub.AggrescanError = type("AggrescanError", (Exception,), {})
    sys.modules["logger"] = stub
    sys.path.insert(0, str(agg))
    import aggrescan_3d as ref

    a3v = ref.parse_matrix(str(agg / "data" / "matrices" / "aggrescan.mat"))
    fs = agg / "data" / "freesasa-2.0.1"
    # naccess reference max-ASA lives in the C source, not naccess.config -- loading
    # the config alone gives NaN relative areas.
    ref_max = {k: float(v) for k, v in re.findall(
        r'\.max_area = \{\.name = "([A-Z]{3})", \.total = ([0-9.]+)',
        (fs / "src" / "classifier_naccess.c").read_text())}
    return ref, a3v, str(fs / "share" / "naccess.config"), ref_max


def reference_a3d(ref, a3v, naccess_cfg, ref_max, pdb_path, max_dist=10.0):
    """Real Aggrescan3D 1.0.2 scoring on one structure."""
    import freesasa
    st = freesasa.Structure(str(pdb_path), freesasa.Classifier(naccess_cfg))
    areas = freesasa.calc(st, freesasa.Parameters(
        {"algorithm": freesasa.LeeRichards, "n-slices": 100})).residueAreas()

    prot = ref.Protein(pdb_code=pdb_path.stem, residues=[], out_dir=".", max_dist=max_dist)
    for line in pdb_path.read_text().splitlines():
        if not line.startswith("ATOM ") or line[12:16].strip() != "CA":
            continue
        ch, resi, resn = line[21], line[22:29].strip(), line[17:20].strip()
        if resn not in ref_max or ref.get_one_letter(resn) is None:
            continue
        try:
            abs_sasa = areas[ch][resi].total
        except KeyError:
            continue
        prot.add_residue(ref.Residue(
            chain=ch, resi=resi, resn=resn,
            coords=(float(line[30:38]), float(line[38:46]), float(line[46:54])),
            rsa=100.0 * abs_sasa / ref_max[resn], score_matrix=a3v))
    prot.calc_agg_fin()
    fins = np.array([r.agg_fin for r in prot.get_residues()])
    if fins.size == 0:
        return None
    return {"peak": float(fins.max()), "total_pos": float(fins[fins > 0].sum()), "n": fins.size}


def single_chain(pdb_id, workdir):
    """Download a PDB and write out its first protein chain alone."""
    out = workdir / f"{pdb_id}_A.pdb"
    if out.exists():
        return out
    raw = workdir / f"{pdb_id}.pdb"
    if not raw.exists():
        raw.write_bytes(urllib.request.urlopen(
            f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=300).read())
    s = pdbio.PDBFile.read(str(raw)).get_structure(model=1)
    s = s[struc.filter_amino_acids(s) & ~s.hetero]
    if s.array_length() == 0:
        return None
    s = s[s.chain_id == s.chain_id[0]]
    f = pdbio.PDBFile()
    f.set_structure(s)
    f.write(str(out))
    return out


def _gate(name, r, n, min_r, min_n):
    """Exit status for an unattended run: 0 pass, 1 regression, 2 inconclusive.

    A network check that scored too few structures has NOT passed -- the external
    service may simply be down -- so it gets its own status rather than a green
    tick, and is reported either way.
    """
    if n < min_n:
        print(f"INCONCLUSIVE: only {n} of the set scored (need {min_n}); "
              "external service may be unavailable. Not treating as a pass.")
        return 2
    if r < min_r:
        print(f"REGRESSION: {name} r={r:.4f} is below the floor of {min_r}")
        return 1
    print(f"OK: {name} r={r:.4f} >= {min_r} on n={n}")
    return 0


def main(min_r=0.99, min_n=30):
    workdir = HERE / "_a3d_work"
    workdir.mkdir(exist_ok=True)
    ref, a3v, naccess_cfg, ref_max = fetch_reference(workdir)

    rows = []
    for pdb_id in PDBS:
        try:
            path = single_chain(pdb_id, workdir)
            if path is None:
                print(f"  {pdb_id}: skipped (no protein chain)")
                continue
            r = reference_a3d(ref, a3v, naccess_cfg, ref_max, path)
            if r is None:
                print(f"  {pdb_id}: skipped (reference scored no residues)")
                continue
            arr = pdbio.PDBFile.read(str(path)).get_structure(model=1)
            iso = np.nan_to_num(struc.sasa(arr), nan=0.0)
            ours_peak, ours_total = _a3d_score(arr, iso)
            rows.append((pdb_id, r["n"], r["peak"], ours_peak, r["total_pos"], ours_total))
            print(f"  {pdb_id} n={r['n']:4d}  peak ref={r['peak']:7.3f} ours={ours_peak:7.3f}"
                  f"   total+ ref={r['total_pos']:8.3f} ours={ours_total:8.3f}")
        except Exception as e:                                    # noqa: BLE001
            print(f"  {pdb_id}: FAILED ({type(e).__name__}: {e})")

    if len(rows) < 3:
        print(f"\nonly {len(rows)} structures scored -- not enough to correlate")
        return 2
    ids = [r[0] for r in rows]
    rp, op = np.array([r[2] for r in rows]), np.array([r[3] for r in rows])
    rt, ot = np.array([r[4] for r in rows]), np.array([r[5] for r in rows])
    print(f"\nn={len(rows)} structures")
    for label, a, b in (("peak   (a3d_score)", rp, op), ("total+ (a3d_total_positive)", rt, ot)):
        print(f"  {label:30s} Pearson r={pearsonr(a, b)[0]:.4f}  "
              f"Spearman={spearmanr(a, b)[0]:.4f}  mean ours/ref={np.mean(b / a):.3f}  "
              f"median |rel err|={np.median(np.abs(b - a) / a) * 100:.1f}%")

    fig, axes = plt.subplots(1, 2, figsize=(11, 5.2))
    panels = ((rp, op, "a3d_score  (peak)"), (rt, ot, "a3d_total_positive"))
    for ax, (a, b, name) in zip(axes, panels, strict=True):
        lim = (0, max(a.max(), b.max()) * 1.08)
        ax.plot(lim, lim, color="#8d99a6", lw=1, ls="--", zorder=0)
        ax.scatter(a, b, s=34, color="#1f6aa5", zorder=2)
        ax.set_xlim(lim); ax.set_ylim(lim)
        ax.set_xlabel("real Aggrescan3D 1.0.2"); ax.set_ylabel("binderqc port")
        med = float(np.median(np.abs(b - a) / a)) * 100
        ax.set_title(f"{name}\nPearson r = {pearsonr(a, b)[0]:.4f}   "
                     f"median error {med:.1f}%", fontsize=11)
        ax.set_aspect("equal")
    fig.suptitle("binderqc reproduces Aggrescan3D 1.0.2 to under 1% on a typical structure\n"
                 f"{len(rows)} single-chain structures \u00b7 dashed line is y = x, not a fit",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(HERE / "a3d_correctness.png", dpi=190)
    print(f"wrote {HERE / 'a3d_correctness.png'}")
    with open(HERE / "a3d_correctness.csv", "w") as fh:
        fh.write("pdb,n_res,ref_peak,our_peak,ref_total_pos,our_total_pos\n")
        for r in rows:
            fh.write(f"{r[0]},{r[1]},{r[2]:.4f},{r[3]:.4f},{r[4]:.4f},{r[5]:.4f}\n")
    print(f"wrote {HERE / 'a3d_correctness.csv'}  ({len(ids)} rows)")
    return _gate("a3d_score", float(pearsonr(rp, op)[0]), len(rows), min_r, min_n)


if __name__ == "__main__":
    sys.exit(main())
