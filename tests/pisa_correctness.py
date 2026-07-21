"""Check our interface area against PISA on a set of well-defined 2-chain
complexes, and draw the agreement figure.

The point is correctness, not binding prediction: does the buried-area number we
get out of biotite match what PISA (the usual reference) reports? For each PDB we
ask PDBePISA for its main (type-1) interface and the two chains it's between, pull
those two chains from RCSB, and recompute the interface area the way PISA defines
it, (dSASA_A + dSASA_B) / 2. We also print the tool's own one-sided binder_bsa.

Everything here is public (RCSB + PDBePISA). It's a standalone script, not part of
the pytest suite (it needs the network), so run it directly:

    pip install -e ".[validation]"
    python tests/pisa_correctness.py

The set is enzyme-inhibitor pairs, a couple of antibody-antigen complexes, and
two de novo minibinders (7JZU, 7JZM). They're all cases with a single clean
biological 2-chain interface; homodimers and multi-copy assemblies (where "the
interface" is ambiguous) are deliberately left out.
"""
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr, linregress
import biotite.structure as struc

from binderqc.core import _load_protein, _chain_lengths

HERE = Path(__file__).parent
PDBS = ["2PTC", "3SGB", "1PPF", "1CHO", "1ACB", "1DFJ", "1EMV", "1FSS", "1CSE",
        "1TEC", "1TGS", "1AVW", "1MAH", "1BRS", "7JZU", "7JZM", "1VFB", "1DVF"]


def pisa_main_interface(pdb):
    """PDBePISA's biggest type-1 interface: (chain_a, chain_b, area). None if it
    doesn't hand back two distinct chains."""
    url = f"https://www.ebi.ac.uk/pdbe/pisa/cgi-bin/interfaces.pisa?{pdb.lower()}"
    root = ET.fromstring(urllib.request.urlopen(url, timeout=60).read())
    best = None
    for iface in root.iter("interface"):
        if iface.findtext("type") != "1":
            continue
        chains = [m.findtext("chain_id") for m in iface.iter("molecule") if m.findtext("chain_id")]
        area = float(iface.findtext("int_area"))
        if len(chains) >= 2 and chains[0] != chains[1] and (best is None or area > best[2]):
            best = (chains[0], chains[1], area)
    return best


def interface_area(path, chain_a, chain_b):
    """Interface area our way, on just these two chains: for each chain, SASA
    alone minus SASA in the pair, averaged. Also return the smaller chain's side
    (that's what the tool reports as binder_bsa)."""
    arr = _load_protein(path)
    pair = arr[np.isin(arr.chain_id, [chain_a, chain_b])]
    lens = _chain_lengths(pair)
    if chain_a not in lens or chain_b not in lens:
        return None
    sasa = np.nan_to_num(struc.sasa(pair), nan=0.0)

    def side(chain):
        m = pair.chain_id == chain
        alone = np.nan_to_num(struc.sasa(pair[m]), nan=0.0)
        return float(np.clip(alone - sasa[m], 0, None).sum())

    a, b = side(chain_a), side(chain_b)
    smaller = a if lens[chain_a] <= lens[chain_b] else b
    return {"our_area": (a + b) / 2, "binder_bsa": smaller}


def main():
    rows = []
    for pdb in PDBS:
        mi = pisa_main_interface(pdb)
        if not mi:
            print(f"  {pdb}: skipped (no clean 2-chain interface from PISA)")
            continue
        ca, cb, pisa_area = mi
        pdb_path = HERE / f"{pdb}.pdb"
        if not pdb_path.exists():
            urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pdb}.pdb", pdb_path)
        res = interface_area(str(pdb_path), ca, cb)
        pdb_path.unlink(missing_ok=True)
        if not res:
            print(f"  {pdb}: skipped (chains {ca}/{cb} not both present)")
            continue
        rows.append({"pdb": pdb, "chains": f"{ca}/{cb}", "pisa": pisa_area, **res})
        print(f"  {pdb} {ca}/{cb}: PISA={pisa_area:7.1f}  ours={res['our_area']:7.1f}")

    pisa = np.array([r["pisa"] for r in rows])
    ours = np.array([r["our_area"] for r in rows])
    pct = 100 * np.abs(ours - pisa) / pisa
    r = pearsonr(ours, pisa)[0]
    fit = linregress(pisa, ours)
    print(f"\nn={len(rows)}  Pearson r={r:.3f}  Spearman={spearmanr(ours, pisa)[0]:.3f}"
          f"  slope={fit.slope:.3f}  median |error|={np.median(pct):.1f}%")

    fig, ax = plt.subplots(figsize=(5.4, 5.2))
    lim = [0, max(pisa.max(), ours.max()) * 1.08]
    ax.fill_between(lim, [0.95 * x for x in lim], [1.05 * x for x in lim],
                    color="#0072B2", alpha=0.08, zorder=0, label="±5%")
    ax.plot(lim, lim, "--", color="#8a8a8a", lw=1.2, zorder=1, label="y = x")
    ax.scatter(pisa, ours, s=55, color="#0072B2", edgecolor="white", lw=0.8, zorder=3,
               label=f"complexes (n={len(rows)})")
    ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect("equal")
    ax.set_xlabel("PISA interface area  (Å²)"); ax.set_ylabel("this tool, interface area  (Å²)")
    ax.set_title("Interface area agrees with PISA", fontsize=13, fontweight="bold", pad=10)
    ax.text(0.04, 0.88, f"Pearson r = {r:.3f}\nslope = {fit.slope:.2f}\nmedian |error| = {np.median(pct):.1f}%",
            transform=ax.transAxes, fontsize=9.5, va="top",
            bbox=dict(boxstyle="round,pad=0.4", fc="#f2f6fa", ec="#8a8a8a", lw=0.6))
    ax.legend(fontsize=8.5, loc="lower right", frameon=False)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.grid(True, color="#ececec", lw=0.6, zorder=0)
    fig.tight_layout()
    fig.savefig(HERE / "pisa_correctness.png", dpi=160, facecolor="white")
    print(f"wrote {HERE/'pisa_correctness.png'}")


if __name__ == "__main__":
    main()
