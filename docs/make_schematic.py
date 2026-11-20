"""Build docs/schematic.png: the real 7JZU render (docs/structure.png, made by
render_structure.py) on the left, and the numbers binderqc actually returns for
it on the right.

The values are not illustrative -- they are this structure's real row, so the
figure shows what you get rather than naming what you get.

Run: python docs/render_structure.py && python docs/make_schematic.py
"""
import os

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

TEAL, INK, SUB, RULE = "#14726a", "#14202a", "#5b6b73", "#dde3e7"
HERE = os.path.dirname(os.path.abspath(__file__))
# The chimenodes render is the hero. render_structure.py writes structure.png as
# the CI-reproducible fallback; swap the name here to use that instead.
STRUCTURE = "structure_chimenodes.png"

# 7JZU chain A (LCB1) against chain B (RBD), from score_structure. Regenerate with:
#   score_structure("tests/data/7JZU_LCB1_RBD.pdb", ["A"], ["B"])[0]
PANEL = [
    ("Interface", [("buried area", "1021 Å²"), ("contacts", "21 res"),
                   ("H-bonds / salt bridges", "17 / 2"), ("packing", "14.5")]),
    ("Pose", [("approach angle", "86°  (across)")]),
    ("Grippability", [("epitope planarity", "3.21 Å")]),
    ("Developability", [("SAP peak / A3D load", "0.33 / 0.90"),
                        ("charge patch (+)", "1.78"),
                        ("net charge · pI", "−8.0 · 4.2")]),
]


def _load_trimmed(path):
    """Composite onto white (the chimenodes render is genuinely transparent) and trim."""
    img = plt.imread(path)
    if img.ndim == 3 and img.shape[2] == 4:
        a = img[..., 3:4]
        img = img[..., :3] * a + 1.0 * (1 - a)
    rgb = img[..., :3] if img.ndim == 3 else np.stack([img] * 3, -1)
    nonwhite = np.any(rgb < 0.96, axis=2)
    ys, xs = np.where(nonwhite)
    if len(xs):
        pad = 6
        img = img[max(ys.min() - pad, 0):ys.max() + pad, max(xs.min() - pad, 0):xs.max() + pad]
    return img


fig = plt.figure(figsize=(13, 5.6), facecolor="white")

# ---- structure ----------------------------------------------------------
ax = fig.add_axes([0.0, 0.05, 0.43, 0.92])
ax.imshow(_load_trimmed(os.path.join(HERE, STRUCTURE)))
ax.axis("off")
ax.text(0.5, -0.012, "7JZU   ·   LCB1 minibinder (blue) + SARS-CoV-2 RBD (grey)",
        transform=ax.transAxes, ha="center", va="top", fontsize=9.5, color=SUB)

# ---- arrow through the package name -------------------------------------
mid = fig.add_axes([0.44, 0.07, 0.13, 0.88]); mid.axis("off")
mid.set_xlim(0, 1); mid.set_ylim(0, 1)
mid.add_patch(FancyBboxPatch((0.08, 0.455), 0.84, 0.10, boxstyle="round,pad=0.022,rounding_size=0.03",
                             fc=TEAL, ec="none", transform=mid.transAxes))
mid.text(0.5, 0.505, "binderqc", transform=mid.transAxes, ha="center", va="center",
         fontsize=15, color="white", fontweight="bold")
for x0, x1 in ((-0.26, 0.05), (0.95, 1.26)):
    mid.add_patch(FancyArrowPatch((x0, 0.505), (x1, 0.505), transform=mid.transAxes,
                                  arrowstyle="-|>", mutation_scale=20, lw=2.0,
                                  color="#9aa6ad", clip_on=False, zorder=0))
mid.text(0.5, 0.40, "CPU only\nno GPU · no network", transform=mid.transAxes,
         ha="center", va="top", fontsize=8.5, color=SUB, linespacing=1.5)

# ---- the actual row it returns ------------------------------------------
right = fig.add_axes([0.575, 0.07, 0.42, 0.88]); right.axis("off")
right.set_xlim(0, 1); right.set_ylim(0, 1)

y = 0.975
for group, items in PANEL:
    right.text(0.0, y, group, fontsize=12, color=INK, fontweight="bold", va="top")
    y -= 0.062
    for label, value in items:
        right.text(0.02, y, label, fontsize=10, color=SUB, va="top")
        right.text(1.0, y, value, fontsize=10.5, color=INK, va="top", ha="right",
                   fontfamily="DejaVu Sans Mono")
        y -= 0.055
    y -= 0.012
    right.plot([0.0, 1.0], [y + 0.028, y + 0.028], color=RULE, lw=0.9,
               transform=right.transAxes, clip_on=False)
    y -= 0.030

right.text(0.0, y - 0.01, "+ tag site, epitope agreement, self-association, sequence liabilities",
           fontsize=9, color=SUB, va="top", style="italic")

fig.savefig(os.path.join(HERE, "schematic.png"), dpi=170,
            bbox_inches="tight", pad_inches=0.12, facecolor="white")
print("wrote docs/schematic.png")
