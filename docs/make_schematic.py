"""Build docs/schematic.png: the real 7JZU render (docs/structure.png, made by
render_structure.py) on the left, the binderqc QC-output flow on the right.
Run: python docs/render_structure.py && python docs/make_schematic.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch

TEAL, GREY, INK, SUB = "#2f6fb0", "#9aa0a6", "#1a2b34", "#5b6b73"
HERE = os.path.dirname(os.path.abspath(__file__))


def _load_trimmed(path):
    img = plt.imread(path)
    rgb = img[..., :3] if img.ndim == 3 else np.stack([img] * 3, -1)
    nonwhite = np.any(rgb < 0.96, axis=2)
    ys, xs = np.where(nonwhite)
    if len(xs):
        pad = 12
        img = img[max(ys.min() - pad, 0):ys.max() + pad, max(xs.min() - pad, 0):xs.max() + pad]
    return img


fig = plt.figure(figsize=(12, 5), facecolor="white")

ax = fig.add_axes([0.01, 0.10, 0.40, 0.84])
ax.imshow(_load_trimmed(os.path.join(HERE, "structure.png")))
ax.axis("off")

fx = fig.add_axes([0, 0, 1, 1]); fx.set_xlim(0, 1); fx.set_ylim(0, 1)
fx.axis("off"); fx.patch.set_alpha(0)
fx.text(0.21, 0.055, "7JZU  ·  LCB1 minibinder (blue) + SARS-CoV-2 RBD (grey)",
        ha="center", fontsize=9.5, color=SUB)

box = FancyBboxPatch((0.45, 0.42), 0.13, 0.16, boxstyle="round,pad=0.012,rounding_size=0.02",
                     linewidth=0, facecolor="#0f7d6b")
fx.add_patch(box)
fx.text(0.515, 0.50, "binderqc", ha="center", va="center", color="white",
        fontsize=15, fontweight="bold")
for x0, x1 in ((0.415, 0.445), (0.585, 0.635)):
    fx.add_patch(FancyArrowPatch((x0, 0.50), (x1, 0.50), arrowstyle="-|>",
                                 mutation_scale=16, color="#b9c2c6", lw=1.6))

rows = [
    ("1", "Interface", "buried area · H-bonds · salt bridges · packing"),
    ("2", "Pose", "approach angle (end-on / across)"),
    ("3", "Grippability", "epitope planarity · aromatic anchors · glyco-occlusion"),
    ("4", "Tag site", "recommended N / C terminus · Cys-SG"),
    ("5", "Developability", "SAP / Aggrescan3D · liabilities · GRAVY · pI · MW"),
]
for (num, title, sub), y in zip(rows, np.linspace(0.86, 0.14, len(rows))):
    fx.add_patch(Circle((0.665, y), 0.022, facecolor="#0f7d6b", linewidth=0))
    fx.text(0.665, y, num, ha="center", va="center", color="white", fontsize=10, fontweight="bold")
    fx.text(0.70, y + 0.03, title, ha="left", va="center", color=INK, fontsize=13, fontweight="bold")
    fx.text(0.70, y - 0.03, sub, ha="left", va="center", color=SUB, fontsize=9.5)

fig.savefig(os.path.join(HERE, "schematic.png"), dpi=150)
plt.close(fig)
print("wrote", os.path.join(HERE, "schematic.png"))
