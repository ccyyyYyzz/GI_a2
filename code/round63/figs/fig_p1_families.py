"""Paper-1 small Results figure -- Study-2 family outcomes (R20 M12).

This holds the six-family empirical acceleration outcomes that were REMOVED
from the paper-1 mechanism figure (`fig_mechanism_p1`, panel c) per R20 M12.
The points are shown as EMPIRICAL OUTCOMES only: measured bucket contrast C_u
versus the number of accelerating instances (out of four). No region labels,
no transition band, no Gamma=1 boundary, and no cosmetic horizontal offsets --
each family sits at its true measured contrast.

Numbers are transcribed verbatim from paper/main_oe.tex (L505-508): family
(measured bucket contrast C_u, accelerating instances / 4):
  maze 0.53 4/4; glyph 0.41 4/4; spokes 0.38 4/4; chirp 0.28 3/4;
  contour 0.27 1/4; microtexture 0.17 0/4.

Run (cwd = repo root D:\\GI_another):
    D:\\Anacondar\\anaconda3\\python.exe code\\round63\\figs\\fig_p1_families.py
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import style_r20 as S

REPO = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
OUT = os.path.join(REPO, "paper", "figs")

# (family, measured C_u, accelerating instances / 4)   -- main_oe.tex L505-508
FAMILIES = [("maze", 0.53, 4), ("glyph", 0.41, 4), ("spokes", 0.38, 4),
            ("chirp", 0.28, 3), ("contour", 0.27, 1), ("microtexture", 0.17, 0)]


def _blend(frac):
    """gray (0/4, no benefit) -> saturated blue (4/4) -- reduced palette."""
    g = np.array([0x6E, 0x6E, 0x6E]) / 255.0
    b = np.array([0x00, 0x72, 0xB2]) / 255.0
    c = (1 - frac) * g + frac * b
    return tuple(c)


def main():
    S.use_style()
    fig, ax = plt.subplots(figsize=(S.COL1_IN, 2.55))
    fig.subplots_adjust(left=0.155, right=0.965, top=0.86, bottom=0.185)

    # per-family label placement (dx, dy, ha) to avoid overlap of the two
    # near-coincident 4/4 families (spokes C_u=0.38, glyph C_u=0.41)
    LAB = {"maze": (0.0, 0.30, "center"), "glyph": (0.028, 0.30, "left"),
           "spokes": (-0.028, 0.30, "right"), "chirp": (0.0, 0.30, "center"),
           "contour": (0.0, -0.42, "center"),
           "microtexture": (0.0, 0.30, "center")}
    for name, cu, npass in FAMILIES:
        col = _blend(npass / 4.0)
        ax.scatter([cu], [npass], s=95, color=col, edgecolor="k",
                   linewidth=0.6, zorder=4)
        dx, dy, ha = LAB[name]
        va = "bottom" if dy > 0 else "top"
        ax.annotate(name, xy=(cu, npass), xytext=(cu + dx, npass + dy),
                    ha=ha, va=va, fontsize=6.8, color=S.INK)

    ax.set_xlim(0.12, 0.58)
    ax.set_ylim(-0.6, 4.7)
    ax.set_yticks([0, 1, 2, 3, 4])
    ax.set_xlabel(r"measured bucket contrast  $C_u$", fontsize=7.8)
    ax.set_ylabel("accelerating instances (of 4)", fontsize=7.6)
    ax.tick_params(labelsize=7.0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    S.panel_title(ax, "acceleration vs. contrast", size=8.0)

    S.save(fig, OUT, "fig_p1_families")
    return 0


if __name__ == "__main__":
    sys.exit(main())
