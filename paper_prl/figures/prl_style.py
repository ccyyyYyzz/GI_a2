"""
Shared PRL figure style for the ROUND63 Rank-Jet Separation Letter.

Publication-grade defaults per R43 §6:
  - vector PDF + PNG output
  - colorblind-safe Okabe-Ito palette
  - fonts >= 7 pt at final PRL width
  - pdf.fonttype = 42 (TrueType, editable/embeddable)

This module contains NO data. All figure scripts read ONLY committed frozen
artifacts from paper_prl/figures/_frozen_sources/ (see each script header).
"""
import matplotlib as mpl

# ---- Okabe-Ito colorblind-safe palette ----
OI = {
    "black":   "#000000",
    "orange":  "#E69F00",
    "skyblue": "#56B4E9",
    "green":   "#009E73",
    "yellow":  "#F0E442",
    "blue":    "#0072B2",
    "vermil":  "#D55E00",
    "purple":  "#CC79A7",
    "grey":    "#999999",
}

# Semantic role colors for this Letter
C_GENERIC = OI["blue"]      # generic direction, m=1 (first-order observable)
C_ORTHO   = OI["vermil"]    # first-order-orthogonal, m=2 (curvature-rescued)
C_MIXED   = OI["green"]     # mixed direction (crossover)
C_ANCHOR  = OI["orange"]    # amplitude anchor / restoration
C_WALL    = OI["grey"]      # mean wall / blind

# PRL column widths (inches)
COL1 = 86.0 / 25.4     # single column  ~3.386 in
COL2 = 180.0 / 25.4    # double column  ~7.087 in


def apply_style():
    mpl.rcParams.update({
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "mathtext.fontset": "dejavusans",
        "font.size": 7.5,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "axes.linewidth": 0.7,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "xtick.minor.width": 0.5,
        "ytick.minor.width": 0.5,
        "lines.linewidth": 1.3,
        "lines.markersize": 4.0,
        "axes.labelpad": 2.0,
        "figure.dpi": 150,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
    })


def save(fig, stem):
    """Write both vector PDF and 600-dpi PNG next to the scripts."""
    fig.savefig(stem + ".pdf")
    fig.savefig(stem + ".png", dpi=600)
    print("wrote", stem + ".pdf", "and", stem + ".png")
