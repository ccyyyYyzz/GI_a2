"""
Shared figure style for PAPER 2 (Optica) — "Engineering What an Optical Sensor
Cannot See" (ROUND63 / R45 §C, architecture FROZEN by TLSG_PASS).

Publication-grade defaults (Optica Publishing Group style):
  - vector PDF + 600-dpi PNG output
  - colorblind-safe Okabe-Ito palette
  - pdf.fonttype = 42 (TrueType, editable/embeddable)
  - fonts >= 7 pt at final column width

This is a verbatim reuse of paper_prl/figures/prl_style.py (paper 1 frozen)
with Optica column widths substituted. It contains NO data: every figure script
reads ONLY committed frozen artifacts from paper2_optica/figures/_frozen_sources/.
Nothing here touches paper_prl/ or any sealed artifact.
"""
import matplotlib as mpl

# ---- Okabe-Ito colorblind-safe palette (identical to paper 1) ----
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

# Semantic role colors for Paper 2 (the three-panel spine)
C_WALL    = OI["blue"]      # exact symmetry wall (score-null orbit)
C_BREAK   = OI["vermil"]    # controlled breaker (finite window / NA clip)
C_CAP     = OI["green"]     # blindness-capacity (certified-blind code family)
C_DPSS    = OI["orange"]    # DPSS concentration route (noncomposable)
C_KNOWN   = OI["blue"]      # Ledger I  — known direction   (T*lambda ~ 1)
C_BLIND   = OI["green"]     # Ledger II — blind existence   (T*lambda ~ sqrt p)
C_EST     = OI["vermil"]    # Ledger III— full estimation   (T*lambda ~ p)
C_KT6     = OI["purple"]    # KT6 operating point (below blind line)
C_JET     = OI["orange"]    # projected jet (restored m=2 channel)
C_GREY    = OI["grey"]

# Optica column widths (inches). Optica single column ~3.5", double ~7.0".
COL1 = 3.5
COL2 = 7.0


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
