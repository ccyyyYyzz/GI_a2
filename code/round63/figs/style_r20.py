"""Shared figure style for the R20 presentation-audit redraw (issue #12).

Implements the R20 M10 / SF7 / SF8 figure rules used by every figure produced
for the R20 ruling:

  * reduced palette -- neutral GRAY = safe/baseline, one saturated colour
    (BLUE) = ridge / high-flux (the protagonist positive result), one secondary
    colour (ORANGE) = jitter / negative / uncertainty / cost / censoring;
  * clean sans-serif axes and labels (one font family);
  * `pdf.fonttype 42` embedded fonts;
  * a 1-3 word in-panel title / result label is allowed (helper `panel_title`);
  * one short causal callout per panel is allowed (helper `callout`).

Final absolute text sizes and single/double-column widths are only APPROXIMATE
here: the exact resize + font-match pass happens AFTER the Optica template port
(R20 M9).  Build widths use the R20 guess of 8.6 cm single / 17.8 cm double
column (`CM`, `COL1_IN`, `COL2_IN`).
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- reduced R20 palette (M10 / SF8) ------------------------------------- #
GRAY   = "#6E6E6E"   # safe / baseline / deterministic reference
GRAY_L = "#B4B4B4"   # light gray (bands, minor reference)
BLUE   = "#0072B2"   # ridge / high-flux / benefit  (the one saturated colour)
BLUE_L = "#7FC1E3"   # light tint of the ridge colour (intermediate series)
ORANGE = "#D55E00"   # jitter / negative / uncertainty / cost / censoring
INK    = "#222222"
BLACK  = "#000000"

# ---- approximate column widths (final sizing deferred to template port) --- #
CM = 1.0 / 2.54
COL1_IN = 8.6 * CM    # single column  ~3.39 in
COL2_IN = 17.8 * CM   # double column  ~7.01 in


def use_style():
    """Install the shared R20 rcParams (sans-serif, embedded fonts)."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "mathtext.fontset": "dejavusans",
        "font.size": 8.0,               # approximate; >=7.5-8 pt final (M9/SF7)
        "axes.linewidth": 0.7,
        "axes.labelsize": 8.0,
        "axes.titlesize": 8.5,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "legend.fontsize": 7.0,
        "lines.linewidth": 1.6,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "figure.dpi": 120,
        "savefig.dpi": 200,
    })


def panel_title(ax, text, x=0.5, y=1.0, ha="center", va="bottom",
                color=INK, size=8.2, weight="bold", pad=2.0):
    """A 1-3 word in-panel title / result label (R20 M10)."""
    ax.annotate(text, xy=(x, y), xycoords="axes fraction", ha=ha, va=va,
                fontsize=size, color=color, fontweight=weight,
                xytext=(0, pad), textcoords="offset points", clip_on=False)


def panel_letter(fig_or_ax, letter, x, y, size=10.0):
    """Consistent upper-left panel letter, embedded in the figure (SF14/M14)."""
    obj = fig_or_ax
    if hasattr(obj, "text") and not hasattr(obj, "transAxes"):
        obj.text(x, y, letter, fontsize=size, fontweight="bold", va="top")
    else:
        obj.text(x, y, letter, transform=obj.transAxes, fontsize=size,
                 fontweight="bold", va="top", ha="left")


def callout(ax, text, xy, xytext, color=ORANGE, size=7.0, ha="left",
            va="center", lw=1.1, rad=-0.2):
    """One short causal callout with a curved arrow (R20 M10)."""
    ax.annotate(text, xy=xy, xytext=xytext, fontsize=size, color=color,
                ha=ha, va=va, fontweight="bold",
                arrowprops=dict(arrowstyle="-|>", color=color, lw=lw,
                                connectionstyle="arc3,rad=%g" % rad))


def save(fig, out_dir, name, also_png=True):
    os.makedirs(out_dir, exist_ok=True)
    pdf = os.path.join(out_dir, name + ".pdf")
    fig.savefig(pdf)
    if also_png:
        fig.savefig(os.path.join(out_dir, name + ".png"), dpi=200)
    plt.close(fig)
    print("[r20] wrote %s (+ .png)" % pdf, flush=True)
    return pdf
