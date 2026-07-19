"""Optics-Express paper-1 MECHANISM / SCHEMATIC figure (Fig. 1).

Paper: paper/main_oe.tex --- "When high flux helps single-pixel imaging:
a contrast--dead-time operating map".

A single three-panel Fig.-1 that a reviewer reads in ~10 s:

  (a) the single-pixel imaging chain  source -> DMD (structured illumination,
      pattern inset) -> object (low- vs high-contrast pair) -> single-pixel
      photon-counting detector (SPD) -> renewal quasi-likelihood (RQL)
      reconstruction;
  (b) the dead-time mechanism:  a photon-arrival timeline with the blind
      interval tau censoring counts after every registered event (lost counts
      x-marked), and the resulting mean detected-vs-incident saturation
      nu_det/nu = rho/(1+rho) that bends away from the ideal at high load,
      with the count-information ridge rho*(nu) = (6 nu)^{1/3} - 2/3 marked;
  (c) the operating map:  bucket contrast C_u versus detector load rho, with
      the photon-limited onset Gamma = 1, and the six Study-2 structural
      families placed by their measured contrast and coloured by the frozen
      acquisition-speed gate outcome (maze 4/4 down to microtexture 0/4).

NUMBER PROVENANCE (every hard-coded value is transcribed from paper/main_oe.tex;
no value is invented, and any quantity not in the manuscript is omitted):
  * conventional load rho <= 0.05 .......................... main_oe.tex L98,167
  * high-flux operating point rho_bar = 0.6 ................ main_oe.tex L368,479
  * exposure T = nu tau, nu = T/tau ....................... main_oe.tex L51,158
  * arrival rate lambda_i = Phi a_i^T x + d ............... main_oe.tex L156-157
  * non-paralyzable per-slot detected fraction rho/(1+rho)  ridge kernel L87
                                    (code/round63/fig_a_ridge_map.py, J_CLT)
  * cube-root ridge rho*(nu) = (6 nu)^{1/3} - 2/3 ......... main_oe.tex Eq.(6),L265
  * contrast-to-count-noise index Gamma = C_u sqrt(nu rho/(1+rho)),
    onset of the photon-limited structural regime at Gamma=1  main_oe.tex L465-466
  * measured family bucket contrast C_u and Study-2 gate outcome
    (out of four instances): maze 0.53 4/4, glyph 0.41 4/4, spokes 0.38 4/4,
    chirp 0.28 3/4, contour 0.27 1/4, microtexture 0.17 0/4  main_oe.tex L505-508
  * dense multiplex-limited contrast C_u ~ 1/sqrt(n) ~ 0.016 (64^2),
    0.040 at the k=512 rung (32^2) ......................... main_oe.tex L437-439,476

This generator is self-contained: it depends only on numpy + matplotlib and the
constants above; it reads no campaign artifact and modifies no other file.

Run (cwd = repo root D:\\GI_another):
    python code\\round63\\figs\\fig_mechanism_p1.py
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))              # code/round63/figs
REPO = os.path.dirname(os.path.dirname(os.path.dirname(HERE))) # repo root
OUT_DIR = os.path.join(REPO, "paper", "figs")
OUT_PDF = os.path.join(OUT_DIR, "fig_mechanism_p1.pdf")
OUT_PNG = os.path.join(OUT_DIR, "fig_mechanism_p1.png")

# ----------------------------------------------------------------------------
# Okabe-Ito / Wong colourblind-safe palette (matched to the paper's other figs)
# ----------------------------------------------------------------------------
BLACK  = "#000000"
BLUE   = "#0072B2"
ORANGE = "#D55E00"
GREEN  = "#009E73"
PINK   = "#CC79A7"
SKY    = "#56B4E9"
YELLOW = "#F0E442"
INK    = "#222222"
BOXFC  = "#EAF2F8"   # light box fill
BOXEC  = "#34495E"   # box edge
GRIDGY = "#BBBBBB"
FAILGY = "#D9D9D9"   # "non-accelerating" pie wedge

# ----------------------------------------------------------------------------
# Manuscript constants (see NUMBER PROVENANCE header)
# ----------------------------------------------------------------------------
RHO_CONV   = 0.05      # conventional photon-counting operating point
RHO_HIFLUX = 0.6       # pre-registered high-flux operating point (rho_bar)
NU_OP      = 2000.0    # high-flux dwell for the Gamma / operating-map panel

# (family, measured bucket contrast C_u, accelerating instances out of 4)
FAMILIES = [
    ("maze",         0.53, 4),
    ("glyph",        0.41, 4),
    ("spokes",       0.38, 4),
    ("chirp",        0.28, 3),
    ("contour",      0.27, 1),
    ("microtexture", 0.17, 0),
]
CU_DENSE = 0.040       # dense k=512 rung (multiplex-limited anchor)


def rho_star(nu):
    """Information-optimal load rho*(nu) = (6 nu)^{1/3} - 2/3  (main_oe.tex Eq.6)."""
    return (6.0 * nu) ** (1.0 / 3.0) - 2.0 / 3.0


def detected_fraction(rho):
    """Mean non-paralyzable detected counts per dead-time slot nu_det/nu = rho/(1+rho)."""
    return rho / (1.0 + rho)


def gamma_boundary_cu(rho, nu=NU_OP):
    """Contrast on the Gamma = 1 photon-limited-onset curve: C_u = 1/sqrt(nu rho/(1+rho))."""
    return 1.0 / np.sqrt(nu * rho / (1.0 + rho))


# ----------------------------------------------------------------------------
# small deterministic icon arrays for panel (a)
# ----------------------------------------------------------------------------
def dmd_pattern(n=8, seed=3):
    """Binary Bernoulli(1/2) illumination pattern a_i (DMD micromirror mask)."""
    return np.random.RandomState(seed).randint(0, 2, size=(n, n)).astype(float)


def object_lowcontrast(n=40):
    """Smooth, weakly modulated object: the multiplex-/contrast-limited case."""
    yy, xx = np.mgrid[0:n, 0:n] / (n - 1.0)
    r = np.sqrt((xx - 0.5) ** 2 + (yy - 0.5) ** 2)
    img = 0.5 + 0.06 * np.exp(-(r / 0.35) ** 2)          # range ~0.50-0.56
    return img


def object_highcontrast(n=64):
    """Bold radial-spoke object: the contrast-adequate case (evokes the spokes family)."""
    yy, xx = np.mgrid[0:n, 0:n] / (n - 1.0)
    th = np.arctan2(yy - 0.5, xx - 0.5)
    rr = np.sqrt((xx - 0.5) ** 2 + (yy - 0.5) ** 2)
    img = 0.5 + 0.5 * np.sign(np.sin(6.0 * th))
    img[rr > 0.48] = 0.0                                  # circular aperture
    return img


# ----------------------------------------------------------------------------
# panel (a) : the single-pixel imaging chain
# ----------------------------------------------------------------------------
def rounded_box(ax, cx, cy, w, h, fc="none", ec=BOXEC, lw=0.8, z=2):
    box = FancyBboxPatch((cx - w / 2, cy - h / 2), w, h,
                         boxstyle="round,pad=0.0,rounding_size=2.2",
                         linewidth=lw, edgecolor=ec, facecolor=fc, zorder=z,
                         mutation_aspect=1.0)
    ax.add_patch(box)
    return box


def light_arrow(ax, x0, x1, y, color=ORANGE, lw=2.0, z=3, ls="-"):
    ar = FancyArrowPatch((x0, y), (x1, y), arrowstyle="-|>", mutation_scale=11,
                         lw=lw, color=color, zorder=z, linestyle=ls,
                         shrinkA=0, shrinkB=0)
    ax.add_patch(ar)


AX_A_W_IN = 0.965 * 7.16      # panel-(a) axes width  (inches) -- see main()
AX_A_H_IN = 0.320 * 5.00      # panel-(a) axes height (inches)


def icon_inset(ax, cx, cy, side_in, ec=BOXEC, lw=0.8):
    """Square image inset at data centre (cx,cy) on the 0..100 schematic axis.

    Sizing corrects for the panel's wide/short pixel aspect so the icon renders
    as a true visual square regardless of the axis limits."""
    wf = side_in / AX_A_W_IN
    hf = side_in / AX_A_H_IN
    iax = ax.inset_axes([cx / 100.0 - wf / 2, cy / 100.0 - hf / 2, wf, hf])
    iax.set_xticks([])
    iax.set_yticks([])
    for s in iax.spines.values():
        s.set_edgecolor(ec)
        s.set_linewidth(lw)
    return iax


def panel_a(ax):
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.set_autoscale_on(False)
    ax.axis("off")

    cy = 58.0
    bw, bh = 15.0, 30.0
    centers = {"src": 9.5, "dmd": 30.0, "obj": 51.0, "spd": 72.0, "rec": 92.5}

    # --- light path (drawn first, under the boxes) ---
    light_arrow(ax, centers["src"] + bw / 2, centers["dmd"] - bw / 2, cy)
    light_arrow(ax, centers["dmd"] + bw / 2, centers["obj"] - bw / 2 - 1, cy)
    light_arrow(ax, centers["obj"] + bw / 2 + 1, centers["spd"] - bw / 2, cy)
    light_arrow(ax, centers["spd"] + bw / 2, centers["rec"] - bw / 2, cy, color=BLUE)

    # --- source ---
    rounded_box(ax, centers["src"], cy, bw, bh)
    ax.plot([centers["src"] - 4.2, centers["src"] + 3.2], [cy + 2, cy + 2],
            color=ORANGE, lw=2.4, solid_capstyle="round", zorder=4)
    ax.add_patch(FancyArrowPatch((centers["src"] + 3.2, cy + 2),
                                 (centers["src"] + 5.4, cy + 2),
                                 arrowstyle="-|>", mutation_scale=8, lw=2.0,
                                 color=ORANGE, zorder=4, shrinkA=0, shrinkB=0))
    ax.text(centers["src"], cy - bh / 2 - 6.5, "light\nsource", ha="center",
            va="top", fontsize=8.0, color=INK)

    # --- DMD with binary pattern inset ---
    rounded_box(ax, centers["dmd"], cy, bw, bh)
    iax = icon_inset(ax, centers["dmd"], cy, 0.36)
    iax.imshow(dmd_pattern(), cmap="gray", vmin=0, vmax=1, aspect="auto",
               interpolation="nearest")
    ax.text(centers["dmd"], cy - bh / 2 - 6.5, "DMD", ha="center", va="top",
            fontsize=8.0, color=INK, fontweight="bold")
    ax.text(centers["dmd"], cy - bh / 2 - 12.0,
            r"pattern $a_i$", ha="center", va="top", fontsize=7.0, color=INK)

    # --- object: low- vs high-contrast pair ---
    ow = 16.5
    rounded_box(ax, centers["obj"], cy, ow, bh)
    iax = icon_inset(ax, centers["obj"] - 3.9, cy + 3.2, 0.30)
    iax.imshow(object_lowcontrast(), cmap="gray", vmin=0, vmax=1, aspect="auto",
               interpolation="bilinear")
    iax = icon_inset(ax, centers["obj"] + 3.9, cy + 3.2, 0.30)
    iax.imshow(object_highcontrast(), cmap="gray", vmin=0, vmax=1, aspect="auto",
               interpolation="nearest")
    # labels clearly BELOW the insets, on the white box background
    ax.text(centers["obj"] - 3.9, cy - 10.6, r"low $C_u$", ha="center",
            va="center", fontsize=7.0, color=INK, zorder=5)
    ax.text(centers["obj"] + 3.9, cy - 10.6, r"high $C_u$", ha="center",
            va="center", fontsize=7.0, color=INK, zorder=5)
    ax.text(centers["obj"], cy - bh / 2 - 6.5, r"object $x$", ha="center",
            va="top", fontsize=8.0, color=INK)

    # --- single-pixel detector (bucket + detector) ---
    rounded_box(ax, centers["spd"], cy, bw, bh)
    # collecting lens funnel -> detector dot
    lx = centers["spd"] - 5.0
    ax.add_patch(plt.Polygon([[lx, cy + 6.5], [lx, cy - 6.5],
                              [lx + 8.5, cy + 2.0], [lx + 8.5, cy - 2.0]],
                             closed=True, facecolor=SKY, edgecolor=BOXEC,
                             lw=0.8, alpha=0.85, zorder=3))
    ax.add_patch(plt.Circle((lx + 10.2, cy), 1.7, facecolor=BLACK,
                            edgecolor=BOXEC, lw=0.8, zorder=4))
    ax.text(centers["spd"], cy - bh / 2 - 6.5, "SPD", ha="center", va="top",
            fontsize=8.0, color=INK, fontweight="bold")
    ax.text(centers["spd"], cy - bh / 2 - 12.0, "photon-\ncounting",
            ha="center", va="top", fontsize=7.0, color=INK)

    # --- RQL reconstruction ---
    rounded_box(ax, centers["rec"], cy, bw, bh, ec=GREEN, lw=0.9)
    iax = icon_inset(ax, centers["rec"], cy + 1.0, 0.34, ec=GREEN)
    iax.imshow(object_highcontrast(48), cmap="gray", vmin=0, vmax=1,
               aspect="auto", interpolation="nearest")
    ax.text(centers["rec"], cy - bh / 2 - 6.5, "RQL", ha="center", va="top",
            fontsize=8.0, color=INK, fontweight="bold")
    ax.text(centers["rec"], cy - bh / 2 - 12.0, r"recon. $\hat{x}$",
            ha="center", va="top", fontsize=7.0, color=INK)

    # --- counts-path label (data-path name; other prose lives in the caption)
    ax.text((centers["spd"] + centers["rec"]) / 2, cy + 3.4,
            r"bucket counts $N_i$", ha="center", va="bottom", fontsize=6.6,
            color=BLUE)


# ----------------------------------------------------------------------------
# panel (b1) : photon-arrival timeline with dead-time censoring
# ----------------------------------------------------------------------------
def panel_b_pulses(ax):
    tau = 1.0
    # deterministic illustrative arrival times (in units of tau)
    arrivals = [0.35, 1.05, 1.30, 2.55, 3.10, 3.35, 4.35, 5.65, 5.95, 7.05, 8.25]
    tmax = 9.0
    # non-paralyzable processing: a registered count blinds the detector for tau
    reg, lost, blocks = [], [], []
    free_at = 0.0
    for t in arrivals:
        if t >= free_at:
            reg.append(t)
            blocks.append((t, t + tau))
            free_at = t + tau
        else:
            lost.append(t)

    ax.set_xlim(-0.2, tmax + 0.2)
    ax.set_ylim(0, 1.62)
    for b0, b1 in blocks:                                   # dead-time blocks
        ax.add_patch(Rectangle((b0, 0), min(b1, tmax) - b0, 1.0,
                               facecolor=ORANGE, alpha=0.16, edgecolor="none",
                               zorder=1))
    ax.axhline(0, color=INK, lw=1.0, zorder=2)
    for t in reg:                                           # registered counts
        ax.plot([t, t], [0, 1.0], color=GREEN, lw=1.8, zorder=4,
                solid_capstyle="round")
        ax.plot([t], [1.0], marker="v", color=GREEN, ms=4.5, zorder=4)
    for t in lost:                                          # censored arrivals
        ax.plot([t], [0.5], marker="x", color=ORANGE, ms=5.5, mew=1.6, zorder=5)
    # tau brace on the first dead-time block (explicit, large dead-time glyph)
    b0, b1 = blocks[0]
    ax.annotate("", xy=(b0, 1.12), xytext=(b1, 1.12),
                arrowprops=dict(arrowstyle="<->", color=INK, lw=1.0))
    ax.text((b0 + b1) / 2, 1.20, r"$\tau$", ha="center", va="bottom",
            fontsize=13.0, color=ORANGE, clip_on=False)
    # (registered/lost/blind symbols and the exposure are glossed in the
    #  caption; no in-figure legend or title)
    ax.set_yticks([])
    ax.set_xticks([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_visible(False)


# ----------------------------------------------------------------------------
# panel (b2) : mean detected-vs-incident saturation + information ridge
# ----------------------------------------------------------------------------
def panel_b_saturation(ax):
    rho = np.geomspace(0.02, 30.0, 400)
    ax.plot(rho, rho, ls="--", lw=1.1, color=GRIDGY, zorder=2,
            label="ideal (no dead time)")
    ax.plot(rho, detected_fraction(rho), lw=2.0, color=BLUE, zorder=4,
            label=r"detected  $\rho/(1{+}\rho)$")

    # deployment band rho <= 1
    ax.axvspan(0.02, 1.0, color=BLUE, alpha=0.06, zorder=0)
    # information ridge rho*(nu) band, nu in [100, 2000]
    r_lo, r_hi = rho_star(100.0), rho_star(2000.0)
    ax.axvspan(r_lo, r_hi, color=ORANGE, alpha=0.12, zorder=0)
    ax.text(np.sqrt(r_lo * r_hi), 0.16, r"ridge $\rho^{*}(\nu)$", ha="center",
            va="bottom", fontsize=6.2, color=ORANGE)

    # operating points (terse value labels; roles named in the caption)
    ax.plot([RHO_CONV], [detected_fraction(RHO_CONV)], marker="o", color=GREEN,
            ms=5.5, mec="k", mew=0.5, zorder=6)
    ax.plot([RHO_HIFLUX], [detected_fraction(RHO_HIFLUX)], marker="o", color=PINK,
            ms=6.0, mec="k", mew=0.5, zorder=6)
    ax.annotate(r"$\rho{=}0.05$", xy=(RHO_CONV, detected_fraction(RHO_CONV)),
                xytext=(0.024, 0.55), fontsize=6.2, color=GREEN, ha="left", va="center",
                arrowprops=dict(arrowstyle="-|>", color=GREEN, lw=0.9))
    ax.annotate(r"$\bar{\rho}{=}0.6$", xy=(RHO_HIFLUX, detected_fraction(RHO_HIFLUX)),
                xytext=(0.9, 0.22), fontsize=6.2, color=PINK, ha="left", va="center",
                arrowprops=dict(arrowstyle="-|>", color=PINK, lw=0.9))

    ax.set_xscale("log")
    ax.set_xlim(0.02, 30.0)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel(r"incident load  $\rho=\lambda\tau$", fontsize=7.6, labelpad=1.5)
    ax.set_ylabel(r"detected / slot  $\nu_{\det}/\nu$", fontsize=7.6, labelpad=2)
    ax.tick_params(labelsize=6.6)
    ax.legend(loc="upper left", fontsize=6.0, frameon=False, handlelength=1.6,
              borderpad=0.2, labelspacing=0.25)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


# ----------------------------------------------------------------------------
# panel (c) : the contrast--dead-time operating map
# ----------------------------------------------------------------------------
def _pie_marker(ax, x, y, frac, base=150.0, ec="k", lw=0.7, z=6):
    """Circular pass-fraction pie: green = accelerating instances, grey = rest."""
    ax.scatter([x], [y], s=base, marker="o", facecolor=FAILGY, edgecolor=ec,
               linewidths=lw, zorder=z)
    if frac > 0:
        ang = np.linspace(np.pi / 2, np.pi / 2 - 2 * np.pi * frac,
                          max(3, int(80 * frac) + 2))
        verts = np.vstack([[0, 0], np.column_stack([np.cos(ang), np.sin(ang)])])
        scale = np.abs(verts).max()
        ax.scatter([x], [y], s=base * scale ** 2, marker=verts, facecolor=GREEN,
                   edgecolor="none", zorder=z + 1)
    ax.scatter([x], [y], s=base, marker="o", facecolor="none", edgecolor=ec,
               linewidths=lw, zorder=z + 2)


def panel_c(ax):
    rho = np.geomspace(0.01, 1.2, 200)
    ymin, ymax = 0.03, 1.75

    # photon-limited onset Gamma = 1 (at nu = 2000)
    ax.plot(rho, gamma_boundary_cu(rho), color=BLACK, lw=1.3, zorder=4)
    # multiplex-limited region: below the Gamma=1 curve
    ax.fill_between(rho, ymin, gamma_boundary_cu(rho), color=GRIDGY, alpha=0.30,
                    zorder=1)

    # "high-flux operation helps" region (above the descriptive transition)
    ax.fill_between(rho, 0.285, ymax, color=GREEN, alpha=0.09, zorder=0)
    # descriptive contrast-controlled transition band (contour 0.27 <-> chirp 0.28)
    ax.axhspan(0.265, 0.285, color=YELLOW, alpha=0.5, zorder=1)

    # conventional low-flux zone
    ax.axvspan(0.01, 0.05, color=SKY, alpha=0.13, zorder=0)
    # high-flux operating column
    ax.axvline(RHO_HIFLUX, color=PINK, lw=1.1, ls="--", zorder=2)

    # marker x-positions: every family sits at the operating load rho_bar = 0.6;
    # the horizontal offsets are purely cosmetic, applied only so the six pies
    # never overlap (chirp/contour are 0.01 apart in C_u) — housekeeping note
    # in the caption. y = true measured bucket contrast C_u, never offset.
    XPOS = {"maze": 0.60, "glyph": 0.74, "spokes": 0.49,
            "chirp": 0.53, "contour": 0.87, "microtexture": 0.60}
    # tidy left-hand label column, evenly spaced in log-C_u; a thin leader runs
    # from each label to the CENTRE of its pie (under the marker) so the
    # label-marker association is unambiguous.
    order = sorted(FAMILIES, key=lambda t: t[1])          # ascending C_u
    label_y = np.geomspace(0.165, 0.575, len(order))
    lab_x = 0.42
    for (name, cu, npass), ly in zip(order, label_y):
        mx = XPOS[name]
        ax.plot([lab_x * 1.05, mx], [ly, cu], color="#888888", lw=0.6,
                zorder=3)
        _pie_marker(ax, mx, cu, npass / 4.0, base=125.0)
        ax.text(lab_x, ly, r"%s  %d/4" % (name, npass), ha="right", va="center",
                fontsize=6.7, color=INK)

    # dense multiplex-limited anchor (Study 1)
    ax.scatter([RHO_HIFLUX], [CU_DENSE], marker="D", s=30, facecolor="white",
               edgecolor=BLACK, linewidths=0.9, zorder=6)
    ax.plot([lab_x * 1.05, RHO_HIFLUX], [CU_DENSE, CU_DENSE],
            color="#888888", lw=0.6, zorder=3)
    ax.text(lab_x, CU_DENSE, "dense (Study 1)", ha="right", va="center",
            fontsize=6.2, color=INK)

    # minimal region / boundary names (full descriptions live in the caption)
    ax.text(0.011, 1.52, "high-flux helps", fontsize=7.2, color=GREEN,
            ha="left", va="center", fontweight="bold")
    ax.text(0.011, 0.037, "multiplex-limited", fontsize=6.4,
            color="#555555", ha="left", va="center")
    ax.text(0.155, gamma_boundary_cu(0.155) * 1.16, r"$\Gamma=1$", fontsize=7.4,
            color=BLACK, ha="left", va="bottom", rotation=-24)
    ax.text(0.011, 0.275, "transition band", fontsize=5.8,
            color="#8a6d00", ha="left", va="center")
    ax.text(RHO_HIFLUX, ymax * 0.96, r"$\bar{\rho}=0.6$", fontsize=6.4,
            color=PINK, ha="center", va="top")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.01, 1.2)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel(r"detector load  $\rho=\lambda\tau$", fontsize=7.8, labelpad=1.5)
    ax.set_ylabel(r"bucket contrast  $C_u$", fontsize=7.8, labelpad=1)
    ax.tick_params(labelsize=6.6)
    ax.set_xticks([0.01, 0.05, 0.1, 0.6, 1.0])
    ax.set_xticklabels(["0.01", "0.05", "0.1", "0.6", "1"])
    ax.set_yticks([0.05, 0.1, 0.2, 0.5, 1.0])
    ax.set_yticklabels(["0.05", "0.1", "0.2", "0.5", "1"])
    ax.minorticks_off()
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


# ----------------------------------------------------------------------------
# figure assembly
# ----------------------------------------------------------------------------
def main():
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 8,
        "axes.linewidth": 0.6,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "mathtext.fontset": "dejavuserif",
    })

    fig = plt.figure(figsize=(7.16, 5.0))
    ax_a = fig.add_axes([0.020, 0.660, 0.965, 0.320])
    ax_b1 = fig.add_axes([0.085, 0.400, 0.375, 0.150])
    ax_b2 = fig.add_axes([0.085, 0.090, 0.375, 0.235])
    ax_c = fig.add_axes([0.605, 0.090, 0.370, 0.470])

    panel_a(ax_a)
    panel_b_pulses(ax_b1)
    panel_b_saturation(ax_b2)
    panel_c(ax_c)

    # panel letters
    fig.text(0.018, 0.965, "(a)", fontsize=10, fontweight="bold", va="top")
    fig.text(0.018, 0.560, "(b)", fontsize=10, fontweight="bold", va="top")
    fig.text(0.545, 0.560, "(c)", fontsize=10, fontweight="bold", va="top")

    os.makedirs(OUT_DIR, exist_ok=True)
    fig.savefig(OUT_PDF)
    fig.savefig(OUT_PNG, dpi=220)
    plt.close(fig)
    print("[fig_mechanism_p1] wrote %s" % OUT_PDF, flush=True)
    print("[fig_mechanism_p1] wrote %s (220 dpi)" % OUT_PNG, flush=True)
    print("[fig_mechanism_p1] rho*(100)=%.2f  rho*(2000)=%.2f" %
          (rho_star(100.0), rho_star(2000.0)), flush=True)
    print("[fig_mechanism_p1] Gamma=1 contrast @rho=0.6,nu=2000: C_u=%.4f" %
          gamma_boundary_cu(0.6), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
