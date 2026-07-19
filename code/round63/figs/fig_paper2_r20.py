"""Paper-2 figures for the R20 presentation ruling (GitHub issue #12).

Builds, into paper2/figs/, the reader-facing figure set mandated by R20:

  fig_mechanism        (M3)  hero Fig.1: (a) simulation-model SPI chain,
                             (b) deterministic vs jittered dead time,
                             (c) HEADLINE J(rho) jitter cap at nu=2000.
  fig_jitter_validation(M4)  (a) peak load vs c_v (log-log, both MC runs +
                             preregistered predictions + fitted slope -0.658),
                             (b) retained information at the deterministic ridge
                             vs jitter (engineering warning).
  fig_speed_results    (M5+M15) positive-results figure: reconstruction rows +
                             fixed-dwell dQ strip + elapsed-time PAVA curves +
                             24-scene S_gate strip + two-resource annotation.
  fig_audit_supp       (M6)  SUPPLEMENT: guard-collapse path (old actiii_a) +
                             certificate status distribution (old actiii_c).
  actiii_b             (M6)  compact single-panel dose-only headroom (restyle).

Sources (all pre-existing; NO new Monte Carlo is run):
  * jitter grids -> results/round63_study2/jitter_sfi_v2_nu2000.log,
    jitter_zoom.log, docs/R14_PREREGISTERED_PREDICTIONS.md  (values transcribed
    as frozen constants below);
  * long-window law J_inf = rho/((1+rho)(1+c^2 rho^2))  (analytic);
  * c=0 finite-window exact count law -> code/round63/physics.fisher_exact;
  * corrected campaign verdicts / per-image tables ->
    results/round63_m1/CORRECTION_2026-07-19/M1_VERDICTS_SPEC_CORRECTED_R19.json;
  * elapsed-time quality curves -> results/round63_m1/shards/*.csv
    (optical_time_s = elapsed T_opt; PSNR_rad quality);
  * reconstructions -> results/round63_m1/ACTIII_D_RECON.npz;
  * certificate cells / guard trajectory -> results/round63_m1/ (as actiii).

Run (cwd = repo root D:\\GI_another):
    D:\\Anacondar\\anaconda3\\python.exe code\\round63\\figs\\fig_paper2_r20.py
Individual figs: pass names, e.g. `... fig_paper2_r20.py mech jitter speed`.
"""
import os
import sys
import glob
import csv
import json
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.transforms as mtransforms
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
R63 = os.path.dirname(HERE)
REPO = os.path.dirname(os.path.dirname(R63))
sys.path.insert(0, HERE)
import style_r20 as S

OUT = os.path.join(REPO, "paper2", "figs")
M1 = os.path.join(REPO, "results", "round63_m1")
SHARDS = os.path.join(M1, "shards")
CORR = os.path.join(M1, "CORRECTION_2026-07-19")
D_NPZ = os.path.join(M1, "ACTIII_D_RECON.npz")

# ========================================================================== #
#  FROZEN JITTER MEASUREMENTS (transcribed; no new Monte Carlo)               #
#  results/round63_study2/jitter_sfi_v2_nu2000.log  (nu=2000, n_mc=100k grid) #
# ========================================================================== #
MC_RHO = np.array([1.200, 1.389, 1.607, 1.860, 2.153, 2.491, 2.883, 3.337,
                   3.862, 4.470, 5.173, 5.986, 6.928, 8.018, 9.280, 10.739,
                   12.429, 14.384, 16.647, 19.266, 22.297, 25.805, 29.864,
                   34.563, 40.000])
MC_J = {  # c_v -> J(rho; nu=2000) measured
    0.00: np.array([0.5513, 0.5885, 0.6242, 0.6555, 0.6848, 0.7176, 0.7424,
                    0.7721, 0.7904, 0.8230, 0.8468, 0.8549, 0.8781, 0.8867,
                    0.9066, 0.9126, 0.9236, 0.9284, 0.9394, 0.9354, 0.9452,
                    0.9385, 0.9379, 0.9295, 0.9198]),
    0.02: np.array([0.5478, 0.5879, 0.6256, 0.6546, 0.6869, 0.7155, 0.7493,
                    0.7749, 0.7946, 0.8164, 0.8310, 0.8445, 0.8595, 0.8662,
                    0.8700, 0.8731, 0.8712, 0.8571, 0.8464, 0.8194, 0.7751,
                    0.7408, 0.6969, 0.6382, 0.5709]),
    0.05: np.array([0.5504, 0.5805, 0.6188, 0.6502, 0.6830, 0.7097, 0.7297,
                    0.7538, 0.7709, 0.7833, 0.7848, 0.7819, 0.7799, 0.7612,
                    0.7345, 0.7099, 0.6646, 0.6129, 0.5559, 0.4908, 0.4262,
                    0.3570, 0.2953, 0.2404, 0.1922]),
}
DET_RIDGE = 22.3          # deterministic optimum rho* (measured peak; (6nu)^1/3-2/3=22.23)
DET_JPEAK = 0.9452        # measured c=0 peak
CV05_OPT = 5.7            # c_v=0.05 optimum (long-window law + Colab 400k = 5.700)
J_AT_RIDGE_CV05 = 0.4262  # measured J(22.3; c_v=0.05)
DEPLOY_RHO = 22.25        # deployed RIDGE-SCAT32 load (paper2 table)

# M4 peak-vs-c_v data (frozen)
CV_GRID = np.array([0.02, 0.05, 0.10, 0.20, 0.30])
PRED_R14 = np.array([10.40, 5.69, 3.50, 2.30, 1.63])          # preregistered exact optimum
CV_ZOOM = np.array([0.02, 0.05, 0.10, 0.20])                  # local 150k continuous peaks
RHO_ZOOM = np.array([9.645, 5.154, 3.752, 2.106])
CV_COLAB = np.array([0.02, 0.05, 0.10, 0.20])                 # Colab 400k continuous peaks
RHO_COLAB = np.array([9.610, 5.700, 3.399, 2.051])
FIT_SLOPE = -0.658                                            # pooled 8-pt log-log fit
# retained information at deterministic ridge rho=22.297 (nu=2000 grid, last-but-few col)
RET_CV = np.array([0.00, 0.02, 0.05, 0.10, 0.20, 0.30])
RET_J = np.array([0.9452, 0.7751, 0.4262, 0.1608, 0.0479, 0.0245])


def Jlw(rho, c):
    """Long-window jitter law J_inf = rho/((1+rho)(1+c^2 rho^2))."""
    return rho / ((1.0 + rho) * (1.0 + c ** 2 * rho ** 2))


def Jfin(rho, nu):
    """c=0 finite-window exact count-information efficiency (physics.fisher_exact)."""
    sys.path.insert(0, R63)
    from physics import fisher_exact
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return np.array([fisher_exact(float(r), float(nu), 1.0) * r ** 2 / nu
                         for r in np.atleast_1d(rho)])


# ========================================================================== #
#  M3 -- hero mechanism figure                                                #
# ========================================================================== #
def _scat32_mask(seed=3232, occ=0.12, n=32):
    rng = np.random.default_rng(seed)
    return (rng.random((n, n)) < occ).astype(float)


def _chain_panel(ax):
    """(a) minimal computational SPI chain -- one SCAT32 mask inset only."""
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 40)
    ax.axis("off")
    y = 22.0
    bh, bw = 15.0, 12.5
    centers = {"src": 8.0, "dmd": 29.0, "obj": 50.0, "spd": 71.0, "rec": 92.0}

    def rbox(cx, ec=S.INK, lw=0.8):
        ax.add_patch(FancyBboxPatch((cx - bw / 2, y - bh / 2), bw, bh,
                     boxstyle="round,pad=0.3,rounding_size=1.4",
                     fc="none", ec=ec, lw=lw, zorder=2))

    def beam(x0, x1, color=S.ORANGE, lw=2.0):
        ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle="-|>",
                     mutation_scale=11, color=color, lw=lw, zorder=3,
                     shrinkA=0, shrinkB=0))

    def below(cx, txt):
        ax.text(cx, y - bh / 2 - 2.0, txt, ha="center", va="top",
                fontsize=7.4, color=S.INK)

    # source + power multiplier (plain arrow-up glyph, NO dial halo/ticks)
    rbox(centers["src"])
    ax.plot([centers["src"] - 3.2, centers["src"] + 2.2], [y + 2.4, y + 2.4],
            color=S.ORANGE, lw=2.2, solid_capstyle="round", zorder=4)
    ax.add_patch(FancyArrowPatch((centers["src"], y + 2.4),
                 (centers["src"], y + 5.6), arrowstyle="-|>", mutation_scale=8,
                 color=S.BLUE, lw=1.8, zorder=4, shrinkA=0, shrinkB=0))
    ax.text(centers["src"], y + 6.4, r"$\Phi$", ha="center", va="bottom",
            fontsize=9.5, color=S.BLUE, fontweight="bold")
    below(centers["src"], "source")
    ax.text(centers["src"], y - bh / 2 - 6.2, "global\nmultiplier", ha="center",
            va="top", fontsize=6.4, color=S.BLUE)

    beam(centers["src"] + bw / 2, centers["dmd"] - bw / 2)
    # DMD with ONE SCAT32 mask inset
    rbox(centers["dmd"])
    iax = ax.inset_axes([(centers["dmd"] - 4.4) / 100.0, (y - 4.4) / 40.0,
                         8.8 / 100.0, 8.8 / 40.0])
    iax.imshow(_scat32_mask(), cmap="gray", vmin=0, vmax=1,
               interpolation="nearest")
    iax.set_xticks([]); iax.set_yticks([])
    for sp in iax.spines.values():
        sp.set_edgecolor(S.BLUE); sp.set_linewidth(0.7)
    below(centers["dmd"], r"DMD, SCAT32 $a_i$")

    beam(centers["dmd"] + bw / 2, centers["obj"] - bw / 2)
    rbox(centers["obj"])
    ax.text(centers["obj"], y, r"$x$", ha="center", va="center", fontsize=11,
            color=S.INK)
    below(centers["obj"], "object")

    beam(centers["obj"] + bw / 2, centers["spd"] - bw / 2)
    rbox(centers["spd"])
    lx = centers["spd"] - 4.4
    ax.fill([lx, lx + 5.6, lx], [y + 3.2, y, y - 3.2], color=S.GRAY_L,
            zorder=3)
    ax.add_patch(plt.Circle((lx + 6.6, y), 1.3, color=S.BLACK, zorder=4))
    below(centers["spd"], "bucket detector")

    beam(centers["spd"] + bw / 2, centers["rec"] - bw / 2, color=S.BLUE)
    ax.text((centers["spd"] + centers["rec"]) / 2, y + 2.2,
            r"counts $N_i$", ha="center", va="bottom", fontsize=6.6,
            color=S.BLUE)
    rbox(centers["rec"], ec=S.BLUE, lw=1.0)
    ax.text(centers["rec"], y, "RQL", ha="center", va="center", fontsize=8.0,
            color=S.BLUE, fontweight="bold")
    below(centers["rec"], r"recon. $\hat{x}$")
    # single "simulation model" marker so it is not read as a bench experiment
    S.panel_title(ax, "simulation model", x=0.5, y=1.0, color=S.GRAY,
                  weight="normal", size=7.6)


def _timeline_panel(ax):
    """(b) deterministic vs jittered dead time -- two stacked short timelines."""
    ax.set_xlim(0, 10.5)
    ax.set_ylim(-0.4, 2.9)
    ax.axis("off")
    arrivals = [0.55, 1.15, 1.55, 2.75, 3.25, 3.55, 4.55, 5.75, 6.05, 7.15,
                8.35, 9.1]
    tau = 0.95

    def render(y0, holds, tag, tag_color):
        ax.plot([0.2, 10.3], [y0, y0], color=S.INK, lw=1.0)
        ready = 0.0
        first = None
        for i, t in enumerate(arrivals):
            w = holds[i]
            if t >= ready:
                ax.add_patch(Rectangle((t, y0), w, 0.95, fc=S.ORANGE,
                             ec="none", alpha=0.16, zorder=1))
                ax.plot([t, t], [y0, y0 + 0.95], color=S.BLUE, lw=1.6,
                        zorder=3, solid_capstyle="round")
                ax.plot([t], [y0 + 0.95], marker="v", color=S.BLUE, ms=3.6,
                        zorder=3)
                if first is None:
                    first = (t, w)
                ready = t + w
            else:
                ax.plot([t], [y0 + 0.45], marker="x", color=S.ORANGE, ms=5.0,
                        mew=1.4, zorder=4)
        ax.text(0.1, y0 + 0.45, tag, ha="right", va="center", fontsize=6.6,
                color=tag_color, rotation=90)
        return first

    # top: fixed hold tau (deterministic)
    f1 = render(1.6, [tau] * len(arrivals), "fixed  " + r"$\tau$", S.GRAY)
    t0, w0 = f1
    ax.annotate("", xy=(t0, 1.6 + 1.18), xytext=(t0 + w0, 1.6 + 1.18),
                arrowprops=dict(arrowstyle="<->", color=S.GRAY, lw=0.9))
    ax.text(t0 + w0 / 2, 1.6 + 1.24, r"$\tau$", ha="center", va="bottom",
            fontsize=9.5, color=S.GRAY)
    # bottom: jittered holds B_i (random)
    rng = np.random.default_rng(7)
    jit = tau * (1.0 + 0.55 * rng.standard_normal(len(arrivals)))
    jit = np.clip(jit, 0.45, 1.9)
    f2 = render(0.0, list(jit), "jittered  " + r"$B_i$", S.ORANGE)
    t0b, w0b = f2
    ax.annotate("", xy=(t0b, 0.0 + 1.18), xytext=(t0b + w0b, 0.0 + 1.18),
                arrowprops=dict(arrowstyle="<->", color=S.ORANGE, lw=0.9))
    ax.text(t0b + w0b / 2, 0.0 + 1.24, r"$B_i$", ha="center", va="bottom",
            fontsize=9.0, color=S.ORANGE)
    S.panel_title(ax, "hidden hold variance", x=0.5, y=1.0, color=S.ORANGE,
                  weight="normal", size=7.4)


def _headline_panel(ax):
    """(c) HEADLINE: J(rho) jitter cap at nu=2000."""
    rho_s = np.geomspace(0.6, 40.0, 240)
    # smooth model curves
    ax.plot(rho_s, Jfin(rho_s, 2000.0), color=S.GRAY, lw=1.8, zorder=4,
            label=r"$c_v=0$")
    ax.plot(rho_s, Jlw(rho_s, 0.02), color=S.ORANGE, lw=1.5, ls="--",
            alpha=0.75, zorder=4, label=r"$c_v=0.02$")
    ax.plot(rho_s, Jlw(rho_s, 0.05), color=S.ORANGE, lw=1.9, zorder=5,
            label=r"$c_v=0.05$")
    # measured Monte-Carlo points (nu=2000)
    for c, col, mk in ((0.00, S.GRAY, "o"), (0.02, S.ORANGE, "^"),
                       (0.05, S.ORANGE, "s")):
        ax.plot(MC_RHO, MC_J[c], mk, color=col, ms=3.0, mfc="white",
                mew=0.8, zorder=6, ls="none")
    # deterministic ridge + jittered optimum markers (result labels)
    ax.axvline(DET_RIDGE, color=S.GRAY, ls=":", lw=1.0, zorder=2)
    ax.plot([DET_RIDGE], [DET_JPEAK], "o", color=S.GRAY, ms=6, mec="k",
            mew=0.5, zorder=7)
    ax.plot([CV05_OPT], [Jlw(CV05_OPT, 0.05)], "s", color=S.ORANGE, ms=6,
            mec="k", mew=0.5, zorder=7)
    ax.annotate(r"det. $\bar\rho^{*}{\approx}22$", xy=(DET_RIDGE, DET_JPEAK),
                xytext=(11.5, 1.005), fontsize=6.6, color=S.GRAY, ha="left")
    ax.annotate(r"$\bar\rho^{*}{\approx}5.7$",
                xy=(CV05_OPT, Jlw(CV05_OPT, 0.05)), xytext=(4.2, 0.585),
                fontsize=6.8, color=S.ORANGE, ha="center",
                arrowprops=dict(arrowstyle="-|>", color=S.ORANGE, lw=0.8))
    # deployed ridge operating point (protagonist, blue)
    ax.axvline(DEPLOY_RHO, color=S.BLUE, ls="-", lw=0.9, alpha=0.30, zorder=1)
    # ~55% loss -- the one causal callout
    ax.annotate("", xy=(DET_RIDGE, J_AT_RIDGE_CV05),
                xytext=(DET_RIDGE, DET_JPEAK),
                arrowprops=dict(arrowstyle="<->", color=S.ORANGE, lw=1.3))
    ax.annotate(r"$\sim$55% lost", xy=(DET_RIDGE, 0.685), xytext=(8.0, 0.47),
                fontsize=7.2, color=S.ORANGE, fontweight="bold", ha="center",
                arrowprops=dict(arrowstyle="-|>", color=S.ORANGE, lw=1.0,
                                connectionstyle="arc3,rad=0.25"))
    ax.text(0.97, 0.05, r"$\nu=2000$", transform=ax.transAxes, fontsize=6.6,
            color=S.GRAY, va="bottom", ha="right")
    ax.set_xscale("log")
    ax.set_xlim(0.6, 40)
    ax.set_ylim(0, 1.08)
    ax.set_xlabel(r"detector load  $\bar\rho=\tau\,\mathbb{E}[\lambda]$")
    ax.set_ylabel(r"count-information efficiency  $J(\bar\rho)$")
    ax.set_xticks([1, 2, 5, 10, 22.3, 40])
    ax.set_xticklabels(["1", "2", "5", "10", "22", "40"])
    ax.legend(loc="upper left", frameon=False, fontsize=6.6, ncol=1,
              handlelength=1.6, borderpad=0.2, labelspacing=0.3,
              bbox_to_anchor=(-0.01, 0.86))
    S.panel_title(ax, "jitter caps the optimal load", x=0.5, y=1.0,
                  color=S.INK, size=8.0)


def fig_mechanism():
    S.use_style()
    fig = plt.figure(figsize=(S.COL2_IN, 4.5))
    gs = gridspec.GridSpec(2, 2, height_ratios=[1.0, 1.28],
                           width_ratios=[1.0, 1.05],
                           left=0.075, right=0.975, top=0.90, bottom=0.115,
                           hspace=0.42, wspace=0.24)
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])
    _chain_panel(ax_a)
    _timeline_panel(ax_b)
    _headline_panel(ax_c)
    fig.text(0.012, 0.965, "(a)", fontsize=10, fontweight="bold", va="top")
    fig.text(0.012, 0.545, "(b)", fontsize=10, fontweight="bold", va="top")
    fig.text(0.505, 0.545, "(c)", fontsize=10, fontweight="bold", va="top")
    S.save(fig, OUT, "fig_mechanism")


# ========================================================================== #
#  M4 -- jitter validation                                                    #
# ========================================================================== #
def fig_jitter_validation():
    S.use_style()
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(S.COL2_IN, 3.0))
    fig.subplots_adjust(left=0.085, right=0.975, top=0.87, bottom=0.17,
                        wspace=0.32)

    # ---- (a) peak load vs c_v (log-log) ---- #
    cc = np.geomspace(0.015, 0.35, 100)
    # -2/3 reference law anchored at the c=0 finite-window ridge cap constant
    ref = (2.0) ** (-1.0 / 3.0) * cc ** (-2.0 / 3.0)
    axa.plot(cc, ref, color=S.GRAY, ls="--", lw=1.3, zorder=2,
             label=r"$-2/3$ law")
    # fitted slope -0.658 (pooled), anchored at c_v=0.05 geometric mean peak
    anchor_c, anchor_r = 0.05, np.sqrt(RHO_ZOOM[1] * RHO_COLAB[1])
    fit = anchor_r * (cc / anchor_c) ** FIT_SLOPE
    axa.plot(cc, fit, color=S.BLUE, ls="-", lw=1.5, zorder=3,
             label=r"fit slope $-0.658$")
    # preregistered exact-optimum prediction (theorem / matched asymptotic)
    axa.plot(CV_GRID, PRED_R14, color=S.INK, ls="none", marker="D", ms=5,
             mfc="none", mew=1.0, zorder=5, label="predicted (R14)")
    # two independent Monte-Carlo runs
    axa.plot(CV_ZOOM, RHO_ZOOM, color=S.ORANGE, ls="none", marker="o", ms=5,
             zorder=6, label="MC run 1 (150k)")
    axa.plot(CV_COLAB, RHO_COLAB, color=S.BLUE, ls="none", marker="^", ms=5,
             zorder=6, label="MC run 2 (400k)")
    axa.set_xscale("log"); axa.set_yscale("log")
    axa.set_xlabel(r"hold-time jitter  $c_v$")
    axa.set_ylabel(r"information-optimal load  $\bar\rho^{*}$")
    axa.set_xlim(0.015, 0.35)
    axa.set_ylim(1.3, 16)
    axa.set_xticks([0.02, 0.05, 0.1, 0.2, 0.3])
    axa.set_xticklabels(["0.02", "0.05", "0.1", "0.2", "0.3"])
    axa.set_yticks([2, 3, 5, 10])
    axa.set_yticklabels(["2", "3", "5", "10"])
    axa.legend(loc="upper right", frameon=False, fontsize=6.3,
               handlelength=1.6, borderpad=0.2, labelspacing=0.3)
    S.panel_title(axa, r"$\bar\rho^{*}\propto c_v^{-2/3}$", color=S.INK,
                  size=8.0)

    # ---- (b) retained information at the deterministic ridge ---- #
    ret = 100.0 * RET_J / RET_J[0]
    axb.plot(RET_CV, ret, color=S.ORANGE, lw=1.8, marker="o", ms=4.5,
             zorder=4)
    axb.axhline(100.0, color=S.GRAY, ls="--", lw=1.0, zorder=1)
    axb.axhline(50.0, color=S.GRAY_L, ls=":", lw=1.0, zorder=1)
    axb.fill_between(RET_CV, 0, ret, color=S.ORANGE, alpha=0.08, zorder=0)
    axb.set_xlabel(r"hold-time jitter  $c_v$")
    axb.set_ylabel(r"info. retained at det. ridge $\bar\rho{\approx}22$  (%)")
    axb.set_xlim(-0.005, 0.31)
    axb.set_ylim(0, 108)
    S.callout(axb, "operating the\nzero-jitter ridge\nwastes >50%",
              xy=(0.05, 100.0 * RET_J[2] / RET_J[0]), xytext=(0.11, 66.0),
              color=S.ORANGE, size=6.6, rad=-0.25)
    S.panel_title(axb, "engineering warning", color=S.ORANGE, weight="normal",
                  size=7.6)

    fig.text(0.012, 0.965, "(a)", fontsize=10, fontweight="bold", va="top")
    fig.text(0.525, 0.965, "(b)", fontsize=10, fontweight="bold", va="top")
    S.save(fig, OUT, "fig_jitter_validation")


# ========================================================================== #
#  Speed-figure data helpers (M5)                                             #
# ========================================================================== #
def _pava(y, w=None):
    """Equal-weight block-collapse PAVA (non-decreasing). Unit: [3,2,1]->[2,2,2].

    Standard pool-adjacent-violators with an explicit stack of blocks, each a
    [mean, weight, size] triple; adjacent violators are merged (block collapse),
    not weight-duplicated (R19 D3)."""
    y = np.asarray(y, float)
    n = y.size
    w = np.ones(n) if w is None else np.asarray(w, float)
    stack = []  # each: [mean, weight, size]
    for k in range(n):
        cur = [y[k], w[k], 1]
        while stack and stack[-1][0] >= cur[0]:
            pm, pw, ps = stack.pop()
            tot = pw + cur[1]
            cur = [(pm * pw + cur[0] * cur[1]) / tot, tot, ps + cur[2]]
        stack.append(cur)
    res = []
    for m, _ww, sz in stack:
        res.extend([m] * sz)
    return np.array(res)


def _speed_curves(img, arm):
    """(T_opt, PSNR_rad) seed-mean per dwell for (arm, img), sorted by time."""
    acc = {}
    for p in sorted(glob.glob(os.path.join(SHARDS, "*.csv"))):
        if os.path.basename(p).startswith("M1_CERT"):
            continue
        with open(p, encoding="utf-8") as f:
            for r in csv.DictReader(f):
                if r["image"] != img:
                    continue
                sid = r["shard_id"]
                if sid[len("M1_"):sid.rindex("_")] != arm:
                    continue
                try:
                    t = float(r["optical_time_s"]); q = float(r["PSNR_rad"])
                except (ValueError, KeyError):
                    continue
                acc.setdefault(round(t, 12), []).append(q)
    ts = np.array(sorted(acc))
    qs = np.array([np.mean(acc[t]) for t in ts])
    return ts, qs


def _cross_logtime(t, q_pava, target):
    """First elapsed time where the monotone PAVA curve reaches target
    (linear interp in log-time). None if never reached."""
    if q_pava[-1] < target:
        return None
    if q_pava[0] >= target:
        return t[0]
    for i in range(1, len(t)):
        if q_pava[i] >= target:
            lo, hi = q_pava[i - 1], q_pava[i]
            f = (target - lo) / (hi - lo) if hi > lo else 0.0
            return float(np.exp(np.log(t[i - 1]) +
                                f * (np.log(t[i]) - np.log(t[i - 1]))))
    return None


# corrected per-image speed table (frozen, R19) -- image -> (S_gate, status, Q90)
def _load_corrected():
    j = json.load(open(os.path.join(
        CORR, "M1_VERDICTS_SPEC_CORRECTED_R19.json")))
    return j


# ========================================================================== #
#  M5 + M15 -- positive-results speed figure                                  #
# ========================================================================== #
FAM_ORDER = ["chirp", "contour", "glyph", "maze", "microtexture", "spokes"]
# per-image dQ (5-seed mean, frozen R19 table)
DQ = {
    "m1_chirp_0": 2.361180, "m1_chirp_1": 1.270420, "m1_chirp_2": 1.005020,
    "m1_chirp_3": 1.538840, "m1_contour_0": -3.612180, "m1_contour_1": -5.794460,
    "m1_contour_2": -4.738300, "m1_contour_3": 4.151300, "m1_glyph_0": 4.734300,
    "m1_glyph_1": 6.005160, "m1_glyph_2": 2.430660, "m1_glyph_3": 2.212580,
    "m1_maze_0": 6.329740, "m1_maze_1": 6.436420, "m1_maze_2": 5.583480,
    "m1_maze_3": 6.408240, "m1_microtexture_0": 1.529520,
    "m1_microtexture_1": 1.252680, "m1_microtexture_2": 1.525840,
    "m1_microtexture_3": 1.448760, "m1_spokes_0": -0.765240,
    "m1_spokes_1": -5.084720, "m1_spokes_2": 4.530260, "m1_spokes_3": 2.195000,
}
# per-image S_gate + status + Q90 target (frozen R19 table)
SG = {
    "m1_chirp_0": (17.800203, "NORMAL", 6.914794),
    "m1_chirp_1": (20.459721, "NORMAL", 7.389890),
    "m1_chirp_2": (18.431468, "NORMAL", 7.261996),
    "m1_chirp_3": (19.261646, "NORMAL", 7.404890),
    "m1_contour_0": (2.940635, "NORMAL", 9.099688),
    "m1_contour_1": (0.0, "FAST_RIGHT_CENSORED", 9.584150),
    "m1_contour_2": (0.0, "FAST_RIGHT_CENSORED", 10.442378),
    "m1_contour_3": (18.784826, "NORMAL", 15.186056),
    "m1_glyph_0": (18.992440, "NORMAL", 17.970960),
    "m1_glyph_1": (21.500155, "NORMAL", 18.454512),
    "m1_glyph_2": (18.312756, "NORMAL", 9.343498),
    "m1_glyph_3": (20.593897, "NORMAL", 9.230494),
    "m1_maze_0": (21.440543, "NORMAL", 13.789082),
    "m1_maze_1": (26.814231, "NORMAL", 18.569958),
    "m1_maze_2": (26.844183, "NORMAL", 15.161556),
    "m1_maze_3": (29.012398, "NORMAL", 15.054804),
    "m1_microtexture_0": (18.142114, "NORMAL", 14.435742),
    "m1_microtexture_1": (1.0, "SAFE_RANGE_UNINFORMATIVE", None),
    "m1_microtexture_2": (18.014590, "NORMAL", 14.874426),
    "m1_microtexture_3": (15.355503, "NORMAL", 15.629622),
    "m1_spokes_0": (46.150824, "NORMAL", 7.989012),
    "m1_spokes_1": (20.384410, "NORMAL", 9.191172),
    "m1_spokes_2": (21.208912, "NORMAL", 11.207580),
    "m1_spokes_3": (22.034111, "NORMAL", 8.462386),
}
S_MEDIAN, S_LB = 19.127043, 18.328492
DQ_MEDIAN = 1.866920
FAM_COL = {"chirp": S.BLUE, "contour": S.ORANGE, "glyph": S.BLUE,
           "maze": S.BLUE, "microtexture": S.GRAY, "spokes": S.BLUE}
# three disclosed speed scenes
SPEED_SCENES = [("m1_chirp_3", "median"), ("m1_contour_1", "fast-censored"),
                ("m1_microtexture_1", "uninformative")]
# two reconstruction rows (median + worst per taste-call 3); best optional
RECON_ROWS = [("median", "m1_chirp_3"), ("worst", "m1_contour_1")]


def _fam_of(img):
    return img.replace("m1_", "").rsplit("_", 1)[0]


def fig_speed_results():
    S.use_style()
    z = np.load(D_NPZ, allow_pickle=True) if os.path.exists(D_NPZ) else None
    fig = plt.figure(figsize=(S.COL2_IN, 6.1))
    outer = gridspec.GridSpec(3, 3, height_ratios=[1.0, 0.94, 0.26],
                              width_ratios=[1.12, 1.0, 1.0],
                              left=0.055, right=0.965, top=0.925, bottom=0.055,
                              wspace=0.34, hspace=0.52)

    # ---- LEFT (rows 0-1): reconstruction rows (GT | SCAT32-060 | RIDGE) ---- #
    gl = gridspec.GridSpecFromSubplotSpec(len(RECON_ROWS), 3, subplot_spec=
                                          outer[0:2, 0], wspace=0.06,
                                          hspace=0.30)
    col_titles = ["truth", r"safe $\bar\rho{=}0.6$", r"ridge $\bar\rho^{*}$"]
    for ri, (role, img) in enumerate(RECON_ROWS):
        if z is None:
            break
        tag = {"median": "median", "worst": "min"}[role]
        gt = z["%s__gt" % tag]
        vmax = float(np.percentile(gt, 99.7)) or 1.0
        cells = [("gt", gt, None),
                 ("060", z["%s__SCAT32-060__recon" % tag],
                  float(z["%s__SCAT32-060__psnr_rad" % tag])),
                 ("ridge", z["%s__RIDGE-SCAT32__recon" % tag],
                  float(z["%s__RIDGE-SCAT32__psnr_rad" % tag]))]
        for ci, (kind, arr, psnr) in enumerate(cells):
            ax = fig.add_subplot(gl[ri, ci])
            ax.imshow(np.clip(arr, 0, vmax), cmap="gray", vmin=0, vmax=vmax,
                      interpolation="nearest")
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_linewidth(0.6)
            if ri == 0:
                ax.set_title(col_titles[ci], fontsize=7.0, pad=2)
            if psnr is not None:
                ax.text(0.05, 0.95, "%.1f" % psnr, transform=ax.transAxes,
                        fontsize=6.6, color="white", va="top", ha="left",
                        bbox=dict(boxstyle="round,pad=0.12", fc="black",
                                  ec="none", alpha=0.55))
            if ci == 0:
                dqi = DQ[img]
                lab = "%s\n%s  %+.1f dB" % (
                    img.replace("m1_", ""),
                    {"median": "median", "worst": "worst"}[role], dqi)
                ax.set_ylabel(lab, fontsize=6.6,
                              color=S.BLUE if dqi >= 0 else S.ORANGE)
    fig.text(0.075, 0.945, "(a)  reconstructions", fontsize=8.0,
             fontweight="bold", va="top")

    # ---- TOP RIGHT (row 0, cols 1-2): fixed-dwell dQ strip (24 scenes) ---- #
    axq = fig.add_subplot(outer[0, 1:3])
    rng = np.random.default_rng(63)
    for fi, fam in enumerate(FAM_ORDER):
        ys = np.array([DQ["m1_%s_%d" % (fam, j)] for j in range(4)])
        xs = fi + (rng.random(4) - 0.5) * 0.42
        cols = [S.BLUE if v >= 0 else S.ORANGE for v in ys]  # positive/negative
        axq.scatter(xs, ys, s=22, color=cols, edgecolor="white",
                    linewidth=0.4, zorder=3)
    axq.axhline(DQ_MEDIAN, color=S.BLUE, lw=1.4, zorder=2)
    axq.axhline(1.0, color=S.GRAY, ls="--", lw=1.0, zorder=2)
    axq.axhline(0.0, color=S.GRAY_L, ls=":", lw=1.0, zorder=1)
    tr = mtransforms.blended_transform_factory(axq.transAxes, axq.transData)
    axq.text(0.02, DQ_MEDIAN, "median +1.87 dB", transform=tr, fontsize=6.6,
             va="bottom", ha="left", color=S.BLUE,
             bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none",
                       alpha=0.8))
    axq.text(0.98, 1.0, "+1 dB bar", transform=tr, fontsize=6.2, va="bottom",
             ha="right", color=S.GRAY)
    axq.set_xticks(range(6))
    axq.set_xticklabels(FAM_ORDER, rotation=32, ha="right", fontsize=6.6)
    axq.set_ylabel(r"fixed-dwell $\Delta$Q (dB)", fontsize=7.4)
    axq.set_xlim(-0.6, 5.6)
    S.panel_title(axq, "fixed-dwell quality: 19/24 positive", size=7.6)
    fig.text(0.40, 0.94, "(b)", fontsize=8.0, fontweight="bold", va="top")

    # ---- (c) PAVA quality-vs-elapsed curves for the 3 disclosed scenes ---- #
    axc = fig.add_subplot(outer[1, 1])
    scene_col = {"m1_chirp_3": S.BLUE, "m1_contour_1": S.ORANGE,
                 "m1_microtexture_1": S.GRAY}
    scene_lab = {"m1_chirp_3": "median", "m1_contour_1": "censored",
                 "m1_microtexture_1": "flat safe"}
    for img, role in SPEED_SCENES:
        col = scene_col[img]
        ts_s, qs_s = _speed_curves(img, "SCAT32-SAFE")
        ts_r, qs_r = _speed_curves(img, "RIDGE-SCAT32")
        pv_s = _pava(qs_s); pv_r = _pava(qs_r)
        axc.plot(ts_s, pv_s, color=col, lw=1.2, ls="--", alpha=0.85, zorder=3)
        axc.plot(ts_r, pv_r, color=col, lw=1.8, ls="-", zorder=4)
        axc.plot(ts_r, qs_r, "o", color=col, ms=2.4, mfc=col, mew=0, zorder=3)
        sgate, status, q90 = SG[img]
        # end-of-curve scene label
        axc.annotate(scene_lab[img], xy=(ts_r[-1], pv_r[-1]),
                     xytext=(2, 0), textcoords="offset points", fontsize=6.0,
                     color=col, va="center", ha="left")
        if q90 is not None:
            axc.axhline(q90, color=col, lw=0.5, ls=":", alpha=0.45, zorder=1)
            tr_c = _cross_logtime(ts_r, pv_r, q90)
            ts_c = _cross_logtime(ts_s, pv_s, q90)
            if tr_c:
                axc.plot([tr_c], [q90], "o", color=col, ms=4.5, mec="k",
                         mew=0.4, zorder=6)
            if ts_c:
                axc.plot([ts_c], [q90], "D", color=col, ms=4.5, mec="k",
                         mew=0.4, zorder=6)
    axc.set_xscale("log")
    axc.set_xlim(1.8e-4, 0.4)
    axc.set_xlabel(r"elapsed optical time  $T_{\rm opt}$ (s)", fontsize=7.2)
    axc.set_ylabel(r"PSNR$_{\rm rad}$ (dB)", fontsize=7.2)
    # style legend: solid=ridge, dashed=safe, o=ridge Q90, D=safe Q90
    handles = [Line2D([0], [0], color=S.INK, lw=1.8, ls="-", label="ridge"),
               Line2D([0], [0], color=S.INK, lw=1.2, ls="--", label="safe"),
               Line2D([0], [0], color=S.INK, marker="o", ls="none", ms=4,
                      label="ridge Q90"),
               Line2D([0], [0], color=S.INK, marker="D", ls="none", ms=4,
                      label="safe Q90")]
    axc.legend(handles=handles, loc="lower right", frameon=False, fontsize=5.6,
               handlelength=1.5, labelspacing=0.2, borderpad=0.2, ncol=2,
               columnspacing=0.8)
    S.panel_title(axc, "ridge reaches quality sooner", size=7.4)
    fig.text(0.395, 0.575, "(c)", fontsize=8.0, fontweight="bold", va="top")

    # ---- (d) 24-scene S_gate strip grouped by family ---- #
    axs = fig.add_subplot(outer[1, 2])
    rng2 = np.random.default_rng(14)
    for fi, fam in enumerate(FAM_ORDER):
        for j in range(4):
            img = "m1_%s_%d" % (fam, j)
            sgate, status, q90 = SG[img]
            xj = fi + (rng2.random() - 0.5) * 0.4
            if status == "NORMAL":
                axs.scatter([xj], [sgate], s=20, color=S.BLUE,
                            edgecolor="white", linewidth=0.3, zorder=3)
            elif status == "FAST_RIGHT_CENSORED":
                axs.scatter([xj], [1.4], s=26, marker="v", color=S.ORANGE,
                            edgecolor="k", linewidth=0.4, zorder=4)
            else:  # SAFE_RANGE_UNINFORMATIVE
                axs.scatter([xj], [1.0], s=26, marker="s", facecolor="none",
                            edgecolor=S.GRAY, linewidth=0.9, zorder=4)
    axs.axhspan(S_LB, S_MEDIAN, color=S.BLUE, alpha=0.12, zorder=0)
    axs.axhline(S_MEDIAN, color=S.BLUE, lw=1.4, zorder=2)
    axs.axhline(3.0, color=S.GRAY, ls="--", lw=1.0, zorder=2)
    tr2 = mtransforms.blended_transform_factory(axs.transAxes, axs.transData)
    axs.text(0.02, S_MEDIAN, "median 19.13", transform=tr2, fontsize=6.6,
             va="bottom", ha="left", color=S.BLUE,
             bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none",
                       alpha=0.85))
    axs.text(0.98, 3.0, "gate S=3", transform=tr2, fontsize=6.2, va="bottom",
             ha="right", color=S.GRAY)
    axs.text(0.02, S_LB, "LB 18.33", transform=tr2, fontsize=6.0, va="top",
             ha="left", color=S.BLUE)
    axs.set_yscale("log")
    axs.set_ylim(0.7, 60)
    axs.set_xticks(range(6))
    axs.set_xticklabels(FAM_ORDER, rotation=32, ha="right", fontsize=6.6)
    axs.set_ylabel(r"elapsed-time speedup  $S_{\rm gate}$", fontsize=7.2)
    axs.set_xlim(-0.6, 5.6)
    # censoring legend
    hl = [Line2D([0], [0], color=S.BLUE, marker="o", ls="none", ms=4,
                 label="normal"),
          Line2D([0], [0], color=S.ORANGE, marker="v", ls="none", ms=4,
                 label="fast-censored"),
          Line2D([0], [0], color=S.GRAY, marker="s", ls="none", mfc="none",
                 ms=4, label="uninformative")]
    axs.legend(handles=hl, loc="upper left", frameon=False, fontsize=5.6,
               handlelength=1.2, labelspacing=0.25, borderpad=0.2)
    S.panel_title(axs, "elapsed speedup: 21/24 > 1", size=7.4)
    fig.text(0.665, 0.575, "(d)", fontsize=8.0, fontweight="bold", va="top")

    # ---- two-resource annotation (M15): benefit vs cost in one field ---- #
    axr = fig.add_subplot(outer[2, :])
    axr.axis("off")
    axr.set_xlim(0, 1); axr.set_ylim(0, 1)
    axr.text(0.005, 0.70, "One global power control", fontsize=7.6,
             fontweight="bold", va="center", color=S.INK)
    axr.text(0.005, 0.22, "fixed-dwell power-for-time,",
             fontsize=6.0, color=S.GRAY, va="center", style="italic")
    axr.text(0.005, 0.02, "not equal-photon", fontsize=6.0, color=S.GRAY,
             va="center", style="italic")
    groups = [("BENEFIT", S.BLUE, [(r"$+1.87$ dB", "fixed-dwell"),
                                   (r"$19.1\times$", "elapsed time")]),
              ("COST", S.ORANGE, [(r"$37.1\times$", "incident dose"),
                                  (r"$2.6\times$", "detected counts")])]
    x0 = 0.40
    for gname, col, items in groups:
        axr.text(x0 - 0.03, 0.95, gname, fontsize=6.0, color=col,
                 fontweight="bold", ha="left", va="center")
        for val, lab in items:
            axr.text(x0, 0.60, val, fontsize=8.6, fontweight="bold", color=col,
                     ha="center", va="center")
            axr.text(x0, 0.12, lab, fontsize=6.2, color=col, ha="center",
                     va="center")
            x0 += 0.15
        x0 += 0.03
    S.save(fig, OUT, "fig_speed_results")


# ========================================================================== #
#  M6 -- supplement audit figure (guard path + certificate status)            #
# ========================================================================== #
def _read_cert_csv(path):
    with open(path, encoding="utf-8") as f:
        hdr = f.readline().strip().split(",")
        rows = []
        for line in f:
            parts = line.rstrip("\n").split(",")
            if len(parts) > len(hdr):
                i0 = hdr.index("solver_status_primary")
                extra = len(parts) - len(hdr)
                parts = (parts[:i0] + [",".join(parts[i0:i0 + 1 + extra])]
                         + parts[i0 + 1 + extra:])
                if len(parts) > len(hdr):
                    i1 = hdr.index("solver_status_retry")
                    extra = len(parts) - len(hdr)
                    parts = (parts[:i1] + [",".join(parts[i1:i1 + 1 + extra])]
                             + parts[i1 + 1 + extra:])
            rows.append(dict(zip(hdr, parts)))
    return rows


def _fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def fig_audit_supp():
    S.use_style()
    # data: guard trajectory (a) + certificate cells (c)
    aj = os.path.join(M1, "ACTIII_A_GUARD_TRAJECTORY.json")
    have_a = os.path.exists(aj)
    cert = []
    for p in sorted(glob.glob(os.path.join(SHARDS, "M1_CERT_*.csv"))):
        cert.extend(_read_cert_csv(p))

    fig, (axg, axb, axe) = plt.subplots(1, 3, figsize=(S.COL2_IN, 2.7),
                                        gridspec_kw={"width_ratios":
                                                     [1.0, 1.0, 0.9]})
    fig.subplots_adjust(left=0.075, right=0.975, top=0.83, bottom=0.32,
                        wspace=0.42)

    # (a) guard-collapse path
    if have_a:
        d = json.load(open(aj))
        alpha = np.array(d["alpha"]); dose = np.array(d["dose_dev"])
        arisk = np.array(d["arisk_ratio"]); mu = np.array(d["mu_min"])
        bars = d["guard_bars"]
        o = np.argsort(alpha)
        alpha, dose, arisk, mu = alpha[o], dose[o], arisk[o], mu[o]
        axg.set_yscale("log")
        axg.axvspan(alpha[alpha > 0].min(), 1.02, color=S.ORANGE, alpha=0.06,
                    zorder=0)
        axg.plot(alpha, dose, color=S.ORANGE, lw=1.6, label="dose dev.")
        axg.plot(alpha, mu, color=S.GRAY, lw=1.6, label=r"floor $\mu_{\min}$")
        axg.plot(alpha, arisk, color=S.BLUE, lw=1.6, label="A-risk")
        for val, col in ((bars["dose_dev"], S.ORANGE),
                         (bars["spectral_mu_min"], S.GRAY),
                         (bars["arisk_ratio"], S.BLUE)):
            axg.axhline(val, color=col, ls="--", lw=0.8, alpha=0.8)
        axg.set_ylim(0.03, 15)
        axg.set_xlim(-0.02, 1.02)
        axg.set_xlabel(r"OED fraction $\alpha=m/972$", fontsize=7.2)
        axg.set_ylabel("guard quantity", fontsize=7.2)
        axg.text(0.97, 0.03, "no feasible point", transform=axg.transAxes,
                 fontsize=6.4, color=S.ORANGE, ha="right", va="bottom",
                 fontweight="bold")
        axg.legend(loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=3,
                   frameon=False, fontsize=6.0, handlelength=1.3,
                   columnspacing=0.9)
    else:
        axg.text(0.5, 0.5, "guard trajectory\ncache absent",
                 transform=axg.transAxes, ha="center", va="center",
                 fontsize=7, color=S.GRAY)
        axg.axis("off")
    S.panel_title(axg, "guard collapse", size=7.6, y=1.0, pad=12)
    axg.annotate("development only", xy=(0.5, 1.0), xycoords="axes fraction",
                 xytext=(0, 2), textcoords="offset points", ha="center",
                 va="bottom", fontsize=6.2, color=S.GRAY, style="italic")

    # (b) certificate status counts by anchor
    anchors = [("nu200_b0.05", 200.0, 0.05), ("nu200_b0.6", 200.0, 0.60),
               ("nu2000_b0.05", 2000.0, 0.05), ("nu2000_b0.6", 2000.0, 0.60)]
    counts = {a[0]: {"COUNTEREXAMPLE": 0, "NUMERICAL_UNRESOLVED": 0,
                     "CERTIFIED": 0} for a in anchors}
    for r in cert:
        k = "nu%g_b%g" % (_fnum(r["nu"]), _fnum(r["b"]))
        if k in counts:
            counts[k][r["status"]] += 1
    x = np.arange(4)
    ce = [counts[a[0]]["COUNTEREXAMPLE"] for a in anchors]
    un = [counts[a[0]]["NUMERICAL_UNRESOLVED"] for a in anchors]
    axb.bar(x, ce, color=S.ORANGE, width=0.64, label="counterexample")
    axb.bar(x, un, bottom=ce, color=S.GRAY, width=0.64, label="unresolved")
    for xi, c_, u_ in zip(x, ce, un):
        axb.text(xi, c_ / 2, str(c_), ha="center", va="center", color="white",
                 fontsize=6.6)
        axb.text(xi, c_ + u_ / 2, str(u_), ha="center", va="center",
                 color="white", fontsize=6.6)
    axb.set_xticks(x)
    axb.set_xticklabels([r"$\nu$200" "\n" r"$b$.05", r"$\nu$200" "\n" r"$b$.6",
                         r"$\nu$2k" "\n" r"$b$.05", r"$\nu$2k" "\n" r"$b$.6"],
                        fontsize=6.2)
    axb.set_ylabel("certificate cells", fontsize=7.2)
    axb.set_ylim(0, 128)
    axb.legend(loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=2,
               frameon=False, fontsize=6.0, handlelength=1.2)
    S.panel_title(axb, "0/480 certified", size=7.6, color=S.ORANGE, y=1.0,
                  pad=12)
    axb.annotate("descriptive", xy=(0.5, 1.0), xycoords="axes fraction",
                 xytext=(0, 2), textcoords="offset points", ha="center",
                 va="bottom", fontsize=6.2, color=S.GRAY, style="italic")

    # (c) ECDF of finite counterexample primal gaps
    gaps = np.array([_fnum(r["primal_gap_lower_per_r"]) for r in cert
                     if r["status"] == "COUNTEREXAMPLE"
                     and np.isfinite(_fnum(r["primal_gap_lower_per_r"]))])
    xs = np.sort(gaps)
    ecdf = np.arange(1, xs.size + 1) / xs.size
    axe.step(xs, ecdf, where="post", color=S.ORANGE, lw=1.6)
    axe.axvline(1e-2, color=S.GRAY, ls="--", lw=1.0)
    axe.text(1.3e-2, 0.5, r"bar $10^{-2}$", fontsize=6.2, color=S.GRAY,
             ha="left", va="center")
    axe.set_xscale("log")
    axe.set_xlim(5e-3, 3.0)
    axe.set_ylim(0, 1.02)
    axe.set_xlabel(r"primal gap / $r$", fontsize=7.2)
    axe.set_ylabel("ECDF", fontsize=7.2)
    S.panel_title(axe, "gap $\\gg$ bar", size=7.6, y=1.0, pad=3)

    fig.text(0.012, 0.97, "(a)", fontsize=9, fontweight="bold", va="top")
    fig.text(0.365, 0.97, "(b)", fontsize=9, fontweight="bold", va="top")
    fig.text(0.70, 0.97, "(c)", fontsize=9, fontweight="bold", va="top")
    S.save(fig, OUT, "fig_audit_supp")


# ========================================================================== #
#  M6 -- compact dose-only headroom (restyle of actiii_b)                     #
# ========================================================================== #
def _read_md_table(path, must_have):
    text = open(path, encoding="utf-8").read()
    lines = text.splitlines()
    i, n = 0, len(lines)
    while i < n:
        ln = lines[i].strip()
        if (ln.startswith("|") and i + 1 < n
                and set(lines[i + 1].strip()) <= set("|-: ")):
            header = [c.strip() for c in ln.strip("|").split("|")]
            rows = []
            j = i + 2
            while j < n and lines[j].strip().startswith("|"):
                rows.append([c.strip()
                             for c in lines[j].strip().strip("|").split("|")])
                j += 1
            if all(m in header for m in must_have):
                return header, rows
            i = j
        else:
            i += 1
    raise SystemExit("no md table with %s in %s" % (must_have, path))


def actiii_b_compact():
    S.use_style()
    path = os.path.join(M1, "R18_GAP_PROBE_REPLICATION.md")
    hdr, rows = _read_md_table(
        path, ["scene", "nu", "b", "DOSE_ONLY_PRIMAL_GAP_PER_R"])
    ig = hdr.index("DOSE_ONLY_PRIMAL_GAP_PER_R")
    isc, inu, ib = hdr.index("scene"), hdr.index("nu"), hdr.index("b")
    data = {}
    for r in rows:
        fam = r[isc].replace("m1_dev_", "")
        nu, b = _fnum(r[inu]), _fnum(r[ib])
        anchor = "A" if (nu == 2000.0 and b == 0.60) else "B"
        data.setdefault(fam, {})[anchor] = _fnum(r[ig])
    order = ["glyph", "chirp", "maze", "spokes", "contour", "microtexture"]
    fams = [f for f in order if f in data]

    fig, ax = plt.subplots(figsize=(S.COL1_IN, 2.5))
    fig.subplots_adjust(left=0.17, right=0.83, top=0.86, bottom=0.26)
    x = np.arange(len(fams))
    gA = [data[f]["A"] for f in fams]
    gB = [data[f]["B"] for f in fams]
    for xi, a, bb in zip(x, gA, gB):
        ax.plot([xi, xi], [min(a, bb), max(a, bb)], color=S.GRAY_L, lw=1.0,
                zorder=1)
    ax.plot(x, gA, "o", color=S.BLUE, ms=6, zorder=3, label=r"A: $\nu$2000")
    ax.plot(x, gB, "D", color=S.ORANGE, ms=5, zorder=3, label=r"B: $\nu$200")
    ax.axhline(0.0, color="black", lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(fams, rotation=30, ha="right", fontsize=6.8)
    ax.set_ylabel(r"dose-only primal gap / $r$", fontsize=7.4)
    ax.set_ylim(0.0, 1.05)
    ax.set_xlim(-0.5, len(fams) - 0.5)
    sec = ax.secondary_yaxis("right", functions=(np.exp, np.log))
    sec.set_ylabel(r"$D$-eff. lower bound", fontsize=7.4)
    sec.set_yticks([1.0, 1.5, 2.0, 2.5])
    ax.legend(loc="lower left", frameon=False, fontsize=6.6, ncol=2,
              handletextpad=0.4, columnspacing=1.0)
    ax.text(0.5, 1.02, "development only; no gate", transform=ax.transAxes,
            ha="center", fontsize=6.4, color=S.GRAY, style="italic")
    S.panel_title(ax, "geometry headroom exists", size=7.6, y=1.06)
    with np.errstate(divide="ignore", invalid="ignore"):
        S.save(fig, OUT, "actiii_b")


# ========================================================================== #
def main(argv=None):
    which = [a.lower() for a in (argv or sys.argv[1:])]
    todo = {
        "mech": fig_mechanism, "m3": fig_mechanism,
        "jitter": fig_jitter_validation, "m4": fig_jitter_validation,
        "speed": fig_speed_results, "m5": fig_speed_results,
        "audit": fig_audit_supp, "m6": fig_audit_supp,
        "actiiib": actiii_b_compact, "b": actiii_b_compact,
    }
    if not which:
        which = ["mech", "jitter", "speed", "audit", "actiiib"]
    for w in which:
        fn = todo.get(w)
        if fn is None:
            print("[r20] unknown figure '%s'" % w, flush=True)
            continue
        try:
            fn()
        except SystemExit as e:
            print("[r20] SKIP %s: %s" % (w, e), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
