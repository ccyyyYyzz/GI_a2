"""
FIG 2 - mean-layer certificate.
One sentence: the mean layer has a reachable optimum and no material residual
count-only headroom in the declared class.

Panels (1x3):
  (a) hidden-state / jitter information ridge with the selected operating point
  (b) the two confirmed cohort statistics (+1.87 dB, 19.13x) with lower bounds
  (c) closure strip: five preregistered interventions below the materiality line,
      harmful arm visually capped

No method names, solver labels, or audit codes.
Sources: T_fig2_ridge_curves.csv, T_fig2_ridge_points.json, T_fig2_cohorts.json,
         T_fig2_closure.csv
"""
import os, sys, collections
sys.path.insert(0, os.path.dirname(__file__))
import common as K
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle

# ---- load sources --------------------------------------------------------
_, ridge_rows = K.load_csv("T_fig2_ridge_curves.csv")
RP = K.load_json("T_fig2_ridge_points.json")
CO = K.load_json("T_fig2_cohorts.json")
_, CL = K.load_csv("T_fig2_closure.csv")

by_nu = collections.defaultdict(list)
for r in ridge_rows:
    by_nu[int(float(r["nu"]))].append((float(r["rho"]), float(r["I_exact"]), float(r["I_clt"])))
for nu in by_nu:
    by_nu[nu].sort()

peaks = {p["nu"]: p for p in RP["ridge_peaks"]}
OP = RP["selected_operating_point"]


def panel_ridge(ax):
    show_nu = [20, 200, 2000]
    shades = {20: K.C["sky"], 200: K.C["blue"], 2000: "#023E5A"}
    for nu in show_nu:
        d = np.array(by_nu[nu])
        ax.plot(d[:, 0], d[:, 1], "-", color=shades[nu], lw=1.6,
                label=r"$\nu=%d$" % nu, zorder=3)
    # ridge locus = peaks
    px = [peaks[nu]["rho_star"] for nu in sorted(peaks)]
    py = [peaks[nu]["I_at_peak"] for nu in sorted(peaks)]
    ax.plot(px, py, "--", color=K.C["gray"], lw=1.1, zorder=2)
    ax.scatter(px, py, s=14, facecolor="white", edgecolor=K.C["gray"],
               linewidth=0.9, zorder=4)
    ax.text(2.15, 0.685, "information ridge  $\\rho^\\star(\\nu)$",
            fontsize=6.8, color=K.C["line"], rotation=19, ha="left", va="bottom")

    # selected operating point (nu=2000)
    ax.plot(OP["rho_star_deterministic"], OP["I_at_peak"], "o", ms=9,
            mfc=K.C["green"], mec="white", mew=1.2, zorder=6)
    ax.annotate("selected operating point\n($\\nu=2000$, on the ridge)",
                (OP["rho_star_deterministic"], OP["I_at_peak"]),
                (3.2, 0.55), fontsize=6.8, color=K.C["green"], fontweight="bold",
                ha="left", va="center",
                arrowprops=dict(arrowstyle="-", color=K.C["green"], lw=0.9,
                                connectionstyle="arc3,rad=0.2"))
    # jitter cap annotation (finite illumination jitter pulls the optimum down)
    jc = RP["jitter_cap"]["example_cv0p05_nu2000"]
    ax.axvline(jc["capped_rho_c"], color=K.C["red"], lw=1.0, ls=(0, (4, 2)), zorder=2)
    ax.text(jc["capped_rho_c"] * 1.05, 0.16,
            "jitter cap\n($c_v{=}0.05$):\n$\\rho^\\star{\\to}%.1f$" % jc["capped_rho_c"],
            fontsize=6.2, color=K.C["red"], ha="left", va="bottom")

    ax.set_xscale("log")
    ax.set_xlim(0.7, 64)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("hidden-state load  $\\rho$")
    ax.set_ylabel("Fisher information per unit time")
    ax.legend(loc="lower right", frameon=False, fontsize=6.6, handlelength=1.3,
              borderpad=0.2, labelspacing=0.25)
    ax.set_title("(a)  the living ridge", fontsize=8.5, loc="left", fontweight="bold")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def _cohort_row(ax, value, lb, thr, xmax, unit, label, npos, ntot, materiality=None):
    ax.add_patch(Rectangle((thr, 0.30), value - thr, 0.34, facecolor=K.C["darkfill"],
                           edgecolor="none", zorder=2))
    # lower-bound whisker
    ax.plot([lb, lb], [0.24, 0.70], color=K.C["ink"], lw=1.4, zorder=4)
    ax.plot(value, 0.47, "o", ms=7, mfc=K.C["green"], mec="white", mew=1.0, zorder=5)
    ax.axvline(thr, color=K.C["line"], lw=0.9, zorder=1)
    if materiality is not None:
        ax.axvline(materiality, color=K.C["gray"], lw=0.9, ls=(0, (3, 2)), zorder=1)
        ax.text(materiality, 0.02, "1 dB bar", ha="center", va="bottom",
                fontsize=5.8, color=K.C["gray"])
    # category label above the bar, left-aligned at threshold (stays in panel)
    ax.text(thr, 0.80, label, ha="left", va="bottom", fontsize=7.4,
            fontweight="bold", color=K.C["ink"])
    ax.text(value, 0.72, "%.2f%s" % (value, unit), ha="center", va="bottom",
            fontsize=8.0, fontweight="bold", color=K.C["green"])
    # n_positive INSIDE the bar (white), clear of the whisker/LB label
    ax.text(thr + 0.04 * (value - thr), 0.47, "%d/%d positive" % (npos, ntot),
            ha="left", va="center", fontsize=6.3, color="white", zorder=6)
    ax.text(lb, 0.20, "LB %.2f" % lb, ha="center", va="top", fontsize=6.2,
            color=K.C["ink"])
    ax.text(thr, 0.04, "parity", ha="left", va="bottom", fontsize=5.8, color=K.C["line"])
    ax.set_xlim(thr - 0.05 * (xmax - thr), xmax + 0.02 * (xmax - thr))
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)


def panel_cohorts(gs, fig):
    inner = gs.subgridspec(2, 1, hspace=0.95)
    a0 = fig.add_subplot(inner[0])
    a1 = fig.add_subplot(inner[1])
    c0, c1 = CO["cohorts"]
    _cohort_row(a0, c0["value"], c0["lower_bound"], 0.0, 2.4, " dB",
                "fixed-dwell quality", c0["n_positive"], c0["n_total"],
                materiality=CO["materiality_bar_dB"])
    a0.set_xlabel("median gain over comparator (dB)", fontsize=6.8, labelpad=1)
    _cohort_row(a1, c1["value"], c1["lower_bound"], 1.0, 22.0, "$\\times$",
                "elapsed speed", c1["n_positive"], c1["n_total"])
    a1.set_xlabel("median elapsed-time speed ratio ($\\times$)", fontsize=6.8, labelpad=1)
    a0.set_title("(b)  certified cohort statistics", fontsize=8.5, loc="left",
                 fontweight="bold")


def panel_closure(ax):
    rows = CL          # 5 interventions
    bar = float(rows[0]["materiality_bar_dB"])
    labels, vals, verds = [], [], []
    for r in rows:
        labels.append(r["mechanism"])
        vals.append(float(r["governing_number"]))
        verds.append(r["verdict"])
    y = np.arange(len(rows))[::-1].astype(float)
    FLOOR = -2.4       # visual cap for the harmful arm
    verd_word = {"FAIL": "fail", "HARM": "harm", "SATURATED": "saturated", "KILL": "kill"}
    for yi, (v, vd, lab) in zip(y, zip(vals, verds, labels)):
        capped = v < FLOOR
        vv = FLOOR + 0.35 if capped else v
        col = K.C["red"] if vd == "HARM" else (K.C["amber"] if vd == "SATURATED" else K.C["gray"])
        ax.barh(yi, vv, height=0.5, color=col, edgecolor="none", zorder=3)
        # mechanism label sits ABOVE its bar, left-aligned (stays inside panel)
        ax.text(FLOOR, yi + 0.34, lab, ha="left", va="bottom", fontsize=6.5,
                color=K.C["ink"])
        if capped:                     # break marks + true value label
            for dx in (0.18, 0.30):
                ax.plot([FLOOR + 0.35 + dx, FLOOR + 0.47 + dx],
                        [yi + 0.16, yi - 0.16], color="white", lw=1.0, zorder=5)
            ax.text(FLOOR + 0.12, yi, "%.1f dB (harm)" % v, ha="left", va="center",
                    fontsize=6.2, color="white", fontweight="bold", zorder=6)
        else:
            tag = ("$\\leq$%.1f dB (%s)" % (v, verd_word[vd])) if vd == "SATURATED" \
                  else ("%.2f dB (%s)" % (v, verd_word[vd]))
            ax.text(max(v, 0) + 0.07, yi, tag, ha="left", va="center", fontsize=6.2,
                    color=K.C["ink"])
    ax.axvline(0, color=K.C["line"], lw=0.8, zorder=2)
    ax.axvline(bar, color=K.C["gray"], lw=1.1, ls=(0, (4, 2)), zorder=2)
    ax.text(bar, y.max() + 0.55, "materiality (1 dB)", ha="right", va="bottom",
            fontsize=6.3, color=K.C["gray"])
    ax.set_xlim(FLOOR - 0.05, 1.75)
    ax.set_ylim(-0.6, y.max() + 0.95)
    ax.set_yticks([])
    ax.set_xticks([-2, -1, 0, 1])
    ax.set_xlabel("deployable image gain (dB)")
    ax.set_title("(c)  closure: five interventions, all sub-materiality", fontsize=8.0,
                 loc="left", fontweight="bold")
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0)


def build():
    K.apply_style()
    fig = plt.figure(figsize=(K.cm(17.8), K.cm(6.6)))
    gs = GridSpec(1, 3, width_ratios=[1.15, 1.15, 1.25], wspace=0.42,
                  left=0.055, right=0.985, top=0.86, bottom=0.16)
    panel_ridge(fig.add_subplot(gs[0]))
    panel_cohorts(gs[1], fig)
    panel_closure(fig.add_subplot(gs[2]))
    fig.suptitle("The mean layer: a reachable optimum with no material residual headroom",
                 x=0.055, y=0.985, ha="left", fontsize=9.5, fontweight="bold")
    return fig


if __name__ == "__main__":
    fig = build()
    K.save(fig, "fig2")
    plt.close(fig)
    print("FIG2 done.")
