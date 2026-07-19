"""Study-1 figure, panel (d) — empirical acquisition-speedup S versus mean load
rho_f, with the count-rate (Fisher-rate) envelope (frozen prereg
docs/ROUND63_FIGURE_PREREG.md panel d).

Source: results/round63_s2_detail/pilot_rows.csv (RQL arm). For every DETAIL-24
image and each fast mean-load rho_f present in the data, the acquisition speedup
S is the Q90 crossing-time ratio of the safe (rho=0.05) vs fast (rho_f) dwell
sweeps, computed with the SAME recipe as
code/round63/dev_pilot_patch.py's endpoint section:

  * per (image, rho): seed-mean PSNR_rad at each nu, sorted by nu;
  * PAVA isotonic regression of the seed-mean PSNR_rad against log nu;
  * safe informativeness floor R = q_safe[-1] - q_safe[0]; if R < 0.5 dB the
    safe range is uninformative  -> S = 1 (censoring bound);
  * Q90 = q_safe[0] + 0.9 R;
  * T_safe = exp(interp(Q90, q_safe, log nu));  if the fast sweep never reaches
    Q90 (q_fast[-1] < Q90) it is right-censored  -> S = 0 (censoring bound);
  * else T_fast = exp(interp(Q90, q_fast, log nu)) and S = T_safe / T_fast.

(Using log nu is identical to using log T_opt for the ratio: T_opt = M nu tau,
so the constant M tau cancels in T_safe / T_fast.)

Plot: per-image S as jittered grey dots on a log y-axis, a median marker per
rho_f, a dashed S = 3 'Study-1 gate bar', and the exact count-rate envelope
ratio(rho_f) = [rho_f/(1+rho_f)] / [0.05/1.05] as a dotted reference. Censored
images (S = 0) are drawn as open markers at a log-axis floor and their counts,
together with the S = 1 uninformative-bound counts, are annotated under the axis.

PROVENANCE: PSNR_rad read verbatim; '' / 'nan' dropped, never imputed.

Run (cwd = repo root D:\\GI_another):
    D:\\Anacondar\\anaconda3\\python.exe code\\round63\\figs\\fig_s1d_svsrho.py
"""
import os
import sys
import csv
import collections

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
R63 = os.path.dirname(HERE)
REPO = os.path.dirname(os.path.dirname(R63))

DATA = os.path.join(REPO, "results", "round63_s2_detail", "pilot_rows.csv")
OUT_DIR = os.path.join(REPO, "paper", "figs")
OUT_PNG = os.path.join(OUT_DIR, "fig_s1d_svsrho.png")
OUT_PDF = os.path.join(OUT_DIR, "fig_s1d_svsrho.pdf")

ARM = "RQL"
RHO_SAFE = 0.05
RHO_FAST_CANDIDATES = [0.3, 0.6, 1.0, 2.0]     # keep those present in the data
R_FLOOR = 0.5                                  # safe informativeness floor (dB)
GATE_BAR = 3.0                                 # Study-1 speedup gate bar
S_FLOOR = 0.1                                  # log-axis placement of S = 0

C_BLUE = "#0072B2"
C_ORANGE = "#D55E00"
C_DOT = "0.55"
JIT_W = 0.035                                  # jitter half-width in rho units


def _to_float(s):
    if s is None:
        return None
    s = s.strip()
    if s == "" or s.lower() == "nan":
        return None
    try:
        return float(s)
    except ValueError:
        return None


def load_rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def pava(y):
    """Pool-adjacent-violators isotonic regression (non-decreasing), verbatim
    from dev_pilot_patch.py."""
    y = list(y)
    n = len(y)
    w = [1.0] * n
    i = 0
    while i < n - 1:
        if y[i] > y[i + 1] + 1e-12:
            ym = (y[i] * w[i] + y[i + 1] * w[i + 1]) / (w[i] + w[i + 1])
            y[i] = y[i + 1] = ym
            w[i] = w[i + 1] = w[i] + w[i + 1]
            i = max(i - 1, 0)
        else:
            i += 1
    return y


def seedmean_curve(rows, image, rho):
    """(log nu array, isotonic seed-mean PSNR_rad array) for one (image, rho)."""
    byn = collections.defaultdict(list)
    for r in rows:
        if r["arm"] != ARM or r["image"] != image:
            continue
        if _to_float(r["rho_bar"]) != rho:
            continue
        pv = _to_float(r["PSNR_rad"])
        if pv is None:
            continue
        byn[_to_float(r["nu"])].append(pv)
    nus = sorted(byn)
    means = [float(np.mean(byn[nu])) for nu in nus]
    return np.log(np.array(nus)), np.array(pava(means))


def s_gate(rows, image, rho_fast):
    """Return (S, tag) following the frozen dev_pilot recipe. tag in
    {'ok', 'safe_uninformative'(S=1), 'fast_censored'(S=0)}."""
    ls, qs = seedmean_curve(rows, image, RHO_SAFE)
    lf, qf = seedmean_curve(rows, image, rho_fast)
    R = qs[-1] - qs[0]
    if R < R_FLOOR:
        return 1.0, "safe_uninformative"
    Q90 = qs[0] + 0.9 * R
    Ts = np.exp(np.interp(Q90, qs, ls))
    if qf[-1] < Q90:
        return 0.0, "fast_censored"
    Tf = np.exp(np.interp(Q90, qf, lf))
    return float(Ts / Tf), "ok"


def envelope(rho_f):
    """Exact count-rate (Fisher-rate) ratio ceiling vs the safe load."""
    return (rho_f / (1.0 + rho_f)) / (RHO_SAFE / (1.0 + RHO_SAFE))


def main():
    rows = load_rows(DATA)
    print("[fig_s1d] loaded %d rows from %s" % (len(rows), DATA), flush=True)

    present = sorted({_to_float(r["rho_bar"]) for r in rows
                      if r["arm"] == ARM})
    rho_fast = [r for r in RHO_FAST_CANDIDATES if r in present]
    print("[fig_s1d] rho_bar present (RQL): %s" % present, flush=True)
    print("[fig_s1d] fast loads used: %s" % rho_fast, flush=True)

    images = sorted({r["image"] for r in rows
                     if r["arm"] == ARM and r["image"].startswith("detail_")})
    print("[fig_s1d] %d DETAIL images" % len(images), flush=True)

    # S[rho_f] = {image: (S, tag)}
    S_by_rho = {}
    for rf in rho_fast:
        d = {img: s_gate(rows, img, rf) for img in images}
        S_by_rho[rf] = d
        vals = np.array([v for v, _t in d.values()])
        n0 = sum(1 for _v, t in d.values() if t == "fast_censored")
        n1 = sum(1 for _v, t in d.values() if t == "safe_uninformative")
        nok = sum(1 for _v, t in d.values() if t == "ok")
        print("[fig_s1d] --- rho_f=%.2f  envelope=%.3f ---" % (rf, envelope(rf)),
              flush=True)
        print("[fig_s1d]   median S = %.4f  (n=%d; ok=%d, S=0 censored=%d, "
              "S=1 uninformative=%d; S>=3: %d)"
              % (float(np.median(vals)), len(d), nok, n0, n1,
                 int(np.sum(vals >= GATE_BAR))), flush=True)
        for img in images:
            v, t = d[img]
            print("[fig_s1d]     %-24s S=%8.4f  %s" % (img, v, t), flush=True)

    # ---------------- figure ----------------
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 8,
        "axes.linewidth": 0.6,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    fig, ax = plt.subplots(figsize=(3.6, 2.9))
    ax.set_yscale("log")

    for i, rf in enumerate(rho_fast):
        d = S_by_rho[rf]
        vals = np.array([v for v, _t in d.values()])
        pos = np.array([v for v, t in d.values() if t != "fast_censored"])
        zero = np.array([v for v, t in d.values() if t == "fast_censored"])
        rng = np.random.RandomState(4000 + i)
        # informative + S=1 uninformative bound -> grey dots
        if pos.size:
            xj = rf + rng.uniform(-JIT_W, JIT_W, size=pos.size)
            ax.scatter(xj, pos, s=9, color=C_DOT, alpha=0.75, linewidths=0,
                       zorder=5)
        # S=0 fast-censored -> open markers at the log-axis floor
        if zero.size:
            xj0 = rf + rng.uniform(-JIT_W, JIT_W, size=zero.size)
            ax.scatter(xj0, np.full(zero.size, S_FLOOR), s=16,
                       facecolors="none", edgecolors=C_DOT, linewidths=0.8,
                       zorder=2)
        med = float(np.median(vals))              # median over ALL images
        med_plot = med if med > 0 else S_FLOOR
        ax.scatter([rf], [med_plot], s=54, marker="D", color=C_ORANGE,
                   edgecolors="black", linewidths=0.7, zorder=4)

    # count-rate (Fisher-rate) envelope, dotted reference
    xr = np.linspace(min(rho_fast) - 0.05, max(rho_fast) + 0.1, 200)
    ax.plot(xr, [envelope(r) for r in xr], ls=":", lw=1.3, color=C_BLUE,
            zorder=3)

    # gate bar S = 3
    ax.axhline(GATE_BAR, ls="--", lw=0.9, color="0.4", zorder=1)
    ax.text(max(rho_fast) + 0.07, GATE_BAR, "Study-1 gate bar (S=3)",
            va="bottom", ha="right", fontsize=6.0, color="0.35", zorder=5)

    ax.set_xlabel(r"fast mean load $\bar\rho_f$")
    ax.set_ylabel(r"acquisition speedup $S$")
    ax.set_xlim(min(rho_fast) - 0.12, max(rho_fast) + 0.18)
    ax.set_ylim(0.07, max(20.0, envelope(max(rho_fast)) * 1.5))
    ax.set_xticks(rho_fast)
    ax.set_xticklabels([("%.1f" % r) for r in rho_fast])

    # censoring-count annotations under the axis
    y_ann = 0.075
    for rf in rho_fast:
        d = S_by_rho[rf]
        n0 = sum(1 for _v, t in d.values() if t == "fast_censored")
        n1 = sum(1 for _v, t in d.values() if t == "safe_uninformative")
        ax.annotate("S0:%d\nS1:%d" % (n0, n1), xy=(rf, y_ann),
                    xycoords=("data", "axes fraction"), ha="center",
                    va="bottom", fontsize=5.4, color="0.4")
    ax.text(0.5, -0.20, "S0 = fast right-censored,  S1 = safe uninformative",
            transform=ax.transAxes, ha="center", va="top", fontsize=5.6,
            color="0.4")

    handles = [
        Line2D([0], [0], color=C_BLUE, ls=":", lw=1.3),
        Line2D([0], [0], marker="D", color=C_ORANGE, ls="none", ms=6,
               markeredgecolor="black", markeredgewidth=0.7),
        Line2D([0], [0], marker="o", color=C_DOT, ls="none", ms=4),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="none",
               markeredgecolor=C_DOT, ls="none", ms=5, markeredgewidth=0.8),
    ]
    labels = ["count-rate ratio", "median $S$", "per-image $S$",
              "censored ($S{=}0$)"]
    ax.legend(handles, labels, loc="lower left", bbox_to_anchor=(0.02, 0.16),
              fontsize=5.8, frameon=True, framealpha=0.9, borderpad=0.4,
              handlelength=1.6, labelspacing=0.3)

    fig.subplots_adjust(left=0.17, right=0.97, top=0.96, bottom=0.24)
    os.makedirs(OUT_DIR, exist_ok=True)
    fig.savefig(OUT_PDF)
    fig.savefig(OUT_PNG, dpi=200)
    plt.close(fig)
    print("[fig_s1d] wrote %s" % OUT_PDF, flush=True)
    print("[fig_s1d] wrote %s (200 dpi)" % OUT_PNG, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
