"""Optics-Express figure, panel (e) — occupancy-ladder gate diagram: fixed-dwell
quality gain and the (frozen, k=16-only) acquisition-speedup gate.

Two stacked single-column subplots sharing the log2 occupancy axis
k in {512, 32, 16, 1} (left -> right; labelled with occupancy percent):

  * TOP    : per-image fixed-dwell quality gain DeltaQ (rho_bar 0.6 vs 0.05 at
    nu = 2000, RQL), as jittered small grey dots (one per detail-32 image) with a
    large median marker per rung. Sources per rung:
        k = 512, 1  -> controls_rows.csv  (per-image seed-mean fast - safe)
        k = 32      -> robustness_rows.csv (per-image seed-mean fast - safe)
        k = 16      -> study2_summary.json per-image DeltaQ_budget (the frozen
                       budget-matched secondary endpoint; read verbatim)
    A dashed line marks 0 dB; a dotted line marks the +1 dB preregistered
    secondary bar (only k=16 is gated).
  * BOTTOM : the acquisition-speedup gate S is carried by k=16 ONLY (frozen
    ruling). Per-image S_gate values (study2_summary.json) are shown as jittered
    dots at the k=16 rung on a log y-axis, with a median marker and a dashed
    S = 3 gate bar. S_gate = 0 (censored) values are drawn as OPEN markers at
    y = 0.1 (log axis cannot render 0; flagged in the printout). The other rungs
    carry no gate statistic (frozen): a light grey band + 'no gate (frozen)'.

PROVENANCE — CSV values are read verbatim; empty / 'nan' cells are dropped from
every mean (never imputed). k=16 statistics are read from the frozen analyzer
output study2_summary.json. Jitter is deterministic (fixed RandomState per rung);
it perturbs only the horizontal position for legibility, never a plotted value.

Run (cwd = repo root D:\\GI_another):
    D:\\Anacondar\\anaconda3\\python.exe code\\round63\\figs\\fig_e_gate.py
"""
import os
import sys
import csv
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))          # code/round63/figs
R63 = os.path.dirname(HERE)                                # code/round63
REPO = os.path.dirname(os.path.dirname(R63))               # repo root (D:\GI_another)

RES = os.path.join(REPO, "results", "round63_study2")
PRIMARY_JSON = os.path.join(RES, "study2_summary.json")
ROBUST_CSV = os.path.join(RES, "robustness_rows.csv")
CONTROLS_CSV = os.path.join(RES, "controls_rows.csv")

OUT_DIR = os.path.join(REPO, "paper", "figs")
OUT_PNG = os.path.join(OUT_DIR, "fig_e_gate.png")
OUT_PDF = os.path.join(OUT_DIR, "fig_e_gate.pdf")

N_PIX = 1024
K_LADDER = [512, 32, 16, 1]        # left -> right (descending occupancy ladder)
K_LABELS = {                       # occupancy-percent labels (match panels b/c)
    512: "512 (50%)",
    32:  "32 (3.1%)",
    16:  "16 (1.6%)",
    1:   "1 (0.1%)",
}
K16 = 16
ARM = "RQL"
NU_FIXED = 2000.0                  # fixed-dwell operating dwell
RHO_SAFE, RHO_FAST = 0.05, 0.6

PREREG_BAR = 1.0                   # +1 dB preregistered secondary bar (k=16 gate)
GATE_BAR = 3.0                     # S = 3 acquisition-speedup gate bar
S_FLOOR = 0.1                      # log-axis placement for censored S_gate = 0

C_BLUE = "#0072B2"                 # colourblind-safe pair
C_ORANGE = "#D55E00"
C_DOT = "0.55"                     # jittered per-image points (grey)
JIT_W = 0.16                       # jitter half-width in log2 units


def _to_float(s):
    """Parse a cell to float; '' / 'nan' -> None (dropped from means, never imputed)."""
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


def deltaq_from_csv(path, k):
    """Per-image fixed-dwell DeltaQ = seed-mean(fast) - seed-mean(safe) at
    (k, nu=NU_FIXED, RQL). Returns {image: delta} over images having both arms."""
    rows = load_rows(path)
    by = {}
    for r in rows:
        if r["arm"] != ARM or _to_float(r["k"]) != float(k):
            continue
        if _to_float(r["nu"]) != NU_FIXED:
            continue
        rho = _to_float(r["rho_bar"])
        pv = _to_float(r["PSNR_rad"])
        if pv is None or rho not in (RHO_SAFE, RHO_FAST):
            continue
        by.setdefault(r["image"], {}).setdefault(rho, []).append(pv)
    out = {}
    for img, d in by.items():
        if RHO_SAFE in d and RHO_FAST in d:
            out[img] = float(np.mean(d[RHO_FAST]) - np.mean(d[RHO_SAFE]))
    return out


def load_k16():
    """k=16 per-image DeltaQ_budget and S_gate from the frozen analyzer output."""
    with open(PRIMARY_JSON) as f:
        j = json.load(f)
    per = j["endpoint_analysis"]["per_image"]
    dq = {p["image"]: float(p["DeltaQ_budget"]) for p in per}
    sg = {p["image"]: float(p["S_gate"]) for p in per}
    return dq, sg


def _jitter(k, n, seed):
    """Deterministic horizontal jitter: k * 2**U(-JIT_W, JIT_W), n points."""
    rng = np.random.RandomState(seed)
    return float(k) * np.power(2.0, rng.uniform(-JIT_W, JIT_W, size=n))


def report(dq_by_k, sg16):
    print("[fig_e] TOP: fixed-dwell DeltaQ (rho 0.6 vs 0.05, nu=%.0f, %s)"
          % (NU_FIXED, ARM), flush=True)
    print("[fig_e] %-6s %-9s %-12s %-8s %-10s"
          % ("k", "source", "median_dB", "n_img", "n_pos(>0)"), flush=True)
    srcs = {512: "controls", 32: "robustness", 16: "summary.json", 1: "controls"}
    for k in K_LADDER:
        vals = np.array(list(dq_by_k[k].values()), dtype=float)
        print("[fig_e] %-6d %-9s %-12.4f %-8d %-10d"
              % (k, srcs[k], float(np.median(vals)), vals.size,
                 int(np.sum(vals > 0.0))), flush=True)
    sg = np.array(list(sg16.values()), dtype=float)
    n_zero = int(np.sum(sg == 0.0))
    print("[fig_e] BOTTOM: k=16 acquisition-speedup S_gate (study2_summary.json)",
          flush=True)
    print("[fig_e]   median S_gate = %.4f  (n=%d images; S>=3: %d; S>1: %d)"
          % (float(np.median(sg)), sg.size, int(np.sum(sg >= GATE_BAR)),
             int(np.sum(sg > 1.0))), flush=True)
    print("[fig_e]   censored S_gate = 0 (drawn OPEN at y=%.1f): %d image(s)"
          % (S_FLOOR, n_zero), flush=True)


def make_figure(dq_by_k, sg16):
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 8,
        "axes.linewidth": 0.6,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    fig, (axt, axb) = plt.subplots(2, 1, sharex=True, figsize=(3.4, 3.2),
                                   gridspec_kw={"height_ratios": [1.3, 1.0]})
    fig.subplots_adjust(left=0.16, right=0.965, top=0.955, bottom=0.17,
                        hspace=0.13)

    # ---------------- TOP: fixed-dwell DeltaQ ---------------- #
    axt.axhline(0.0, ls="--", lw=0.8, color="0.4", zorder=1)
    axt.axhline(PREREG_BAR, ls=":", lw=1.0, color=C_ORANGE, zorder=1)
    axt.text(0.025, 0.955, "preregistered bar (k=16 gate)",
             transform=axt.transAxes, va="top", ha="left",
             fontsize=6.0, color=C_ORANGE,
             bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none",
                       alpha=0.75))
    for i, k in enumerate(K_LADDER):
        vals = np.array(list(dq_by_k[k].values()), dtype=float)
        xj = _jitter(k, vals.size, seed=1000 + i)
        axt.scatter(xj, vals, s=7, color=C_DOT, alpha=0.7, linewidths=0,
                    zorder=2)
        med = float(np.median(vals))
        axt.scatter([k], [med], s=52, marker="D", color=C_BLUE,
                    edgecolors="black", linewidths=0.7, zorder=4)
    axt.set_ylabel(r"fixed-dwell $\Delta$PSNR$_\mathrm{rad}$ (dB)", fontsize=7,
                   labelpad=2)
    axt.tick_params(axis="y", labelsize=7)

    # ---------------- BOTTOM: k=16-only gate S ---------------- #
    axb.set_yscale("log")
    axb.axhline(GATE_BAR, ls="--", lw=0.8, color="0.4", zorder=1)
    axb.text(0.70, GATE_BAR, "gate bar (S=3)",
             transform=axb.get_yaxis_transform(), va="bottom", ha="center",
             fontsize=6.0, color="0.35", zorder=3)
    # non-gate rungs: light grey band + 'no gate (frozen)'
    for k in K_LADDER:
        if k == K16:
            continue
        lo, hi = k * 2.0 ** (-JIT_W - 0.06), k * 2.0 ** (JIT_W + 0.06)
        axb.axvspan(lo, hi, color="0.90", zorder=0)
        axb.text(k, 0.13, "no gate (frozen)", rotation=90, ha="center",
                 va="bottom", fontsize=5.8, color="0.5", zorder=1)
    # k=16 gate points
    sg = np.array(list(sg16.values()), dtype=float)
    pos = sg[sg > 0.0]
    zero = sg[sg == 0.0]
    xj = _jitter(K16, pos.size, seed=2016)
    axb.scatter(xj, pos, s=8, color=C_DOT, alpha=0.75, linewidths=0, zorder=2)
    if zero.size:
        xj0 = _jitter(K16, zero.size, seed=3016)
        axb.scatter(xj0, np.full(zero.size, S_FLOOR), s=14, facecolors="none",
                    edgecolors=C_DOT, linewidths=0.8, zorder=2)
    med_s = float(np.median(sg))     # median over ALL 24 images (censored incl.)
    axb.scatter([K16], [med_s], s=52, marker="D", color=C_ORANGE,
                edgecolors="black", linewidths=0.7, zorder=4)
    axb.set_ylabel(r"acquisition speedup $S$", fontsize=7, labelpad=2)
    axb.set_ylim(0.08, 60.0)
    axb.tick_params(axis="y", labelsize=7)

    # ---------------- shared log2 occupancy x-axis ---------------- #
    axb.set_xscale("log", base=2)
    axb.set_xticks(K_LADDER)
    axb.set_xticklabels([K_LABELS[k] for k in K_LADDER], rotation=30,
                        rotation_mode="anchor", ha="right", va="top")
    axb.set_xlim(700.0, 0.72)        # inverted: 512 (50%) at left, 1 (0.1%) right
    axb.tick_params(axis="x", labelsize=6.0, pad=2)
    axb.minorticks_off()
    axb.set_xlabel(r"pattern occupancy $k$ ($k/%d$)" % N_PIX)

    os.makedirs(OUT_DIR, exist_ok=True)
    fig.savefig(OUT_PDF)
    fig.savefig(OUT_PNG, dpi=200)
    plt.close(fig)
    print("[fig_e] wrote %s" % OUT_PDF, flush=True)
    print("[fig_e] wrote %s (200 dpi)" % OUT_PNG, flush=True)


def main():
    dq_by_k = {
        512: deltaq_from_csv(CONTROLS_CSV, 512),
        32:  deltaq_from_csv(ROBUST_CSV, 32),
        1:   deltaq_from_csv(CONTROLS_CSV, 1),
    }
    dq16, sg16 = load_k16()
    dq_by_k[16] = dq16
    report(dq_by_k, sg16)
    make_figure(dq_by_k, sg16)
    return 0


if __name__ == "__main__":
    sys.exit(main())
