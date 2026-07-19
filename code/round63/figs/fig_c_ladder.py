"""Optics-Express figure, panel (c) — actual bucket-energy contrast C_u and the
photon-limited flux figure Gamma versus pattern occupancy k.

Reads the frozen Study-2 fluxmap table (results/round63_study2/fluxmap_rows.csv)
and plots, on a shared log2 occupancy axis k in {1, 16, 32, 512}:

  * LEFT y-axis (log)  : measured bucket-energy contrast C_u, averaged over ALL
    fluxmap rows at each k. C_u is a property of the pattern realisation
    (image, seed, k) and is invariant to (rho_bar, nu), so the row-mean equals
    the mean over the 18 (image, seed) pairs. Error bars are +/-1 sd over those
    (image, seed) pairs.
  * RIGHT y-axis (log) : the photon-limited flux figure Gamma at the operating
    point (rho_bar = 0.6, nu = 2000), same (image, seed) treatment. A dashed
    line marks the photon-limited threshold Gamma = 1.

PROVENANCE — values are read verbatim from the committed campaign table; no
imputation is performed. Empty / 'nan' cells are dropped from every mean and sd
(never imputed); the fluxmap table has none for the plotted columns.

Run (cwd = repo root D:\\GI_another):
    D:\\Anacondar\\anaconda3\\python.exe code\\round63\\figs\\fig_c_ladder.py
"""
import os
import sys
import csv

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))          # code/round63/figs
R63 = os.path.dirname(HERE)                                # code/round63
REPO = os.path.dirname(os.path.dirname(R63))               # repo root (D:\GI_another)

DATA = os.path.join(REPO, "results", "round63_study2", "fluxmap_rows.csv")
OUT_DIR = os.path.join(REPO, "paper", "figs")
OUT_PNG = os.path.join(OUT_DIR, "fig_c_ladder.png")
OUT_PDF = os.path.join(OUT_DIR, "fig_c_ladder.pdf")

N_PIX = 1024                       # scene pixels n = 32*32 (occupancy = k / n)
K_LADDER = [1, 16, 32, 512]        # log2 occupancy ladder (ascending)
K_LABELS = {                       # exactly as ruled in the manuscript
    1:   "1 (0.1%)",
    16:  "16 (1.6%)",
    32:  "32 (3.1%)",
    512: "512 (50%)",
}
GAMMA_RHO = 0.6                    # operating point for the Gamma series
GAMMA_NU = 2000.0

C_BLUE = "#0072B2"                 # colourblind-safe pair
C_ORANGE = "#D55E00"


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


def _sd(values):
    """Sample standard deviation (ddof=1); nan for < 2 values."""
    a = np.asarray(values, dtype=float)
    return float(np.std(a, ddof=1)) if a.size >= 2 else float("nan")


def compute_cu(rows):
    """C_u central value = mean over ALL rows at k; error = sd over (image,seed) pairs."""
    out = {}
    for k in K_LADDER:
        allvals = []
        pair = {}                                  # (image,seed) -> [C_u ...]
        for r in rows:
            if _to_float(r["k"]) != float(k):
                continue
            v = _to_float(r["C_u"])
            if v is None:
                continue
            allvals.append(v)
            pair.setdefault((r["image"], r["seed"]), []).append(v)
        pair_means = [float(np.mean(v)) for v in pair.values()]
        out[k] = dict(mean=float(np.mean(allvals)), sd=_sd(pair_means),
                      n_rows=len(allvals), n_pairs=len(pair_means),
                      pair_mean=float(np.mean(pair_means)))
    return out


def compute_gamma(rows):
    """Gamma at (rho_bar=GAMMA_RHO, nu=GAMMA_NU); mean and sd over (image,seed) pairs."""
    out = {}
    for k in K_LADDER:
        pair = {}
        for r in rows:
            if _to_float(r["k"]) != float(k):
                continue
            if _to_float(r["rho_bar"]) != GAMMA_RHO or _to_float(r["nu"]) != GAMMA_NU:
                continue
            v = _to_float(r["Gamma"])
            if v is None:
                continue
            pair.setdefault((r["image"], r["seed"]), []).append(v)
        pair_means = [float(np.mean(v)) for v in pair.values()]
        out[k] = dict(mean=float(np.mean(pair_means)), sd=_sd(pair_means),
                      n_pairs=len(pair_means))
    return out


def report(cu, gm):
    print("[fig_c] data = %s" % DATA, flush=True)
    print("[fig_c] LEFT axis : bucket-energy contrast C_u (mean over ALL fluxmap "
          "rows at k; +/-1 sd over image,seed pairs)", flush=True)
    print("[fig_c] %-6s %-12s %-12s %-8s %-8s %-12s"
          % ("k", "C_u_mean", "C_u_sd", "n_rows", "n_pairs", "pair_mean(chk)"), flush=True)
    for k in K_LADDER:
        c = cu[k]
        print("[fig_c] %-6d %-12.6f %-12.6f %-8d %-8d %-12.6f"
              % (k, c["mean"], c["sd"], c["n_rows"], c["n_pairs"], c["pair_mean"]),
              flush=True)
    print("[fig_c] RIGHT axis: Gamma at rho_bar=%.2f, nu=%.0f (mean & +/-1 sd over "
          "image,seed pairs)" % (GAMMA_RHO, GAMMA_NU), flush=True)
    print("[fig_c] %-6s %-12s %-12s %-8s" % ("k", "Gamma_mean", "Gamma_sd", "n_pairs"),
          flush=True)
    for k in K_LADDER:
        g = gm[k]
        print("[fig_c] %-6d %-12.6f %-12.6f %-8d"
              % (k, g["mean"], g["sd"], g["n_pairs"]), flush=True)


def make_figure(cu, gm):
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "mathtext.fontset": "dejavusans",
        "font.size": 8,
        "axes.linewidth": 0.6,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    x = np.array(K_LADDER, dtype=float)
    cu_m = np.array([cu[k]["mean"] for k in K_LADDER])
    cu_s = np.array([cu[k]["sd"] for k in K_LADDER])
    gm_m = np.array([gm[k]["mean"] for k in K_LADDER])
    gm_s = np.array([gm[k]["sd"] for k in K_LADDER])

    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    ax2 = ax.twinx()

    l1 = ax.errorbar(x, cu_m, yerr=cu_s, fmt="o-", color=C_BLUE, ms=4, lw=1.0,
                     capsize=2, elinewidth=0.6, markeredgewidth=0.6,
                     label=r"$C_u$ (all $\bar{\rho},\nu$)")
    l2 = ax2.errorbar(x, gm_m, yerr=gm_s, fmt="s-", color=C_ORANGE, ms=4, lw=1.0,
                      capsize=2, elinewidth=0.6, markeredgewidth=0.6,
                      label=r"$\Gamma$ ($\bar{\rho}{=}0.6,\ \nu{=}2000$)")

    # photon-limited threshold on the Gamma (right) axis
    ax2.axhline(1.0, ls="--", lw=0.8, color=C_ORANGE, alpha=0.7, zorder=0)
    ax2.text(0.03, 1.0, r"photon-limited threshold $\Gamma=1$",
             transform=ax2.get_yaxis_transform(), va="bottom", ha="left",
             fontsize=6.5, color=C_ORANGE)

    # x-axis: log2 occupancy ladder with exact ticks + occupancy-percent labels
    ax.set_xscale("log", base=2)
    ax.set_xticks(K_LADDER)
    ax.set_xticklabels([K_LABELS[k] for k in K_LADDER], rotation=30,
                       rotation_mode="anchor", ha="right", va="top")
    ax.set_xlim(0.7, 730)
    ax.tick_params(axis="x", labelsize=6.5, pad=2)
    ax.minorticks_off()            # suppress log2 minor ticks between the 4 rungs

    ax.set_yscale("log")
    ax2.set_yscale("log")
    ax.set_ylim(0.02, 3.0)
    ax2.set_ylim(0.6, 60.0)

    ax.set_xlabel(r"pattern occupancy $k$ ($k/1024$)")
    ax.set_ylabel(r"measured contrast $C_u$", color=C_BLUE)
    ax2.set_ylabel(r"flux figure $\Gamma$", color=C_ORANGE)
    ax.tick_params(axis="y", labelcolor=C_BLUE)
    ax2.tick_params(axis="y", labelcolor=C_ORANGE)

    handles = [l1, l2]
    labels = [h.get_label() for h in handles]
    ax.legend(handles, labels, loc="upper right", fontsize=6.5, frameon=True,
              framealpha=0.9, borderpad=0.4, handlelength=1.8)

    fig.tight_layout()
    os.makedirs(OUT_DIR, exist_ok=True)
    fig.savefig(OUT_PDF)
    fig.savefig(OUT_PNG, dpi=200)
    plt.close(fig)
    print("[fig_c] wrote %s" % OUT_PDF, flush=True)
    print("[fig_c] wrote %s (200 dpi)" % OUT_PNG, flush=True)


def main():
    rows = load_rows(DATA)
    print("[fig_c] loaded %d fluxmap rows" % len(rows), flush=True)
    cu = compute_cu(rows)
    gm = compute_gamma(rows)
    report(cu, gm)
    make_figure(cu, gm)
    return 0


if __name__ == "__main__":
    sys.exit(main())
