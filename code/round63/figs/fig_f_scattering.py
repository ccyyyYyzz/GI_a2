"""Optics-Express figure, panel (f) — dynamic-scattering result reported in the
Hao-compatible detected-photon flux unit S_det.

Reads the frozen Study-2 dynamic-scattering table
(results/round63_study2/scattering_rows.csv; k = 16, CV(alpha) = 0.20) and plots
radiometric PSNR against the MEASURED per-pattern-per-pixel detected-photon flux
S_det (the Hao 2025 convention), for two mean-load arms:

  * safe       : rho_bar = 0.05
  * high-flux  : rho_bar = 0.6

Each plotted point is one (rho_bar, nu) cell: x = mean S_det over the cell,
y = mean PSNR_rad over the cell, with +/-1 sd error bars over the 18 (image,seed)
pairs in the cell. Points are annotated with their nu value (20, 200, 2000). A
light vertical band marks the Hao et al. reported experimental range
S_det in [0.096, 0.44].

PROVENANCE — values are read verbatim from the committed campaign table; no
imputation is performed. Empty / 'nan' cells are dropped from every mean and sd
(never imputed); the scattering table has none for the plotted columns.

Run (cwd = repo root D:\\GI_another):
    D:\\Anacondar\\anaconda3\\python.exe code\\round63\\figs\\fig_f_scattering.py
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

DATA = os.path.join(REPO, "results", "round63_study2", "scattering_rows.csv")
OUT_DIR = os.path.join(REPO, "paper", "figs")
OUT_PNG = os.path.join(OUT_DIR, "fig_f_scattering.png")
OUT_PDF = os.path.join(OUT_DIR, "fig_f_scattering.pdf")

NU_LADDER = [20.0, 200.0, 2000.0]
HAO_LO, HAO_HI = 0.096, 0.44               # Hao et al. reported experimental range

# (rho_bar, legend label, colour, marker, annotation offset for nu labels)
SERIES = [
    (0.05, r"$\bar{\rho}=0.05$ (safe)",      "#0072B2", "o", (0, 8)),
    (0.6,  r"$\bar{\rho}=0.6$ (high-flux)",  "#D55E00", "s", (0, -13)),
]


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


def compute_series(rows, rho):
    """Per (rho, nu) cell: mean S_det, mean PSNR_rad, +/-1 sd PSNR_rad over pairs."""
    out = {}
    for nu in NU_LADDER:
        sdet, psnr = [], []
        for r in rows:
            if _to_float(r["rho_bar"]) != rho or _to_float(r["nu"]) != nu:
                continue
            sv, pv = _to_float(r["S_det"]), _to_float(r["PSNR_rad"])
            if sv is not None:
                sdet.append(sv)
            if pv is not None:
                psnr.append(pv)
        out[nu] = dict(s_det=float(np.mean(sdet)), psnr=float(np.mean(psnr)),
                       psnr_sd=_sd(psnr), n=len(psnr))
    return out


def report(data_by_rho):
    print("[fig_f] data = %s" % DATA, flush=True)
    print("[fig_f] x = mean S_det per (rho,nu) cell (Hao detected-photon unit); "
          "y = mean PSNR_rad; +/-1 sd over image,seed pairs", flush=True)
    print("[fig_f] Hao et al. experimental range: S_det in [%.3f, %.3f]"
          % (HAO_LO, HAO_HI), flush=True)
    print("[fig_f] %-8s %-8s %-12s %-12s %-12s %-6s"
          % ("rho_bar", "nu", "S_det_mean", "PSNR_mean", "PSNR_sd", "n"), flush=True)
    for rho, _lab, _c, _m, _o in SERIES:
        for nu in NU_LADDER:
            d = data_by_rho[rho][nu]
            print("[fig_f] %-8.2f %-8.0f %-12.6f %-12.4f %-12.4f %-6d"
                  % (rho, nu, d["s_det"], d["psnr"], d["psnr_sd"], d["n"]), flush=True)


def make_figure(data_by_rho):
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 8,
        "axes.linewidth": 0.6,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    fig, ax = plt.subplots(figsize=(3.4, 2.6))

    # Hao experimental range band (drawn first, behind the data)
    ax.axvspan(HAO_LO, HAO_HI, color="0.82", alpha=0.6, lw=0, zorder=0)
    ax.text(np.sqrt(HAO_LO * HAO_HI), 0.30, "Hao et al. exp. range",
            transform=ax.get_xaxis_transform(), rotation=90, va="center",
            ha="center", fontsize=6.0, color="0.35", zorder=1)

    for rho, label, color, marker, (dx, dy) in SERIES:
        d = data_by_rho[rho]
        xs = np.array([d[nu]["s_det"] for nu in NU_LADDER])
        ys = np.array([d[nu]["psnr"] for nu in NU_LADDER])
        es = np.array([d[nu]["psnr_sd"] for nu in NU_LADDER])
        ax.errorbar(xs, ys, yerr=es, fmt=marker + "-", color=color, ms=4, lw=1.0,
                    capsize=2, elinewidth=0.6, markeredgewidth=0.6, label=label,
                    zorder=3)
        for nu, xv, yv in zip(NU_LADDER, xs, ys):
            ax.annotate(r"$\nu{=}%d$" % int(nu), (xv, yv),
                        textcoords="offset points", xytext=(dx, dy),
                        ha="center", fontsize=6.0, color=color, zorder=4)

    ax.set_xscale("log")
    ax.set_xlim(5e-4, 1.3)
    ax.set_ylim(3.0, 17.0)
    ax.set_xlabel(r"measured $S_\mathrm{det}$ (photons pattern$^{-1}$ pixel$^{-1}$)")
    ax.set_ylabel(r"PSNR$_\mathrm{rad}$ (dB)")

    ax.legend(loc="upper left", fontsize=6.5, frameon=True, framealpha=0.9,
              borderpad=0.4, handlelength=1.8)
    ax.text(0.03, 0.04, r"dynamic scattering, CV($\alpha$)$=0.20$",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=6.5,
            color="0.25",
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.7", lw=0.5))

    fig.tight_layout()
    os.makedirs(OUT_DIR, exist_ok=True)
    fig.savefig(OUT_PDF)
    fig.savefig(OUT_PNG, dpi=200)
    plt.close(fig)
    print("[fig_f] wrote %s" % OUT_PDF, flush=True)
    print("[fig_f] wrote %s (200 dpi)" % OUT_PNG, flush=True)


def main():
    rows = load_rows(DATA)
    print("[fig_f] loaded %d scattering rows" % len(rows), flush=True)
    data_by_rho = {rho: compute_series(rows, rho) for rho, _l, _c, _m, _o in SERIES}
    report(data_by_rho)
    make_figure(data_by_rho)
    return 0


if __name__ == "__main__":
    sys.exit(main())
