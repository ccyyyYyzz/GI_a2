"""Optics-Express figure, panel (d) — safe/fast radiometric PSNR and CNR versus
dwell budget nu for the frozen resolution-target subject.

Reads the frozen Study-2 primary table (results/round63_study2/primary_rows.csv;
k = 16, RQL arm) and plots, for the single prereg-frozen subject
`detail32_spokes_0` (round-8 addendum, docs/ROUND63_FIGURE_PREREG_ADDENDUM_LADDER.md
- the radial resolution target chosen a priori by index convention), against the
per-pattern dwell budget nu on a log axis:

  * LEFT y-axis  : seed-mean radiometric PSNR_rad (solid lines + circles), with
    +/-1 sd error bars over the five confirmatory measurement seeds.
  * RIGHT twin y : seed-mean CNR (dashed lines + squares, alpha 0.7), same colours.

Two mean-load arms are shown:
  * safe : rho_bar = 0.05 (blue  #0072B2)
  * fast : rho_bar = 0.6  (orange #D55E00)

NO peak-location claim is drawn: per round-8 §12 / the addendum, the panel must NOT
assert the image optimum coincides with the scalar Fisher count-information ridge,
so no rho* / ridge marker is placed.

PROVENANCE — values are read verbatim from the committed campaign table; no
imputation is performed. Empty / 'nan' CNR cells are dropped from the mean and sd
(never imputed), exactly as the addendum rules ("CNR cells recorded NA are omitted
from the mean without substitution").

Run (cwd = repo root D:\\GI_another):
    D:\\Anacondar\\anaconda3\\python.exe code\\round63\\figs\\fig_d_dwell.py
"""
import os
import sys
import csv

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))          # code/round63/figs
R63 = os.path.dirname(HERE)                                # code/round63
REPO = os.path.dirname(os.path.dirname(R63))               # repo root (D:\GI_another)

DATA = os.path.join(REPO, "results", "round63_study2", "primary_rows.csv")
OUT_DIR = os.path.join(REPO, "paper", "figs")
OUT_PNG = os.path.join(OUT_DIR, "fig_d_dwell.png")
OUT_PDF = os.path.join(OUT_DIR, "fig_d_dwell.pdf")

SUBJECT = "detail32_spokes_0"       # frozen prereg subject (addendum, panel d)
ARM = "RQL"
K_PRIMARY = 16
NU_LADDER = [5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0, 1000.0, 2000.0]

C_BLUE = "#0072B2"                  # colourblind-safe pair (safe arm)
C_ORANGE = "#D55E00"               # fast arm
CNR_ALPHA = 0.7                    # dashed CNR series drawn lighter

# (rho_bar, tag, colour)
SERIES = [
    (0.05, "safe", C_BLUE),
    (0.6,  "fast", C_ORANGE),
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
    """Per nu cell: seed-mean PSNR_rad (+/-1 sd) and seed-mean CNR (+/-1 sd).

    PSNR uses all seed rows; CNR drops empty/'nan' cells before averaging.
    """
    out = {}
    for nu in NU_LADDER:
        psnr, cnr = [], []
        for r in rows:
            if r["image"] != SUBJECT or r["arm"] != ARM:
                continue
            if _to_float(r["k"]) != float(K_PRIMARY):
                continue
            if _to_float(r["rho_bar"]) != rho or _to_float(r["nu"]) != nu:
                continue
            pv = _to_float(r["PSNR_rad"])
            cv = _to_float(r["cnr"])
            if pv is not None:
                psnr.append(pv)
            if cv is not None:
                cnr.append(cv)
        out[nu] = dict(
            psnr=float(np.mean(psnr)) if psnr else float("nan"),
            psnr_sd=_sd(psnr), n_psnr=len(psnr),
            cnr=float(np.mean(cnr)) if cnr else float("nan"),
            cnr_sd=_sd(cnr), n_cnr=len(cnr),
        )
    return out


def report(data_by_rho):
    print("[fig_d] data = %s" % DATA, flush=True)
    print("[fig_d] subject = %s (frozen addendum), arm = %s, k = %d"
          % (SUBJECT, ARM, K_PRIMARY), flush=True)
    print("[fig_d] LEFT: seed-mean PSNR_rad (+/-1 sd over 5 seeds); "
          "RIGHT: seed-mean CNR (+/-1 sd; NA cells dropped, never imputed)",
          flush=True)
    for rho, tag, _c in SERIES:
        print("[fig_d] --- rho_bar = %.2f (%s) ---" % (rho, tag), flush=True)
        print("[fig_d] %-8s %-6s %-11s %-11s %-6s %-11s %-11s"
              % ("nu", "nP", "PSNR_mean", "PSNR_sd", "nCNR", "CNR_mean", "CNR_sd"),
              flush=True)
        for nu in NU_LADDER:
            d = data_by_rho[rho][nu]
            print("[fig_d] %-8.0f %-6d %-11.4f %-11.4f %-6d %-11.4f %-11.4f"
                  % (nu, d["n_psnr"], d["psnr"], d["psnr_sd"],
                     d["n_cnr"], d["cnr"], d["cnr_sd"]), flush=True)


def make_figure(data_by_rho):
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "mathtext.fontset": "dejavusans",
        "font.size": 8,
        "axes.linewidth": 0.6,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    x = np.array(NU_LADDER, dtype=float)
    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    ax2 = ax.twinx()

    for rho, tag, color in SERIES:
        d = data_by_rho[rho]
        ps = np.array([d[nu]["psnr"] for nu in NU_LADDER])
        pe = np.array([d[nu]["psnr_sd"] for nu in NU_LADDER])
        cs = np.array([d[nu]["cnr"] for nu in NU_LADDER])
        ce = np.array([d[nu]["cnr_sd"] for nu in NU_LADDER])
        # LEFT axis: PSNR_rad, solid + circles
        ax.errorbar(x, ps, yerr=pe, fmt="o-", color=color, ms=4, lw=1.0,
                    capsize=2, elinewidth=0.6, markeredgewidth=0.6, zorder=3)
        # RIGHT axis: CNR, dashed + squares, lighter
        ax2.errorbar(x, cs, yerr=ce, fmt="s--", color=color, ms=3.5, lw=1.0,
                     alpha=CNR_ALPHA, capsize=2, elinewidth=0.6,
                     markeredgewidth=0.6, zorder=2)

    ax.set_xscale("log")
    ax.set_xlim(4.0, 2600.0)
    ax.set_xlabel(r"dwell budget $\nu$ (patterns)")
    ax.set_ylabel(r"PSNR$_\mathrm{rad}$ (dB)")
    ax2.set_ylabel("CNR")

    # headroom so the compact 4-entry legend clears the data
    y0, y1 = ax.get_ylim()
    ax.set_ylim(y0, y1 + 0.30 * (y1 - y0))

    # compact 4-series legend (two axes -> proxy handles), 2 columns
    handles = [
        Line2D([0], [0], color=C_BLUE, marker="o", ls="-", ms=4, lw=1.0),
        Line2D([0], [0], color=C_ORANGE, marker="o", ls="-", ms=4, lw=1.0),
        Line2D([0], [0], color=C_BLUE, marker="s", ls="--", ms=3.5, lw=1.0,
               alpha=CNR_ALPHA),
        Line2D([0], [0], color=C_ORANGE, marker="s", ls="--", ms=3.5, lw=1.0,
               alpha=CNR_ALPHA),
    ]
    labels = ["PSNR safe", "PSNR fast", "CNR safe", "CNR fast"]
    ax.legend(handles, labels, loc="upper left", fontsize=6.0, ncol=2,
              frameon=True, framealpha=0.9, borderpad=0.4, handlelength=1.7,
              columnspacing=1.0, handletextpad=0.5)

    fig.tight_layout()
    os.makedirs(OUT_DIR, exist_ok=True)
    fig.savefig(OUT_PDF)
    fig.savefig(OUT_PNG, dpi=200)
    plt.close(fig)
    print("[fig_d] wrote %s" % OUT_PDF, flush=True)
    print("[fig_d] wrote %s (200 dpi)" % OUT_PNG, flush=True)


def main():
    rows = load_rows(DATA)
    print("[fig_d] loaded %d primary rows" % len(rows), flush=True)
    data_by_rho = {rho: compute_series(rows, rho) for rho, _t, _c in SERIES}
    report(data_by_rho)
    make_figure(data_by_rho)
    return 0


if __name__ == "__main__":
    sys.exit(main())
