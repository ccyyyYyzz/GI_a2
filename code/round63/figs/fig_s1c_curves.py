"""Study-1 figure, panel (c) — RQL radiometric PSNR versus optical acquisition
time, DETAIL-24 cohort (frozen prereg docs/ROUND63_FIGURE_PREREG.md panel c:
"the DETAIL-24 cohort seed-mean curves (RQL, both rho) ... per-family thin lines
+ cohort median band. No subject selection freedom remains (all 24 shown).")

Source: results/round63_s2_detail/pilot_rows.csv (S2A_DETAIL merged rows).
For the RQL arm at both mean-load points (safe rho=0.05 blue, fast rho=0.6
orange), each of the 24 DETAIL images gets a seed-mean PSNR_rad curve versus the
OPTICAL TIME  T_opt = M * nu * tau  (tau = 50 ns, M = 4096) on a log x-axis.
Thin per-image lines (alpha 0.25) + a thick cohort-median line with a shaded
inter-quartile band per rho.

PROVENANCE: PSNR_rad read verbatim; '' / 'nan' cells dropped, never imputed.

Run (cwd = repo root D:\\GI_another):
    D:\\Anacondar\\anaconda3\\python.exe code\\round63\\figs\\fig_s1c_curves.py
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
OUT_PNG = os.path.join(OUT_DIR, "fig_s1c_curves.png")
OUT_PDF = os.path.join(OUT_DIR, "fig_s1c_curves.pdf")

ARM = "RQL"
TAU = 50e-9
M = 4096

C_BLUE = "#0072B2"      # safe rho = 0.05
C_ORANGE = "#D55E00"    # fast rho = 0.6
SERIES = [(0.05, "safe", C_BLUE), (0.6, "fast", C_ORANGE)]


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


def build(rows, rho):
    """Return (nus_sorted, {image: {nu: seed_mean_PSNR}}) for the RQL arm at rho.
    Only DETAIL images (name starts 'detail_') are used; all 24 are kept."""
    acc = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in rows:
        if r["arm"] != ARM:
            continue
        if not r["image"].startswith("detail_"):
            continue
        if _to_float(r["rho_bar"]) != rho:
            continue
        pv = _to_float(r["PSNR_rad"])
        if pv is None:
            continue
        acc[r["image"]][_to_float(r["nu"])].append(pv)
    nus = sorted({nu for img in acc.values() for nu in img})
    per_img = {}
    for img, byn in acc.items():
        per_img[img] = {nu: float(np.mean(byn[nu])) for nu in byn}
    return nus, per_img


def main():
    rows = load_rows(DATA)
    print("[fig_s1c] loaded %d rows from %s" % (len(rows), DATA), flush=True)

    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 8,
        "axes.linewidth": 0.6,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    fig, ax = plt.subplots(figsize=(3.4, 2.7))

    for rho, tag, color in SERIES:
        nus, per_img = build(rows, rho)
        t_opt = np.array([M * nu * TAU for nu in nus])
        print("[fig_s1c] rho=%.2f (%s): %d images, nu grid %s"
              % (rho, tag, len(per_img), nus), flush=True)
        # per-image thin lines
        mat = np.full((len(per_img), len(nus)), np.nan)
        for i, (img, byn) in enumerate(sorted(per_img.items())):
            y = np.array([byn.get(nu, np.nan) for nu in nus])
            mat[i] = y
            ax.plot(t_opt, y, color=color, lw=0.5, alpha=0.25, zorder=1)
        # cohort median + IQR band
        med = np.nanmedian(mat, axis=0)
        q1 = np.nanpercentile(mat, 25, axis=0)
        q3 = np.nanpercentile(mat, 75, axis=0)
        ax.fill_between(t_opt, q1, q3, color=color, alpha=0.18, lw=0, zorder=2)
        ax.plot(t_opt, med, color=color, lw=2.0, zorder=4)
        print("[fig_s1c]   cohort median PSNR_rad by T_opt:", flush=True)
        for nu, to, mv in zip(nus, t_opt, med):
            print("[fig_s1c]     nu=%-6.0f T_opt=%.4g s  median=%6.3f dB"
                  % (nu, to, mv), flush=True)

    ax.set_xscale("log")
    ax.set_xlabel(r"optical time $T_\mathrm{opt}=M\,\nu\,\tau$ (s)")
    ax.set_ylabel(r"PSNR$_\mathrm{rad}$ (dB)")

    handles = [
        Line2D([0], [0], color=C_BLUE, lw=2.0),
        Line2D([0], [0], color=C_ORANGE, lw=2.0),
        Line2D([0], [0], color="0.4", lw=0.6, alpha=0.5),
    ]
    labels = [r"safe $\bar\rho=0.05$ (median)",
              r"fast $\bar\rho=0.6$ (median)",
              "per-image (24)"]
    ax.legend(handles, labels, loc="lower right", fontsize=6.0, frameon=True,
              framealpha=0.9, borderpad=0.4, handlelength=1.8)

    fig.tight_layout()
    os.makedirs(OUT_DIR, exist_ok=True)
    fig.savefig(OUT_PDF)
    fig.savefig(OUT_PNG, dpi=200)
    plt.close(fig)
    print("[fig_s1c] wrote %s" % OUT_PDF, flush=True)
    print("[fig_s1c] wrote %s (200 dpi)" % OUT_PNG, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
