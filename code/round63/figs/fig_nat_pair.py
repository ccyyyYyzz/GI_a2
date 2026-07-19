"""Study-1 figure — natural-subject pair (stl_02, the NATURAL-24 regime-boundary
subject). A 1x3 row: truth, safe/full RQL (rho=0.05, nu=2000), and fast RQL
(rho=0.6, nu=2000).

Both reconstructions go through the FROZEN production path exactly as
code/round63/figs/render_recons.py does: campaign.run_cell verbatim with
campaign.MET.main_metrics monkey-patched to capture the raw reconstruction.
Display range: vmin 0, vmax = 1.2 * truth.max() (frozen rule). Each cell is
titled with its arm + PSNR_rad (radiometric, non-rescaled;
10 log10(truth.max()^2 / mean((clip(xhat,0)-truth)^2))).

Run (cwd = repo root D:\\GI_another):
    D:\\Anacondar\\anaconda3\\python.exe code\\round63\\figs\\fig_nat_pair.py
"""
import os
import sys
import time

import numpy as np

os.chdir(r"D:\GI_another")
sys.path.insert(0, r"D:\GI_another\code\round63")
import campaign  # noqa: E402

CAP = []
_orig = campaign.MET.main_metrics


def _capture(xh, x, side, **kw):
    CAP.append((np.array(xh, copy=True), np.array(x, copy=True)))
    return _orig(xh, x, side, **kw)


campaign.MET.main_metrics = _capture

SIDE = 64
M = 4096
PATTERN = "bern50"
SEED = 0
IMAGE = "stl_02"
IMAGESET = "conf"

# (arm, rho, nu, label)
COLUMNS = [
    ("RQL", 0.05, 2000, "safe/full"),
    ("RQL", 0.6,  2000, "RQL fast"),
]

OUT_DIR = os.path.join(r"D:\GI_another", "paper", "figs")
OUT_PNG = os.path.join(OUT_DIR, "fig_nat_pair.png")
OUT_PDF = os.path.join(OUT_DIR, "fig_nat_pair.pdf")


def recon(arm, rho, nu):
    CAP.clear()
    cell = {"side": SIDE, "nu": float(nu), "rho_bar": float(rho), "seed": SEED,
            "M": M, "pattern": PATTERN, "arms": [arm],
            "imageset": IMAGESET, "images": [IMAGE],
            "cell_id": f"fignat_{IMAGE}_{arm}_{rho}_{nu}", "audit": False}
    campaign.run_cell(cell)
    xh, x = CAP[-1]
    return np.maximum(xh, 0.0).reshape(SIDE, SIDE), x.reshape(SIDE, SIDE)


def psnr_rad(xhat, truth):
    mse = float(np.mean((xhat - truth) ** 2))
    if mse <= 0:
        return float("inf")
    return 10.0 * np.log10(float(truth.max()) ** 2 / mse)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 8,
        "axes.linewidth": 0.6,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    truth = None
    cells = []
    t0 = time.time()
    for i, (arm, rho, nu, label) in enumerate(COLUMNS, start=1):
        ts = time.time()
        xhat, x = recon(arm, rho, nu)
        truth = x
        p = psnr_rad(xhat, truth)
        cells.append((label, xhat, p))
        print("[fig_nat] %d/%d %-9s rho=%.2f nu=%-4d PSNR_rad=%6.2f dB (%.1fs)"
              % (i, len(COLUMNS), arm, rho, nu, p, time.time() - ts), flush=True)
    print("[fig_nat] recons done in %.1fs" % (time.time() - t0), flush=True)

    vmax = 1.2 * float(truth.max())
    # 4th panel: intensity profile along the frozen center row (y=32,
    # stated a priori) — makes the RADIOMETRIC difference visible where the
    # grayscale rendering cannot (prior-limited cohort; structure recovery
    # is not claimed on naturals, absolute flux accuracy is).
    fig, axes = plt.subplots(1, 4, figsize=(6.4, 1.85),
                             gridspec_kw={"width_ratios": [1, 1, 1, 1.25]})

    def show(ax, im, title):
        ax.imshow(im, cmap="gray", vmin=0.0, vmax=vmax, interpolation="nearest")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(title, fontsize=8, pad=3)

    show(axes[0], truth, "truth")
    for c, (label, xhat, p) in enumerate(cells, start=1):
        show(axes[c], xhat, "%s\n%.1f dB" % (label, p))

    ROW = 32
    side = int(np.sqrt(truth.size))
    ax = axes[3]
    ax.plot(truth.reshape(side, side)[ROW], color="k", lw=0.9, label="truth")
    ax.plot(np.maximum(cells[0][1], 0.0).reshape(side, side)[ROW],
            color="#0072B2", lw=0.9, label="safe/full")
    ax.plot(np.maximum(cells[1][1], 0.0).reshape(side, side)[ROW],
            color="#D55E00", lw=0.9, label="RQL fast")
    ax.set_title("row %d profile" % ROW, fontsize=8, pad=3)
    ax.set_xticks([])
    ax.yaxis.tick_right()
    ax.tick_params(labelsize=6, length=2)
    leg = ax.legend(fontsize=5.5, loc="lower left", framealpha=0.85)
    leg.get_frame().set_linewidth(0.0)

    fig.subplots_adjust(left=0.01, right=0.99, top=0.86, bottom=0.06,
                        wspace=0.12)
    os.makedirs(OUT_DIR, exist_ok=True)
    fig.savefig(OUT_PDF)
    fig.savefig(OUT_PNG, dpi=200)
    plt.close(fig)
    print("[fig_nat] wrote %s" % OUT_PDF, flush=True)
    print("[fig_nat] wrote %s (200 dpi)" % OUT_PNG, flush=True)

    print("\n[fig_nat] === PSNR_rad (dB), stl_02, vmax=%.3g ===" % vmax,
          flush=True)
    print("[fig_nat]     %-10s %7.3f dB" % ("truth", float("nan")), flush=True)
    for (label, _xh, p) in cells:
        print("[fig_nat]     %-10s %7.3f dB" % (label, p), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
