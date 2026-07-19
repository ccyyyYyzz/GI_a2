"""Study-1 figure, panel (b) — reconstruction image array (frozen prereg
docs/ROUND63_FIGURE_PREREG.md, main figure panel b).

Two subject rows, six columns per row (truth + five operating points). Every
reconstruction goes through the FROZEN production path exactly as
code/round63/figs/render_recons.py does: campaign.run_cell verbatim, with
campaign.MET.main_metrics monkey-patched to capture each arm's raw
reconstruction. Nothing frozen is written; audit is disabled (display-only,
does not alter eta*/reconstruction — matches render_recons).

Subjects (frozen prereg):
  * detail_glyph_0  (imageset 'detail24', side 64) PLUS a zoom crop
    x in [2,42), y in [24,40) (the middle glyph row 'FBPEW'), shown as a
    dedicated crop row aligned under the full-frame glyph row.
  * stl_02          (imageset 'conf', side 64), full frame, no crop.

Columns (frozen prereg / task): truth; safe/full = (rho=0.05, nu=2000);
safe/short = (rho=0.05, nu=50); POISSON-LIN, PRECORRECT, RQL all at the
high-flux point (rho=0.6, nu=50). Seed 0, M=4096, pattern bern50.
The two safe columns are the flagship RQL arm at the safe operating point
(RQL is the production reconstructor; "safe" names the operating point, not a
different arm) — resolved toward the prereg column list, which names the
high-flux columns by arm and the safe columns by operating point.

Display range per row: vmin 0, vmax = 1.2 * truth.max() (frozen).
Each reconstruction cell is titled with its arm + PSNR_rad, where
PSNR_rad = 10 log10(truth.max()^2 / mean((clip(xhat,0)-truth)^2)) — the
radiometric, non-rescaled PSNR (spec D2 s4 primary metric), identical to the
campaign PSNR_rad definition.

Run (cwd = repo root D:\\GI_another):
    D:\\Anacondar\\anaconda3\\python.exe code\\round63\\figs\\fig_s1b_images.py
"""
import os
import sys
import time

import numpy as np

os.chdir(r"D:\GI_another")
sys.path.insert(0, r"D:\GI_another\code\round63")
import campaign  # noqa: E402

# ---- FROZEN production-path capture (verbatim render_recons.py pattern) ----
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

# crop for the glyph row (frozen): x in [2,42), y in [24,40)  -> img[24:40, 2:42]
CROP_Y = (24, 40)
CROP_X = (2, 42)

# (arm, rho, nu, column-label)
COLUMNS = [
    ("RQL",         0.05, 2000, "safe/full"),
    ("RQL",         0.05,   50, "safe/short"),
    ("POISSON-LIN", 0.6,    50, "POISSON-LIN"),
    ("PRECORRECT",  0.6,    50, "PRECORRECT"),
    ("RQL",         0.6,    50, "RQL"),
]

# (subject image, imageset, has_crop)
SUBJECTS = [
    ("detail_glyph_0", "detail24", True),
    ("stl_02",         "conf",     False),
]

OUT_DIR = os.path.join(r"D:\GI_another", "paper", "figs")
OUT_PNG = os.path.join(OUT_DIR, "fig_s1b_images.pdf".replace(".pdf", ".png"))
OUT_PDF = os.path.join(OUT_DIR, "fig_s1b_images.pdf")


def recon(image, imageset, arm, rho, nu):
    """Run one arm at one operating point through campaign.run_cell (frozen
    production path) and capture its raw reconstruction. Returns
    (xhat_clipped_64x64, truth_64x64)."""
    CAP.clear()
    cell = {"side": SIDE, "nu": float(nu), "rho_bar": float(rho), "seed": SEED,
            "M": M, "pattern": PATTERN, "arms": [arm],
            "imageset": imageset, "images": [image],
            "cell_id": f"figs1b_{image}_{arm}_{rho}_{nu}", "audit": False}
    campaign.run_cell(cell)
    xh, x = CAP[-1]
    return np.maximum(xh, 0.0).reshape(SIDE, SIDE), x.reshape(SIDE, SIDE)


def psnr_rad(xhat, truth):
    """Radiometric PSNR: data range = truth.max(), MSE on clipped recon."""
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

    ncol = 1 + len(COLUMNS)                     # truth + 5 arms
    t_start = time.time()
    n_recon = sum(len(COLUMNS) for _ in SUBJECTS)
    done = 0

    # ---- run every reconstruction, store per subject ----
    # results[subject] = dict(truth=.., cells=[(label, xhat, psnr), ...])
    results = {}
    for image, imageset, _crop in SUBJECTS:
        cells = []
        truth = None
        for (arm, rho, nu, label) in COLUMNS:
            t0 = time.time()
            xhat, x = recon(image, imageset, arm, rho, nu)
            truth = x
            p = psnr_rad(xhat, truth)
            cells.append((label, xhat, p))
            done += 1
            print("[fig_s1b] %2d/%2d  %-14s %-11s rho=%.2f nu=%-4d "
                  "PSNR_rad=%6.2f dB  (%.1fs)"
                  % (done, n_recon, image, label, rho, nu, p, time.time() - t0),
                  flush=True)
        results[image] = dict(truth=truth, cells=cells)

    print("[fig_s1b] all %d reconstructions done in %.1fs"
          % (n_recon, time.time() - t_start), flush=True)

    # ---- layout: 3 rows (glyph full, glyph crop, stl_02 full) x ncol ----
    glyph = "detail_glyph_0"
    stl = "stl_02"
    fig = plt.figure(figsize=(1.35 * ncol, 4.35))
    gs = fig.add_gridspec(3, ncol, height_ratios=[1.0, 0.46, 1.0],
                          hspace=0.28, wspace=0.06)

    def show(ax, im, vmax, title=None, title_color="black"):
        ax.imshow(im, cmap="gray", vmin=0.0, vmax=vmax,
                  interpolation="nearest")
        ax.set_xticks([])
        ax.set_yticks([])
        if title is not None:
            ax.set_title(title, fontsize=7.5, pad=2, color=title_color)

    # row 0: glyph full frame
    gtruth = results[glyph]["truth"]
    gvmax = 1.2 * float(gtruth.max())
    ax = fig.add_subplot(gs[0, 0])
    show(ax, gtruth, gvmax, title="truth")
    for c, (label, xhat, p) in enumerate(results[glyph]["cells"], start=1):
        ax = fig.add_subplot(gs[0, c])
        show(ax, xhat, gvmax, title="%s\n%.1f dB" % (label, p))

    # row 1: glyph crop  x in [2,42), y in [24,40)
    def crop(im):
        return im[CROP_Y[0]:CROP_Y[1], CROP_X[0]:CROP_X[1]]

    ax = fig.add_subplot(gs[1, 0])
    show(ax, crop(gtruth), gvmax)
    for c, (label, xhat, p) in enumerate(results[glyph]["cells"], start=1):
        ax = fig.add_subplot(gs[1, c])
        show(ax, crop(xhat), gvmax)

    # row 2: stl_02 full frame
    struth = results[stl]["truth"]
    svmax = 1.2 * float(struth.max())
    ax = fig.add_subplot(gs[2, 0])
    show(ax, struth, svmax, title="truth")
    for c, (label, xhat, p) in enumerate(results[stl]["cells"], start=1):
        ax = fig.add_subplot(gs[2, c])
        show(ax, xhat, svmax, title="%s\n%.1f dB" % (label, p))

    # row labels (left margin)
    fig.text(0.012, 0.80, "detail_glyph_0", rotation=90, va="center",
             ha="left", fontsize=7.5)
    fig.text(0.012, 0.545, "glyph crop\nx[2,42) y[24,40)", rotation=90,
             va="center", ha="left", fontsize=6.0)
    fig.text(0.012, 0.235, "stl_02", rotation=90, va="center", ha="left",
             fontsize=7.5)

    fig.subplots_adjust(left=0.075, right=0.99, top=0.93, bottom=0.01)

    os.makedirs(OUT_DIR, exist_ok=True)
    fig.savefig(OUT_PDF)
    fig.savefig(OUT_PNG, dpi=200)
    plt.close(fig)
    print("[fig_s1b] wrote %s" % OUT_PDF, flush=True)
    print("[fig_s1b] wrote %s (200 dpi)" % OUT_PNG, flush=True)

    # ---- final PSNR summary ----
    print("\n[fig_s1b] === PSNR_rad summary (dB) ===", flush=True)
    for image, _iset, _crop in SUBJECTS:
        print("[fig_s1b] %s (vmax=%.3g):" % (image, 1.2 * float(
            results[image]["truth"].max())), flush=True)
        for (label, _xh, p) in results[image]["cells"]:
            print("[fig_s1b]     %-12s %7.3f dB" % (label, p), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
