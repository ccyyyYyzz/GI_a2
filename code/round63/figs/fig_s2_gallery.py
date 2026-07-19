"""Study-2 reconstruction gallery (k=16, RQL, frozen production path).

Rows and their DISCLOSED selection rules:
  1. detail32_spokes_0    — the preregistered panel-(d) display subject
                            (docs/ROUND63_FIGURE_PREREG_ADDENDUM_LADDER.md);
  2. detail32_maze_2      — the largest per-scene fixed-dwell gain in the
                            frozen verdict table (+9.70 dB; a disclosed
                            BEST case, labeled as such);
  3. detail32_microtexture_1 — failing-family representative (honest
                            counterpoint; ties to the C_u mechanism).
Columns: truth | safe rho=0.05 nu=2000 | fast rho=0.6 nu=200 (1/10 dwell) |
fast rho=0.6 nu=2000. Seed 0. Display range per row: [0, 1.2*truth.max].
Captures xhat by wrapping MET.main_metrics (render_recons.py pattern).
"""
import os
import sys

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

ROWS = [
    ("detail32_spokes_0", "spokes (prereg subject)"),
    ("detail32_maze_2", "maze (largest gain, best case)"),
    ("detail32_microtexture_1", "microtexture (failing family)"),
]
COLS = [(0.05, 2000, "safe $\\nu$=2000"),
        (0.6, 200, "fast $\\nu$=200 (1/10 dwell)"),
        (0.6, 2000, "fast $\\nu$=2000")]


def recon(image, rho, nu):
    CAP.clear()
    cell = {"side": 32, "nu": float(nu), "rho_bar": float(rho), "seed": 0,
            "M": 1024, "pattern": "sparsek", "k": 16, "arms": ["RQL"],
            "imageset": "detail32", "images": [image],
            "cell_id": f"gal_{image}_{rho}_{nu}", "audit": False}
    campaign.run_cell(cell)
    xh, x = CAP[-1]
    return np.maximum(xh, 0.0).reshape(32, 32), x.reshape(32, 32)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "serif", "font.size": 8,
                         "axes.linewidth": 0.6, "pdf.fonttype": 42})

    fig, axes = plt.subplots(3, 4, figsize=(6.4, 5.1))
    for r, (image, rowlabel) in enumerate(ROWS):
        panels = []
        truth = None
        for (rho, nu, clabel) in COLS:
            xh, x = recon(image, rho, nu)
            truth = x
            mse = float(np.mean((xh - x) ** 2))
            p = 10 * np.log10(float(x.max()) ** 2 / mse)
            panels.append((xh, clabel, p))
            print(f"[gal] {image} rho={rho} nu={nu} PSNR_rad={p:.2f}",
                  flush=True)
        vmax = 1.2 * truth.max()

        def show(ax, im, title):
            ax.imshow(im, cmap="gray", vmin=0, vmax=vmax,
                      interpolation="nearest")
            ax.set_xticks([]); ax.set_yticks([])
            if title:
                ax.set_title(title, fontsize=7.5, pad=3)

        show(axes[r][0], truth, "truth" if r == 0 else None)
        axes[r][0].set_ylabel(rowlabel, fontsize=7)
        for c, (im, clabel, p) in enumerate(panels, start=1):
            ttl = f"{clabel}\n{p:.1f} dB" if r == 0 else f"{p:.1f} dB"
            show(axes[r][c], im, ttl)

    fig.subplots_adjust(left=0.06, right=0.99, top=0.93, bottom=0.01,
                        wspace=0.06, hspace=0.12)
    out = r"D:\GI_another\paper\figs\fig_s2_gallery"
    fig.savefig(out + ".pdf")
    fig.savefig(out + ".png", dpi=200)
    print("[gal] saved", out, flush=True)


if __name__ == "__main__":
    main()
