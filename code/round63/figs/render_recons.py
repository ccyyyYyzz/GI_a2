"""Render actual reconstructions through the FROZEN production path.

Display-only: runs campaign.run_cell verbatim (k=16 primary geometry,
confirmatory seeds) and captures each arm's reconstruction by wrapping
MET.main_metrics. Nothing frozen is modified; no metric row is written.
"""
import sys, os
import numpy as np

os.chdir(r"D:\GI_another")
sys.path.insert(0, r"D:\GI_another\code\round63")
import campaign

CAP = []
_orig = campaign.MET.main_metrics
def capture(xh, x, side, **kw):
    CAP.append((np.array(xh, copy=True), np.array(x, copy=True)))
    return _orig(xh, x, side, **kw)
campaign.MET.main_metrics = capture

def one(image, rho, nu):
    CAP.clear()
    cell = {"side": 32, "nu": nu, "rho_bar": rho, "seed": 0, "M": 1024,
            "pattern": "sparsek", "k": 16, "arms": ["RQL"],
            "imageset": "detail32", "images": [image],
            "cell_id": f"render_{image}_{rho}_{nu}", "audit": False}
    rows = campaign.run_cell(cell)
    xh, x = CAP[-1]
    psnr = [r for r in rows][0]
    return np.maximum(xh, 0.0).reshape(32, 32), x.reshape(32, 32)

jobs = [
    ("detail32_spokes_0", [(0.05, 2000, "safe nu=2000"), (0.6, 200, "fast nu=200 (1/10 time)"), (0.6, 2000, "fast nu=2000")]),
    ("detail32_microtexture_1", [(0.05, 2000, "safe nu=2000"), (0.6, 2000, "fast nu=2000")]),
]

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ncol = 4
fig, axes = plt.subplots(2, ncol, figsize=(9.6, 5.1))
for r, (image, pts) in enumerate(jobs):
    panels = []
    truth = None
    for (rho, nu, label) in pts:
        xh, x = one(image, rho, nu)
        truth = x
        mse = float(np.mean((xh - x) ** 2))
        psnr = 10 * np.log10(float(x.max()) ** 2 / mse)
        panels.append((xh, f"{label}\nPSNR_rad {psnr:.1f} dB"))
        print(image, label, f"PSNR_rad={psnr:.2f}")
    vmax = 1.2 * truth.max()
    axes[r][0].imshow(truth, cmap="gray", vmin=0, vmax=vmax, interpolation="nearest")
    axes[r][0].set_title(f"{image}\ntruth", fontsize=8)
    for c, (im, ttl) in enumerate(panels, start=1):
        axes[r][c].imshow(im, cmap="gray", vmin=0, vmax=vmax, interpolation="nearest")
        axes[r][c].set_title(ttl, fontsize=8)
    for c in range(len(panels) + 1, ncol):
        axes[r][c].axis("off")
for ax in axes.ravel():
    ax.set_xticks([]); ax.set_yticks([])
fig.tight_layout()
out = r"D:\GI_another\results\round63_study2\recon_showcase.png"
fig.savefig(out, dpi=180)
print("saved", out)
