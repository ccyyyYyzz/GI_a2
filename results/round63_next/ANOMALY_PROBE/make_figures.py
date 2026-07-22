"""3-panel visuals (True planted / A2 unconstrained-DL / A3 information-guided)
for the 4 most illustrative resolvable-anomaly cases, clean regime (VQAE prior).
Anomaly bounding box marked in red; a zoomed inset shows the anomaly region."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

OUT = Path(r"D:/GI_another/results/round63_next/ANOMALY_PROBE")
S = 64

def make_delta_bbox(kind, cy, cx):
    """Reproduce anomaly footprint bbox (+2px) for marking. Mirrors anomaly_probe.make_delta."""
    yy, xx = np.mgrid[0:S, 0:S]
    d = np.zeros((S, S))
    if kind == "blob5_bright":
        d[((yy-cy)**2+(xx-cx)**2) <= 2.5**2] = 1
    elif kind == "blob6_dark":
        d[((yy-cy)**2+(xx-cx)**2) <= 3.0**2] = 1
    elif kind == "bar_3x8_bright":
        d[cy:cy+3, cx:cx+8] = 1
    elif kind == "ring_r5_bright":
        d[(((yy-cy)**2+(xx-cx)**2) <= 2.5**2) & (((yy-cy)**2+(xx-cx)**2) > 1.2**2)] = 1
    elif kind == "point_1px":
        d[cy, cx] = 1
    elif kind == "point_2px":
        d[cy:cy+2, cx:cx+2] = 1
    ys, xs = np.nonzero(d)
    y0 = max(0, ys.min()-2); y1 = min(S, ys.max()+1+2)
    x0 = max(0, xs.min()-2); x1 = min(S, xs.max()+1+2)
    return y0, y1, x0, x1

def main():
    res = json.loads((OUT/"anomaly_probe_results.json").read_text(encoding="utf-8"))
    recs = res["regimes"]["clean"]["records"]
    npz = np.load(OUT/"_recon_cache.npz")
    plans = json.loads((OUT/"_plans.json").read_text(encoding="utf-8"))["plans"]

    # pick 4 illustrative resolvable cases (seed 101): highest witnessed fraction,
    # distinct images, distinct shapes if possible
    seed = 101
    cand = [r for r in recs if r["anomaly_class"] == "resolvable" and r["seed"] == seed and r["arm"] == "A3"]
    cand.sort(key=lambda r: -r["witnessed_fraction"])
    picked = []
    seen_shapes = set()
    for r in cand:
        if r["anomaly"] not in seen_shapes:
            picked.append(r); seen_shapes.add(r["anomaly"])
        if len(picked) == 4:
            break
    while len(picked) < 4 and len(picked) < len(cand):
        for r in cand:
            if r not in picked:
                picked.append(r); break

    fig, axes = plt.subplots(4, 3, figsize=(7.2, 9.6))
    col_titles = ["True (planted)", "A2  unconstrained DL\n(prior rewrites all)", "A3  information-guided\n(range locked, prior in null)"]
    for row, r in enumerate(picked):
        img = r["image"]; kind = r["anomaly"]
        # find the resolvable anomaly's plotted position from the plan
        pl = plans[str(img)]
        cy, cx = None, None
        for (k, yy, xx) in pl:
            if k == kind:
                cy, cx = yy, xx
        y0, y1, x0, x1 = make_delta_bbox(kind, cy, cx)
        xt = npz[f"{img}_{seed}_xtrue"].reshape(S, S)
        a2 = npz[f"{img}_{seed}_A2"].reshape(S, S)
        a3 = npz[f"{img}_{seed}_A3"].reshape(S, S)
        # bbox metrics for annotation
        def bbpsnr(x):
            m = np.zeros((S, S), bool); m[y0:y1, x0:x1] = True
            mse = np.mean((x[m]-xt[m])**2); return 10*np.log10(1/max(mse, 1e-12))
        panels = [(xt, col_titles[0], None),
                  (a2, col_titles[1], bbpsnr(a2)),
                  (a3, col_titles[2], bbpsnr(a3))]
        for col, (im, title, psnr) in enumerate(panels):
            ax = axes[row, col]
            ax.imshow(im, cmap="gray", vmin=0, vmax=1, interpolation="nearest")
            ax.add_patch(Rectangle((x0-0.5, y0-0.5), x1-x0, y1-y0, fill=False, edgecolor="red", lw=1.4))
            ax.set_xticks([]); ax.set_yticks([])
            if row == 0:
                ax.set_title(title, fontsize=8.5)
            lab = f"{kind}\nwf={r['witnessed_fraction']:.2f}" if col == 0 else (f"bbox PSNR {psnr:.1f} dB" if psnr is not None else "")
            ax.set_xlabel(lab, fontsize=7.5)
        # zoom insets on A2 and A3 (anomaly region)
        for col, im in [(1, a2), (2, a3)]:
            ax = axes[row, col]
            axins = ax.inset_axes([0.62, 0.62, 0.36, 0.36])
            pad = 4
            zy0, zy1 = max(0, y0-pad), min(S, y1+pad); zx0, zx1 = max(0, x0-pad), min(S, x1+pad)
            axins.imshow(im[zy0:zy1, zx0:zx1], cmap="gray", vmin=0, vmax=1, interpolation="nearest")
            axins.set_xticks([]); axins.set_yticks([])
            for sp in axins.spines.values():
                sp.set_edgecolor("yellow"); sp.set_linewidth(1.0)
    fig.suptitle("Anomaly fidelity: unconstrained DL prior (A2) erases resolvable anomalies;\n"
                 "information-guided reconstruction (A3, range locked) preserves them  "
                 "[clean regime, real lane0 rate05 operator + VQAE prior]", fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(OUT/"anomaly_fidelity_3panel.png", dpi=150)
    print("wrote", OUT/"anomaly_fidelity_3panel.png")
    print("cases:", [(r["image"], r["anomaly"], round(r["witnessed_fraction"], 2)) for r in picked])

if __name__ == "__main__":
    main()
