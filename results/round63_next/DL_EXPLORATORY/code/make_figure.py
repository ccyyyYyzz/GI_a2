#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""EXPLORATORY / DEV-ONLY summary figure for the SCGI vs DLGI head-to-head."""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "D:/GI_another/results/round63_next/DL_EXPLORATORY"
hh = json.load(open(os.path.join(OUT, "head_to_head_results.json")))
og = json.load(open(os.path.join(OUT, "off_grid_results.json")))

cells = hh["cells"]
labels = [f"tc={c['tc']:.0f}\ncv={c['cv']:.2f}" for c in cells]
arms = ["arm0", "L_exp", "L_mix", "D"]
arm_names = {"arm0": "Arm 0 (no corr)", "L_exp": "Arm L SCGI-exp (OOD)",
             "L_mix": "Arm L SCGI-mix", "D": "Arm D DLGI"}
colors = {"arm0": "#b0b0b0", "L_exp": "#e07b39", "L_mix": "#3b7dd8", "D": "#2ca24c"}

fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))

# Panel A: on-grid scene PSNR
x = np.arange(len(cells)); w = 0.2
for i, a in enumerate(arms):
    ax[0].bar(x + (i - 1.5) * w, [c["psnr"][a] for c in cells], w,
              label=arm_names[a], color=colors[a])
for j, c in enumerate(cells):
    ax[0].hlines(c["psnr"]["oracle"], x[j] - 2 * w, x[j] + 2 * w,
                 color="k", ls="--", lw=1.2, label="oracle (true gain)" if j == 0 else None)
ax[0].set_xticks(x); ax[0].set_xticklabels(labels)
ax[0].set_ylabel("scene PSNR (dB)")
ax[0].set_title("(A) On-grid scene PSNR\n6 DEV bridge scenes x 3 seeds")
ax[0].legend(fontsize=7.5, loc="upper right")
ax[0].grid(axis="y", alpha=0.3)

# Panel B: medium t_c and CV rel-err (mean over on-grid cells)
def mean_re(kind, arm):
    return float(np.mean([c[kind][arm] for c in cells]))
groups = ["D\n(product)", "L_mix\n(trash,in-fam)", "L_exp\n(trash,OOD)", "oracle\nmonitor"]
gk = ["D", "L_mix", "L_exp", "oracle"]
tc_re = [mean_re("tc_relerr", k) for k in gk]
cv_re = [mean_re("cv_relerr", k) for k in gk]
xg = np.arange(len(groups))
ax[1].bar(xg - 0.2, tc_re, 0.4, label="t_c rel-err", color="#7a4fbf")
ax[1].bar(xg + 0.2, cv_re, 0.4, label="CV rel-err", color="#d94f8a")
ax[1].axhline(0.5, color="r", ls=":", lw=1, label="0.5 (uninformative)")
ax[1].set_xticks(xg); ax[1].set_xticklabels(groups, fontsize=8)
ax[1].set_ylabel("median relative error")
ax[1].set_title("(B) Medium products: 'their trash vs our product'\nmean over on-grid OU cells")
ax[1].legend(fontsize=8); ax[1].grid(axis="y", alpha=0.3)

# Panel C: OOD / family-shift collapse (PSNR gain over floor)
def gain_over_floor(cell_list, arm):
    return [c["psnr"][arm] - c["psnr"]["arm0"] for c in cell_list]
allc = cells
xL = np.arange(len(allc)); w2 = 0.28
ax[2].bar(xL - w2, gain_over_floor(allc, "L_exp"), w2, label="SCGI-exp (assumed exp family)", color=colors["L_exp"])
ax[2].bar(xL, gain_over_floor(allc, "L_mix"), w2, label="SCGI-mix (family matched)", color=colors["L_mix"])
ax[2].bar(xL + w2, gain_over_floor(allc, "D"), w2, label="DLGI (model-based, no training)", color=colors["D"])
ax[2].axhline(0, color="k", lw=0.8)
ax[2].set_xticks(xL); ax[2].set_xticklabels(labels)
ax[2].set_ylabel("PSNR gain over no-correction floor (dB)")
ax[2].set_title("(C) Assumed-corruption-family fragility\nSCGI-exp tested out-of-family on OU")
ax[2].legend(fontsize=7.5); ax[2].grid(axis="y", alpha=0.3)

fig.suptitle("EXPLORATORY / DEV-ONLY  -  SCGI-style learned correction vs DLGI on frozen OU-drift protocol "
             "(faithful-in-spirit, not preregistered)", fontsize=10, y=1.02)
fig.tight_layout()
fn = os.path.join(OUT, "head_to_head_figure.png")
fig.savefig(fn, dpi=130, bbox_inches="tight")
print("saved", fn)
