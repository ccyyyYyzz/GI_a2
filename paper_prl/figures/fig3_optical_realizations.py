"""
FIGURE 3 (R43 §6.4) — Two optical realizations and the honest boundary.

Full-width, three compact parts, all from COMMITTED FROZEN artifacts:
  (1) Thin screen sealed 2% power   -> CONFIRMATORY_RESULTS.json @ b37c841 (D2, D4)
  (2) Calibration boundary          -> CONFIRMATORY_RESULTS.json @ b37c841 (D3, D5)
  (3) Complete-scrambling endpoint  -> SCRAMBLE_RESULTS.json      @ ed7a1e0
extracted verbatim to paper_prl/figures/_frozen_sources/.
No data generation; every number is read from the committed JSON.

Single caption sentence (R43 §6.4):
 "The order law survives both weak and complete scrambling; practical power
  survives, while attribution calibration has a measured domain."
"""
import json, os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LogNorm
from prl_style import apply_style, save, COL2, C_GENERIC, C_ORTHO, C_ANCHOR, C_WALL, OI

HERE = os.path.dirname(os.path.abspath(__file__))
FS = os.path.join(HERE, "_frozen_sources")
conf = json.load(open(os.path.join(FS, "CONFIRMATORY_RESULTS.committed.json")))
scr = json.load(open(os.path.join(FS, "SCRAMBLE_RESULTS.committed.json")))
apply_style()

fig = plt.figure(figsize=(COL2, 2.55))
# Panel order leads with the physics (complete scrambling: image dies, detection
# lives), then thin-screen power, then the calibration boundary (register fix, F4).
gs = GridSpec(1, 3, figure=fig, width_ratios=[0.85, 1.20, 1.05], wspace=0.50,
              left=0.065, right=0.988, top=0.80, bottom=0.16)
ax_scr = fig.add_subplot(gs[0, 0])    # PART 3: complete scrambling (physics lead)
ax_thin = fig.add_subplot(gs[0, 1])   # PART 1: thin-screen 2% power map
ax_cal = fig.add_subplot(gs[0, 2])    # PART 2: calibration boundary

# ============ PART 1: thin-screen 2% power operating map (77/81) =============
mp = conf["D2"]["analytic_map_2pct"]
sfs = [0.3, 0.6, 1.0]; shapes = ["flat", "k^-1", "k^-2"]
kwfs = [1, 2, 4]; claims = [1.25, 1.5, 1.8]
rows = [(sf, sh) for sf in sfs for sh in shapes]      # 9
cols = [(k, c) for k in kwfs for c in claims]          # 9
M = np.full((9, 9), np.nan); fail = np.zeros((9, 9), bool)
for e in mp:
    r = rows.index((e["sf"], e["shape"])); c = cols.index((e["kwf"], e["claim"]))
    M[r, c] = e["T_det"]; fail[r, c] = not e["pass_4096"]
im = ax_thin.imshow(M, cmap="viridis_r", norm=LogNorm(vmin=400, vmax=9000),
                aspect="auto", origin="upper")
for r in range(9):
    for c in range(9):
        if fail[r, c]:
            ax_thin.plot(c, r, "x", color=OI["vermil"], ms=5.5, mew=1.5, zorder=5)
# 3x3 super-block separators (each block = 3 shapes x 3 claims)
for p in (2.5, 5.5):
    ax_thin.axhline(p, color="white", lw=1.1); ax_thin.axvline(p, color="white", lw=1.1)
# 7pt group tick labels at block centers
ax_thin.set_xticks([1, 4, 7]); ax_thin.set_yticks([1, 4, 7])
ax_thin.set_xticklabels(["1", "2", "4"], fontsize=7)
ax_thin.set_yticklabels(["0.3", "0.6", "1.0"], fontsize=7)
ax_thin.set_xlabel("code weight  (shape inner)", fontsize=7, labelpad=2)
ax_thin.set_ylabel("static factor  (claim inner)", fontsize=7, labelpad=2)
ax_thin.tick_params(length=2.0)
cb = fig.colorbar(im, ax=ax_thin, fraction=0.046, pad=0.03)
cb.set_label(r"$T_{\det}$ (banks)", fontsize=7); cb.ax.tick_params(labelsize=6.5, length=1.5)
ax_thin.set_title("thin screen: 2% beyond-band power", fontsize=7.4, pad=8)
ax_thin.text(0.5, 1.13,
         r"$\mathbf{77/81}$ cells detectable  |  best $\mathbf{453}$ banks  |  MC LCB $\mathbf{0.990}$",
         transform=ax_thin.transAxes, ha="center", va="bottom", fontsize=5.9)
# fixed-vs-fresh callout (D4 best cell)
d4 = conf["D4"]["rows"][0]
ax_thin.text(0.02, -0.30,
         rf"fixed vs fresh (best cell): $\mathbf{{{int(d4['fixed_latency'])}}}$ vs "
         rf"$\mathbf{{{int(d4['fresh_latency'])}}}$ banks  ($\times$ marks = 4 slow cells)",
         transform=ax_thin.transAxes, ha="left", va="top", fontsize=6.3)

# ============ PART 2: calibration boundary (honest threshold graphic) ========
d3 = conf["D3"]; d5 = conf["D5"]["rows"]
def d5fa(axis, lvl):
    return [r for r in d5 if r["axis"] == axis and r["level"] == lvl][0]["nontarget_fa"]
labels = ["in-band", "amplitude", "timescale", "rot 10%", "rot 20%", "slope -1"]
vals = [d3["fa"]["inband"], d3["fa"]["amplitude"], d3["fa"]["timescale"],
        d5fa("basis_rotation", 0.1), d5fa("basis_rotation", 0.2), d5fa("spectral_slope", -1)]
bar_bar = 0.05
colors = [OI["green"] if v <= bar_bar else C_ORTHO for v in vals]
xs = np.arange(len(labels))
ax_cal.bar(xs, vals, width=0.66, color=colors, edgecolor="black", lw=0.5, zorder=3)
ax_cal.axhline(bar_bar, color="black", ls="--", lw=0.9, zorder=4)
ax_cal.text(-0.35, bar_bar + 0.003, "frozen 5% bar", ha="left", va="bottom", fontsize=6.2)
for x, v in zip(xs, vals):
    ax_cal.text(x, v + 0.003, f"{v:.3f}", ha="center", va="bottom", fontsize=6.0)
ax_cal.set_xticks(xs); ax_cal.set_xticklabels(labels, fontsize=7, rotation=25, ha="right")
ax_cal.set_ylabel("non-target false-alarm rate", fontsize=7)
ax_cal.set_ylim(0, 0.12)
ax_cal.set_title("calibration boundary (1% threshold)", fontsize=7.4, pad=8)
ax_cal.text(0.5, 1.13,
         r"TPR $\mathbf{0.988}$ | bal. acc. $\mathbf{0.916}$ | $\leq\mathbf{10\%}$ rot.",
         transform=ax_cal.transAxes, ha="center", va="bottom", fontsize=6.2)
ax_cal.tick_params(direction="in", right=True, length=2.2)
# pass/fail legend (upper-right, clear of tallest bar)
from matplotlib.patches import Patch
ax_cal.legend(handles=[Patch(fc=OI["green"], ec="black", lw=0.4, label=r"$\leq$ 5% bar"),
                    Patch(fc=C_ORTHO, ec="black", lw=0.4, label="> 5% bar")],
           loc="upper right", frameon=False, fontsize=6.2, handlelength=1.0,
           handleheight=1.0, labelspacing=0.2, borderaxespad=0.2)

# ============ PART 3: complete-scrambling endpoint ===========================
mean_coh = scr["valC_mean_detector"]["coherent"]["auc"]
mean_ort = scr["valC_mean_detector"]["ortho"]["auc"]
duel = scr["detection"]["coherent"]["eps_0.05"]["8192"]
duel_auc = duel["auc"]; duel_dp = duel["d_prime"]
rank1 = scr["valB_rank1"]["frac_cov_energy_in_OhO_direction"]
mean_null = scr["valA_mean_null"]["delta_mean_rel"]

bars = [mean_coh, mean_ort, duel_auc]
blabels = ["mean\n(generic)", "mean\n(orthog.)", "covariance\nduel 5%"]
bcolors = [C_WALL, C_WALL, C_GENERIC]
xs3 = np.arange(3)
ax_scr.bar(xs3, bars, width=0.62, color=bcolors, edgecolor="black", lw=0.5, zorder=3)
ax_scr.axhline(0.5, color="black", ls="--", lw=0.9, zorder=4)
ax_scr.text(2.45, 0.515, "chance", ha="right", va="bottom", fontsize=5.8)
ax_scr.set_ylim(0.45, 1.03)
ax_scr.set_xticks(xs3); ax_scr.set_xticklabels(blabels, fontsize=7)
ax_scr.set_ylabel("detection AUC", fontsize=7)
ax_scr.set_title("complete scrambling: image dies,\ndetection lives", fontsize=7.4, pad=8)
ax_scr.text(2, duel_auc + 0.012, rf"$d'={duel_dp:.1f}$", ha="center", va="bottom",
         fontsize=6.0, color=C_GENERIC)
ax_scr.text(0.5, 0.545, "mean channel\nblind (only DC)", ha="center", va="bottom",
         fontsize=5.4, color="black")
ax_scr.text(0.0, -0.30,
         rf"rank-one covariance $\mathbf{{99.99\%}}$;"
         "\n" rf"mean null $\Delta=\mathbf{{6.6\times10^{{-19}}}}$",
         transform=ax_scr.transAxes, ha="left", va="top", fontsize=6.3)
ax_scr.tick_params(direction="in", right=True, length=2.2)

save(fig, os.path.join(HERE, "fig3_optical_realizations"))
