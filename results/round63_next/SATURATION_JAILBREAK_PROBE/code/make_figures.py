"""Figures for the MOLT jailbreak probe."""
import json
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import saturation_core as sc
import gi_operator as op

OUT = os.path.dirname(HERE)
FIG = os.path.join(OUT, "figs")
DATA = "data/r63_bridge_scenes"
os.makedirs(FIG, exist_ok=True)
rng = np.random.default_rng(7)

# ---- FIG 1: P0 mechanism -- product separates, coherent collapses ----------
x = op.load_scene(os.path.join(DATA, "bridge_contour_0.npz"))
mA = (np.random.default_rng(11).random(op.P) < 0.15).astype(float)
mB = (np.random.default_rng(12).random(op.P) < 0.05).astype(float)
vA = (mA * x); vA = vA[vA > 0]
vB = (mB * x); vB = vB[vB > 0]; vB *= vA.sum() / vB.sum()      # equal p1
t = np.logspace(-1.3, 0.75, 60); a = t / vA.sum()
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
ax[0].plot(t, sc.Q_exp(a, vA), 'C0-', lw=2, label="scene A (low $p_2$)")
ax[0].plot(t, sc.Q_exp(a, vB), 'C3-', lw=2, label="scene B (high $p_2$, same $p_1$)")
ax[0].plot(t, sc.Q_gamma(a, vA, 10), 'C0--', lw=1, alpha=.7, label="A, Gamma $k$=10")
ax[0].plot(t, sc.Q_gamma(a, vB, 10), 'C3--', lw=1, alpha=.7, label="B, Gamma $k$=10")
ax[0].set_title("Intensity-additive (product): $p_2$ channel LIVES")
ax[0].set_xlabel("load $t=\\rho p_1$"); ax[0].set_ylabel("no-fire fraction $r$")
ax[0].set_xscale("log"); ax[0].legend(fontsize=8)
ax[1].plot(t, sc.coherent_survival_exact(a, vA), 'C0-', lw=2, label="scene A")
ax[1].plot(t, sc.coherent_survival_exact(a, vB), 'C3--', lw=2, label="scene B (identical)")
ax[1].set_title("Coherent collapse: $r=1/(1+t)$, only $p_1$ (DEAD)")
ax[1].set_xlabel("load $t=\\rho p_1$"); ax[1].set_ylabel("no-fire fraction $r$")
ax[1].set_xscale("log"); ax[1].legend(fontsize=8)
fig.suptitle("Gate P0: the mechanism is life-or-death on intensity vs coherent addition")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_p0_mechanism.png"), dpi=130)
plt.close(fig)

# ---- FIG 2: T3 null-pair scenes x, x', delta ------------------------------
d = np.load(os.path.join(OUT, "t3_scenes.npz"))
xg, xpg, dg = d["x"], d["xp"], d["delta"]
fig, ax = plt.subplots(1, 3, figsize=(11, 3.8))
im0 = ax[0].imshow(xg, cmap="gray", vmin=0, vmax=1); ax[0].set_title("$x$")
im1 = ax[1].imshow(xpg, cmap="gray", vmin=0, vmax=1); ax[1].set_title("$x' = x+\\delta$")
im2 = ax[2].imshow(dg, cmap="RdBu_r", vmin=-np.abs(dg).max(), vmax=np.abs(dg).max())
ax[2].set_title("$\\delta \\in \\ker M$  ($\\|\\delta\\|$=%.2f)" % np.linalg.norm(dg))
for a_, im in zip(ax, (im0, im1, im2)):
    a_.axis("off"); fig.colorbar(im, ax=a_, fraction=0.046)
fig.suptitle("T3: null-space pair -- all 51 linear buckets identical to machine precision")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_t3_scenes.png"), dpi=130)
plt.close(fig)
# also individual PNGs
for arr, nm, cmap in [(xg, "t3_x", "gray"), (xpg, "t3_xprime", "gray")]:
    plt.imsave(os.path.join(FIG, nm + ".png"), arr, cmap=cmap, vmin=0, vmax=1)
plt.imsave(os.path.join(FIG, "t3_delta.png"), dg, cmap="RdBu_r",
           vmin=-np.abs(dg).max(), vmax=np.abs(dg).max())

# ---- FIG 3: d'-vs-budget (from t3 json) with P2/P3 markers -----------------
tj = json.load(open(os.path.join(OUT, "t3_nullspace_break.json")))
bud = np.array(tj["dprime_curve"]["budget"]); dfl = np.array(tj["dprime_curve"]["flat"])
dop = np.array(tj["dprime_curve"]["optimal"])
camp = tj["linear_campaign"]
fig, ax = plt.subplots(figsize=(6.4, 4.4))
ax.loglog(bud / camp, dfl, 'C0-', lw=2, label="flat allocation (textured pair)")
ax.loglog(bud / camp, dop, 'C1--', lw=2, label="optimal single-(mask,level)")
ax.axhline(3, color='C3', ls=':', lw=1.5, label="$d'$=3 bar")
ax.axvline(1, color='k', ls='--', lw=1, alpha=.6, label="linear campaign")
ax.set_xlabel("budget / linear campaign"); ax.set_ylabel("discriminability $d'$")
ax.set_title("T3/P2: null-space break needs $\\sim$10$^4$x the linear campaign for $d'$=3")
ax.legend(fontsize=8); ax.grid(alpha=.3, which="both")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_dprime_budget.png"), dpi=130)
plt.close(fig)

# ---- FIG 4: P3 rank doubling singular spectrum -----------------------------
M = op.build_operator()
LO, HI = 0.08, 0.92
xr = np.clip(op.load_scene(os.path.join(DATA, "bridge_microtex_0.npz")), LO, HI)
svM = np.linalg.svd(M, compute_uv=False)
J = np.vstack([M, 2 * M * xr[None, :]])
svJ = np.linalg.svd(J, compute_uv=False)
fig, ax = plt.subplots(figsize=(6.4, 4.4))
ax.semilogy(np.arange(1, len(svM) + 1), svM, 'C0o-', ms=4, label="$M$ (rank 51)")
ax.semilogy(np.arange(1, len(svJ) + 1), svJ, 'C3s-', ms=3,
            label="$[M;\\,2M\\,\\mathrm{diag}(x)]$ (rank 102)")
ax.axvline(51, color='C0', ls=':', alpha=.6); ax.axvline(102, color='C3', ls=':', alpha=.6)
ax.set_xlabel("singular-value index"); ax.set_ylabel("singular value")
ax.set_title("P3: the $x^2$ channel doubles measurement rank ($51\\to102$)")
ax.legend(fontsize=9); ax.grid(alpha=.3, which="both"); ax.set_xlim(0, 115)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_p3_rank_doubling.png"), dpi=130)
plt.close(fig)

print("wrote figures to", FIG)
for f in sorted(os.listdir(FIG)):
    print("  ", f)
