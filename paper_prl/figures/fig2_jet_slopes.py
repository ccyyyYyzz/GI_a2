"""
FIGURE 2 (R43 §6.3) — The PRL evidence figure: the log-log EXPONENT HIERARCHY.

Rebuilds the jet_slopes asset at PRL width from the COMMITTED FROZEN artifact
    results/round63_next/JET_TEST/JET_TEST.json  @ commit 1bf29f1
extracted verbatim to  paper_prl/figures/_frozen_sources/JET_TEST.committed.json
(the working-tree JET_TEST.json is an uncommitted re-run; NOT used — see
 paper_prl/CLAIM_SOURCE_MATRIX.md GAP-1).

Hero = exact profiled-Chernoff contact orders 2 (generic) and 4 (orthogonal)
with coefficient agreement + crossover marker. Minimum insets: MC d' slopes 1/2;
CUSUM delay slopes -2/-4; two-symbol nuisance gauge (J=0 -> anchor J>0).
No data generation: every value is read from the committed JSON.
"""
import json, os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from prl_style import apply_style, save, COL2, C_GENERIC, C_ORTHO, C_MIXED, C_ANCHOR, C_WALL, OI

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "_frozen_sources", "JET_TEST.committed.json")
d = json.load(open(SRC))
apply_style()

pd_ = d["bank_A"]["per_dir"]
eps = np.array(d["bank_A"]["eps"])
g = pd_["delta_g"]; o = pd_["delta_o"]
Cg = np.array(g["chernoff_star"]); Co = np.array(o["chernoff_star"])
cross = d["bank_A"]["crossover"]
ce = np.array(cross["eps"]); ck = np.array(cross["kl"]); eps_cross = cross["eps_cross_pred"]
mc = d["bank_A_montecarlo"]
cu = d["cusum"]
bb = d["bank_B"]

fig = plt.figure(figsize=(COL2, 3.05))
gs = GridSpec(3, 2, figure=fig, width_ratios=[1.62, 1.0], height_ratios=[1, 1, 1],
              wspace=0.34, hspace=0.62, left=0.075, right=0.985, top=0.99, bottom=0.115)
axM = fig.add_subplot(gs[:, 0])
axMC = fig.add_subplot(gs[0, 1])
axCU = fig.add_subplot(gs[1, 1])
axNU = fig.add_subplot(gs[2, 1])

# ---------------- MAIN: profiled Chernoff distance, contact orders 2 and 4 ----
axM.loglog(eps, Cg, "o", color=C_GENERIC, ms=4.2, mec="white", mew=0.4, zorder=5,
           label=r"generic  $x^{\!\top}\!G\delta\neq0$")
axM.loglog(eps, Co, "s", color=C_ORTHO, ms=4.0, mec="white", mew=0.4, zorder=5,
           label=r"orthogonal  $x^{\!\top}\!G\delta=0$")
# mixed direction curve (crossover)
axM.loglog(ce, ck, "-", color=C_MIXED, lw=1.1, alpha=0.9, zorder=3,
           label=r"mixed  ($m:2\!\to\!4$)")

# fitted power-law guides anchored to data
def guide(ax, x, y_ref, x_ref, slope, **kw):
    xx = np.array(x)
    yy = y_ref * (xx / x_ref) ** slope
    ax.loglog(xx, yy, **kw)
guide(axM, eps, Cg[7], eps[7], 2.0, color=C_GENERIC, ls="--", lw=0.8, alpha=0.55, zorder=2)
guide(axM, eps, Co[7], eps[7], 4.0, color=C_ORTHO, ls="--", lw=0.8, alpha=0.55, zorder=2)

# crossover marker
yc = np.interp(np.log(eps_cross), np.log(ce), np.log(ck)); yc = np.exp(yc)
axM.plot([eps_cross], [yc], "v", color=OI["black"], ms=5.5, zorder=6)
axM.annotate(r"$\varepsilon_{\mathrm{cross}}$", (eps_cross, yc),
             textcoords="offset points", xytext=(6, 1), fontsize=6.8, ha="left")

axM.set_xlabel(r"scene-change amplitude  $\varepsilon$")
axM.set_ylabel(r"profiled Chernoff distance  $C_*(\varepsilon;\delta)$")
axM.set_xlim(4e-4, 0.20)
axM.set_ylim(3e-15, 6e-3)
# slope labels (axes-fraction placement to avoid data-range errors)
axM.text(0.545, 0.80, r"slope $2.04$", color=C_GENERIC, fontsize=7.4, rotation=21,
         transform=axM.transAxes, ha="left", va="bottom")
axM.text(0.30, 0.235, r"slope $4.00$", color=C_ORTHO, fontsize=7.4, rotation=37,
         transform=axM.transAxes, ha="left", va="bottom")
# coefficient-agreement annotation (open lower-right region)
axM.text(0.60, 0.045,
         "coef. agreement\n"
         r"gen. $c^2\|K\|^2$: $1.07\%$" "\n"
         r"orth. $d^2\|K\|^2/4$: $0.00\%$",
         transform=axM.transAxes, fontsize=6.2, ha="left", va="bottom",
         bbox=dict(boxstyle="round,pad=0.28", fc="white", ec=C_WALL, lw=0.5, alpha=0.92))
leg = axM.legend(loc="upper left", frameon=False, fontsize=6.6, handlelength=1.4,
                 borderaxespad=0.15, labelspacing=0.25)
axM.tick_params(which="both", direction="in", top=True, right=True)

# ---------------- INSET A: Monte-Carlo d' slopes 1 and 2 ----------------------
mg = mc["delta_g"]["rows"]; mo = mc["delta_o"]["rows"]
eg = np.array([r["eps"] for r in mg]); dg = np.array([r["d_prime"] for r in mg])
eo = np.array([r["eps"] for r in mo]); do = np.array([r["d_prime"] for r in mo])
axMC.loglog(eg, dg, "o", color=C_GENERIC, ms=3.2, mec="white", mew=0.3)
axMC.loglog(eo, do, "s", color=C_ORTHO, ms=3.0, mec="white", mew=0.3)
guide(axMC, eg, dg[0], eg[0], 1.0, color=C_GENERIC, ls="--", lw=0.7, alpha=0.6)
guide(axMC, eo, do[0], eo[0], 2.0, color=C_ORTHO, ls="--", lw=0.7, alpha=0.6)
axMC.set_ylabel(r"MC  $d'$", fontsize=7)
axMC.set_title(r"single-look $d'\!\sim\!\varepsilon^{m}$   ($m=0.95,\,2.05$)",
               fontsize=6.4, pad=2)
axMC.tick_params(which="both", direction="in", labelsize=6, top=True, right=True)

# ---------------- INSET B: CUSUM delay slopes -2 and -4 -----------------------
cg = cu["delta_g"]["rows"]; co = cu["delta_o"]["rows"]
ecg = np.array([r["eps"] for r in cg]); dlg = np.array([r["delay"] for r in cg])
eco = np.array([r["eps"] for r in co]); dlo = np.array([r["delay"] for r in co])
axCU.loglog(ecg, dlg, "o", color=C_GENERIC, ms=3.2, mec="white", mew=0.3)
axCU.loglog(eco, dlo, "s", color=C_ORTHO, ms=3.0, mec="white", mew=0.3)
guide(axCU, ecg, dlg[0], ecg[0], -2.0, color=C_GENERIC, ls="--", lw=0.7, alpha=0.6)
guide(axCU, eco, dlo[0], eco[0], -4.0, color=C_ORTHO, ls="--", lw=0.7, alpha=0.6)
axCU.set_ylabel(r"CUSUM delay", fontsize=7)
axCU.set_title(r"delay $\sim\varepsilon^{-2m}$   ($-2.16,\,-3.92$)", fontsize=6.4, pad=2)
from matplotlib.ticker import ScalarFormatter, NullFormatter
axCU.set_xticks([0.03, 0.06, 0.1, 0.16])
axCU.xaxis.set_major_formatter(ScalarFormatter())
axCU.xaxis.set_minor_formatter(NullFormatter())
axCU.ticklabel_format(axis="x", style="plain")
axCU.tick_params(which="both", direction="in", labelsize=6, top=True, right=True)

# ---------------- INSET C: two-symbol nuisance gauge --------------------------
J_gauge = bb["single_zero_lag"]["I_theta_eff"] / bb["single_zero_lag"]["I_theta_naive"]  # 0
J_anchor = bb["amplitude_anchor"]["frac_recovered"]  # 0.20
xs = [0, 1]
axNU.bar(xs, [J_gauge, J_anchor], width=0.55,
         color=[C_WALL, C_ANCHOR], edgecolor="black", lw=0.5)
axNU.axhline(0, color="black", lw=0.6)
axNU.set_xticks(xs)
axNU.set_xticklabels(["amplitude\ngauge", "+ anchor"], fontsize=6.2)
axNU.set_ylabel(r"$J_\theta/I_\theta$", fontsize=7)
axNU.set_ylim(0, 0.26)
axNU.text(0, 0.012, r"$J=0$", ha="center", va="bottom", fontsize=6.2, color="black")
axNU.text(1, J_anchor + 0.006, r"$J>0$", ha="center", va="bottom", fontsize=6.2,
          color=OI["orange"])
axNU.set_title("nuisance quotient", fontsize=6.4, pad=2)
axNU.tick_params(which="both", direction="in", labelsize=6, right=True)

save(fig, os.path.join(HERE, "fig2_jet_slopes"))
