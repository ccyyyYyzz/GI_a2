"""
PAPER 2 (Optica) — FIGURE 2: Jet transmutation (R45 §C3 fig2 spec).

KT2 measured normal form  C(eps) = A eps^2 + B eps^4  (a first-order mean leak on
top of a near-contact curvature-rescued covariance direction). Projecting out the
measured mean-leak tangent zeroes A and moves the strict local observability class
from m=1 (slope 2) to m=2 (slope 4) WITHOUT changing the scene perturbation.

Panel L: reconstructed C(eps) for the four KT2 cells (from committed fitted A,B):
  - unprojected total  A eps^2 + B eps^4  -> slope 2 at small eps, bending to 4
  - projected          B eps^4            -> slope 4 everywhere
  with slope-2 / slope-4 guides and the crossover eps_cross = sqrt(A/B).
Panel R: the exact logistic jet-flow (R45 §B2 eqs 2.4/2.6) m_eff(eps)=(1+2y)/(1+y),
  y=(eps/eps_cross)^2, flowing between the integer contact-order fixed points 1 and 2,
  with the four measured projected slopes (all 4.00) and unprojected low-eps slopes.

DATA: KT2_jet_transmutation.json @05f1029 (committed), extracted to _frozen_sources/.
Curves are analytic evaluations of the COMMITTED fitted coefficients (A,B) — no data
generation. Cross-anchor: JET_TEST @1bf29f1 integer KL orders 2.038/4.000.
"""
import json
import os
import numpy as np
import matplotlib.pyplot as plt
from o2_style import apply_style, save, COL2, OI, C_BREAK, C_JET, C_GREY

HERE = os.path.dirname(os.path.abspath(__file__))
FS = os.path.join(HERE, "_frozen_sources")
kt2 = json.load(open(os.path.join(FS, "KT2_jet_transmutation.committed.json")))

cells = kt2["cells"]
eps = np.logspace(np.log10(kt2["params"]["eps_range"][0]),
                  np.log10(kt2["params"]["eps_range"][1]), 200)
apply_style()
fig, (axL, axR) = plt.subplots(1, 2, figsize=(COL2, 3.0))

# ---- Panel L: C(eps) unprojected vs projected ----
rep = cells[0]  # z1=10, z2=1 representative
A = rep["directions"]["orth_firstorder_null"]["A_mean_eps2_coeff"]
B = rep["directions"]["orth_firstorder_null"]["B_curv_eps4_coeff"]
eps_cross = rep["directions"]["orth_firstorder_null"]["eps_cross"]
C_unproj = A * eps**2 + B * eps**4
C_proj = B * eps**4
axL.loglog(eps, C_unproj, color=C_BREAK, lw=1.8,
           label=r"unprojected  $A\varepsilon^2+B\varepsilon^4$")
axL.loglog(eps, C_proj, color=C_JET, lw=1.8, ls="--",
           label=r"projected  $B\varepsilon^4$  ($m{=}2$)")
# slope guides
xg = np.array([4e-3, 2e-2])
axL.loglog(xg, 3e-4 * (xg / xg[0])**2, color="0.55", lw=0.9)
axL.text(xg[1] * 1.15, 3e-4 * (xg[1] / xg[0])**2, "slope 2", fontsize=6, color="0.4")
xg2 = np.array([0.12, 0.5])
axL.loglog(xg2, 5e-2 * (xg2 / xg2[0])**4, color="0.55", lw=0.9)
axL.text(xg2[1] * 1.05, 5e-2 * (xg2[1] / xg2[0])**4, "slope 4", fontsize=6, color="0.4")
axL.axvline(eps_cross, color=C_GREY, lw=0.8, ls=(0, (2, 2)))
axL.text(eps_cross * 1.06, 3e-9, r"$\varepsilon_{\rm cross}$", fontsize=6.5,
         color="0.35", rotation=90, va="bottom")
axL.set_xlabel(r"perturbation amplitude  $\varepsilon$")
axL.set_ylabel(r"Chernoff information  $C(\varepsilon)$")
axL.set_title("Mean-leak projection restores the curvature class", fontsize=7.2, pad=3)
axL.legend(loc="upper left", fontsize=6.0, handlelength=1.6)
axL.text(0.03, 0.03,
         f"cell $z_1{{=}}{rep['z1_mm']:.0f}$, $z_2{{=}}{rep['z2_mm']:.0f}$ mm\n"
         r"$\bot$-residual $\leq5\times10^{-18}$",
         transform=axL.transAxes, fontsize=5.8, color="0.35", va="bottom")

# ---- Panel R: logistic jet flow m_eff(eps) ----
for i, c in enumerate(cells):
    o = c["directions"]["orth_firstorder_null"]
    ec = o["eps_cross"]
    y = (eps / ec)**2
    m_eff = (1.0 + 2.0 * y) / (1.0 + y)
    lab = rf"$z_1{{=}}{c['z1_mm']:.0f},z_2{{=}}{c['z2_mm']:.0f}$"
    axR.semilogx(eps, m_eff, lw=1.4, alpha=0.9, label=lab)
axR.axhline(1.0, color="0.6", lw=0.8, ls=":")
axR.axhline(2.0, color="0.6", lw=0.8, ls=":")
axR.text(3.5e-3, 1.03, "$m{=}1$ leak-dominated fixed pt", fontsize=5.8, color="0.4")
axR.text(3.5e-3, 1.90, "$m{=}2$ curvature fixed pt", fontsize=5.8, color="0.4")
axR.axhline(1.5, color=C_GREY, lw=0.6, ls=(0, (2, 2)))
axR.set_ylim(0.9, 2.1)
axR.set_xlabel(r"perturbation amplitude  $\varepsilon$")
axR.set_ylabel(r"effective contact order  $m_{\rm eff}$")
axR.set_title(r"Exact jet flow  $m_{\rm eff}=\dfrac{1+2y}{1+y}$", fontsize=7.2, pad=3)
axR.legend(loc="center right", fontsize=5.6, handlelength=1.4, labelspacing=0.25)
# projected-slope callout (all four cells = 4.00)
proj = kt2["verdict"]["projected_cov_slopes"]
unp = kt2["verdict"]["unprojected_lowEps_slopes"]
axR.text(0.03, 0.045,
         f"projected cov slope = {proj[0]:.2f} (all 4 cells)\n"
         f"unprojected low-$\\varepsilon$ slope = {min(unp):.2f}–{max(unp):.2f}",
         transform=axR.transAxes, fontsize=5.8, color=C_JET, va="bottom", fontweight="bold")

save(fig, os.path.join(HERE, "fig2_jet_transmutation"))
