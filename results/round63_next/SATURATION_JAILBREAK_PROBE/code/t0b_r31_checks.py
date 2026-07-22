"""
R31-ruling unit checks (fold into T0): exact finite-C covariance (S2.1), per-gate
variance decomposition (S2.2), per-photon information ridge t* table (S4.1/S4.4),
and the annealed FIM closed form (S3.2) vs the ideal I_{beta2beta2}=Ct^4/(4(e^t-1)).
"""
import json
import os
import sys
import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import saturation_core as sc

OUT = os.path.dirname(HERE)
rng = np.random.default_rng(650_031)
rep = {}

# --- (1) exact finite-C covariance  Cov[r_C(ai),r_C(aj)]=(r(ai+aj)-r ri rj)/C --
print("(1) exact finite-C covariance (ruling S2.1) vs MC:")
v = np.array([0.3, 0.7, 0.5, 0.9, 0.2])
C = 4000
a1, a2 = 0.8, 1.3
NR = 4000
r1 = np.empty(NR); r2 = np.empty(NR)
for i in range(NR):
    W = sc.draw_W(C, len(v), rng)
    r1[i] = sc.empirical_survival(np.array([a1]), v, W)[0]
    r2[i] = sc.empirical_survival(np.array([a2]), v, W)[0]
cov_mc = np.cov(r1, r2)[0, 1]
cov_ana = sc.r_cov_exact(a1, a2, v, C)
print("   Cov[r(a1),r(a2)]  MC=%.4e  analytic=%.4e  ratio=%.3f"
      % (cov_mc, cov_ana, cov_mc / cov_ana))
# relative sd form
var_mc = r1.var(); rr = sc.Q_exp(np.array([a1]), v)[0]; rr2 = sc.Q_exp(np.array([2*a1]), v)[0]
relsd_ana = np.sqrt(rr2 / rr**2 - 1) / np.sqrt(C)
print("   sd[r(a1)]/r  MC=%.4e  analytic=%.4e" % (np.sqrt(var_mc)/rr, relsd_ana))
rep["covariance_check"] = dict(cov_mc=float(cov_mc), cov_analytic=float(cov_ana),
    ratio=float(cov_mc/cov_ana))

# --- (2) per-gate variance decomposition (S2.2): Var(rhat)= (r(2a)-r^2)/C + (r-r(2a))/(CG)
print("\n(2) per-gate variance decomposition (ruling S2.2) vs MC (fixed W):")
C = 3600; G = 30; a = np.array([1.0])
NR = 3000
W = sc.draw_W(C, len(v), rng)                      # FIXED W
f = sc.cell_fire_prob(a, v, W)[0]
rhats = np.empty(NR)
for i in range(NR):
    S = (rng.random((G, C)) < f[None, :]).sum(1)
    rhats[i] = 1 - S.mean() / C
var_mc = rhats.var()
r = sc.Q_exp(a, v)[0]; r2 = sc.Q_exp(2*a, v)[0]
# NOTE: this fixed-W variance is the per-gate (binomial) term only; the quenched
# term is 0 here because W is fixed and we look at fluctuation across gates.
var_binom = (r - r2) / (C * G) + 0*(r2 - r**2)/C
# the labeled formula's binomial piece for THIS fixed W:
var_pred = (f*(1-f)).sum() / C**2 / G
print("   Var(rhat|fixedW) MC=%.4e  predicted (Poisson-binomial/G)=%.4e"
      % (var_mc, var_pred))
rep["pergate_var_check"] = dict(var_mc=float(var_mc), var_pred=float(var_pred))

# --- (3) per-photon information ridge t* table (S4.1/S4.4) -------------------
print("\n(3) per-photon information ridge  t=(2j-1)(1-e^-t):")
ridge = {}
for j in (2, 3, 4):
    tj = brentq(lambda t: t - (2*j-1)*(1-np.exp(-t)), 1e-6, 40)
    ridge[f"p{j}"] = dict(t_star=float(tj), Cr_at_3600=float(3600*np.exp(-tj)))
    print("   p%d: t*=%.4f  Cr@C=3600=%.1f" % (j, tj, 3600*np.exp(-tj)))
assert abs(ridge["p2"]["t_star"] - 2.821439) < 1e-4
rep["information_ridge"] = ridge

# --- (4) annealed FIM (S3.2) vs ideal I_{beta2beta2}=Ct^4/(4(e^t-1)) ---------
print("\n(4) annealed FIM (S3.2) vs ideal I_beta2 = C t^4/(4(e^t-1)) at t=3.54:")
# near-flat scene so p2/p1^2 small; compare I_{22} (p2, known p1) to ideal via chain rule
n = 400; vflat = np.full(n, 0.5); p1 = vflat.sum()
t = 3.54; a_lv = np.array([t / p1]); G = np.array([1.0])
I = sc.annealed_fim(a_lv, G, vflat, C, order=2)
# I_beta2 = p1^4 * I_p2 (beta2=p2/p1^2 -> dp2 = p1^2 dbeta2), known-p1 diagonal
I_beta2_fim = p1**4 * I[1, 1]
I_beta2_ideal = C * t**4 / (4*(np.exp(t) - 1))
print("   I_beta2 (from FIM)=%.1f   ideal=%.1f   ratio=%.3f"
      % (I_beta2_fim, I_beta2_ideal, I_beta2_fim / I_beta2_ideal))
rep["fim_vs_ideal"] = dict(I_beta2_fim=float(I_beta2_fim),
    I_beta2_ideal=float(I_beta2_ideal), ratio=float(I_beta2_fim/I_beta2_ideal))

with open(os.path.join(OUT, "t0b_r31_checks.json"), "w") as fh:
    json.dump(rep, fh, indent=2, default=float)
print("\nR31 unit checks done -> t0b_r31_checks.json")
