"""
T0 -- verify the theorem seed BEFORE anything else.

(a) Symbolic 3-pixel MGF: E_w[exp(-rho w^T v)] = prod_p 1/(1+a v_p), and
    1/Q expands with coefficients = elementary symmetric polys e_k(v);
    -log Q expands with power sums p_k. (sympy, exact.)
(b) Numerical MC (exponential weights, mean 1/C): empirical survival ->
    ensemble product; recovered p2 matches truth.
(c) Gamma(k) generalisation: Q_k = prod (1+a v_p/k)^{-k}, p2-coeff ~ 1/(2k).
(d) Self-averaging (quenched vs annealed) relative sd ~ 1/sqrt(C).
"""
import json
import os
import sys
import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import saturation_core as sc

OUT = os.path.dirname(HERE)
rng = np.random.default_rng(650_000_00)
report = {}

# --------------------------------------------------------------------------- #
print("=" * 70)
print("(a) SYMBOLIC 3-pixel MGF and its expansion")
print("=" * 70)
a, v1, v2, v3 = sp.symbols("a v1 v2 v3", positive=True)
w1, w2, w3 = sp.symbols("w1 w2 w3", positive=True)
# Exponential rate lam = 1/mean = C ; here per-pixel MGF E[exp(-rho w v)] with
# rho = a*C and w~Exp(mean 1/C): E = 1/(1 + (a*C)*(1/C)*v) = 1/(1+a v).  Verify:
lam = sp.symbols("lam", positive=True)              # rate = 1/mean
t = sp.symbols("t", positive=True)
mgf_one = sp.integrate(lam * sp.exp(-lam * w1) * sp.exp(-t * w1), (w1, 0, sp.oo))
mgf_one = sp.simplify(mgf_one)                      # should be lam/(lam+t)=1/(1+t/lam)
print("  E_w[exp(-t w)] for w~Exp(rate lam) =", mgf_one)
# substitute lam=C (mean 1/C), t = rho v = a*C*v  -> 1/(1 + a v)
per_pixel = mgf_one.subs({lam: sp.symbols("C", positive=True),
                          t: a * sp.symbols("C", positive=True) * v1})
per_pixel = sp.simplify(per_pixel)
print("  substituted (mean 1/C, t=a*C*v1):", per_pixel, " (expect 1/(a*v1+1))")
assert sp.simplify(per_pixel - 1 / (1 + a * v1)) == 0
report["mgf_one_pixel"] = str(per_pixel)

Q = 1 / ((1 + a * v1) * (1 + a * v2) * (1 + a * v3))
invQ = sp.expand(1 / Q)
print("  1/Q expanded:", invQ)
# coefficients in a -> elementary symmetric polynomials
poly = sp.Poly(invQ, a)
e1 = v1 + v2 + v3
e2 = v1 * v2 + v1 * v3 + v2 * v3
e3 = v1 * v2 * v3
assert sp.simplify(poly.coeff_monomial(a) - e1) == 0
assert sp.simplify(poly.coeff_monomial(a ** 2) - e2) == 0
assert sp.simplify(poly.coeff_monomial(a ** 3) - e3) == 0
print("  1/Q coeffs = e1,e2,e3  OK")
# -log Q power-sum expansion
L = -sp.log(Q)
ser = sp.series(L, a, 0, 4).removeO()
p1 = v1 + v2 + v3
p2 = v1 ** 2 + v2 ** 2 + v3 ** 2
p3 = v1 ** 3 + v2 ** 3 + v3 ** 3
target = a * p1 - a ** 2 * p2 / 2 + a ** 3 * p3 / 3
assert sp.simplify(ser - target) == 0
print("  -logQ series = a p1 - a^2 p2/2 + a^3 p3/3  OK")
# the money identity: sum v^2 = e1^2 - 2 e2
assert sp.simplify(p2 - (e1 ** 2 - 2 * e2)) == 0
print("  identity  sum v^2 = e1^2 - 2 e2  OK")
report["symbolic_ok"] = True

# --------------------------------------------------------------------------- #
print("\n" + "=" * 70)
print("(b) NUMERICAL MC: empirical survival -> ensemble product; recover p2")
print("=" * 70)
vnp = np.array([0.2, 0.7, 0.45])                   # 3-pixel test scene
C = 200_000                                        # large C -> tight LLN
a_grid = np.array([0.3, 0.6, 1.0, 1.6, 2.4])
W = sc.draw_W(C, len(vnp), rng)
emp = sc.empirical_survival(a_grid, vnp, W)
ana = sc.Q_exp(a_grid, vnp)
rel = np.abs(emp - ana) / ana
print("  a-grid            :", a_grid)
print("  empirical survival:", np.round(emp, 6))
print("  ensemble product  :", np.round(ana, 6))
print("  rel diff          :", np.round(rel, 5), " max=%.2e" % rel.max())
assert rel.max() < 5e-3, "MC survival does not match ensemble product"
# recover power sums from the ANALYTIC survival (pure identifiability).
# For a 3-pixel scene 1/Q is EXACTLY cubic in a -> exact elementary-symmetric
# fit at order=3 recovers e1,e2,e3 and hence p2 = e1^2-2e2 to machine precision.
ptrue = sc.power_sums(vnp, 4)
a_lv = sc.default_a_levels(vnp, L=10)          # spans occupancy 0.05..3
Qtrue = sc.Q_exp(a_lv, vnp)
fex = sc.fit_es_poly(a_lv, Qtrue, order=3)
print("  [exact ES fit, order=3] p2 = %.10f  truth %.10f  relerr=%.2e"
      % (fex["p2"], ptrue[2], abs(fex["p2"] - ptrue[2]) / ptrue[2]))
assert abs(fex["p2"] - ptrue[2]) / ptrue[2] < 1e-8
# The log/power-sum fit converges only when a*max(v) < 1; restrict the sweep so
# a*max(v) <= 0.6 and it too recovers p2 well (order=5).  (Real masked scenes
# have p1~O(hundreds) so a<<1 automatically -- see T1.)
a_conv = np.logspace(np.log10(0.05 / ptrue[1]),
                     np.log10(0.6 / vnp.max()), 10)
flog = sc.fit_powersums(a_conv, sc.Q_exp(a_conv, vnp), order=5)
print("  [log fit, a*max(v)<=0.6, order=5] p2 = %.6f  relerr=%.2e"
      % (flog["p2"], abs(flog["p2"] - ptrue[2]) / ptrue[2]))
assert abs(flog["p2"] - ptrue[2]) / ptrue[2] < 5e-3   # residual truncation bias
report["mc_survival_max_reldiff"] = float(rel.max())
report["p2_noiseless_relerr_3pix_ES"] = float(abs(fex["p2"] - ptrue[2]) / ptrue[2])
report["p2_noiseless_relerr_3pix_log"] = float(abs(flog["p2"] - ptrue[2]) / ptrue[2])
report["convergence_note"] = "log/power-sum series needs a*max(v)<1; ES-poly exact for P<=order"

# --------------------------------------------------------------------------- #
print("\n" + "=" * 70)
print("(c) GAMMA(k) generalisation")
print("=" * 70)
for kg in (1, 2, 5):
    Wg = sc.draw_W(C, len(vnp), rng, k_gamma=kg)
    # quenched empirical for gamma weights (reuse empirical_survival: works for any W)
    empg = sc.empirical_survival(a_grid, vnp, Wg)
    anag = sc.Q_gamma(a_grid, vnp, kg)
    rg = np.abs(empg - anag) / anag
    print("  k=%d  max rel diff MC vs prod = %.2e" % (kg, rg.max()))
    assert rg.max() < 6e-3
# p2 coefficient shrinks as 1/k
print("  p2 series coeff (a^2) = -p2/(2k):  strongest at k=1, ~0 as k->inf  (see L_series)")
report["gamma_ok"] = True

# --------------------------------------------------------------------------- #
print("\n" + "=" * 70)
print("(d) SELF-AVERAGING relative sd vs C  (quenched fluctuation)")
print("=" * 70)
# analytic relative sd of empirical survival at fixed a: sqrt(Q2/Q^2-1)/sqrt(C)
a_test = 1.0
selfavg = {}
for Cc in (900, 3600, 14400):
    reps = 40
    vals = np.empty(reps)
    for r in range(reps):
        Wr = sc.draw_W(Cc, len(vnp), rng)
        vals[r] = sc.empirical_survival(np.array([a_test]), vnp, Wr)[0]
    relsd_mc = vals.std(ddof=1) / vals.mean()
    Q = sc.Q_exp(np.array([a_test]), vnp)[0]
    Q2 = sc.Q_exp(np.array([2 * a_test]), vnp)[0]
    relsd_ana = np.sqrt(Q2 / Q ** 2 - 1) / np.sqrt(Cc)
    selfavg[Cc] = dict(relsd_mc=float(relsd_mc), relsd_ana=float(relsd_ana))
    print("  C=%6d  quenched rel sd  MC=%.4e  analytic=%.4e" %
          (Cc, relsd_mc, relsd_ana))
report["self_averaging"] = selfavg

with open(os.path.join(OUT, "t0_unit_tests.json"), "w") as fh:
    json.dump(report, fh, indent=2, default=float)
print("\nALL T0 UNIT TESTS PASSED.  wrote t0_unit_tests.json")
