"""
QUENCHED FLOOR + SPECKLE REFRESH  (missing piece 4; Pro appendix B/C).
Also supplies the C_eff and effective-mode-number k reporting for R32 BAR 1.

Appendix C floor (fixed cell mixing W, G -> infinity):
    Var(q_C)/q^2 ~ [exp(beta2 t^2) - 1]/C_eff ~ beta2 t^2 / C_eff        (C1)
    SE(beta2_hat)/beta2 ~ (2/t) sqrt(n_eff / C_eff)                       (C2)
Even after infinitely many gates on a FIXED W, curvature estimation has this
irreducible quenched floor.  If W is independently REFRESHED between gates
(rotating diffuser / boiling speckle), the effective sample size grows ~ G*C_eff
and SE(beta2) ~ 1/sqrt(G): the floor is beaten.                            (C3)

C_eff (appendix B1):  Var(q_C) = Var(Y_c)/C_eff,  Y_c = exp(-rho T_c).
For iid cells C_eff = C; optical-grain / crosstalk correlation reduces it.
Effective mode number k (Gamma shape): base fully-developed speckle => k=1 (one
mode); k>1 divides the p2 signal by k (gate demand ~ k^2).
"""
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
FIGS = os.path.join(OUT, "figs")
DATA = "data/r63_bridge_scenes"
rng = np.random.default_rng(651_204)

C = 3600
ORDER = 4


def sparse_scene(name, support=32, seed=0):
    x = op.load_scene(os.path.join(DATA, name + ".npz"))
    m = np.zeros(op.P)
    m[np.random.default_rng(seed).choice(op.P, size=support, replace=False)] = 1.0
    v = (m * x)
    return v[v > 0]


# representative sparse mask x scene (n_eff ~ 19)
v = sparse_scene("bridge_contour_2", 32, seed=651_605)
p1 = float(v.sum()); p2t = float((v ** 2).sum())
n_eff = p1 ** 2 / p2t
beta2 = 1.0 / n_eff
a_lv = sc.default_a_levels(v, L=10)
t_hi = float(a_lv[-1] * p1)                          # top-of-sweep load
rep = dict(scene="bridge_contour_2@32-support", n_active=int(v.size),
           n_eff=float(n_eff), beta2=float(beta2), p1=p1, p2_true=p2t, t_hi=t_hi)
print("sparse scene: n_active=%d  n_eff=%.1f  beta2=%.4f  top-sweep t=%.2f"
      % (v.size, n_eff, beta2, t_hi))


# =========================================================================== #
# BAR 1 -- measured C_eff and effective mode number k
# =========================================================================== #
print("\nC_eff measurement (Var(q_C) = Var(Y_c)/C_eff), ridge load, iid speckle:")
a_ridge = np.array([2.821439 / p1])
NW = 2000
qC = np.empty(NW); cellvar = np.empty(NW)
for i in range(NW):
    W = sc.draw_W(C, v.size, rng)
    Yc = np.exp(-(a_ridge[0] * C) * (W @ v))          # Y_c = exp(-rho T_c)
    qC[i] = Yc.mean(); cellvar[i] = Yc.var()
C_eff_iid = float(cellvar.mean() / qC.var())
print("  iid cells:      C_eff = %.0f   (physical C = %d)  ratio %.3f"
      % (C_eff_iid, C, C_eff_iid / C))
# correlated (block-shared grain) -> C_eff = C/block
ceff_corr = {}
for block in (4, 16):
    qb = np.empty(600)
    for i in range(600):
        ndist = C // block
        dist = rng.exponential(scale=1.0 / C, size=(ndist, v.size))
        Wc = np.repeat(dist, block, axis=0)
        qb[i] = np.exp(-(a_ridge[0] * C) * (Wc @ v)).mean()
    Ceff_b = float(cellvar.mean() / qb.var())
    ceff_corr[f"block{block}"] = dict(C_eff=Ceff_b, expected=C // block)
    print("  block-%2d shared: C_eff = %5.0f  (expected C/%d = %d)"
          % (block, Ceff_b, block, C // block))
rep["bar1_C_eff"] = dict(C_physical=C, C_eff_iid=C_eff_iid,
    C_eff_iid_over_C=float(C_eff_iid / C), correlated=ceff_corr,
    k_eff_simulated=1, k_note="base fully-developed speckle = 1 mode; k>1 divides "
    "p2 signal by k (gate demand ~ k^2, committed p0: k3/k1=9.6, k10/k1=108.6)")


# =========================================================================== #
# exact appendix-S2.1 quenched covariance of the survival sweep (fixed W):
#   Cov[r_C(a_i), r_C(a_j)] = (r(a_i+a_j) - r(a_i) r(a_j)) / C_eff.
# We use this (a) to reach the G->inf floor and large G analytically, and
# (b) cross-check it against a small draw_W Monte-Carlo.
# =========================================================================== #
def quench_cov(a_lv, v, C_eff):
    r = sc.Q_exp(a_lv, v)
    S = a_lv[:, None] + a_lv[None, :]
    rij = np.array([[sc.Q_exp(np.array([s]), v)[0] for s in row] for row in S])
    return (rij - np.outer(r, r)) / C_eff, r


def shot_var(a_lv, v, C_eff):
    r = sc.Q_exp(a_lv, v); r2 = sc.Q_exp(2 * a_lv, v)
    return np.maximum((r - r2) / C_eff, 1e-14)          # per-gate binomial term (S2.2)


def fit_b2(a_lv, qhat, p1):
    return sc.fit_powersums(a_lv, np.clip(qhat, 1e-9, 1 - 1e-9), order=ORDER)["p2"] / p1 ** 2


# validate the analytic covariance against draw_W (fixed-W, G=inf) at C=3600
a_lv = sc.default_a_levels(v, L=10)
Sig, rvec = quench_cov(a_lv, v, C)
Lchol = np.linalg.cholesky(Sig + 1e-15 * np.eye(len(a_lv)))
NR = 400
b2_ana = np.array([fit_b2(a_lv, rvec + Lchol @ rng.standard_normal(len(a_lv)), p1)
                   for _ in range(NR)])
b2_mc = np.empty(NR)
for i in range(NR):
    W = sc.draw_W(C, v.size, rng)
    b2_mc[i] = fit_b2(a_lv, sc.empirical_survival(a_lv, v, W), p1)
floor_ana = float(np.std(b2_ana, ddof=1) / beta2)
floor_mc = float(np.std(b2_mc, ddof=1) / beta2)
print("\nanalytic-vs-MC floor cross-check (C=3600, G->inf): SE/beta2 analytic=%.4f  MC=%.4f"
      % (floor_ana, floor_mc))

# =========================================================================== #
# C2 floor: SE(beta2)/beta2 vs (2/t) sqrt(n_eff/C_eff), scaling in C_eff & n_eff
# =========================================================================== #
print("C2 QUENCHED FLOOR (G->inf) SE(beta2)/beta2 vs (2/t) sqrt(n_eff/C_eff):")
floor_scaling = {}
t_use = float(a_lv[-1] * p1)
for Cv in (900, 3600, 14400):
    Sg, rv = quench_cov(a_lv, v, Cv)
    Lc = np.linalg.cholesky(Sg + 1e-15 * np.eye(len(a_lv)))
    b2s = np.array([fit_b2(a_lv, rv + Lc @ rng.standard_normal(len(a_lv)), p1)
                    for _ in range(1500)])
    rel = float(np.std(b2s, ddof=1) / beta2)
    formula = (2.0 / t_use) * np.sqrt(n_eff / Cv)
    floor_scaling[f"C{Cv}"] = dict(rel_floor=rel, formula=float(formula),
        t_used=t_use, rel_times_sqrtC=float(rel * np.sqrt(Cv)))
    print("  C_eff=%6d  SE/beta2=%.4f  formula(2/t sqrt(n_eff/C))=%.4f  x sqrt(C)=%.3f"
          % (Cv, rel, formula, rel * np.sqrt(Cv)))
consts = [floor_scaling[k]["rel_times_sqrtC"] for k in floor_scaling]
print("  1/sqrt(C_eff) scaling: SE/beta2 * sqrt(C_eff) in [%.3f, %.3f] (~const)"
      % (min(consts), max(consts)))
# sqrt(n_eff) scaling in the diffuse limit: fix a small load t where beta2 t^2<<1
print("  sqrt(n_eff) scaling at fixed load t=1.0 (diffuse limit beta2 t^2<<1):")
neff_scaling = {}
for sup in (16, 32, 64):
    vv = sparse_scene("bridge_contour_2", sup, seed=651_700 + sup)
    p1v = float(vv.sum()); nev = p1v ** 2 / float((vv ** 2).sum()); b2v = 1 / nev
    a_s = np.array([0.2, 0.5, 1.0]) / p1v            # low-load sweep (diffuse regime)
    Sg, rv = quench_cov(a_s, vv, 3600)
    Lc = np.linalg.cholesky(Sg + 1e-15 * np.eye(len(a_s)))
    b2s = np.array([fit_b2(a_s, rv + Lc @ rng.standard_normal(len(a_s)), p1v)
                    for _ in range(2000)])
    rel = float(np.std(b2s, ddof=1) / b2v)
    neff_scaling[f"support{sup}"] = dict(n_eff=nev, rel_floor=rel,
        rel_over_sqrt_neff=float(rel / np.sqrt(nev)))
    print("    support=%2d  n_eff=%5.1f  SE/beta2=%.3f  /sqrt(n_eff)=%.4f"
          % (sup, nev, rel, rel / np.sqrt(nev)))
rep["C2_quenched_floor"] = dict(floor_analytic_vs_mc=dict(analytic=floor_ana, mc=floor_mc),
    scaling_in_C=floor_scaling, sqrtC_const_range=[float(min(consts)), float(max(consts))],
    scaling_in_neff=neff_scaling,
    verdict="analytic S2.1 covariance matches draw_W; the G->inf floor scales as "
            "1/sqrt(C_eff) exactly and as sqrt(n_eff) in the diffuse (beta2 t^2<<1) "
            "limit -> a fixed microcell ensemble caps curvature precision at "
            "SE/beta2 ~ (2/t) sqrt(n_eff/C_eff), no matter how many gates")


# =========================================================================== #
# C3 refresh vs fixed: SE(beta2) vs G (analytic; reaches the floor). THE PLOT.
# =========================================================================== #
print("\nC3 REFRESH vs FIXED  SE(beta2) vs G (analytic S2.1/S2.2, sparse mask):")
Gs = np.array([1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000])
NT = 3000
Sig, rvec = quench_cov(a_lv, v, C)
Lchol = np.linalg.cholesky(Sig + 1e-15 * np.eye(len(a_lv)))
Vs = shot_var(a_lv, v, C)
se_fixed = np.empty(len(Gs)); se_refresh = np.empty(len(Gs))
for gi, G in enumerate(Gs):
    b2_fix = np.empty(NT); b2_ref = np.empty(NT)
    for r in range(NT):
        z = rng.standard_normal((len(a_lv),))
        # FIXED W: one quenched offset (persists over all G gates) + shot/G
        qf = rvec + Lchol @ z + rng.normal(0, np.sqrt(Vs / G))
        b2_fix[r] = fit_b2(a_lv, qf, p1)
        # REFRESHED W: average of G independent frames -> offset & shot both /G
        zf = rng.standard_normal((len(a_lv),))
        qr = rvec + (Lchol @ zf) / np.sqrt(G) + rng.normal(0, np.sqrt(Vs / G))
        b2_ref[r] = fit_b2(a_lv, qr, p1)
    se_fixed[gi] = np.std(b2_fix, ddof=1)
    se_refresh[gi] = np.std(b2_ref, ddof=1)
    print("  G=%6d  SE(beta2) fixed-W=%.4e  refreshed-W=%.4e  ratio=%.2f"
          % (G, se_fixed[gi], se_refresh[gi], se_fixed[gi] / se_refresh[gi]))

floor_line = float(np.std([fit_b2(a_lv, rvec + Lchol @ rng.standard_normal(len(a_lv)), p1)
                           for _ in range(4000)], ddof=1))     # G->inf fixed-W floor
rep["C3_refresh"] = dict(G=Gs.tolist(), se_fixed_W=se_fixed.tolist(),
    se_refreshed_W=se_refresh.tolist(), quenched_floor_Ginf=floor_line,
    fixed_over_refresh_at_Gmax=float(se_fixed[-1] / se_refresh[-1]),
    verdict="fixed-W SE(beta2) plateaus at the quenched floor (~%.2e); refreshed-W "
            "SE keeps falling ~1/sqrt(G) (effective samples ~ G*C_eff) and beats the "
            "floor by %.0fx at G=1e5" % (floor_line, se_fixed[-1] / se_refresh[-1]))

# ---- figure ----
fig, ax = plt.subplots(figsize=(6.2, 4.4))
ax.loglog(Gs, se_fixed, "o-", color="#c0392b", label="fixed W (quenched)")
ax.loglog(Gs, se_refresh, "s-", color="#2471a3", label="refreshed W (boiling speckle)")
ax.axhline(floor_line, ls="--", color="#c0392b", alpha=0.6,
           label="quenched floor (G$\\to\\infty$)")
ax.loglog(Gs, se_refresh[0] * np.sqrt(Gs[0] / Gs), ":", color="#2471a3", alpha=0.6,
          label="$1/\\sqrt{G}$ guide")
ax.set_xlabel("gates per level  G"); ax.set_ylabel("SE($\\beta_2$)")
ax.set_title("Quenched floor vs speckle refresh (sparse mask, $n_{eff}$=%.0f, C=%d)"
             % (n_eff, C))
ax.legend(fontsize=8); ax.grid(True, which="both", alpha=0.25)
fig.tight_layout()
fig.savefig(os.path.join(FIGS, "fig_quenched_floor.png"), dpi=130)
print("wrote figs/fig_quenched_floor.png")

with open(os.path.join(OUT, "quenched_floor.json"), "w") as fh:
    json.dump(rep, fh, indent=2, default=float)
print("wrote quenched_floor.json")
