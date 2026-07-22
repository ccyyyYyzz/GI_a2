"""
T2 -- FISHER / NOISE  (the effect size in photons).

Correct minimum-variance answer = Cramer-Rao bound (MLE, nuisances p1,p3,p4
marginalised), NOT an OLS plug-in.

Observation model (self-averaged / annealed): at sweep level l we observe
hatQ_l = 1 - Sbar_l/C with
    E[hatQ_l] = Q_l(theta),   Var(hatQ_l) = (Q_l - Q2_l)/(C M_l),
where theta = (p1,p2,p3,p4) are the power sums, Q_l = exp(-L_l),
L_l = sum_k (-1)^{k+1} a_l^k p_k /k, and Q2_l = Q(2 a_l).
Sensitivity: dQ_l/dp_k = Q_l * (-1)^k a_l^k / k.
Fisher (Gaussian): I_{jk} = sum_l M_l dQ_l/dp_j dQ_l/dp_k / sigma_l^2 (sigma at M=1).
CRB: Cov(theta_hat) >= I^{-1};  se(p2) = sqrt((I^{-1})_{22}).

RESOURCE (budget): total INCIDENT photons per pattern over the sweep,
    B_sat = sum_l M_l * rho_l * p1 = M * C * sum_l occ_l    (flat M; occ_l=a_l p1).
LINEAR ARM reference: ordinary bucket estimates p1 with rel se 1/sqrt(B) in the
linear (non-saturating) regime -> B_LIN photons/pattern for rel se sqrt-scale.

HEADLINE LAW: nonlinear/linear signal fraction in L(a) is
    [a^2 p2/2] / [a p1] = (occ/2) * p2/p1^2 = (occ/2)/n_eff,
    n_eff = (sum v)^2 / sum v^2  (participation ratio).
=> Fisher(p2) per photon ~ 1/n_eff^2 => B_sat(target) ~ n_eff^2.  The nonlinear
functional's effect size VANISHES precisely when p2 is nontrivial (many pixels).

BAR: rel se(p2) <= 10% at per-pattern budget <= 4x the linear arm's.
"""
import json
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import saturation_core as sc
import gi_operator as op

OUT = os.path.dirname(HERE)
DATA = "data/r63_bridge_scenes"
rng = np.random.default_rng(651_303)

C = 3600
L = 12
ORDER = 4
TARGET = 0.10                    # rel se(p2) target
B_LIN = 1.0e4                    # linear-arm incident photons/pattern (1% p1)
SCENES = ["bridge_contour_0", "bridge_twopop_1", "bridge_microtex_2",
          "bridge_control_0", "bridge_contour_2", "bridge_twopop_3"]


def crb_se_p2(a_lv, v, C, M_per_level, order=ORDER):
    """CRB se(p2) with flat M per level. Returns (se_p2, Q, dico)."""
    v = np.asarray(v, float); v = v[v > 0]
    Q = sc.Q_exp(a_lv, v)
    Q2 = sc.Q_exp(2 * a_lv, v)
    sig2 = (Q - Q2) / C                                # per-window var of hatQ
    dQ = np.stack([Q * ((-1.0) ** k) * a_lv ** k / k
                   for k in range(1, order + 1)], axis=1)   # (L, order)
    I = np.zeros((order, order))
    Ml = np.atleast_1d(M_per_level) * np.ones(len(a_lv))
    for l in range(len(a_lv)):
        I += Ml[l] * np.outer(dQ[l], dQ[l]) / sig2[l]
    cov = np.linalg.inv(I)
    return float(np.sqrt(cov[1, 1])), Q, dict(cond=float(np.linalg.cond(I)))


report = {"C": C, "L": L, "order": ORDER, "target_rel_se": TARGET,
          "B_LIN": B_LIN, "scenes": {}}
mults = []
print("scene                dens  n_eff   p1      p2     relse@M1   M_need(10%)"
      "   B_sat        x_linear")

# primary design: 50%-density masks (canonical GI patterns) + density sweep
for name in SCENES:
    x = op.load_scene(os.path.join(DATA, name + ".npz"))
    rng_m = np.random.default_rng(651_400)
    m = (rng_m.random(op.P) < 0.5).astype(float)          # canonical 50% mask
    v = m * x
    vp = v[v > 0]
    p = sc.power_sums(vp, 4); p1, p2 = float(p[1]), float(p[2])
    n_eff = p1 ** 2 / p2
    a_lv = sc.default_a_levels(vp, L=L)
    se1, Q, dico = crb_se_p2(a_lv, vp, C, 1.0)
    rel1 = se1 / p2
    M_need = (rel1 / TARGET) ** 2
    sum_occ = float(np.sum(a_lv * p1))
    B_sat = M_need * C * sum_occ
    mult = B_sat / B_LIN
    mults.append(mult)
    report["scenes"][name] = dict(density=0.5, n_active=int(vp.size),
        n_eff=float(n_eff), p1=p1, p2=p2, relse_at_M1=float(rel1),
        M_need=float(M_need), B_sat=float(B_sat), budget_multiple=float(mult),
        fisher_cond=dico["cond"])
    print("%-20s 0.50 %6.1f %7.1f %7.1f  %8.2e  %9.2e  %9.2e  %8.2e"
          % (name, n_eff, p1, p2, rel1, M_need, B_sat, mult))

# density-scaling law (headline): B_sat ~ n_eff^2 on one scene
print("\nDENSITY / n_eff^2 SCALING LAW (bridge_contour_2):")
x = op.load_scene(os.path.join(DATA, "bridge_contour_2.npz"))
law = {}
for d in (0.5, 0.2, 0.1, 0.05, 0.02):
    rng_m = np.random.default_rng(651_500)
    m = (rng_m.random(op.P) < d).astype(float)
    v = (m * x); vp = v[v > 0]
    if vp.size < 5:
        continue
    p = sc.power_sums(vp, 4); p1, p2 = float(p[1]), float(p[2])
    n_eff = p1 ** 2 / p2
    a_lv = sc.default_a_levels(vp, L=L)
    se1, _, _ = crb_se_p2(a_lv, vp, C, 1.0)
    M_need = (se1 / p2 / TARGET) ** 2
    B_sat = M_need * C * float(np.sum(a_lv * p1))
    law[f"d{d}"] = dict(n_active=int(vp.size), n_eff=float(n_eff),
        B_sat=float(B_sat), mult=float(B_sat / B_LIN),
        B_over_neff2=float(B_sat / n_eff ** 2))
    print("  d=%.2f  n_active=%4d  n_eff=%6.1f  B_sat=%.3e (%.2e x lin)  "
          "B_sat/n_eff^2=%.3e" % (d, vp.size, n_eff, B_sat, B_sat / B_LIN,
                                  B_sat / n_eff ** 2))
report["neff2_scaling_law"] = law
consts = [law[k]["B_over_neff2"] for k in law]
report["neff2_law_const_range"] = [float(min(consts)), float(max(consts))]
print("  B_sat/n_eff^2 in [%.2e, %.2e]  (approx constant => B_sat ~ n_eff^2)"
      % (min(consts), max(consts)))

# MC validation of the CRB at an ACHIEVABLE M on the sparsest design
print("\nMC VALIDATION of CRB (d=0.02, achievable M):")
rng_m = np.random.default_rng(651_500)
m = (rng_m.random(op.P) < 0.02).astype(float)
v = (m * x); vp = v[v > 0]
p2t = float(sc.power_sums(vp, 2)[2]); p1t = float(sc.power_sums(vp, 2)[1])
a_lv = sc.default_a_levels(vp, L=L)
M_mc = 20000
se_crb, _, _ = crb_se_p2(a_lv, vp, C, float(M_mc))
Qm = sc.Q_exp(a_lv, vp); Q2m = sc.Q_exp(2 * a_lv, vp)
wgt = 1.0 / ((Qm - Q2m) / (C * M_mc) / Qm ** 2)          # GLS weights -> CRB
NMC = 400
p2hat = np.empty(NMC)
for i in range(NMC):
    W = sc.draw_W(C, len(vp), rng)                       # fresh W (annealed)
    S, f = sc.simulate_bucket(a_lv, vp, W, M_mc, rng, fast=True)
    Qh = np.clip(1 - S.mean(axis=1) / C, 1e-9, 1 - 1e-9)
    p2hat[i] = sc.fit_powersums(a_lv, Qh, order=ORDER, weights=wgt)["p2"]
se_mc = p2hat.std(ddof=1)
print("  CRB se(p2)=%.4f  MC se(p2)=%.4f  ratio MC/CRB=%.3f  (rel bias=%.3f)"
      % (se_crb, se_mc, se_mc / se_crb, (p2hat.mean() - p2t) / p2t))
# tiny EXACT-Bernoulli cross-check confirms the fast Gaussian draw is faithful
NX, MX = 120, 800
sxg = np.empty(NX); sxb = np.empty(NX)
for i in range(NX):
    W = sc.draw_W(C, len(vp), rng)
    Sg, _ = sc.simulate_bucket(a_lv, vp, W, MX, rng, fast=True)
    Sb, _ = sc.simulate_bucket(a_lv, vp, W, MX, rng, fast=False)
    sxg[i] = np.clip(1 - Sg.mean(1) / C, 1e-9, 1 - 1e-9)[-1]
    sxb[i] = np.clip(1 - Sb.mean(1) / C, 1e-9, 1 - 1e-9)[-1]
print("  exact-vs-fast top-level survival sd: exact=%.4e fast=%.4e (ratio %.3f)"
      % (sxb.std(), sxg.std(), sxg.std() / sxb.std()))
report["mc_validation"] = dict(density=0.02, M=M_mc, se_crb=float(se_crb),
    se_mc=float(se_mc), ratio=float(se_mc / se_crb),
    rel_bias=float((p2hat.mean() - p2t) / p2t), NMC=NMC,
    exact_vs_fast_sd_ratio=float(sxg.std() / sxb.std()))

report["summary"] = dict(
    primary_density=0.5,
    max_budget_multiple_50pct=float(np.max(mults)),
    min_budget_multiple_50pct=float(np.min(mults)),
    best_case_sparse_multiple=float(min(law[k]["mult"] for k in law)),
    BAR="rel se(p2) <= 10% at per-pattern budget <= 4x linear arm",
    PASS=bool(np.min(mults) <= 4.0),
    verdict_note=("per-pattern p2 recovery costs ~n_eff^2 x the linear arm; "
                  "even sparsest sane design is >1e3 x -> BAR NOT MET"))
print("\n50%%-mask budget multiple: min=%.2e max=%.2e   sparsest(d=0.02)=%.2e"
      "   BAR<=4x -> %s"
      % (report["summary"]["min_budget_multiple_50pct"],
         report["summary"]["max_budget_multiple_50pct"],
         report["summary"]["best_case_sparse_multiple"],
         "PASS" if report["summary"]["PASS"] else "FAIL"))

with open(os.path.join(OUT, "t2_fisher_noise.json"), "w") as fh:
    json.dump(report, fh, indent=2, default=float)
print("wrote t2_fisher_noise.json")
