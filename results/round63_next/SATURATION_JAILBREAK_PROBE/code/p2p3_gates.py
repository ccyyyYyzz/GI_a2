"""
GATE P2 (Fisher effect size) + GATE P3 (inverse-space teeth).  Ruling S5, S8.

Key fact for null pairs: x' = x + h with M h = 0 => Delta p1_k = M_k h = 0 for
every mask, so the linear buckets are IDENTICAL and the ONLY changed parameter
is p2_k = M_k(x^{o2}), with Delta p2_k = M_k(2 x o h + h^{o2}).  The
discriminability is therefore purely in the p2 direction, nuisance-profiled
against the (unknown) p1 per mask.

P2 target (ruling eq 5.5, order-2, nuisance-profiled): per null pair
    d'^2 = sum_k (Delta p2_k)^2 * I^{prof}_k(p2),
I^{prof}_k(p2) = Schur complement over p1 of the annealed order-2 FIM at the
three-level design (budget-allocated gates), i.e. p1 (== global power alpha)
profiled out.  Gate: median d' >= 3 over the frozen bank AND >=80% of pairs
>= 2, at the matched (linear-campaign) incident-photon budget.

P3: rank + Fisher-weighted singular spectrum of J=[M; 2M diag(x)] (expect
min(N,2K)=102 rank doubling); ADVERSARIAL worst-case joint-fiber pair (h in
ker[M; M diag(x)] -> first-order Delta p2 = 0) vs random-pair d'.
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
rng = np.random.default_rng(651_909)
C = 3600
LO, HI = 0.08, 0.92
B_LIN = 1.0e4
LIN_CAMPAIGN = op.K * B_LIN                       # matched incident-photon budget

M = op.build_operator()                           # 51 x 1024
side = op.SIDE

# base scene (interior-clamped so null-space adds stay in [0,1] with no clip)
xbase = np.clip(op.load_scene(os.path.join(DATA, "bridge_microtex_0.npz")), LO, HI)

# --- null-space basis of M (for random h in ker M) --------------------------
U, S, Vt = np.linalg.svd(M, full_matrices=True)
kerM = Vt[op.K:].T                                # (N, N-K) orthonormal basis of ker M
print("rank(M) = %d   dim ker(M) = %d" % (np.linalg.matrix_rank(M), kerM.shape[1]))


def valid_pair(h_dir):
    """Scale h_dir (in ker M) so x' = x + a*h stays in [0,1] (no clip); return
    (x, x', h) with a 0.9x-safety amplitude and unit-ish visible scale."""
    h = h_dir / np.abs(h_dir).max()
    pos = h > 1e-12; neg = h < -1e-12
    amax = np.inf
    if pos.any(): amax = min(amax, np.min((HI + 0.08 - xbase[pos]) / h[pos]) if False else np.min((1 - xbase[pos]) / h[pos]))
    if neg.any(): amax = min(amax, np.min((0.0 - xbase[neg]) / h[neg]))
    a = 0.9 * amax
    return xbase, xbase + a * h, a * h


def build_gates(p1):
    """Three-level design at matched per-mask budget = B_LIN incident photons."""
    a_lv, G = sc.three_level_design(p1, budget_incident=B_LIN, C=C)
    return a_lv, G


def dprime_pair(x, xp, budget_per_mask=B_LIN, profile_p1=True):
    """Nuisance-profiled null-pair d' (order-2, three-level design)."""
    d2 = 0.0
    for k in range(op.K):
        m = M[k]
        vx = (m * x)
        vxa = vx[vx > 0]
        p1 = float(vxa.sum())
        if p1 <= 0:
            continue
        dp2 = float(np.sum(m * (xp ** 2 - x ** 2)))       # Delta p2_k (Delta p1=0)
        a_lv, G = sc.three_level_design(p1, budget_incident=budget_per_mask, C=C)
        I = sc.annealed_fim(a_lv, G, vxa, C, order=2)      # 2x2 (p1,p2)
        if profile_p1:
            I_p2 = I[1, 1] - I[1, 0] * I[0, 0] ** -1 * I[0, 1]  # Schur over p1
        else:
            I_p2 = I[1, 1]
        d2 += dp2 ** 2 * I_p2
    return float(np.sqrt(max(d2, 0.0)))


# =========================================================================== #
# GATE P2 -- frozen null-pair bank
# =========================================================================== #
print("\n=== GATE P2: null-pair bank d' (matched linear-campaign budget) ===")
NBANK = 40
dprimes = []
dp2_norms = []
for b in range(NBANK):
    coeff = rng.standard_normal(kerM.shape[1])
    h_dir = kerM @ coeff
    x, xp, h = valid_pair(h_dir)
    dpr = dprime_pair(x, xp)
    dprimes.append(dpr)
    dp2_norms.append(float(np.linalg.norm([np.sum(M[k] * (xp ** 2 - x ** 2))
                                           for k in range(op.K)])))
dprimes = np.array(dprimes)
median_d = float(np.median(dprimes))
frac2 = float(np.mean(dprimes >= 2)); frac3 = float(np.mean(dprimes >= 3))
print("  bank size %d  median d'=%.3f  frac>=2 = %.0f%%  frac>=3 = %.0f%%  "
      "(range %.3f..%.3f)" % (NBANK, median_d, 100 * frac2, 100 * frac3,
                              dprimes.min(), dprimes.max()))
# budget multiple needed for median d'=3 (d' ~ sqrt(budget))
budget_mult_for_d3 = (3.0 / median_d) ** 2
print("  budget multiple for median d'=3: %.1fx the linear campaign" % budget_mult_for_d3)

report = dict(C=C, K=op.K, linear_campaign=LIN_CAMPAIGN, bank_size=NBANK,
    P2=dict(median_dprime=median_d, frac_ge2=frac2, frac_ge3=frac3,
            dprime_min=float(dprimes.min()), dprime_max=float(dprimes.max()),
            budget_mult_for_median_d3=float(budget_mult_for_d3),
            BAR="median d'>=3 AND >=80% pairs >=2 at matched budget",
            PASS=bool(median_d >= 3 and frac2 >= 0.8)))

# Gamma-k degradation of the bank median (p2 signal /k -> d' /k)
gk = {}
for k in (3, 10):
    gk[f"k{k}"] = dict(median_dprime=float(median_d / k),
                       budget_mult_for_d3=float((3.0 / (median_d / k)) ** 2))
report["P2"]["gamma_k_degradation"] = gk
print("  Gamma-k: median d' -> %.3f (k=3), %.3f (k=10)  [signal/k]"
      % (median_d / 3, median_d / 10))

# =========================================================================== #
# GATE P3 -- inverse-space teeth
# =========================================================================== #
print("\n=== GATE P3: rank doubling + adversarial fiber ===")
xr = xbase
J = np.vstack([M, 2 * M * xr[None, :]])            # [M; 2 M diag(x)]  (2K x N)
rankM = int(np.linalg.matrix_rank(M))
rankJ = int(np.linalg.matrix_rank(J))
svJ = np.linalg.svd(J, compute_uv=False)
print("  rank(M)=%d   rank[M;2M diag(x)]=%d  (expect min(N,2K)=%d)"
      % (rankM, rankJ, min(op.P, 2 * op.K)))
print("  J singular spectrum: max=%.3e  min(nonzero,#=%d)=%.3e  cond=%.2e"
      % (svJ[0], rankJ, svJ[rankJ - 1], svJ[0] / svJ[rankJ - 1]))

# adversarial: h in ker([M; M diag(x)]) -> first-order Delta p2 = 0 -> worst case
Jlin = np.vstack([M, M * xr[None, :]])             # [M; M diag(x)]  (2K x N)
Ul, Sl, Vtl = np.linalg.svd(Jlin, full_matrices=True)
rank_lin = int(np.sum(Sl > 1e-9))
ker_joint = Vtl[rank_lin:].T                       # dim N - 2K
print("  dim ker[M; M diag(x)] = %d  (joint first-order fiber = %d%% of N)"
      % (ker_joint.shape[1], round(100 * ker_joint.shape[1] / op.P)))
# worst-case adversarial pair (residual only from h^2 term)
adv_d = []
for b in range(20):
    coeff = rng.standard_normal(ker_joint.shape[1])
    h_dir = ker_joint @ coeff
    x, xp, h = valid_pair(h_dir)
    adv_d.append(dprime_pair(x, xp))
adv_d = np.array(adv_d)
print("  ADVERSARIAL worst-case d' (h in first-order joint fiber): median=%.4f "
      "max=%.4f   vs random-pair median d'=%.3f"
      % (np.median(adv_d), adv_d.max(), median_d))

# restricted-prior injectivity: on a d-dim random support/subspace with 2K>=d,
# is the joint Jacobian full rank d?
print("  restricted-prior injectivity (2K>=d closes the fiber):")
prior_rank = {}
for d in (40, 90, 102, 150):
    B = np.linalg.qr(rng.standard_normal((op.P, d)))[0]    # random d-dim tangent
    Jr = J @ B                                              # (2K x d)
    rk = int(np.linalg.matrix_rank(Jr))
    prior_rank[f"d{d}"] = dict(tangent_dim=d, joint_rank=rk, injective=bool(rk >= d))
    print("    d=%3d  joint rank on tangent=%3d  injective=%s (2K=%d)"
          % (d, rk, rk >= d, 2 * op.K))

report["P3"] = dict(rank_M=rankM, rank_J=rankJ, expected_rank=min(op.P, 2 * op.K),
    rank_doubling=bool(rankJ == min(op.P, 2 * op.K)),
    J_sv_max=float(svJ[0]), J_sv_min_nonzero=float(svJ[rankJ - 1]),
    dim_joint_first_order_fiber=int(ker_joint.shape[1]),
    adversarial_median_dprime=float(np.median(adv_d)),
    adversarial_max_dprime=float(adv_d.max()),
    random_pair_median_dprime=median_d,
    restricted_prior_injectivity=prior_rank,
    honest_claim="fiber HALVING (N-K -> N-2K) + closes low-dim prior tangents "
                 "(2K>=d); NOT dense-scene null-space jailbreak")

with open(os.path.join(OUT, "p2p3_gates.json"), "w") as fh:
    json.dump(report, fh, indent=2, default=float)
print("\nP2 PASS=%s   P3 rank-doubling=%s" % (report["P2"]["PASS"],
                                              report["P3"]["rank_doubling"]))
print("wrote p2p3_gates.json")
