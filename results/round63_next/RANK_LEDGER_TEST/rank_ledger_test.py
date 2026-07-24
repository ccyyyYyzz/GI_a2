#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""RANK_LEDGER_TEST == gate CFG-A (Commutant-Fiber Gate, exact Schur ledger).

ROUND63 divergence synthesis rows #4/#1 (results/round63_next/FIVE_FRONT_DIVERGENCE/
SYNTHESIS.md) EXTENDED by the GPT R47 architecture ruling
(docs/ROUND63_GPT_ROUND47_RULING_RAW.md Q2.2 / Q4.1 CFG-A). CPU, numpy/scipy.
PREREGISTERED: the frozen block (predicted integers + derivation) is written BEFORE compute.

The stable object is the pair (symmetry rep rho, self-adjoint commutant A_G). The paper
carries a RANK LEDGER WITH DISTINCT ROWS, not one scalar 13:
  r_alg  = dim A_G                         (independent admissible quadratic kernels)
  d_base = rank D Phi_G(x) generic         (locally independent quotient coordinates)
  d_fiber= dim V - d_base                  (continuously hidden state fiber = commutant gauge orbit)
Phi_G(x) = (X_alpha X_alpha^dagger)_alpha  (block-Gram map, eq 1.6).

DERIVATION (independent; matches GPT's closed forms exactly).  Z_n regular real rep over m
channels: 2 real isotypic blocks (k=0, k=n/2) + (n/2-1) complex-conjugate blocks, each with
multiplicity m.  Self-adjoint commutant A_G = 2*Herm_m(R) + (n/2-1)*Herm_m(C):
  r_alg (n,m) = 2*m(m+1)/2 + (n/2-1) m^2 = (n/2) m^2 + m.
On the generic nonzero stratum each real block gives a rank-1 real Gram (m base dims) and each
complex block a rank-1 Hermitian Gram (2m-1 base dims), the lost datum being the block's global
complex phase:
  d_base(n,m) = 2m + (n/2-1)(2m-1) = nm - n/2 + 1,
  d_fiber(n,m)= dim V - d_base       = n/2 - 1     (one Fourier phase per complex block).
Free R^24 quadratic cone: trivial group -> full Sym(R^24): r_alg=300, d_base=24, d_fiber=0.

CFG-A BARS (exact integers; only floating-point rank certification tolerance):
  1. r_alg  ranks  == 13 / 26 / 39 / 300  EXACTLY.
  2. generic base-Jacobian ranks d_base == 13 / 19 / 21 / 24 EXACTLY (rank of [G_l x]).
  3. continuous hidden-fiber dims d_fiber == 11 / 5 / 3 / 0: CONSTRUCT equal-base states
     x -> exp(theta J_k) x (per-complex-block Fourier phase), verify EVERY admissible response
     preserved to numerical floor, and rank{J_k x} == d - d_base.
  4. a perturbation transverse to each fiber changes at least one admissible response.
  5. every integer stable across a declared SVD-threshold interval AND independent seeds.
A miss on any exact integer is a REPORTED KILL of the strong complementarity wording (fallback
architecture pre-defined) -- no retuning.  R2 (independent-reseed L-sweep knee at L=r_alg) is a
soft companion bar from synthesis row #4.
"""
import os, json, time
import numpy as np
from scipy.linalg import expm

HERE = os.path.dirname(os.path.abspath(__file__))
D = 24
L_PROBE = 700               # > 300 so the free cone is fully sampled
SV_TOL = 1e-8               # nominal rank cutoff (rank = #{sv > SV_TOL * s_max})
TOL_INTERVAL = [1e-6, 1e-7, 1e-8, 1e-9, 1e-10]   # declared stability interval
SEEDS = [20260724, 11, 2029]                     # independent reseeds
FLOOR = 1e-9               # response-preservation numerical floor (relative to |q| scale)

t0 = time.time()

# ------------------------------------------------------------------ representations
def shift(n):
    P = np.zeros((n, n))
    for i in range(n):
        P[i, (i - 1) % n] = 1.0
    return P

def cyclic_reps(n):
    P = shift(n)
    return [np.linalg.matrix_power(P, g) for g in range(n)]

def tensor_trivial(reps_n, m):
    Im = np.eye(m)
    return [np.kron(r, Im) for r in reps_n]

def reynolds_mats(reps, n_samples, rng):
    """random symmetric matrices projected onto the symmetric commutant (Reynolds)."""
    d = reps[0].shape[0]; G = len(reps)
    mats = []
    for _ in range(n_samples):
        S = rng.standard_normal((d, d)); S = 0.5 * (S + S.T)
        P = np.zeros((d, d))
        for r in reps:
            P += r @ S @ r.T
        P /= G
        mats.append(0.5 * (P + P.T))
    return mats

def rank_at_tol(sv, tol):
    return int((sv > tol * sv.max()).sum())

def vec_rank(mats, tol=SV_TOL):
    iu = np.triu_indices(mats[0].shape[0])
    rows = np.array([Gm[iu] for Gm in mats])
    sv = np.linalg.svd(rows, compute_uv=False)
    return rank_at_tol(sv, tol), sv

def jacobian_rank(mats, x, tol=SV_TOL):
    """rank of D Phi_G(x) = rank span{ G_l x } (rows 2 G_l x of the response Jacobian)."""
    J = np.array([Gm @ x for Gm in mats])          # (L_PROBE, d)
    sv = np.linalg.svd(J, compute_uv=False)
    return rank_at_tol(sv, tol), sv

# ---- Fourier-phase (commutant gauge) generators J_k = (c_k s_k^T - s_k c_k^T) (x) I_m ----
def phase_generators(n, m):
    """skew-symmetric generators of the per-complex-block Fourier phase rotation on R^{n*m}."""
    js = np.arange(n)
    gens = []
    for k in range(1, n // 2):                     # complex blocks k=1..n/2-1
        c = np.sqrt(2.0 / n) * np.cos(2 * np.pi * k * js / n)
        s = np.sqrt(2.0 / n) * np.sin(2 * np.pi * k * js / n)
        Jn = np.outer(c, s) - np.outer(s, c)       # rotate the mode-k (cos,sin) plane
        gens.append(np.kron(Jn, np.eye(m)))
    return gens

# ------------------------------------------------------------------ FROZEN BLOCK (pre-compute)
PRED = {
    "Z24_on_R24":   dict(n=24, m=1, r_alg=13, d_base=13, d_fiber=11),
    "Z12xI2_on_R24": dict(n=12, m=2, r_alg=26, d_base=19, d_fiber=5),
    "Z8xI3_on_R24":  dict(n=8,  m=3, r_alg=39, d_base=21, d_fiber=3),
    "free_on_R24":   dict(n=None, m=None, r_alg=300, d_base=24, d_fiber=0),
}
frozen = {
    "written_before_compute": True,
    "gate": "CFG-A (Commutant-Fiber Gate) per docs/ROUND63_GPT_ROUND47_RULING_RAW.md Q2.2/Q4.1",
    "object": "(rho, self-adjoint commutant A_G); rank ledger with distinct rows r_alg/d_base/d_fiber",
    "formulas": {"r_alg": "(n/2)m^2 + m", "d_base": "nm - n/2 + 1", "d_fiber": "n/2 - 1"},
    "derivation_agrees_with_GPT_closed_form": True,
    "predicted": PRED,
    "bars": {
        "CFG_A_1_r_alg": "13/26/39/300 exact",
        "CFG_A_2_d_base": "13/19/21/24 exact (rank D Phi_G at generic scene)",
        "CFG_A_3_d_fiber": "11/5/3/0 exact; constructed phase orbit preserves all responses to floor",
        "CFG_A_4_transverse": "off-fiber perturbation changes >=1 admissible response",
        "CFG_A_5_stability": "all integers invariant across SVD-tol interval %s and seeds %s" %
                             (TOL_INTERVAL, SEEDS),
        "R2_L_sweep": "independent-reseed accuracy knee at L = r_alg(Z24)=13 (soft companion)",
    },
    "L_PROBE": L_PROBE, "SV_TOL": SV_TOL, "TOL_INTERVAL": TOL_INTERVAL, "SEEDS": SEEDS,
    "d": D, "response_floor": FLOOR,
    "synthesis_reference": {"r_alg": [13, 26, 39, 300], "d_base": [13, 19, 21, 24],
                            "d_fiber": [11, 5, 3, 0]},
}
rep = {"test": "RANK_LEDGER_TEST", "gate": "CFG-A",
       "ref": "SYNTHESIS.md rows #4/#1 + GPT R47 ruling Q2.2/Q4.1", "frozen": frozen, "results": {}}

def family_reps(name):
    p = PRED[name]
    if name == "free_on_R24":
        return [np.eye(D)]
    return tensor_trivial(cyclic_reps(p["n"]), p["m"])

# ================================================================== CFG-A items 1-5
per_family = {}
all_r_alg = all_d_base = all_d_fiber = all_transverse = all_stable = True

for name, pr in PRED.items():
    reps = family_reps(name)
    # --- multi-seed r_alg + threshold-interval stability ---
    r_alg_by_seed = {}
    r_alg_tol_stable = True
    sv_vec_ref = None
    for sd in SEEDS:
        mats = reynolds_mats(reps, L_PROBE, np.random.default_rng(sd))
        rk, sv = vec_rank(mats)
        r_alg_by_seed[sd] = rk
        if sd == SEEDS[0]:
            mats0, sv_vec_ref = mats, sv
        ranks_over_tol = [rank_at_tol(sv, tt) for tt in TOL_INTERVAL]
        if len(set(ranks_over_tol)) != 1 or ranks_over_tol[0] != pr["r_alg"]:
            r_alg_tol_stable = False
    r_alg_ok = (len(set(r_alg_by_seed.values())) == 1 and
                list(r_alg_by_seed.values())[0] == pr["r_alg"] and r_alg_tol_stable)

    # --- d_base: generic base-Jacobian rank, multi-seed x + threshold stability ---
    d_base_by_seed = {}
    d_base_tol_stable = True
    for sd in SEEDS:
        rx = np.random.default_rng(sd + 999).standard_normal(D)     # generic scene
        rk, svJ = jacobian_rank(mats0, rx)
        d_base_by_seed[sd] = rk
        ranks_over_tol = [rank_at_tol(svJ, tt) for tt in TOL_INTERVAL]
        if len(set(ranks_over_tol)) != 1 or ranks_over_tol[0] != pr["d_base"]:
            d_base_tol_stable = False
    d_base_meas = d_base_by_seed[SEEDS[0]]
    d_base_ok = (len(set(d_base_by_seed.values())) == 1 and
                 d_base_meas == pr["d_base"] and d_base_tol_stable)

    # --- d_fiber: (a) nullity d - d_base ; (b) explicit Fourier-phase construction ---
    d_fiber_nullity = D - d_base_meas
    xg = np.random.default_rng(SEEDS[0] + 999).standard_normal(D)
    base_scale = float(np.median([abs(xg @ Gm @ xg) for Gm in mats0]) + 1e-30)
    if name == "free_on_R24":
        gens = []
        fiber_tan_rank = 0
        max_resp_deriv = 0.0
        max_resp_finite = 0.0
    else:
        gens = phase_generators(pr["n"], pr["m"])
        fiber_tan = np.array([Jk @ xg for Jk in gens])              # (n/2-1, d)
        svF = np.linalg.svd(fiber_tan, compute_uv=False)
        fiber_tan_rank = rank_at_tol(svF, SV_TOL)
        # first-order response preservation: d/dtheta q_l = 2 x^T G_l J_k x  ~ 0
        max_resp_deriv = max(abs(2.0 * xg @ Gm @ (Jk @ xg)) / base_scale
                             for Gm in mats0 for Jk in gens)
        # finite-angle equal-base state x' = exp(theta J_k) x preserves ALL responses
        max_resp_finite = 0.0
        for Jk in gens:
            xp = expm(0.7 * Jk) @ xg
            df = max(abs(xp @ Gm @ xp - xg @ Gm @ xg) / base_scale for Gm in mats0)
            max_resp_finite = max(max_resp_finite, df)
    d_fiber_ok = (d_fiber_nullity == pr["d_fiber"] and fiber_tan_rank == pr["d_fiber"]
                  and max_resp_deriv < FLOOR and max_resp_finite < FLOOR)

    # --- transverse perturbation changes >=1 response ---
    rp = np.random.default_rng(SEEDS[0] + 55).standard_normal(D)
    if gens:
        Ft = np.array([Jk @ xg for Jk in gens])
        # project rp onto the complement of the fiber tangent
        Qb, _ = np.linalg.qr(Ft.T)                                  # ortho basis of fiber tangent
        rp_perp = rp - Qb @ (Qb.T @ rp)
    else:
        rp_perp = rp
    rp_perp = rp_perp / (np.linalg.norm(rp_perp) + 1e-30)
    max_transverse_change = max(abs(2.0 * xg @ Gm @ rp_perp) / base_scale for Gm in mats0)
    transverse_ok = max_transverse_change > 1e-3

    all_r_alg &= r_alg_ok; all_d_base &= d_base_ok; all_d_fiber &= d_fiber_ok
    all_transverse &= transverse_ok
    all_stable &= (r_alg_tol_stable and d_base_tol_stable)
    per_family[name] = dict(
        predicted=pr,
        r_alg=dict(measured=r_alg_by_seed, tol_stable=bool(r_alg_tol_stable), exact=bool(r_alg_ok)),
        d_base=dict(measured=d_base_by_seed, tol_stable=bool(d_base_tol_stable), exact=bool(d_base_ok)),
        d_fiber=dict(nullity=d_fiber_nullity, constructed_tangent_rank=int(fiber_tan_rank),
                     max_response_deriv=float(max_resp_deriv),
                     max_response_finite_angle=float(max_resp_finite), exact=bool(d_fiber_ok)),
        transverse=dict(max_response_change=float(max_transverse_change), ok=bool(transverse_ok)),
    )
    print("[CFG-A] %-15s r_alg=%s/%d d_base=%s/%d d_fiber(null=%d,tan=%d)/%d "
          "presv(d=%.1e,f=%.1e) transv=%.2e  %s" %
          (name, list(set(r_alg_by_seed.values())), pr["r_alg"], list(set(d_base_by_seed.values())),
           pr["d_base"], d_fiber_nullity, fiber_tan_rank, pr["d_fiber"], max_resp_deriv,
           max_resp_finite, max_transverse_change,
           "OK" if (r_alg_ok and d_base_ok and d_fiber_ok and transverse_ok) else "MISS"))

rep["results"]["CFG_A"] = dict(
    per_family=per_family,
    item1_r_alg_all_exact=bool(all_r_alg),
    item2_d_base_all_exact=bool(all_d_base),
    item3_d_fiber_all_exact=bool(all_d_fiber),
    item4_transverse_all_ok=bool(all_transverse),
    item5_stability_all_ok=bool(all_stable),
    verdict="PASS" if (all_r_alg and all_d_base and all_d_fiber and all_transverse and all_stable)
            else "FAIL")

# ================================================================== R2: L-sweep knee at r_alg
def sym_circulant(rng, n):
    S = rng.standard_normal((n, n)); S = 0.5 * (S + S.T)
    reps = cyclic_reps(n)
    return 0.5 * (sum(r @ S @ r.T for r in reps) / n + (sum(r @ S @ r.T for r in reps) / n).T)

def ridge_acc(Feat, y, Feat_te, y_te, ridge=1e-2):
    F = Feat - Feat.mean(0, keepdims=True); Fs = F.std(0, keepdims=True) + 1e-9
    F = F / Fs
    Fte = (Feat_te - Feat.mean(0, keepdims=True)) / Fs
    A = F.T @ F + ridge * np.eye(F.shape[1])
    w = np.linalg.solve(A, F.T @ (2.0 * y - 1.0))
    return float(np.mean((Fte @ w > 0).astype(int) == y_te))

def feats(X, banks):
    return np.stack([np.einsum('ni,ij,nj->n', X, G, X) for G in banks], axis=1)

rng_task = np.random.default_rng(77)
Q = sym_circulant(rng_task, D)
Xtr = rng_task.standard_normal((4000, D)); Xte = rng_task.standard_normal((4000, D))
vtr = np.einsum('ni,ij,nj->n', Xtr, Q, Xtr); thr = np.median(vtr)
ytr = (vtr > thr).astype(int)
vte = np.einsum('ni,ij,nj->n', Xte, Q, Xte); yte = (vte > thr).astype(int)
L_grid = list(range(1, 21)); n_reseed = 6; acc_curve = []
for L in L_grid:
    accs = []
    for rs in range(n_reseed):
        rb = np.random.default_rng(1000 + 137 * L + rs)
        banks = [sym_circulant(rb, D) for _ in range(L)]
        accs.append(ridge_acc(feats(Xtr, banks), ytr, feats(Xte, banks), yte))
    acc_curve.append(float(np.mean(accs)))
acc = np.array(acc_curve); eps = 0.003
plateau = acc[12:].mean(); knee = next((L for i, L in enumerate(L_grid) if acc[i] >= plateau - eps), None)
rep["results"]["R2_L_sweep_knee"] = dict(
    L_grid=L_grid, acc_curve=acc_curve, plateau_acc_Lge13=float(plateau), knee_L=knee,
    predicted_knee=13, acc_at_L12=float(acc[11]), acc_at_L13=float(acc[12]), acc_at_L20=float(acc[-1]),
    verdict="PASS" if knee == 13 else ("SOFT_PASS" if (knee and abs(knee - 13) <= 1) else "FAIL"))
print("[R2] knee_L=%s (pred 13) plateau=%.4f acc[12,13,20]=%.4f/%.4f/%.4f" %
      (knee, plateau, acc[11], acc[12], acc[-1]))

# ================================================================== exploratory defect (no bar)
reps24 = cyclic_reps(24)
base_mats = reynolds_mats(reps24, 400, np.random.default_rng(20260724))
iu = np.triu_indices(D)
rows_base = [Gm[iu] for Gm in base_mats]
incr = []
for k in range(5):
    extra = []
    for d0 in range(k):
        Gd = np.zeros((D, D)); Gd[d0, d0] = 1.0; extra.append(Gd[iu])
    rowsk = np.array(rows_base + extra)
    sv = np.linalg.svd(rowsk, compute_uv=False)
    incr.append(rank_at_tol(sv, SV_TOL))
rep["results"]["exploratory_defect_masked_increment"] = dict(
    n_defects=list(range(5)), rank=incr, increment=[incr[i + 1] - incr[i] for i in range(4)],
    note="each masked single-site diagonal-gain defect adds +1 admissible dof (breaks Z_24 translation)")
print("[EXP] defect ranks k=0..4: %s incr %s" % (incr, [incr[i + 1] - incr[i] for i in range(4)]))

# ================================================================== overall
cfg = rep["results"]["CFG_A"]["verdict"]
r2 = rep["results"]["R2_L_sweep_knee"]["verdict"]
overall = "PASS" if (cfg == "PASS" and r2 in ("PASS", "SOFT_PASS")) else \
          ("PARTIAL_r_alg_only" if rep["results"]["CFG_A"]["item1_r_alg_all_exact"] else "FAIL")
rep["verdict"] = dict(CFG_A=cfg, R2=r2, overall=overall,
                      note="strong Noether complementarity holds iff CFG_A==PASS; if only "
                           "item1 lands use narrower 'commutant-rank selection rule' wording")
rep["runtime_sec"] = round(time.time() - t0, 1)
with open(os.path.join(HERE, "RANK_LEDGER_TEST.json"), "w") as f:
    json.dump(rep, f, indent=2)
print("\n=== RANK_LEDGER_TEST / CFG-A overall=%s (CFG_A=%s R2=%s) %.1fs ===" %
      (overall, cfg, r2, rep["runtime_sec"]))
