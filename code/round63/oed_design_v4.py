"""Dead-time-aware OED v4 -- R13 freeze-audit compliance rebuild.

Implements docs/ROUND63_GPT_ROUND13_RULING_RAW.md on top of v1/v2/v3
(ADD-only; v3 kernel primitives are reused, its solver is superseded):

Sec 2   balanced 52-row pre-scan: 32 paired-4x4-block rows (support 32,
        amplitude n/32=32), all 16 8x8-block rows (amplitude 16), 4 quadrant
        rows (amplitude 4).  Identities held to machine precision:
        (1/52) sum_m a_m = 1 (all-ones row average) and per-pixel cumulative
        dose exactly equal (52 amplitude units per pixel).
        V0 from ACTUAL row loads: rho_m = rho_bar * a_m.xhat,
        q_m = a_m/(a_m.xhat), V0 = sum_m nu J_exact(rho_m,nu) q_m q_m^T.
        No xhat-dependent acquisition scaling (the v2/v3 oracle row
        normalization is retired).
A1      trust diagnostics replaced: J_mean = (d_{log rho} E N)^2/(nu Var N),
        mean_info_efficiency = J_mean/J_exact >= 0.90, and the finite-window
        RQL location bias b = |log(rho_dagger/rho)| <= 0.05 with
        rho_dagger = E N/(nu - E N).  Field names J_mean /
        mean_info_efficiency / rql_logload_bias.
A2      two-part peak guard: G3-SHAPE max w_j <= 4 with the ACTUAL support
        size k >= 16, and G3-PHYSICAL max_ij a_ij <= 1536 (uniform k=16 atom
        at the absolute load cap rho=24: 64*24).  Realized peak reported.
A3      every arm returns exactly 52 + 972 physical rows (ridge_fixed_972).
A5      A5_ROW_FEASIBLE_LOAD_BOUND computed per realized cell:
        rho_r <= 1.05 * M * dbar / max_j g_rj.
Sec 3   full R10 minimum dictionary, variable support size k: scattered k=16
        and k=32, solid 4x4, Lblob6x6, bars 1x16/2x8/4x4/8x2/16x1, compact
        32-pixel square (6x6 minus corners) + 4x8 + 8x4, the 16-pixel
        annulus; every cyclic translation, p in {0,1}; SHA256 manifest.
Sec 6   declared task subspace r=200 (top eigenvectors of the FIXED*
        information matrix under identical pre-scan/dwell/resources), fixed
        numerical ridge eps0 = 1e-9 tr(B^T V_FIXED* B)/r frozen before the
        solve.  Additive relaxed certificate
          G_rel = max_a(d_a - theta c_a) - sum_a xi_a(d_a - theta c_a)
                  + theta(b - c^T xi),
        minimized over theta >= 0, field RELAXED_KW_UPPER_BOUND, hard gate
        G_rel/r <= 1e-3.  Exact-design gates Eff_D >= 0.95,
        A-risk <= 1.05 * FIXED*, V >= 0.5 V_FIXED* on B; fail states
        DICTIONARY_RANK_FAIL / CONTINUOUS_CERTIFICATE_FAIL /
        SPECTRAL_GUARD_FAIL / A_RISK_GUARD_FAIL / DOSE_GUARD_FAIL /
        ROUNDING_FAIL.
Sec 5   any dose fallback is MATERIALIZED: FIXED_DOSE is an exact 972-row
        multiset (dose-balanced by construction) and the mixture is realized
        through actual row counts before any guard is quoted.

Solver: active-set pairwise Frank-Wolfe on the projected r x r problem with
CLOSED-FORM rank-2 exact line search (det(V + g(h_j h_j^T - h_i h_i^T)) is
quadratic in g) and Woodbury inverse updates; outer loop re-sweeps the full
dictionary and evaluates the theta-minimized additive certificate.
"""

import hashlib
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import oed_design as v1                    # noqa: E402
import oed_design_v3 as v3                 # noqa: E402
from oed_design_v3 import (_log_fact_table, _log_sf_pois)  # noqa: E402

import oed_design_v2 as v2                 # noqa: E402

RHO_CAP = 24.0
PEAK_PHYSICAL = 1536.0                     # A2 G3-PHYSICAL: 64 * 24

# fast row-level J lookup: the committed exact table (cross-checked to 5e-9
# against fisher_ridge.csv). Used ONLY for per-row information assembly
# (V0, candidate/exact info matrices) where hundreds of distinct continuous
# loads occur; every FROZEN node (palette level, ridge target, admission
# decision) stays on the exact kernel_eval4 kernel.
_JT = None


def _J(rho, nu):
    global _JT
    if _JT is None:
        _JT = v2.get_J_source()
    return float(_JT(float(rho), float(nu)))


_LF_CACHE = {}


def p_ceil_exact(rho, nu):
    """Exact A4 ceiling probability P(Pois(rho) >= nu), closed-form tail."""
    nu = int(nu)
    rho = float(rho)
    if rho <= 0:
        return 0.0
    key = nu
    if key not in _LF_CACHE:
        _LF_CACHE[key] = _log_fact_table(
            int(max(nu + 2, rho + 12 * math.sqrt(rho + 1))) + 400)
    lf = _LF_CACHE[key]
    ls = _log_sf_pois(nu, rho, lf)
    return float(np.exp(ls)) if ls > -700 else 0.0
CEIL_TARGET = 0.01
CEIL_ATOM = 0.05
EFF_MIN = 0.90                             # A1.1 mean-info efficiency
BIAS_MAX = 0.05                            # A1.2 RQL log-load bias
R_SUBSPACE = 200
CERT_TOL_PER_R = 1e-3
W_CLIP = (0.25, 4.0)


# ========================================================================== #
#  A1-amended exact kernel                                                    #
# ========================================================================== #
_K4 = {}


def kernel_eval4(rho, nu):
    """v3 exact kernel + A1 fields: J_mean, mean_info_efficiency,
    rql_logload_bias (all from the same exact PMF pass; dEN/dtheta via the
    telescoping identity E N = sum_{m>=1} G_m => dEN/dtheta = sum Gdot_m)."""
    key = (round(float(rho), 14), int(nu))
    if key in _K4:
        return _K4[key]
    base = dict(v3.kernel_eval(rho, nu))
    rho = float(rho)
    nu = int(nu)
    m = np.arange(0, nu + 2)
    t = nu - (m - 1.0)
    z = rho * np.maximum(t, 0.0)
    zmax = rho * nu
    lf = _log_fact_table(int(max(nu + 2, zmax + 12 * math.sqrt(zmax + 1))) + 200)
    with np.errstate(invalid="ignore", divide="ignore"):
        lgd = np.where((m >= 1) & (z > 0),
                       np.log(np.maximum(m, 1)) +
                       m * np.log(np.maximum(z, 1e-300)) - z - lf[m],
                       -np.inf)
    dEN = float(np.exp(lgd[np.isfinite(lgd)]).sum())     # dE[N]/dtheta
    EN = base["EN"]
    p, mm = base["p"], base["m"]
    VarN = float((p * (mm - EN) ** 2).sum())
    J_mean = dEN ** 2 / (nu * VarN) if VarN > 0 else 0.0
    eff = J_mean / base["J_exact"] if base["J_exact"] > 0 else 0.0
    rho_dag = EN / (nu - EN) if nu > EN else float("inf")
    bias = abs(math.log(rho_dag / rho)) if 0 < rho_dag < np.inf else float("inf")
    base.update({"J_mean": float(J_mean),
                 "mean_info_efficiency": float(eff),
                 "rql_logload_bias": float(bias)})
    _K4[key] = base
    return base


def kernel_certify4(rho, nu, h=1e-5):
    """G7 numerical certification + A1 guard verdicts at one (rho, nu)."""
    c = dict(v3.kernel_certify(rho, nu, h=h))
    k = kernel_eval4(rho, nu)
    c.update({"J_mean": k["J_mean"],
              "mean_info_efficiency": k["mean_info_efficiency"],
              "rql_logload_bias": k["rql_logload_bias"],
              "eff_ok": bool(k["mean_info_efficiency"] >= EFF_MIN),
              "bias_ok": bool(k["rql_logload_bias"] <= BIAS_MAX),
              "p_ceil": k["p_ceil"]})
    c["pass"] = bool(c["pass"] and 0 <= k["J_mean"] <= k["J_exact"] * (1 + 1e-9))
    return c


def _admissible4(rho, nu, ceil_cap):
    k = kernel_eval4(rho, nu)
    return (0 < rho <= RHO_CAP and k["p_ceil"] <= ceil_cap
            and k["mean_info_efficiency"] >= EFF_MIN
            and k["rql_logload_bias"] <= BIAS_MAX)


def ridge_target4(nu):
    """R11 ridge target under the A1-amended admissibility."""
    nodes = np.geomspace(0.05, 64.0, 49)
    Jn = [kernel_eval4(r, nu)["J_exact"] for r in nodes]
    i_u = int(np.argmax(Jn))
    rho_unc = v3._golden_max(lambda r: kernel_eval4(r, nu)["J_exact"],
                             nodes[max(i_u - 1, 0)],
                             nodes[min(i_u + 1, nodes.size - 1)])
    adm_nodes = nodes[nodes <= RHO_CAP + 1e-12]
    adm = np.array([_admissible4(r, nu, CEIL_TARGET) for r in adm_nodes])
    if not adm.any():
        return {"nu": int(nu), "verdict": "RIDGE_TARGET_FAIL"}
    Ja = np.where(adm, [kernel_eval4(r, nu)["J_exact"] for r in adm_nodes],
                  -np.inf)
    i_a = int(np.argmax(Ja))
    rho_R = v3._golden_max(
        lambda r: (kernel_eval4(r, nu)["J_exact"]
                   if _admissible4(r, nu, CEIL_TARGET) else -np.inf),
        adm_nodes[max(i_a - 1, 0)],
        min(adm_nodes[min(i_a + 1, adm_nodes.size - 1)], RHO_CAP))
    if not _admissible4(rho_R, nu, CEIL_TARGET):
        rho_R = float(adm_nodes[i_a])
    k = kernel_eval4(rho_R, nu)
    if abs(math.log(rho_R / rho_unc)) < 1e-4:
        reason = "NONE"
    elif rho_unc > RHO_CAP and abs(rho_R - RHO_CAP) / RHO_CAP < 1e-3:
        reason = "RHO_CAP"
    else:
        up = kernel_eval4(min(rho_R * 1.01, 64.0), nu)
        reason = ("CEILING" if up["p_ceil"] > CEIL_TARGET else "RQL_TRUST")
    return {"nu": int(nu),
            "rho_ridge_exact_unconstrained": float(rho_unc),
            "rho_R_production": float(rho_R),
            "ridge_clip_reason": reason,
            "J_exact_at_target": k["J_exact"],
            "J_mean_at_target": k["J_mean"],
            "mean_info_efficiency": k["mean_info_efficiency"],
            "rql_logload_bias": k["rql_logload_bias"],
            "p_ceiling_scalar": k["p_ceil"], "verdict": "OK"}


def palette4(nu, fast=True):
    """R11 fast palette (with rho_R) or R13 Sec-4 safe palette."""
    if not fast:
        raw = [0.025, 0.05, 0.10]
        ridge = None
    else:
        ridge = ridge_target4(nu)
        rR = ridge["rho_R_production"]
        raw = sorted({0.30, 0.60, 1.00, min(3.00, rR), 0.5 * rR, rR})
    lev, recs = [], []
    for r in raw:
        if lev and abs(r - lev[-1]) <= 1e-10 * max(r, lev[-1]):
            continue
        ok = _admissible4(r, nu, CEIL_ATOM)
        cert = kernel_certify4(r, nu)
        recs.append({"rho": float(r), "admitted": bool(ok and cert["pass"]),
                     "cert": cert["verdict"], **{k_: cert[k_] for k_ in
                     ("J_exact", "J_mean", "mean_info_efficiency",
                      "rql_logload_bias", "p_ceil")}})
        if ok and cert["pass"]:
            lev.append(float(r))
    return {"nu": int(nu), "levels": lev, "records": recs, "ridge": ridge}


# ========================================================================== #
#  Sec 2: balanced 52-row pre-scan + non-oracle V0                            #
# ========================================================================== #
def balanced_prescan_52(side=32):
    """(52 x n) matrix with mean row == all-ones and exactly equal per-pixel
    cumulative dose (each pixel receives 32+16+4 = 52 amplitude units)."""
    n = side * side
    P = np.zeros((52, n))
    r = 0
    # 32 fine rows: block q (row-major on the 8x8 grid of 4x4 blocks)
    # paired with q+32; support 64 px? no: two 4x4 blocks = 32 px, amp n/32.
    for q in range(32):
        for blk in (q, q + 32):
            by, bx = divmod(blk, 8)
            ys, xs = 4 * by, 4 * bx
            for dy in range(4):
                for dx in range(4):
                    P[r, (ys + dy) * side + (xs + dx)] = n / 32.0
        r += 1
    # 16 mid rows: all 8x8 blocks, amp n/64
    for blk in range(16):
        by, bx = divmod(blk, 4)
        ys, xs = 8 * by, 8 * bx
        for dy in range(8):
            for dx in range(8):
                P[r, (ys + dy) * side + (xs + dx)] = n / 64.0
        r += 1
    # 4 coarse rows: quadrants, amp n/256
    for blk in range(4):
        by, bx = divmod(blk, 2)
        ys, xs = 16 * by, 16 * bx
        for dy in range(16):
            for dx in range(16):
                P[r, (ys + dy) * side + (xs + dx)] = n / 256.0
        r += 1
    assert r == 52
    return P


def V0_prescan(P, xhat, nu, rho_bar, B=None):
    """Non-oracle V0: rho_m = rho_bar * (a_m.xhat) with the FROZEN physical
    rows a_m (no rescaling); q_m = a_m/(a_m.xhat); per-row exact kernels."""
    loads = P @ xhat
    r = B.shape[1] if B is not None else P.shape[1]
    V0 = np.zeros((r, r))
    for m in range(P.shape[0]):
        rho_m = rho_bar * float(loads[m])
        q = P[m] / float(loads[m])
        if B is not None:
            q = B.T @ q
        V0 += nu * _J(rho_m, nu) * np.outer(q, q)
    return V0


# ========================================================================== #
#  Sec 3: full variable-k dictionary                                          #
# ========================================================================== #
def _rect(h, w):
    return [(i, j) for i in range(h) for j in range(w)]


def dictionary_families(side=32):
    """R10 minimum dictionary: name -> offset list (variable k)."""
    fams = {}
    rs = np.random.RandomState(63)
    idx16 = rs.choice(side * side, size=16, replace=False)
    fams["scat16"] = [(int(i // side), int(i % side)) for i in idx16]
    rs32 = np.random.RandomState(6332)
    idx32 = rs32.choice(side * side, size=32, replace=False)
    fams["scat32"] = [(int(i // side), int(i % side)) for i in idx32]
    fams["solid4x4"] = _rect(4, 4)
    fams["lblob6x6"] = list(v1._default_shapes(side)["Lblob6x6"])
    fams["bar1x16"] = _rect(1, 16)
    fams["bar2x8"] = _rect(2, 8)
    fams["rect4x4"] = _rect(4, 4)          # R10 names 4x4 in the bar family too
    fams["bar8x2"] = _rect(8, 2)
    fams["bar16x1"] = _rect(16, 1)
    sq = [(i, j) for i in range(6) for j in range(6)
          if (i, j) not in ((0, 0), (0, 5), (5, 0), (5, 5))]
    fams["sq32c"] = sq                      # compact 32-px square variant
    fams["rect4x8"] = _rect(4, 8)
    fams["rect8x4"] = _rect(8, 4)
    fams["ring16"] = list(v1._default_shapes(side)["ring5"])
    for name, offs in fams.items():
        assert len(set(offs)) == len(offs), name
        assert len(offs) >= 16, name
    return fams


def dictionary_manifest(side=32):
    """SHA256 manifest of every family's full translation support table."""
    fams = dictionary_families(side)
    man = {}
    h_all = hashlib.sha256()
    for name in sorted(fams):
        sidx = v1._shape_support(fams[name], side)
        sha = hashlib.sha256(
            np.ascontiguousarray(sidx, dtype=np.int64).tobytes()).hexdigest()
        man[name] = {"k": len(fams[name]), "n_translations": sidx.shape[0],
                     "powers": [0, 1], "sha256": sha}
        h_all.update(sha.encode())
    return {"families": man, "sha256_global": h_all.hexdigest(),
            "side": side}


def build_blocks(xhat, side=32):
    """Per-(family, p) blocks: IDX (n,k), GVAL (n,k) load-normalized
    directions (g.xhat = 1), gmax, w-stat guards. Variable k per block."""
    fams = dictionary_families(side)
    blocks = []
    for name in sorted(fams):
        offs = fams[name]
        k = len(offs)
        sidx = v1._shape_support(offs, side)
        xh_on = xhat[sidx]
        for p in (0, 1):
            w0 = (xh_on / xh_on.mean(axis=1, keepdims=True)) ** float(p)
            w1 = np.clip(w0, W_CLIP[0], W_CLIP[1])
            w = w1 / w1.mean(axis=1, keepdims=True)
            shape_ok = (w.max(axis=1) <= 4.0 + 1e-9)         # G3-SHAPE, real k
            load = (w * xh_on).sum(axis=1)
            g = w / load[:, None]
            blocks.append({"family": name, "p": p, "k": k, "IDX": sidx,
                           "GVAL": g, "gmax": g.max(axis=1),
                           "shape_ok": shape_ok})
    return blocks


# ========================================================================== #
#  Sec 6: task subspace + fixed ridge                                         #
# ========================================================================== #
def info_matrix_full(rows_A, xhat, nu, rho_bar, P=None):
    """Full 1024-dim information of physical rows (+ optional pre-scan):
    V = sum_rows nu J(rho_bar * a.xhat) (a/(a.xhat)) (.)^T."""
    n = rows_A.shape[1]
    V = np.zeros((n, n))
    stack = [rows_A] if P is None else [P, rows_A]
    for A in stack:
        loads = A @ xhat
        for s in range(A.shape[0]):
            a = A[s]
            idx = np.nonzero(a)[0]
            q = a[idx] / float(loads[s])
            J = _J(rho_bar * float(loads[s]), nu)
            V[np.ix_(idx, idx)] += nu * J * np.outer(q, q)
    return V


def subspace_from_fixedstar(V_fix_full, r=R_SUBSPACE):
    """B = top-r eigenvectors of the FIXED* information matrix; fixed ridge
    eps0 = 1e-9 tr(B^T V_FIXED* B)/r (frozen before any OED solve)."""
    evals, evecs = np.linalg.eigh(V_fix_full)
    B = evecs[:, -r:]
    trBVB = float(evals[-r:].sum())
    eps0 = 1e-9 * trBVB / r
    return B, eps0, trBVB


# ========================================================================== #
#  FIXED_DOSE: exact 972-row dose-balanced fallback (Sec 5)                   #
# ========================================================================== #
def fixed_dose_972(side=32, target_mean_load=0.60, band=0.045, n_iter=3000):
    """Exact 972-row dose-balanced multiset (R13 Sec 5 FIXED_DOSE).

    Construction: the bar1x16 translate family (horizontal 16-px strips,
    amplitude n/16) minus a STRUCTURED 52-translate drop chosen so no pixel
    loses more than one covering strip (drop x-offset 0 for all 32 y-rows
    and x-offset 16 for y-rows 0..19), followed by deterministic per-row
    multiplicative (Sinkhorn-style) amplitude balancing to the +/-band
    per-pixel dose target. The lblob16-based drop cannot satisfy the band
    (double-decremented pixels); the strip family can. Returns
    (rows, dose_dev)."""
    n = side * side
    drop = {(y, 0) for y in range(side)} | {(y, 16) for y in range(20)}
    rows_list = []
    for y in range(side):
        for x0 in range(side):
            if (y, x0) in drop:
                continue
            row = np.zeros(n)
            for j in range(16):
                row[y * side + (x0 + j) % side] = n / 16.0
            rows_list.append(row)
    base = np.array(rows_list)
    assert base.shape[0] == 972, base.shape
    c = np.ones(972)
    supp = [np.nonzero(base[s])[0] for s in range(972)]
    for _ in range(n_iter):
        dose = (c[:, None] * base).sum(axis=0)
        rel = dose / dose.mean()
        if np.abs(rel - 1.0).max() <= band - 0.005:
            break
        for s in range(972):
            c[s] /= rel[supp[s]].mean() ** 0.7
    rows = c[:, None] * base
    rows *= 972.0 / c.sum()
    dose = rows.sum(axis=0)
    dev = float(np.abs(dose / dose.mean() - 1.0).max())
    return rows, dev


def ridge_fixed_972(xhat, nu, side=32):
    """A3-compliant RIDGE-FIXED: the SAME frozen 972-row lblob16 multiset as
    the fixed baselines (first 972 translates), ONE global multiplier
    calibrated to mean pre-scan-predicted load rho_R(nu), then the R11
    bisection safety clip under the A1-amended guards."""
    n = side * side
    from patterns import make_patterns
    A972 = make_patterns("lblob16", n, n, 0)["A"][:972]
    base_load = A972 @ xhat
    ridge = ridge_target4(nu)
    rho_R = ridge["rho_R_production"]
    mult = rho_R / float(base_load.mean())
    grid = np.geomspace(max(1e-3, 1e-3 * rho_R), max(1.05 * mult
                        * float(base_load.max()), 1.0), 61)
    kv = [kernel_eval4(r, int(nu)) for r in grid]
    lg = np.log(grid)
    eff_g = np.array([k["mean_info_efficiency"] for k in kv])
    bias_g = np.array([k["rql_logload_bias"] for k in kv])
    pc_g = np.array([k["p_ceil"] for k in kv])

    def ok(m_):
        ll = np.log(np.clip(m_ * base_load, grid[0], grid[-1]))
        pc = np.interp(ll, lg, pc_g)
        return (pc.mean() <= CEIL_TARGET and pc.max() <= CEIL_ATOM
                and np.interp(ll, lg, eff_g).min() >= EFF_MIN
                and np.interp(ll, lg, bias_g).max() <= BIAS_MAX)

    clipped = False
    if not ok(mult):
        clipped = True
        lo, hi = 0.0, mult
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if ok(mid):
                lo = mid
            else:
                hi = mid
        mult = lo
    loads = mult * base_load
    achieved = float(loads.mean())
    return {"A_main": mult * A972, "loads": loads, "ridge": ridge,
            "requested_mean_load": float(rho_R),
            "achieved_mean_load": achieved,
            "global_multiplier": float(mult),
            "bisection_clip_applied": bool(clipped),
            "RIDGE_GUARD_CLIPPED": bool(achieved < 0.90 * rho_R),
            "n_main_rows": int(A972.shape[0])}


# ========================================================================== #
#  the v4 design solver                                                       #
# ========================================================================== #
def _cert_grel(d, c, xi_d, xi_c, dbar, budget):
    """theta-minimized additive dual bound G_rel (Sec 6.1):
    G(theta) = max_a(d_a - theta c_a) - sum_a xi_a(d_a - theta c_a)
               + theta(b - c^T xi)
             = max_a(d_a - theta c_a) - <xi, d> + theta * b,
    valid for every theta >= 0; minimized by deterministic ternary search."""
    def f(theta):
        return float(np.max(d - theta * c) - dbar + theta * budget)
    # piecewise-linear in theta: golden search on [0, theta_hi]
    th_hi = max(1.0, float(np.max(d) / max(np.min(c), 1e-9)))
    lo, hi = 0.0, th_hi
    for _ in range(80):
        m1_ = lo + (hi - lo) / 3
        m2_ = hi - (hi - lo) / 3
        if f(m1_) <= f(m2_):
            hi = m2_
        else:
            lo = m1_
    th = 0.5 * (lo + hi)
    cands = [0.0, th]
    vals = [f(t) for t in cands]
    i = int(np.argmin(vals))
    return float(vals[i]), float(cands[i])


def toy_cert_check():
    """Coordinator-mandated convention check: unconstrained D-opt on the two
    axis atoms in R^2 has the known optimum xi = (1/2, 1/2) with
    d_a = x_a^T M^-1 x_a = 2 = r and G_rel = 0; a dominated third atom must
    not disturb it. Returns (G_rel_at_optimum, d_values)."""
    X = np.array([[1.0, 0.0], [0.0, 1.0], [0.5, 0.5]])
    xi = np.array([0.5, 0.5, 0.0])
    M = sum(w * np.outer(x, x) for w, x in zip(xi, X))
    Mi = np.linalg.inv(M)
    d = np.array([x @ Mi @ x for x in X])
    dbar = float((xi * d).sum())
    c = np.ones(3)
    g, th = _cert_grel(d, c, None, 1.0, dbar, 1.0)
    return float(g), d.tolist()


def design_v4(xhat, nu, budget_mean, palette, B, eps0, V_fix_B,
              rho_prescan, M_rows=972, max_outer=60, max_inner=5000,
              dose_band=0.05, verbose=False, side=32,
              fixed_dose_rows=None):
    """R13-compliant budgeted design on the projected r=200 subspace.

    Returns dict with rows (M_rows x n RELATIVE-irradiance load-profile
    rows whose predicted loads are the atom levels), atoms, certificate
    (RELAXED_KW_UPPER_BOUND), exact-design guards, verdict/fail state, and
    the full disclosure block."""
    n = side * side
    r = B.shape[1]
    levels = np.asarray(palette["levels"], dtype=float)
    L = levels.size
    nuJ = np.array([nu * kernel_eval4(l_, int(nu))["J_exact"]
                    for l_ in levels])
    blocks = build_blocks(xhat, side)
    # per-geometry projected directions QB and admission masks
    QBs, metas = [], []
    for b_ in blocks:
        QB = np.einsum('tk,tkr->tr', b_["GVAL"], B[b_["IDX"]])
        QBs.append(QB)
        metas.append(b_)
    QB = np.vstack(QBs)                                     # (G_tot, r)
    G_tot = QB.shape[0]
    fam_of = np.concatenate([[i] * m["IDX"].shape[0]
                             for i, m in enumerate(metas)])
    t_of = np.concatenate([np.arange(m["IDX"].shape[0]) for m in metas])
    gmax = np.concatenate([m["gmax"] for m in metas])
    shape_ok = np.concatenate([m["shape_ok"] for m in metas])
    # A2 G3-PHYSICAL per (geometry, level): rho * gmax <= 1536
    ALLOW = shape_ok[:, None] & (gmax[:, None] * levels[None, :]
                                 <= PEAK_PHYSICAL + 1e-9)
    if not ALLOW.any():
        return {"verdict": "DICTIONARY_RANK_FAIL", "reason": "no admissible atoms"}

    V0 = V0_prescan(balanced_prescan_52(side), xhat, nu, rho_prescan, B)
    I_r = np.eye(r)

    def assemble(xi):
        u_geo = M_rows * (xi * nuJ[None, :]).sum(axis=1)
        nz = np.nonzero(u_geo)[0]
        M = (QB[nz] * u_geo[nz, None]).T @ QB[nz]
        return V0 + M + eps0 * I_r

    dose_of = None                      # full-pixel dose of a measure xi

    def _dose(xi):
        w_geo = (xi * levels[None, :]).sum(axis=1)
        d = np.zeros(n)
        off = 0
        for mb in metas:
            g_ = mb["IDX"].shape[0]
            w = w_geo[off:off + g_]
            nz = np.nonzero(w)[0]
            if nz.size:
                d += np.bincount(mb["IDX"][nz].ravel(),
                                 weights=(w[nz, None] * mb["GVAL"][nz]).ravel(),
                                 minlength=n)
            off += g_
        return d
    dose_of = _dose

    def full_cert(xi, Vinv):
        quad = np.einsum('gr,rs,gs->g', QB, Vinv, QB)
        D = M_rows * quad[:, None] * nuJ[None, :]
        D = np.where(ALLOW, D, -np.inf)
        dbar = float((xi * np.where(ALLOW, D, 0.0)).sum())
        xi_c = float((xi @ levels).sum())
        g_, th_ = _cert_grel(D[ALLOW],
                             np.repeat(levels[None, :], G_tot, 0)[ALLOW],
                             None, xi_c, dbar, budget_mean)
        return g_, th_, D, dbar, xi_c

    # ================= PHASE 1: relaxed (no-dose) pairwise FW ============== #
    # vertex decomposition of the budget polytope: singles (c<=b) and
    # budget-boundary level pairs; pairwise moves stay feasible by convexity.
    lo_lvl = int(np.argmax(ALLOW.any(axis=0)))
    xi = np.zeros((G_tot, L))
    a0 = (int(np.argmax(ALLOW[:, lo_lvl])), lo_lvl)
    xi[a0] = 1.0
    weights = {("s", a0): 1.0}

    def vertex_atoms(v):
        if v[0] == "s":
            return [(v[1], 1.0)]
        (glo, llo), (ghi, lhi) = v[1], v[2]
        wlo = (levels[lhi] - budget_mean) / (levels[lhi] - levels[llo])
        return [((glo, llo), wlo), ((ghi, lhi), 1.0 - wlo)]

    V = assemble(xi)
    Vinv = np.linalg.inv(V)
    traj = []
    cert, theta_star = np.inf, 0.0
    n_inner_tot = 0
    stall = 0
    best_cert = np.inf
    for outer in range(max_outer):
        cert, theta_star, D, dbar, xi_c = full_cert(xi, Vinv)
        traj.append((outer, float(cert / r)))
        if verbose:
            print("    [v4] outer=%2d G_rel/r=%.3e logdet=%.4f load=%.4f"
                  % (outer, cert / r, np.linalg.slogdet(V)[1], xi_c),
                  flush=True)
        if cert / r <= CERT_TOL_PER_R:
            break
        if cert < best_cert * 0.995:
            best_cert, stall = cert, 0
        else:
            stall += 1
            if stall >= 3 or n_inner_tot >= max_inner:
                break                                   # plateau: report it
        # candidate pool: support atoms + top sensitivities
        s_adj = D - theta_star * levels[None, :]
        top = np.column_stack(np.unravel_index(
            np.argsort(s_adj, axis=None)[::-1][:300], s_adj.shape))
        pool = {tuple(t) for t in top.tolist()}
        pool |= {tuple(t) for t in
                 np.column_stack(np.nonzero(xi > 0)).tolist()}
        pool = [p_ for p_ in pool if ALLOW[p_[0], p_[1]]]
        pg = np.array([p_[0] for p_ in pool])
        pl = np.array([p_[1] for p_ in pool])
        Qp = QB[pg]
        sclp = M_rows * nuJ[pl]
        cp = levels[pl]
        for _ in range(max(max_inner // max_outer, 150)):
            n_inner_tot += 1
            dp = sclp * np.einsum('ar,rs,as->a', Qp, Vinv, Qp)
            # LMO over the polytope restricted to the pool
            feas = cp <= budget_mean + 1e-12
            best_v, best_val = None, -np.inf
            if feas.any():
                i_s = int(np.argmax(np.where(feas, dp, -np.inf)))
                best_v, best_val = ("s", pool[i_s]), float(dp[i_s])
            for llo in range(L):
                if levels[llo] > budget_mean:
                    continue
                m_lo = np.where(pl == llo, dp, -np.inf)
                if not np.isfinite(m_lo).any():
                    continue
                i_lo = int(np.argmax(m_lo))
                for lhi in range(L):
                    if levels[lhi] <= budget_mean:
                        continue
                    m_hi = np.where(pl == lhi, dp, -np.inf)
                    if not np.isfinite(m_hi).any():
                        continue
                    i_hi = int(np.argmax(m_hi))
                    wlo = (levels[lhi] - budget_mean) / (levels[lhi]
                                                         - levels[llo])
                    val = wlo * dp[i_lo] + (1 - wlo) * dp[i_hi]
                    if val > best_val:
                        best_val = float(val)
                        best_v = ("p", pool[i_lo], pool[i_hi])
            # away vertex: worst active vertex by gradient
            dmap = {p_: float(dp[k]) for k, p_ in enumerate(pool)}
            away_v, away_val = None, np.inf
            for v_, w_ in weights.items():
                if w_ <= 1e-14:
                    continue
                val = sum(co * dmap.get(a_, 0.0)
                          for a_, co in vertex_atoms(v_))
                if val < away_val:
                    away_val, away_v = val, v_
            if best_v is None or away_v is None or best_v == away_v:
                break
            if best_val - away_val <= 1e-9 * max(abs(dbar), 1.0):
                break
            # move direction: gamma * (v_fw - v_away), rank <= 4
            coef = {}
            for a_, co in vertex_atoms(best_v):
                coef[a_] = coef.get(a_, 0.0) + co
            for a_, co in vertex_atoms(away_v):
                coef[a_] = coef.get(a_, 0.0) - co
            coef = {a_: co for a_, co in coef.items() if abs(co) > 1e-15}
            if not coef:
                break
            Ucols, svec = [], []
            for (g_, l_), co in coef.items():
                Ucols.append(math.sqrt(M_rows * nuJ[l_]) * QB[g_])
                svec.append(co)
            U = np.column_stack(Ucols)
            svec = np.array(svec)
            K = U.T @ Vinv @ U
            gmax_w = weights[away_v]

            def phi(gam):
                Mm = np.eye(len(svec)) + gam * (svec[:, None] * K)
                s_, ldm = np.linalg.slogdet(Mm)
                return ldm if s_ > 0 else -np.inf
            glo_, ghi_ = 0.0, gmax_w
            for _ in range(45):
                m1_ = glo_ + (ghi_ - glo_) / 3
                m2_ = ghi_ - (ghi_ - glo_) / 3
                if phi(m1_) < phi(m2_):
                    glo_ = m1_
                else:
                    ghi_ = m2_
            g_star = 0.5 * (glo_ + ghi_)
            if phi(g_star) <= 1e-12:
                break
            weights[away_v] = weights.get(away_v, 0.0) - g_star
            if weights[away_v] <= 1e-14:
                weights.pop(away_v, None)
            weights[best_v] = weights.get(best_v, 0.0) + g_star
            for a_, co in coef.items():
                xi[a_] = max(xi[a_] + g_star * co, 0.0)
            C = g_star * np.diag(svec)
            Minv2 = np.linalg.inv(np.linalg.inv(C) + K)
            VU = Vinv @ U
            Vinv = Vinv - VU @ Minv2 @ VU.T
        V = assemble(xi)
        Vinv = np.linalg.inv(V)
    cert_relaxed = cert / r
    traj_relaxed = traj

    # ============ PHASE 2: dose-feasible design (G5 in-optimizer) ========== #
    # multiplicative ascent + IPF dose projection + exact budget repair on
    # the same projected objective; its G_rel is reported as the dose-
    # feasible certificate (R13 Sec 6.1 wording) with the plateau trajectory.
    xi2 = xi.copy()
    traj_dose = []
    for it2 in range(150):
        V2 = assemble(xi2)
        Vi2 = np.linalg.inv(V2)
        if it2 % 10 == 0 or it2 == 149:
            g_, th_, D2, dbar2, xc2 = full_cert(xi2, Vi2)
            dv = dose_of(xi2)
            traj_dose.append((it2, float(g_ / r),
                              float(np.abs(dv / dv.mean() - 1).max())))
        else:
            quad = np.einsum('gr,rs,gs->g', QB, Vi2, QB)
            D2 = M_rows * quad[:, None] * nuJ[None, :]
            D2 = np.where(ALLOW, D2, -np.inf)
        pos = np.where(ALLOW, np.maximum(D2, 0.0), 0.0)
        zb = float((xi2 * pos).sum())
        if zb <= 0:
            break
        xi2 = xi2 * np.sqrt(pos / zb)
        xi2[~ALLOW] = 0.0
        xi2 /= xi2.sum()
        # IPF dose projection toward the inner band
        for _ in range(8):
            dv = dose_of(xi2)
            rel = dv / dv.mean()
            if np.abs(rel - 1.0).max() <= dose_band - 0.015:
                break
            rc = np.clip(rel, 0.5, 2.0)
            off = 0
            fac_all = np.ones(G_tot)
            for mb in metas:
                g_ = mb["IDX"].shape[0]
                fac_all[off:off + g_] = 1.0 / (
                    (mb["GVAL"] * rc[mb["IDX"]]).sum(axis=1)
                    / mb["GVAL"].sum(axis=1))
                off += g_
            xi2 = xi2 * fac_all[:, None]
            xi2[~ALLOW] = 0.0
            xi2 /= xi2.sum()
        # exact budget repair: mix toward each geometry's lowest level
        cur = float((xi2 @ levels).sum())
        if cur > budget_mean + 1e-15:
            low = np.zeros_like(xi2)
            lo_idx = np.argmax(ALLOW, axis=1)
            low[np.arange(G_tot), lo_idx] = xi2.sum(axis=1)
            lmin = float((low @ levels).sum())
            t_ = 1.0 if lmin >= budget_mean else \
                (cur - budget_mean) / (cur - lmin)
            xi2 = (1 - t_) * xi2 + t_ * low
        xi2[xi2 < xi2.max() * 1e-16] = 0.0
        xi2 /= xi2.sum()
    # pure-IPF feasibility tail: drive the continuous dose INTO the band
    # (no further ascent; the certificate is recomputed afterwards and any
    # increase is the honest price of the binding dose constraint)
    for _ in range(800):
        dv = dose_of(xi2)
        rel = dv / dv.mean()
        if np.abs(rel - 1.0).max() <= dose_band - 0.015:
            break
        rc = np.clip(rel, 0.5, 2.0)
        off = 0
        fac_all = np.ones(G_tot)
        for mb in metas:
            g_ = mb["IDX"].shape[0]
            fac_all[off:off + g_] = 1.0 / (
                (mb["GVAL"] * rc[mb["IDX"]]).sum(axis=1)
                / mb["GVAL"].sum(axis=1))
            off += g_
        xi2 = xi2 * fac_all[:, None]
        xi2[~ALLOW] = 0.0
        xi2 /= xi2.sum()
        cur = float((xi2 @ levels).sum())
        if cur > budget_mean + 1e-15:
            low = np.zeros_like(xi2)
            lo_idx = np.argmax(ALLOW, axis=1)
            low[np.arange(G_tot), lo_idx] = xi2.sum(axis=1)
            lmin = float((low @ levels).sum())
            t_ = 1.0 if lmin >= budget_mean else \
                (cur - budget_mean) / (cur - lmin)
            xi2 = (1 - t_) * xi2 + t_ * low
    V = assemble(xi2)
    Vinv = np.linalg.inv(V)
    cert_d, theta_star, _D, _db, _xc = full_cert(xi2, Vinv)
    cert_dose = cert_d / r
    xi = xi2                                     # rounding uses the dose-
    relaxed_kw = cert_dose                       # feasible design (Sec 6.1)
    cert_ok = relaxed_kw <= CERT_TOL_PER_R
    sign, ld_cont = np.linalg.slogdet(V)

    # ---- rounding to M_rows exact rows ------------------------------------ #
    flat = xi.ravel()
    supp = np.where(flat > flat.max() * 1e-9)[0]
    w_s = flat[supp] / flat[supp].sum()
    quota = w_s * M_rows
    cnt = np.floor(quota).astype(np.int64)
    rem = int(M_rows - cnt.sum())
    order = np.lexsort((supp, -(quota - cnt)))
    for i2 in range(rem):
        cnt[order[i2]] += 1
    atoms = []                                # (family, t, p, rho, count)
    for fa, c_ in zip(supp, cnt):
        if c_ == 0:
            continue
        g_, l_ = divmod(int(fa), L)
        mb = metas[int(fam_of[g_])]
        atoms.append((mb["family"], int(t_of[g_]), int(mb["p"]),
                      float(levels[l_]), int(c_)))

    def rows_of(atom_list):
        rows = np.zeros((sum(a[4] for a in atom_list), n))
        rr = 0
        for (fam, t_, p_, rho_, c_) in atom_list:
            bi = next(i for i, m in enumerate(metas)
                      if m["family"] == fam and m["p"] == p_)
            row = np.zeros(n)
            row[metas[bi]["IDX"][t_]] = rho_ * metas[bi]["GVAL"][t_]
            rows[rr:rr + c_] = row
            rr += c_
        return rows

    rows = rows_of(atoms)
    # budget demotion (deterministic) on counts
    def total_load(al):
        return float(sum(a[3] * a[4] for a in al))
    guard_flags = []
    while total_load(atoms) > budget_mean * M_rows + 1e-9:
        atoms.sort(key=lambda a: -a[3])
        fam, t_, p_, rho_, c_ = atoms[0]
        lo = float(levels[0])
        if rho_ <= lo + 1e-12:
            break
        atoms[0] = (fam, t_, p_, rho_, c_ - 1)
        atoms = [a for a in atoms if a[4] > 0]
        merged = False
        for i2, a in enumerate(atoms):
            if a[:3] == (fam, t_, p_) and abs(a[3] - lo) < 1e-12:
                atoms[i2] = (fam, t_, p_, lo, a[4] + 1)
                merged = True
                break
        if not merged:
            atoms.append((fam, t_, p_, lo, 1))
        rows = rows_of(atoms)

    # ---- post-rounding dose polish (strict descent, same-level swaps,
    # budget/row neutral) before any mixture is considered ------------------ #
    def _single_row(fam, t_, p_, rho_):
        bi = next(i for i, m in enumerate(metas)
                  if m["family"] == fam and m["p"] == p_)
        row = np.zeros(n)
        row[metas[bi]["IDX"][t_]] = rho_ * metas[bi]["GVAL"][t_]
        return row

    dosev = rows.sum(axis=0) / rows.shape[0]
    band_in = dose_band - 0.004

    def _pen(dv):
        rl = np.abs(dv / dv.mean() - 1.0)
        return float((np.maximum(rl - band_in, 0.0) ** 2).sum())

    Ppen = _pen(dosev)
    n_polish = 0
    while n_polish < 1200 and Ppen > 0.0:
        rel = dosev / dosev.mean() - 1.0
        jo = int(np.argmax(rel))
        ju = int(np.argmin(rel))
        best = (0.0, None)
        for ai, (fam, t_, p_, rho_, c_) in enumerate(atoms):
            bi = next(i for i, m in enumerate(metas)
                      if m["family"] == fam and m["p"] == p_)
            if jo not in metas[bi]["IDX"][t_]:
                continue
            r_rem = _single_row(fam, t_, p_, rho_) / rows.shape[0]
            # candidates: same family+p, translations covering ju
            cov = np.nonzero((metas[bi]["IDX"] == ju).any(axis=1))[0][:48]
            for t2 in cov:
                if t2 == t_:
                    continue
                r_add = np.zeros(n)
                r_add[metas[bi]["IDX"][t2]] = rho_ * metas[bi]["GVAL"][t2]
                r_add /= rows.shape[0]
                dP = _pen(dosev - r_rem + r_add) - Ppen
                if dP < best[0] - 1e-15:
                    best = (dP, (ai, int(t2), r_rem, r_add))
        if best[1] is None:
            break
        ai, t2, r_rem, r_add = best[1]
        fam, t_, p_, rho_, c_ = atoms[ai]
        atoms[ai] = (fam, t_, p_, rho_, c_ - 1)
        atoms = [a for a in atoms if a[4] > 0]
        merged = False
        for k2, a in enumerate(atoms):
            if a[:4] == (fam, t2, p_, rho_):
                atoms[k2] = (fam, t2, p_, rho_, a[4] + 1)
                merged = True
                break
        if not merged:
            atoms.append((fam, t2, p_, rho_, 1))
        dosev = dosev - r_rem + r_add
        Ppen = _pen(dosev)
        n_polish += 1
    rows = rows_of(atoms)

    # ---- guards on the exact returned rows -------------------------------- #
    def guard_suite(rows_, label):
        loads_ = rows_ @ xhat
        dose = rows_.sum(axis=0) / rows_.shape[0]
        dmean = dose.mean()
        dose_dev = float(np.abs(dose / dmean - 1.0).max())
        peak = float(rows_.max())
        V_ex = V0.copy()
        for s in range(rows_.shape[0]):
            idx = np.nonzero(rows_[s])[0]
            q = B[idx].T @ (rows_[s, idx] / float(loads_[s]))
            J = _J(float(loads_[s]), nu)
            V_ex += nu * J * np.outer(q, q)
        V_ex += eps0 * I_r
        _s, ld_ex = np.linalg.slogdet(V_ex)
        eff_D = math.exp((ld_ex - ld_cont) / r)
        arisk = float(np.trace(np.linalg.inv(V_ex)))
        arisk_fix = float(np.trace(np.linalg.inv(V_fix_B + eps0 * I_r)))
        lam_f = np.linalg.eigvalsh(V_fix_B + eps0 * I_r)
        Lf = np.linalg.cholesky(V_fix_B + eps0 * I_r)
        C = np.linalg.solve(Lf, np.linalg.solve(Lf, V_ex.T).T)
        mu_min = float(np.linalg.eigvalsh(C).min())
        p95_load = float(np.percentile(loads_, 95))
        a5_bound = 1.05 * rows_.shape[0] * dmean / max(
            [rows_[s].max() / max(loads_[s], 1e-300)
             for s in range(rows_.shape[0])])
        ceil_w = float(np.mean([p_ceil_exact(float(l_), int(nu))
                                for l_ in loads_]))
        return {"label": label,
                "incident_sum_rho": float(loads_.sum()),
                "budget_pass": bool(loads_.sum()
                                    <= budget_mean * rows_.shape[0] + 1e-9),
                "detected_pred": float(sum(
                    nu * l_ / (1.0 + l_) for l_ in loads_)),
                "dose_dev": dose_dev,
                "dose_pass": bool(dose_dev <= dose_band),
                "peak_physical": peak,
                "peak_pass": bool(peak <= PEAK_PHYSICAL + 1e-9),
                "peak_ratio": peak / PEAK_PHYSICAL,
                "eff_D": eff_D, "effD_pass": bool(eff_D >= 0.95),
                "a_risk": arisk, "a_risk_fix": arisk_fix,
                "a_risk_pass": bool(arisk <= 1.05 * arisk_fix),
                "mu_min_vs_fix": mu_min,
                "spectral_pass": bool(mu_min >= 0.5 - 1e-9),
                "ceiling_weighted": ceil_w,
                "ceiling_pass": bool(ceil_w <= 0.01),
                "rho_5": float(np.percentile(loads_, 5)),
                "rho_50": float(np.percentile(loads_, 50)),
                "rho_95": p95_load, "rho_max": float(loads_.max()),
                "A5_ROW_FEASIBLE_LOAD_BOUND": float(a5_bound)}

    gsuite = guard_suite(rows, "rounded")
    mixture_m = None
    if not gsuite["dose_pass"] and fixed_dose_rows is not None:
        # Sec 5: MATERIALIZED exact mixture with FIXED_DOSE rows
        fd = fixed_dose_rows * budget_mean          # relative rows -> loads
        for m_keep in range(M_rows, -1, -max(M_rows // 60, 1)):
            oed_atoms = []
            kept = 0
            for (fam, t_, p_, rho_, c_) in atoms:
                take = min(c_, m_keep - kept)
                if take > 0:
                    oed_atoms.append((fam, t_, p_, rho_, take))
                    kept += take
            n_fd = M_rows - kept
            stride_rows = fd[np.linspace(0, fd.shape[0] - 1, n_fd)
                             .astype(int)] if n_fd else np.zeros((0, n))
            rows_mix = np.vstack([rows_of(oed_atoms), stride_rows]) \
                if oed_atoms else stride_rows
            gm = guard_suite(rows_mix, "mixture_m%d" % kept)
            if gm["dose_pass"] and gm["budget_pass"] and gm["peak_pass"]:
                rows, gsuite, mixture_m = rows_mix, gm, kept
                atoms = oed_atoms + [("FIXED_DOSE", -1, 0, float(l_), 1)
                                     for l_ in (stride_rows @ xhat)]
                break

    # ---- verdict ----------------------------------------------------------- #
    if not cert_ok:
        verdict = "CONTINUOUS_CERTIFICATE_FAIL"
    elif not gsuite["dose_pass"]:
        verdict = "DOSE_GUARD_FAIL"
    elif not gsuite["effD_pass"]:
        verdict = "ROUNDING_FAIL"
    elif not gsuite["spectral_pass"]:
        verdict = "SPECTRAL_GUARD_FAIL"
    elif not gsuite["a_risk_pass"]:
        verdict = "A_RISK_GUARD_FAIL"
    else:
        verdict = "OK"

    return {"rows": rows, "atoms": atoms, "verdict": verdict,
            "RELAXED_KW_UPPER_BOUND": float(relaxed_kw),
            "cert_relaxed_solver": float(cert_relaxed),
            "traj_relaxed": traj_relaxed, "traj_dose": traj_dose,
            "theta_star": float(theta_star),
            "cert_pass": bool(cert_ok), "logdet_cont": float(ld_cont),
            "guards": gsuite, "mixture_rows_kept": mixture_m,
            "levels": [float(l_) for l_ in levels],
            "r": r, "eps0": float(eps0),
            "n_atoms_admissible": int(ALLOW.sum()),
            "dict_manifest_sha": dictionary_manifest(side)["sha256_global"]}
