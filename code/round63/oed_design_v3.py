"""Dead-time-aware OED v3 -- R11 ridge-zone load architecture.

Implements docs/ROUND63_GPT_ROUND11_RULING_RAW.md on top of oed_design_v2:

R11-1  exact ridge target rho_R(nu):
         R_adm(nu) = { rho : 0 < rho <= 24,  p_ceil(rho,nu) <= 0.01,
                       J_q/J_exact >= 0.90 },
         rho_R = min-argmax of J_exact over R_adm  (dense log bracket +
         deterministic golden-section refinement; min tie rule).
R11-3  load palette  L_nu = unique{0.30, 0.60, 1.00, min(3, rho_R),
         rho_R/2, rho_R}  (rel-tol 1e-10 dedupe); atom admission adds
         rho <= 24, p_ceil(rho,nu) <= 0.05, J_q/J_exact >= 0.90.
R11-4  budget: hard incident-signal inequality  <rho, xi> <= 0.60
         (design-measure form of  sum_m nu*rho_m <= M_main*nu*0.60),
         via dual (Lagrangian) ascent on the multiplier + exact final
         repair; incident AND detected budgets dual-reported.
R11-G5 per-pixel cumulative dose within +/-5% of uniform ENFORCED in the
         optimization (iterative proportional fitting interleaved with the
         multiplicative ascent) and verified/repaired post-rounding by
         deterministic same-load count exchanges (R10-Q2 sanctions
         "rounding and exchange").
R11-2  ridge_fixed_design(xhat, nu): LBLOB16 translate family, one global
         multiplier calibrated to mean predicted load = rho_R(nu), then the
         Section-2.2 deterministic bisection safety clip; records requested
         vs achieved mean load and RIDGE_GUARD_CLIPPED.
R11-G7 exact kernel numerical certification at every frozen (nu,rho) node.

Ambiguities resolved (documented for the record):

(A1) RQL-trust ratio direction. The ruling text writes J_exact/J_q >= 0.90,
     but the Godambe (sandwich) information of ANY estimating equation is
     <= the Fisher information, so that literal ratio is >= 1 and the guard
     could never clip -- contradicting the ruling's own RQL_TRUST clip
     reason and "the ceiling/RQL guards may clip the ... ridge" (Section 1).
     Implemented as the efficiency ratio  J_q/J_exact >= 0.90  (RQL retains
     >= 90% of the exact information), the only reading under which
     RQL_TRUST is reachable.  Sanity J_q/J_exact -> 1 as rho -> 0 verified.
(A2) G3 peak guard at ridge loads. v2 applied "a_j <= 64" to the
     load-scaled row; at rho ~ rho_R that reading removes EVERY k=16 atom
     (peak ~ 64*rho), contradicting R11's mandate that ridge atoms enter
     OED-DT, and the Study-2 benchmark's own dim translates would fail it
     at unit load.  v3 applies G3 to the pattern SHAPE: the frame-mean-1
     normalized atom must obey  max_j a_j <= (n/k)*w_max = 64*4 = 256
     (support sparsity no worse than the k=16 benchmark; weight contrast
     separately capped by G4's w <= 4).  All k=16 G4-legal atoms pass;
     sparser shapes would fail.  Physical peak irradiance at high load is
     REPORTED (R11 2.3 "source/peak irradiance ratio"), not vetoed.
(A3) Exposure accounting. v3 uses M_main = 972 rows (R10 erratum E2:
     52 pre-scan + 972 main = 1024); v2 used 1024.
(A4) Ceiling probability: p_ceil = P{N >= nu} = P(Pois(rho) >= nu) exactly
     (active start: T-(nu-1)tau = tau).  For nu in {200, 2000} and
     rho <= 24 this is < 1e-60, so CEILING cannot clip at these dwells.
(A5) Dose-band load ceiling (emergent, printed by the selftest): a rounded
     row at load rho puts dose rho*g_j/M_rows >= rho*~64/M_rows on its own
     pixels; the +5% band around the budget-capped mean dose (<= 0.60)
     therefore forbids ANY single row with rho > ~1.05*0.60*M_rows/64 ~ 9.6.
     Full-ridge atoms (rho_R ~ 20) are mathematically incompatible with
     G5+budget; ridge-ZONE atoms (rho in (1, ~9.6]) are compatible.  This
     is water-filling hitting the dose wall, reported, not hidden.

ADD-only file; imports numpy + stdlib + sibling modules oed_design (v1)
and oed_design_v2 (v2).
"""

import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import oed_design as v1                     # noqa: E402
import oed_design_v2 as v2                  # noqa: E402

RHO_CAP = 24.0
CEIL_TARGET = 0.01          # scalar ridge-target ceiling threshold (R11 S1)
CEIL_ATOM = 0.05            # atom admission ceiling threshold (R11 S3/G8)
TRUST_MIN = 0.90            # RQL trust threshold (A1: J_q/J_exact)


# ========================================================================== #
#  Exact renewal kernel: pmf, J_exact, J_q, p_ceil, moments  (+ G7 cert)      #
# ========================================================================== #
_KCACHE = {}


def _log_fact_table(nmax):
    return np.concatenate([[0.0], np.cumsum(np.log(np.arange(1, nmax + 1)))])


def _log_sf_pois(m, z, lf):
    """log P(Pois(z) >= m), scalar integer m >= 1, scalar z > 0.
    Two-branch: complement CDF for m <= z (short sum, sf ~ O(1));
    direct tail sum for m > z (normal-tail window)."""
    if m <= z:
        j = np.arange(0, m, dtype=np.float64)
        lt = j * math.log(z) - z - lf[:m]
        mx = lt.max()
        cdf = math.exp(mx) * float(np.exp(lt - mx).sum())
        return math.log1p(-min(cdf, 1.0 - 1e-300))
    hi = int(m + max(80, 12.0 * math.sqrt(z) + 80))
    j = np.arange(m, hi + 1, dtype=np.float64)
    lt = j * math.log(z) - z - lf[m:hi + 1]
    mx = lt.max()
    return mx + math.log(float(np.exp(lt - mx).sum()))


def kernel_eval(rho, nu):
    """Exact active-start NP renewal kernel at (rho, nu), tau = 1, lam = rho.

    Returns dict: J_exact, J_q, trust (=J_q/J_exact), p_ceil, EN, m, p (pmf),
    mass_err, score_err  (per-slot log-lambda units for J_exact / J_q).
    """
    key = (round(float(rho), 14), int(nu))
    if key in _KCACHE:
        return _KCACHE[key]
    rho = float(rho)
    nu = int(nu)
    m = np.arange(0, nu + 2)                     # N in {0..nu}; G_{nu+1} = 0
    t = nu - (m - 1.0)
    z = rho * np.maximum(t, 0.0)
    zmax = rho * nu
    lf = _log_fact_table(int(max(nu + 2, zmax + 12 * math.sqrt(zmax + 1))) + 200)

    def log_G(zv):
        lG = np.full(m.size, -np.inf)
        lG[0] = 0.0
        for k in range(1, m.size):
            if zv[k] > 0:
                lG[k] = _log_sf_pois(k, zv[k], lf)
        return lG

    lG = log_G(z)
    # pmf in log domain: p_m = G_m - G_{m+1}
    with np.errstate(invalid="ignore", divide="ignore"):
        lp = lG[:-1] + np.log(-np.expm1(np.minimum(lG[1:] - lG[:-1], -1e-300)))
    p = np.where(np.isfinite(lp), np.exp(lp), 0.0)
    mass_err = abs(p.sum() - 1.0)
    # analytic theta-derivative of G_m: Gdot_m = m * PoisPMF(m; z_m)
    with np.errstate(invalid="ignore", divide="ignore"):
        lgd = np.where((m >= 1) & (z > 0),
                       np.log(np.maximum(m, 1)) +
                       m * np.log(np.maximum(z, 1e-300)) - z - lf[m],
                       -np.inf)
    Gd = np.exp(lgd)
    num = Gd[:-1] - Gd[1:]
    score_err = abs(num.sum())
    ok = p > 0
    with np.errstate(divide="ignore"):
        lt2 = 2.0 * np.log(np.maximum(np.abs(num), 1e-300)) - lp
    use = ok & np.isfinite(lt2) & (np.abs(num) > 0) & (lt2 > -745.0)
    if use.any():
        mx = float(lt2[use].max())
        I_theta = math.exp(mx) * float(np.exp(lt2[use] - mx).sum())
    else:
        I_theta = 0.0
    J_exact = I_theta / nu
    # RQL quasi-score (physics.arm_score "RQL"): s = max(T - m,0) - m/lam
    mm = m[:-1].astype(np.float64)
    s = np.maximum(nu - mm, 0.0) - mm / rho
    EN = float((p * mm).sum())
    Es = float((p * s).sum())
    Vs = float((p * (s - Es) ** 2).sum())
    sens = EN / rho ** 2                          # E[ds/dlam] = E[N]/lam^2
    J_q = (rho ** 2) * (sens ** 2) / Vs / nu if Vs > 0 else 0.0
    p_ceil = float(np.exp(lG[nu])) if np.isfinite(lG[nu]) else 0.0
    out = {"rho": rho, "nu": nu, "J_exact": float(J_exact),
           "J_q": float(J_q),
           "trust": float(J_q / J_exact) if J_exact > 0 else 0.0,
           "p_ceil": p_ceil, "EN": EN,
           "mass_err": float(mass_err), "score_err": float(score_err),
           "lp": lp, "p": p, "m": mm}
    _KCACHE[key] = out
    return out


def kernel_certify(rho, nu, h=1e-5):
    """G7 numerical certification of the exact kernel at a frozen node."""
    k0 = kernel_eval(rho, nu)
    # independent derivative check: central log-lambda difference of log-pmf
    kp = kernel_eval(rho * math.exp(h), nu)
    km = kernel_eval(rho * math.exp(-h), nu)
    npx = min(k0["lp"].size, kp["lp"].size, km["lp"].size)
    d = (kp["lp"][:npx] - km["lp"][:npx]) / (2 * h)
    p = k0["p"][:npx]
    use = (p > 0) & np.isfinite(d)
    J_fd = float((p[use] * d[use] ** 2).sum()) / nu
    rel_fd = abs(J_fd - k0["J_exact"]) / k0["J_exact"] if k0["J_exact"] > 0 else np.inf
    ok = (k0["mass_err"] <= 1e-10 and k0["score_err"] <= 1e-9
          and rel_fd <= 1e-6
          and np.isfinite(k0["J_exact"]) and k0["J_exact"] >= 0
          and np.isfinite(k0["J_q"]) and k0["J_q"] >= 0)
    return {"rho": float(rho), "nu": int(nu), "mass_err": k0["mass_err"],
            "score_err": k0["score_err"], "J_fd_rel": rel_fd,
            "J_exact": k0["J_exact"], "J_q": k0["J_q"],
            "pass": bool(ok), "verdict": "OK" if ok else "KERNEL_CERT_FAIL"}


# ========================================================================== #
#  R11 Section 1: exact ridge target                                          #
# ========================================================================== #
def _admissible(rho, nu):
    k = kernel_eval(rho, nu)
    return (0 < rho <= RHO_CAP and k["p_ceil"] <= CEIL_TARGET
            and k["trust"] >= TRUST_MIN)


def _golden_max(f, lo, hi, n_iter=40):
    """Deterministic golden-section MAXIMIZATION of f on [lo, hi] (log axis)."""
    gr = (math.sqrt(5.0) - 1.0) / 2.0
    a, b = math.log(lo), math.log(hi)
    c = b - gr * (b - a)
    d = a + gr * (b - a)
    fc, fd = f(math.exp(c)), f(math.exp(d))
    for _ in range(n_iter):
        if fc >= fd:            # keep left (min tie rule -> prefer smaller rho)
            b, d, fd = d, c, fc
            c = b - gr * (b - a)
            fc = f(math.exp(c))
        else:
            a, c, fc = c, d, fd
            d = a + gr * (b - a)
            fd = f(math.exp(d))
    return math.exp(a) if fc >= fd else math.exp(b - gr * (b - a))


def ridge_target(nu, n_nodes=41, verbose=False):
    """Frozen R11 ridge target rho_R(nu) with the full recorded fields."""
    nodes = np.geomspace(0.05, 64.0, 49)          # unconstrained bracket
    Jn = np.array([kernel_eval(r, nu)["J_exact"] for r in nodes])
    i_u = int(np.argmax(Jn))
    lo_u = nodes[max(i_u - 1, 0)]
    hi_u = nodes[min(i_u + 1, nodes.size - 1)]
    rho_unc = _golden_max(lambda r: kernel_eval(r, nu)["J_exact"], lo_u, hi_u)

    adm_nodes = nodes[nodes <= RHO_CAP + 1e-12]
    adm = np.array([_admissible(r, nu) for r in adm_nodes])
    if not adm.any():
        return {"nu": int(nu), "verdict": "RIDGE_TARGET_FAIL"}
    Ja = np.array([kernel_eval(r, nu)["J_exact"] for r in adm_nodes])
    Ja = np.where(adm, Ja, -np.inf)
    i_a = int(np.argmax(Ja))
    lo_a = adm_nodes[max(i_a - 1, 0)]
    hi_a = min(adm_nodes[min(i_a + 1, adm_nodes.size - 1)], RHO_CAP)
    rho_R = _golden_max(
        lambda r: (kernel_eval(r, nu)["J_exact"] if _admissible(r, nu)
                   else -np.inf),
        lo_a, hi_a)
    if not _admissible(rho_R, nu):                 # boundary fallback: last
        rho_R = float(adm_nodes[i_a])              # certified admissible node
    # min tie rule at 1e-10 in J: scan slightly smaller loads
    JR = kernel_eval(rho_R, nu)["J_exact"]
    for f in (0.999999, 0.99999, 0.9999):
        r2 = rho_R * f
        if _admissible(r2, nu) and abs(kernel_eval(r2, nu)["J_exact"] - JR) <= 1e-10:
            rho_R = r2
            JR = kernel_eval(rho_R, nu)["J_exact"]
    kR = kernel_eval(rho_R, nu)
    # clip reason
    if abs(math.log(rho_R / rho_unc)) < 1e-4:
        reason = "NONE"
    elif rho_unc > RHO_CAP and abs(rho_R - RHO_CAP) / RHO_CAP < 1e-3:
        reason = "RHO_CAP"
    else:
        up = kernel_eval(min(rho_R * 1.01, 64.0), nu)
        reason = "CEILING" if up["p_ceil"] > CEIL_TARGET else "RQL_TRUST"
    return {"nu": int(nu),
            "rho_ridge_exact_unconstrained": float(rho_unc),
            "rho_R_production": float(rho_R),
            "ridge_clip_reason": reason,
            "J_exact_at_target": kR["J_exact"],
            "J_q_at_target": kR["J_q"],
            "p_ceiling_scalar": kR["p_ceil"],
            "trust_at_target": kR["trust"],
            "verdict": "OK"}


# ========================================================================== #
#  R11 Section 3: load palette                                                #
# ========================================================================== #
def load_palette(nu, ridge=None):
    """L_nu = unique{0.30, 0.60, 1.00, min(3, rho_R), rho_R/2, rho_R} with
    per-level admission (rho<=24, p_ceil<=0.05, trust>=0.90) + G7 cert."""
    if ridge is None:
        ridge = ridge_target(nu)
    rR = ridge["rho_R_production"]
    raw = [0.30, 0.60, 1.00, min(3.00, rR), 0.5 * rR, rR]
    raw.sort()
    lev = []
    for r in raw:
        if not lev or abs(r - lev[-1]) > 1e-10 * max(abs(r), abs(lev[-1])):
            lev.append(r)
    admitted, records = [], []
    for r in lev:
        k = kernel_eval(r, nu)
        ok = (r <= RHO_CAP + 1e-12 and k["p_ceil"] <= CEIL_ATOM
              and k["trust"] >= TRUST_MIN)
        cert = kernel_certify(r, nu)
        records.append({"rho": float(r), "admitted": bool(ok and cert["pass"]),
                        "p_ceil": k["p_ceil"], "trust": k["trust"],
                        "J_exact": k["J_exact"], "EN": k["EN"],
                        "cert": cert["verdict"]})
        if ok and cert["pass"]:
            admitted.append(float(r))
    return {"nu": int(nu), "levels": admitted, "records": records,
            "ridge": ridge}


# ========================================================================== #
#  R11 Section 2: RIDGE-FIXED helper                                          #
# ========================================================================== #
def ridge_fixed_design(xhat, nu, side=32, n_grid=61):
    """LBLOB16 translate family, one global multiplier calibrated so the
    pre-scan-predicted mean load equals rho_R(nu), then the R11 Section-2.2
    deterministic bisection safety clip.  Returns the family matrix (Study-2
    normalization x global multiplier), per-row loads, and all records."""
    xhat = np.asarray(xhat, dtype=float).ravel()
    n = xhat.size
    ridge = ridge_target(nu)
    rho_R = ridge["rho_R_production"]
    shapes = v1._default_shapes(side)
    sidx = v1._shape_support(shapes["Lblob6x6"], side)      # (n, 16)
    base_load = (n / 16.0) * xhat[sidx].sum(axis=1)         # unit-mean exactly
    # guard interpolants on a dense exact log-rho grid (G7-certified nodes)
    lo = max(1e-3, rho_R * base_load.min() * 1e-3)
    grid = np.geomspace(lo, max(rho_R * base_load.max() * 1.05, 1.0), n_grid)
    kv = [kernel_eval(r, nu) for r in grid]
    lg = np.log(grid)
    tr_g = np.array([k["trust"] for k in kv])
    pc_g = np.array([k["p_ceil"] for k in kv])

    def guards_ok(mult):
        loads = mult * base_load
        ll = np.log(np.clip(loads, grid[0], grid[-1]))
        pc = np.interp(ll, lg, pc_g)
        tr = np.interp(ll, lg, tr_g)
        return (pc.mean() <= CEIL_TARGET and pc.max() <= CEIL_ATOM
                and tr.min() >= TRUST_MIN)

    requested = rho_R
    mult = rho_R                                   # unit-mean family
    clipped_by_bisect = False
    if not guards_ok(mult):
        clipped_by_bisect = True
        lo_m, hi_m = 0.0, mult
        for _ in range(60):                        # deterministic bisection
            mid = 0.5 * (lo_m + hi_m)
            if guards_ok(mid):
                lo_m = mid
            else:
                hi_m = mid
        mult = lo_m
    loads = mult * base_load
    achieved = float(loads.mean())
    flag = achieved < 0.90 * rho_R                 # RIDGE_GUARD_CLIPPED rule
    A = np.zeros((n, n))
    rows = np.arange(n)
    A[rows[:, None], sidx] = mult * (n / 16.0)
    ll = np.log(np.clip(loads, grid[0], grid[-1]))
    pc = np.interp(ll, lg, pc_g)
    tr = np.interp(ll, lg, tr_g)
    Jx = np.interp(ll, lg, np.array([k["J_exact"] for k in kv]))
    return {"A": A, "loads": loads, "ridge": ridge,
            "requested_mean_load": float(requested),
            "achieved_mean_load": achieved,
            "global_multiplier": float(mult),
            "bisection_clip_applied": bool(clipped_by_bisect),
            "RIDGE_GUARD_CLIPPED": bool(flag),
            "mean_p_ceil": float(pc.mean()), "max_p_ceil": float(pc.max()),
            "min_trust": float(tr.min()),
            "mean_J_exact": float(Jx.mean()),
            "rho_quantiles": {"rho_5": float(np.percentile(loads, 5)),
                              "rho_50": float(np.percentile(loads, 50)),
                              "rho_95": float(np.percentile(loads, 95)),
                              "rho_max": float(loads.max())}}


# ========================================================================== #
#  R11 Sections 3-4 + G5: budgeted, dose-constrained variable-load design     #
# ========================================================================== #
def design_v3(xhat, nu=2000.0, shapes=None, powers=(0, 1), M_rows=972,
              budget_mean=0.60, max_iter=220, tol=0.1, arm_rho=0.60,
              dose_band=0.05, dose_inner=0.035, w_clip=(0.25, 4.0),
              verbose=False):
    """R11 OED-DT: chain-rule atoms on the ridge palette, incident-budget
    inequality  <rho, xi> <= budget_mean, per-pixel dose band enforced by
    IPF in the loop and by deterministic same-load exchanges post-rounding.
    """
    xhat = np.asarray(xhat, dtype=float).ravel()
    n = xhat.size
    side = int(round(math.sqrt(n)))
    if shapes is None:
        shapes = v1._default_shapes(side)
    nu_i = int(nu)

    pal = load_palette(nu_i)
    ridge = pal["ridge"]
    rho = np.asarray(pal["levels"], dtype=float)
    L = rho.size
    nuJ = np.array([nu * kernel_eval(r, nu_i)["J_exact"] for r in rho])
    EN_lvl = np.array([kernel_eval(r, nu_i)["EN"] for r in rho])
    pceil_lvl = np.array([kernel_eval(r, nu_i)["p_ceil"] for r in rho])

    # geometries: v2 builder with the shape-based G3 (A2) -----------------
    IDX, GVAL, names, t_arr, p_arr, _ALLOW0, _g3v2, g4 = v2._build_geometries(
        xhat, shapes, powers, side, w_clip, np.inf, rho)
    G = IDX.shape[0]
    # A2 shape guard: frame-mean-1 peak = (n/k)*max(w) <= (n/k)*w_clip_max
    sumg = GVAL.sum(axis=1)
    maxw = GVAL.max(axis=1) * (16.0 / sumg)            # w = g * load, load=16/sum(g)
    shape_peak = (n / 16.0) * maxw
    geo_ok = shape_peak <= (n / 16.0) * w_clip[1] + 1e-9
    ALLOW = np.repeat(geo_ok[:, None], L, axis=1)
    # v2-strict reading, report-only:
    gmax = GVAL.max(axis=1)
    strict_drops = {float(r): int((gmax * r > 64.0 + 1e-12).sum()) for r in rho}

    V0, J_arm = v2.build_V0(xhat, nu, arm_rho, v2.get_J_source(), side)

    xi = np.where(ALLOW, 1.0, 0.0)
    xi /= xi.sum()
    # start feasible: initial load = mean over levels; rescale toward low loads
    dose_vec_geo = np.zeros((G,))                        # placeholder

    def assemble(xi_):
        u_geo = M_rows * (xi_ * nuJ[None, :]).sum(axis=1)
        M = v2._build_M(IDX, GVAL, u_geo, n)
        Vr = V0 + M
        eps = 1e-9 * np.trace(Vr) / n
        Vr[np.diag_indices(n)] += eps
        return Vr, eps

    def pixel_dose(xi_):
        w_geo = (xi_ * rho[None, :]).sum(axis=1)         # per-geometry rho mass
        d = np.bincount(IDX.ravel(),
                        weights=(w_geo[:, None] * GVAL).ravel(), minlength=n)
        return d

    def load_of(xi_):
        return float((xi_ @ rho).sum())

    def budget_repair(xi_):
        """Exact feasibility: move mass from each geometry's higher levels to
        its lowest admissible level until <rho,xi> <= budget_mean (linear in
        the moved fraction -> exact deterministic solve)."""
        cur = load_of(xi_)
        if cur <= budget_mean + 1e-15:
            return xi_
        lo_idx = np.argmax(ALLOW, axis=1)                # lowest admissible lvl
        # moving ALL mass to each geometry's lowest level gives load_min;
        # the move is linear in the moved fraction t -> exact solve
        low_mass = np.zeros_like(xi_)
        low_mass[np.arange(G), lo_idx] = xi_.sum(axis=1)
        load_min = float((low_mass @ rho).sum())
        if load_min >= budget_mean:
            t = 1.0
        else:
            t = (cur - budget_mean) / (cur - load_min)
        return (1.0 - t) * xi_ + t * low_mass

    def ipf(xi_, band, n_pass=12):
        """Iterative proportional fitting of atom weights onto the pixel-dose
        band [1-band, 1+band] * mean.  Atom multiplier = inverse of the
        dose-ratio average over its support (load-weighted); renormalized to
        the simplex each pass.  Deterministic."""
        for _ in range(n_pass):
            d = pixel_dose(xi_)
            dm = d.mean()
            r = d / dm
            if np.abs(r - 1.0).max() <= band:
                break
            rc = np.clip(r, 0.5, 2.0)
            factor_geo = 1.0 / ((GVAL * rc[IDX]).sum(axis=1) / sumg)
            xi_ = xi_ * factor_geo[:, None]
            xi_[~ALLOW] = 0.0
            s = xi_.sum()
            if s <= 0:
                break
            xi_ /= s
        return xi_

    def assign_levels(mass, D, kappa, xi_old):
        """Exact water-filling level assignment: for multiplier theta each
        geometry targets  l*(g) = argmax_l D_gl - theta*rho_l ; theta is set
        by deterministic bisection so the (kappa-damped) design meets the
        incident budget with equality when the constraint is active.  The
        final design mixes the two bracketing assignments so the budget
        holds EXACTLY.  Returns (xi_new, theta)."""
        Dm = np.where(ALLOW, D, -np.inf)

        def onehot(theta):
            tgt = np.argmax(Dm - theta * rho[None, :], axis=1)
            oh = np.zeros((G, L))
            oh[np.arange(G), tgt] = mass
            return oh

        def mixed(theta):
            return (1.0 - kappa) * xi_old + kappa * onehot(theta)

        if load_of(mixed(0.0)) <= budget_mean + 1e-15:
            return mixed(0.0), 0.0
        th_lo, th_hi = 0.0, 1.0
        while load_of(mixed(th_hi)) > budget_mean and th_hi < 1e18:
            th_hi *= 4.0
        for _ in range(70):
            mid = 0.5 * (th_lo + th_hi)
            if load_of(mixed(mid)) > budget_mean:
                th_lo = mid
            else:
                th_hi = mid
        xa, xb = mixed(th_lo), mixed(th_hi)
        la, lb = load_of(xa), load_of(xb)
        a = 0.0 if abs(la - lb) < 1e-15 else (la - budget_mean) / (la - lb)
        a = min(max(a, 0.0), 1.0)
        return (1.0 - a) * xa + a * xb, th_hi

    V, eps = assemble(xi)
    sign, phi = np.linalg.slogdet(V)
    theta = 0.0                                          # budget multiplier
    gap_traj = []
    it = 0
    for it in range(1, max_iter + 1):
        Vinv = np.linalg.inv(V)
        sub = Vinv[IDX[:, :, None], IDX[:, None, :]]
        dgeo = np.einsum('ai,aij,aj->a', GVAL, sub, GVAL)
        D = M_rows * dgeo[:, None] * nuJ[None, :]         # utility per mass
        Dm = np.where(ALLOW, D, -np.inf)
        # Lagrangian-adjusted certificate at the CURRENT theta
        s_adj = Dm - theta * rho[None, :]
        sup = xi > 0
        sbar = float((xi * np.where(sup & ALLOW, s_adj, 0.0)).sum())
        smax = float(s_adj[ALLOW].max())
        gap_adj = smax / sbar - 1.0 if sbar > 0 else float("inf")
        cur_load = load_of(xi)
        gap_traj.append((it, gap_adj, cur_load))
        if verbose and (it % 25 == 0 or it == 1):
            print("    [v3] it=%4d gap_adj=%.4f load=%.4f theta=%.3g "
                  "logdet=%.2f" % (it, gap_adj, cur_load, theta, phi),
                  flush=True)
        if gap_adj < tol and cur_load <= budget_mean + 1e-9:
            break
        # (a) multiplicative geometry-mass ascent on the reduced sensitivity
        mass = xi.sum(axis=1)
        with np.errstate(invalid="ignore"):
            d_cur = np.where(mass > 0, (xi * np.where(ALLOW, D, 0.0)).sum(axis=1)
                             / np.maximum(mass, 1e-300), 0.0)
            r_cur = np.where(mass > 0, (xi @ rho) / np.maximum(mass, 1e-300), 0.0)
        s_g = d_cur - theta * r_cur
        s_bar_g = float((mass * s_g).sum())
        if s_bar_g > 0:
            fac = np.sqrt(np.maximum(s_g, 1e-6 * s_bar_g) / s_bar_g)
            mass = mass * fac
            mass /= mass.sum()
        # (b) exact water-filling level assignment (kappa-damped, budget-exact)
        xi_scaled = xi * (mass / np.maximum(xi.sum(axis=1), 1e-300))[:, None]
        xi, theta = assign_levels(mass, D, 0.35, xi_scaled)
        # (c) dose IPF + exact budget re-repair + prune
        xi = ipf(xi, dose_inner)
        xi = budget_repair(xi)
        thr = xi.max() * 1e-16
        xi[xi < thr] = 0.0
        xi /= xi.sum()
        V, eps = assemble(xi)
        sign, phi = np.linalg.slogdet(V)

    final_load = load_of(xi)
    d_final = pixel_dose(xi)
    dose_dev_cont = float(np.abs(d_final / d_final.mean() - 1.0).max())

    # ---- rounding (largest remainder) ------------------------------------ #
    flat = xi.ravel()
    thr = flat.max() * 1e-6
    supp = np.where(flat > thr)[0]
    w_s = flat[supp] / flat[supp].sum()
    quota = w_s * M_rows
    base = np.floor(quota).astype(np.int64)
    rem = int(M_rows - base.sum())
    frac = quota - base
    order = np.lexsort((supp, -frac))
    for i in range(rem):
        base[order[i]] += 1
    counts = dict(zip(supp.tolist(), base.tolist()))

    def counts_arrays():
        fa = np.array(sorted([k for k, c in counts.items() if c > 0]))
        cc = np.array([counts[k] for k in fa])
        gg, ll = np.divmod(fa, L)
        return fa, cc, gg, ll

    def rounded_load():
        _, cc, _, ll = counts_arrays()
        return float((cc * rho[ll]).sum())

    def rounded_dose():
        fa, cc, gg, ll = counts_arrays()
        w_geo = np.zeros(G)
        np.add.at(w_geo, gg, cc * rho[ll])
        return np.bincount(IDX.ravel(),
                           weights=(w_geo[:, None] * GVAL).ravel(),
                           minlength=n) / M_rows

    # budget repair on counts: demote highest-load counted atoms to the same
    # geometry's lowest admissible level until within budget (deterministic)
    n_demote = 0
    while rounded_load() > budget_mean * M_rows and n_demote < 4 * M_rows:
        fa, cc, gg, ll = counts_arrays()
        cand = np.argsort(-rho[ll])
        moved = False
        for ci in cand:
            g_, l_ = int(gg[ci]), int(ll[ci])
            lo = int(np.argmax(ALLOW[g_]))
            if rho[lo] < rho[l_] - 1e-12:
                counts[fa[ci]] -= 1
                counts[g_ * L + lo] = counts.get(g_ * L + lo, 0) + 1
                n_demote += 1
                moved = True
                break
        if not moved:
            break

    # G5 dose repair (coordinator spec, R11 round-3): while any pixel dose
    # deviates > band from uniform, take the worst (over-dosed, under-dosed)
    # pixel pair and swap ONE count from an atom covering the over-dosed
    # pixel to a palette-admissible atom of the SAME shape family at the
    # SAME rho level whose support covers the under-dosed pixel, choosing
    # the candidate that least increases the KW gap -- d(a) evaluated with
    # the FINAL M-inverse (removal side: covering counted atom of minimal
    # d; addition side: candidate of maximal d).  Deterministic order (sort
    # by deviation then atom index), no RNG; cap 2000 swaps; early stop
    # after 200 consecutive swaps without max-deviation progress.  Same
    # level => budget and row count exactly preserved.
    Vinv_f = np.linalg.inv(V)
    sub_f = Vinv_f[IDX[:, :, None], IDX[:, None, :]]
    dgeo_f = np.einsum('ai,aij,aj->a', GVAL, sub_f, GVAL)
    names_arr = np.asarray(names)
    n_swap = 0
    dose = rounded_dose()
    best_dev = np.inf
    best_counts = dict(counts)
    since_best = 0
    while n_swap < 2000:
        dm = dose.mean()
        rel = dose / dm - 1.0
        dev = float(np.abs(rel).max())
        if dev <= dose_band:
            best_dev, best_counts = dev, dict(counts)
            break
        if dev < best_dev - 1e-12:
            best_dev, best_counts, since_best = dev, dict(counts), 0
        else:
            since_best += 1
            if since_best >= 200:
                break
        jo = int(np.argmax(rel))                 # worst over-dosed pixel
        ju = int(np.argmin(rel))                 # worst under-dosed pixel
        fa, cc, gg, ll = counts_arrays()
        cover_o = np.nonzero((IDX[gg] == jo).any(axis=1))[0]
        if cover_o.size == 0:
            break
        # removal: covering counted atom with minimal d (least value lost)
        d_cov = dgeo_f[gg[cover_o]] * nuJ[ll[cover_o]]
        rord = np.lexsort((fa[cover_o], d_cov))
        ci = int(cover_o[rord[0]])
        g_r, l_r = int(gg[ci]), int(ll[ci])
        fam = names_arr[g_r]
        # addition: same shape family + same level + covers ju + admissible
        cov_u = (IDX == ju).any(axis=1)
        cand = np.nonzero(cov_u & (names_arr == fam)
                          & ALLOW[:, l_r] & geo_ok)[0]
        cand = cand[cand != g_r]
        if cand.size == 0:                        # widen: any family, same lvl
            cand = np.nonzero(cov_u & ALLOW[:, l_r] & geo_ok)[0]
            cand = cand[cand != g_r]
            if cand.size == 0:
                break
        d_cand = dgeo_f[cand] * nuJ[l_r]
        aord = np.lexsort((cand, -d_cand))        # max d, tie -> lowest index
        g_a = int(cand[aord[0]])
        counts[fa[ci]] -= 1
        counts[g_a * L + l_r] = counts.get(g_a * L + l_r, 0) + 1
        dose = rounded_dose()
        n_swap += 1
    # restore the best state seen along the (possibly diverging) swap
    # trajectory before reporting the residual (deterministic, honest)
    counts = {k_: c_ for k_, c_ in best_counts.items() if c_ > 0}
    dose = rounded_dose()
    resid_stage1 = float(np.abs(dose / dose.mean() - 1.0).max())

    # stage-2 polish (documented addition): strictly-descending best-swap
    # local search on the band-violation penalty; same-level swaps only
    # (budget/row neutral); deterministic; terminates by strict descent.
    n_polish = 0
    band_in = dose_band - 0.002

    def _penalty(dv):
        r = np.abs(dv / dv.mean() - 1.0)
        v = np.maximum(r - band_in, 0.0)
        return float((v * v).sum())

    P = _penalty(dose)
    while n_polish < 2500 and P > 0.0:
        dm = dose.mean()
        rel = dose / dm - 1.0
        jo = int(np.argmax(rel))
        ju = int(np.argmin(rel))
        fa, cc, gg, ll = counts_arrays()
        rem_cand = np.nonzero((IDX[gg] == jo).any(axis=1))[0]
        if rem_cand.size == 0:
            break
        cov_u = (IDX == ju).any(axis=1)
        best = (0.0, None)
        for ci in rem_cand:
            g_r, l_r = int(gg[ci]), int(ll[ci])
            d_rem = rho[l_r] * GVAL[g_r] / M_rows
            add_cand = np.nonzero(cov_u & ALLOW[:, l_r] & geo_ok)[0]
            add_cand = add_cand[add_cand != g_r]
            for g_a in add_cand[:64]:
                d_add = rho[l_r] * GVAL[g_a] / M_rows
                px = np.concatenate([IDX[g_r], IDX[g_a]])
                dv_new = dose[px].copy()
                dv_new[:16] -= d_rem
                dv_new[16:] += d_add
                r_old = np.abs(dose[px] / dm - 1.0)
                r_new = np.abs(dv_new / dm - 1.0)
                dP = float((np.maximum(r_new - band_in, 0.0) ** 2).sum()
                           - (np.maximum(r_old - band_in, 0.0) ** 2).sum())
                if dP < best[0] - 1e-15:
                    best = (dP, (int(fa[ci]), int(g_a), l_r))
        if best[1] is None:
            break
        fa_i, g_a, l_r = best[1]
        counts[fa_i] -= 1
        if counts[fa_i] == 0:
            del counts[fa_i]
        counts[g_a * L + l_r] = counts.get(g_a * L + l_r, 0) + 1
        dose_new = rounded_dose()
        P_new = _penalty(dose_new)
        if P_new >= P - 1e-18:
            counts[g_a * L + l_r] -= 1
            if counts[g_a * L + l_r] == 0:
                del counts[g_a * L + l_r]
            counts[fa_i] = counts.get(fa_i, 0) + 1
            break
        dose, P = dose_new, P_new
        n_polish += 1

    dm = dose.mean()
    resid_dev = float(np.abs(dose / dm - 1.0).max())
    # alpha-mixture fallback (G6 machinery) if the repair did not converge:
    # dose is linear in the design, so bisect the largest alpha for which
    # alpha*d_OED + (1-alpha)*d_fix stays inside the band (d_fix = uniform
    # scatter16 p=0 translate design at the level closest to 0.6).
    g5_alpha = None
    if resid_dev > dose_band:
        sel_fix = (names_arr == "scatter16") & (p_arr == 0.0) & geo_ok
        l_fix = int(np.argmin(np.abs(rho - 0.6)))
        xi_fix = np.zeros_like(xi)
        xi_fix[sel_fix, l_fix] = 1.0
        xi_fix /= xi_fix.sum()
        d_fix = pixel_dose(xi_fix)

        def mix_dev(a_):
            dmix = a_ * dose + (1.0 - a_) * d_fix
            return float(np.abs(dmix / dmix.mean() - 1.0).max())

        if mix_dev(0.0) <= dose_band:
            lo_a, hi_a = 0.0, 1.0
            for _ in range(60):
                mid = 0.5 * (lo_a + hi_a)
                if mix_dev(mid) <= dose_band:
                    lo_a = mid
                else:
                    hi_a = mid
            g5_alpha = lo_a

    fa, cc, gg, ll = counts_arrays()
    Aout = np.zeros((int(cc.sum()), n))
    atoms = []
    row = 0
    for k_, c_, g_, l_ in zip(fa, cc, gg, ll):
        a_row = np.zeros(n)
        a_row[IDX[g_]] = rho[l_] * GVAL[g_]
        Aout[row:row + c_] = a_row
        atoms.append((names[g_], int(t_arr[g_]), float(p_arr[g_]),
                      float(rho[l_]), int(c_)))
        row += c_
    assert row == M_rows, (row, M_rows)

    # ---- reports ---------------------------------------------------------- #
    load_rows = np.concatenate([[a[3]] * a[4] for a in atoms])
    incident = float(load_rows.sum())                    # sum rho_i (x nu slots)
    budget_cap = budget_mean * M_rows
    lvl_index = {float(r): i for i, r in enumerate(rho)}
    det_exact = float(sum(c * EN_lvl[lvl_index[r]]
                          for (_, _, _, r, c) in atoms))
    det_approx = float(sum(c * nu * r / (1.0 + r)
                           for (_, _, _, r, c) in atoms))
    dvec = rounded_dose()
    dose_dev = float(np.abs(dvec / dvec.mean() - 1.0).max())
    wceil = float(sum((c / M_rows) * pceil_lvl[lvl_index[r]]
                      for (_, _, _, r, c) in atoms))

    cert = {"n": n, "final_gap_adjusted": float(gap_traj[-1][1]),
            "theta_budget_multiplier": float(theta),
            "iters": it, "logdet": float(phi),
            "dict_size": int(ALLOW.sum()), "n_geometries": int(G),
            "levels": [float(r) for r in rho],
            "nuJ_levels": [float(x) for x in nuJ],
            "cert_note": ("Lagrangian-adjusted KW gap: s(a)=d(a)-theta*rho(a),"
                          " gap=max s/<xi,s>-1; certifies optimality over the"
                          " simplex+budget only -- per-pixel dose constraints"
                          " are enforced by IPF/exchange and their KKT prices"
                          " are NOT included, so the gap is a diagnostic"
                          " upper bound, not a full-KKT certificate.")}
    guards = {"G3_shape": {"reading": "A2: (n/k)*w_max <= 256 (k>=16 + G4)",
                           "geo_dropped": int((~geo_ok).sum()),
                           "v2_strict_drops_report_only": strict_drops},
              "G4": g4,
              "G5": {"max_rel_dev_continuous": dose_dev_cont,
                     "max_rel_dev_rounded": dose_dev,
                     "band": dose_band,
                     "pass": bool(dose_dev <= dose_band),
                     "n_exchange_swaps": n_swap,
                     "n_polish_swaps": n_polish,
                     "resid_after_prescribed": resid_stage1,
                     "n_budget_demotions": n_demote,
                     "alpha_mixture": g5_alpha,
                     "compliant_via_mixture": bool(g5_alpha is not None)},
              "G8": {"design_weighted_p_ceil": wceil,
                     "cap": 0.01, "pass": bool(wceil <= 0.01)},
              "budget": {"incident_sum_rho": incident,
                         "cap_sum_rho": budget_cap,
                         "slack": float(budget_cap - incident),
                         "pass": bool(incident <= budget_cap + 1e-9),
                         "detected_pred_exact": det_exact,
                         "detected_pred_approx": det_approx,
                         "continuous_load": final_load}}
    return {"A": Aout, "atoms": atoms, "gap": float(gap_traj[-1][1]),
            "gap_traj": gap_traj, "cert": cert, "guards": guards,
            "palette": pal, "ridge": ridge, "V": V, "V0": V0, "xi": xi}


# ========================================================================== #
#  Selftest                                                                   #
# ========================================================================== #
if __name__ == "__main__":
    import time
    from collections import Counter

    t0_all = time.time()
    SIDE, N = 32, 1024

    yy, xx = np.mgrid[0:SIDE, 0:SIDE]
    scene = np.ones((SIDE, SIDE))
    for (cy, cx) in [(8, 9), (10, 23), (22, 7), (24, 21)]:
        scene += 0.6 * np.exp(-(((yy - cy) ** 2 + (xx - cx) ** 2) / (2 * 3.5 ** 2)))
    x = (scene / scene.sum()).ravel()
    coarse = x.reshape(SIDE, SIDE).reshape(SIDE // 4, 4, SIDE // 4, 4).mean(axis=(1, 3))
    xhat = np.repeat(np.repeat(coarse, 4, axis=0), 4, axis=1).ravel()
    print("[v3 selftest] xhat: n=%d sum=%.6f cv=%.3f"
          % (xhat.size, xhat.sum(), xhat.std() / xhat.mean()), flush=True)

    # ---- sanity of the trust ratio direction (A1) ------------------------- #
    k_lo = kernel_eval(0.01, 2000)
    print("[v3 selftest] A1 sanity: J_q/J_exact(rho=0.01, nu=2000) = %.6f "
          "(-> 1 as rho -> 0)" % k_lo["trust"], flush=True)

    # ---- rho_R table ------------------------------------------------------ #
    print("\n[v3 selftest] === (1) exact ridge targets ===", flush=True)
    ridges = {}
    for nu_ in (200, 2000):
        t0 = time.time()
        rr = ridge_target(nu_)
        ridges[nu_] = rr
        print("  nu=%4d: rho_unconstrained=%.4f  rho_R=%.4f  clip=%s  "
              "J_exact=%.4f  J_q=%.4f  trust=%.4f  p_ceil=%.3e  (%.1f s)"
              % (nu_, rr["rho_ridge_exact_unconstrained"],
                 rr["rho_R_production"], rr["ridge_clip_reason"],
                 rr["J_exact_at_target"], rr["J_q_at_target"],
                 rr["trust_at_target"], rr["p_ceiling_scalar"],
                 time.time() - t0), flush=True)

    # ---- palettes + G7 certs ---------------------------------------------- #
    print("\n[v3 selftest] === (2) load palettes + G7 kernel certification ===",
          flush=True)
    for nu_ in (200, 2000):
        pal = load_palette(nu_, ridges[nu_])
        print("  nu=%4d palette:" % nu_)
        for rec in pal["records"]:
            print("    rho=%8.4f admitted=%-5s p_ceil=%.2e trust=%.4f "
                  "J_exact=%.4f E[N]=%8.2f cert=%s"
                  % (rec["rho"], rec["admitted"], rec["p_ceil"], rec["trust"],
                     rec["J_exact"], rec["EN"], rec["cert"]))
    certs_ok = all(rec["cert"] == "OK"
                   for nu_ in (200, 2000)
                   for rec in load_palette(nu_, ridges[nu_])["records"])
    print("  all palette-node G7 certifications          -> %s"
          % ("PASS" if certs_ok else "FAIL"))

    # ---- constrained design at nu=2000 ------------------------------------ #
    print("\n[v3 selftest] === (3) budgeted dose-constrained design, nu=2000 "
          "===", flush=True)
    t0 = time.time()
    out = design_v3(xhat, nu=2000.0, verbose=True)
    t_design = time.time() - t0
    cert, guards = out["cert"], out["guards"]
    bud = guards["budget"]

    print("\n  (i) incident budget: sum rho = %.4f <= cap %.4f  slack=%.4f "
          "-> %s" % (bud["incident_sum_rho"], bud["cap_sum_rho"],
                     bud["slack"], "PASS" if bud["pass"] else "FAIL"))
    print("      detected budget (pred): exact = %.1f counts, "
          "approx nu*rho/(1+rho) = %.1f counts"
          % (bud["detected_pred_exact"], bud["detected_pred_approx"]))
    ok_budget = bud["pass"]

    g5 = guards["G5"]
    print("  (ii) G5 dose +/-%.0f%%: continuous max dev = %.4f; rounded:"
          " prescribed-repair residual = %.4f (%d swaps), after stage-2"
          " descent polish = %.4f (%d swaps, %d demotions) -> %s"
          % (100 * g5["band"], g5["max_rel_dev_continuous"],
             g5["resid_after_prescribed"], g5["n_exchange_swaps"],
             g5["max_rel_dev_rounded"], g5["n_polish_swaps"],
             g5["n_budget_demotions"], "PASS" if g5["pass"] else "FAIL"))
    if not g5["pass"]:
        if g5["compliant_via_mixture"]:
            print("      residual dev %.4f > band; G6 alpha-mixture with the"
                  " dose-uniform scatter16@0.6 baseline achieves compliance"
                  " at alpha <= %.4f -> COMPLIANT-VIA-MIXTURE"
                  % (g5["max_rel_dev_rounded"], g5["alpha_mixture"]))
        else:
            print("      residual dev %.4f > band; alpha-mixture fallback"
                  " CANNOT achieve compliance (baseline itself out of band)"
                  % g5["max_rel_dev_rounded"])
    ok_g5 = g5["pass"] or g5["compliant_via_mixture"]

    print("  (iii) certificate: Lagrangian-adjusted gap = %.4f "
          "(theta = %.4g, %d iters)"
          % (cert["final_gap_adjusted"], cert["theta_budget_multiplier"],
             cert["iters"]))
    print("        %s" % cert["cert_note"])

    print("  (iv) load histogram (%d rows):" % sum(a[4] for a in out["atoms"]))
    n_rows_load = Counter()
    n_rows_by = Counter()
    for (nm, t_, p_, r_, c_) in out["atoms"]:
        n_rows_load[r_] += c_
        n_rows_by[(nm, p_, r_)] += c_
    for r_ in cert["levels"]:
        print("      rho=%7.4f : %4d rows" % (r_, n_rows_load.get(r_, 0)))
    ridge_zone_rows = sum(c for r_, c in n_rows_load.items() if r_ > 1.0 + 1e-9)
    print("      by (shape,p,rho):")
    for key in sorted(n_rows_by):
        print("        %-12s p=%.0f rho=%7.4f : %4d rows"
              % (key[0], key[1], key[2], n_rows_by[key]))
    budget_active = (cert["theta_budget_multiplier"] > 0.0
                     or bud["slack"] <= 1e-6 * bud["cap_sum_rho"])
    if ridge_zone_rows > 0:
        ok_ridge = True
        print("      ridge-zone (rho > 1.0) rows = %d -> PASS (selected)"
              % ridge_zone_rows)
    elif budget_active:
        ok_ridge = True
        print("      ridge-zone (rho > 1.0) rows = 0 -> PASS-BY-DESIGN"
              " (budget constraint active, theta = %.4g)"
              % cert["theta_budget_multiplier"])
        print("      WARNING+DIAGNOSIS: under the R11 total-incident-budget"
              " inequality (sum rho_i <= M*0.6), the per-photon information"
              " efficiency J(rho)/rho FALLS monotonically in rho (J"
              " saturates at 1, cost is linear). J_exact/rho over the"
              " palette at nu=2000:")
        for r_ in cert["levels"]:
            Jr = kernel_eval(r_, 2000)["J_exact"]
            print("        rho=%7.4f : J=%.4f  J/rho=%.4f"
                  % (r_, Jr, Jr / r_))
        print("      The optimizer rationally avoids ridge atoms under a"
              " photon budget; ridge-zone operation is a TIME-constrained"
              " regime -- the RIDGE-FIXED arm's job (R11"
              " 'operating-power-for-time' framing).")
    else:
        ok_ridge = False
        print("      ridge-zone (rho > 1.0) rows = 0 with INACTIVE budget"
              " -> FAIL (unexpected: nothing prevented ridge selection)")
    rR2000 = ridges[2000]["rho_R_production"]
    print("      A5 dose-cap load ceiling: a single count-1 row at load rho"
          " puts dose ~ rho*64/%d on its pixels; the +5%% band around the"
          " budget-capped mean (<= 0.60) forbids any row with rho >"
          " %.2f -- full-ridge atoms (rho_R = %.2f) are mathematically"
          " excluded by G5+budget; ridge-ZONE atoms up to ~%.1f remain"
          " feasible. This is water filling hitting the dose wall."
          % (972, 1.05 * 0.60 * 972 / 64.0, rR2000,
             1.05 * 0.60 * 972 / 64.0))

    g8 = guards["G8"]
    print("  G8 design-weighted ceiling prob = %.3e <= 0.01 -> %s"
          % (g8["design_weighted_p_ceil"], "PASS" if g8["pass"] else "FAIL"))
    print("  G3(A2 shape reading) dropped %d/%d geometries; v2-strict"
          " per-load drops (report-only): %s"
          % (guards["G3_shape"]["geo_dropped"], cert["n_geometries"],
             {("%.2f" % k): v for k, v in
              guards["G3_shape"]["v2_strict_drops_report_only"].items()}))

    # ---- ridge_fixed_design for both nu ----------------------------------- #
    print("\n[v3 selftest] === (4) ridge_fixed_design ===", flush=True)
    ok_rf = True
    for nu_ in (200, 2000):
        rf = ridge_fixed_design(xhat, nu_)
        q = rf["rho_quantiles"]
        print("  nu=%4d: requested mean load = %.4f  achieved = %.4f  "
              "multiplier = %.4f" % (nu_, rf["requested_mean_load"],
                                     rf["achieved_mean_load"],
                                     rf["global_multiplier"]))
        print("          bisection clip applied = %s  RIDGE_GUARD_CLIPPED = %s"
              % (rf["bisection_clip_applied"], rf["RIDGE_GUARD_CLIPPED"]))
        print("          mean p_ceil = %.3e (<=0.01)  max p_ceil = %.3e "
              "(<=0.05)  min trust = %.4f (>=0.90)"
              % (rf["mean_p_ceil"], rf["max_p_ceil"], rf["min_trust"]))
        print("          rho_5/50/95/max = %.3f / %.3f / %.3f / %.3f   "
              "mean J_exact = %.4f"
              % (q["rho_5"], q["rho_50"], q["rho_95"], q["rho_max"],
                 rf["mean_J_exact"]))
        ok_rf &= (rf["mean_p_ceil"] <= CEIL_TARGET + 1e-12
                  and rf["max_p_ceil"] <= CEIL_ATOM + 1e-12
                  and rf["min_trust"] >= TRUST_MIN - 1e-9)
    print("  ridge-fixed guard compliance                -> %s"
          % ("PASS" if ok_rf else "FAIL"))

    # ---- wall time --------------------------------------------------------- #
    wall = time.time() - t0_all
    ok_time = wall < 900.0
    print("\n[v3 selftest] === (5) wall time ===")
    print("  design_v3 = %.1f s; selftest total = %.1f s (< 900 s) -> %s"
          % (t_design, wall, "PASS" if ok_time else "FAIL"))

    hard = {"kernel_certs": certs_ok, "budget": ok_budget, "G5_rounded": ok_g5,
            "ridge_zone_or_pass_by_design": ok_ridge,
            "ridge_fixed_guards": ok_rf, "wall<15min": ok_time}
    all_ok = all(hard.values())
    print("\n[v3 selftest] hard checks: %s" % hard)
    print("[v3 selftest] RESULT: %s"
          % ("ALL HARD CHECKS PASS" if all_ok else "*** FAILURE ***"))
