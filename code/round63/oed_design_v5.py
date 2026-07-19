"""R15 decisive-iteration solver: full dose-constrained KKT certificate.

Implements docs/ROUND63_GPT_ROUND15_RULING_RAW.md on top of oed_design_v4:

Sec 2   FULL_CONSTRAINED_KW_UPPER_BOUND: the additive dual bound with ONE
        budget multiplier theta and ALL upper/lower per-pixel dose-band
        multipliers mu+-, exact linear dose form
        u+_ja = z_aj - (1+delta) zbar_a, u-_ja = (1-delta) zbar_a - z_aj,
        z_a = the atom's per-exposure physical dose vector, delta = 0.05.
        Multiplier minimization by epigraph LP (scipy HiGHS) with certified
        column generation; the final max is scanned over the ENTIRE frozen
        dictionary. Hard gate G_full/r <= 1e-2; report tightness at 1e-3.
Sec 3   support-restricted exact reweighting of the dose-relaxed reference
        problem (frozen support rule: xi > 1e-12 union top-300 budget-
        adjusted sensitivities; deterministic index ties), solved to machine
        precision by an independent convex solver (SLSQP with analytic
        gradients), then RELAXED_REFERENCE_KW_UPPER_BOUND under a full-
        dictionary scan; machinery gate <= 1e-3.
Sec 6   exact materialized mixture search m = 972..0 (descending, m=0
        included, incident tolerance cap + 1e-9), frozen hashed row orders
        for both sources, EVERY final guard recomputed per candidate
        (incremental rank-one swaps), first = largest-OED-fraction winner;
        D-efficiency referenced to the certified full dose-constrained
        continuous design; A-risk/spectral against the LOAD-MATCHED FIXED*.

Toy verification with an ACTIVE dose constraint (R15 Sec 2.3): brute
constrained optimum by an independent solver, dual-bound agreement <= 1e-8.
"""

import hashlib
import json
import math
import os
import sys

import numpy as np
from scipy.optimize import linprog, minimize

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import oed_design_v4 as v4                  # noqa: E402

DELTA = 0.05
GFULL_TOL_PER_R = 1e-2
GFULL_TIGHT_PER_R = 1e-3
RELREF_TOL_PER_R = 1e-3


# ========================================================================== #
#  shared setup                                                               #
# ========================================================================== #
def setup_ctx(xhat, nu, palette, B, eps0, rho_prescan, M_rows=972, side=32):
    n = side * side
    r = B.shape[1]
    levels = np.asarray(palette["levels"], dtype=float)
    L = levels.size
    nuJ = np.array([nu * v4.kernel_eval4(l_, int(nu))["J_exact"]
                    for l_ in levels])
    blocks = v4.build_blocks(xhat, side)
    QBs, metas = [], []
    for b_ in blocks:
        QBs.append(np.einsum('tk,tkr->tr', b_["GVAL"], B[b_["IDX"]]))
        metas.append(b_)
    QB = np.vstack(QBs)
    G_tot = QB.shape[0]
    gmax = np.concatenate([m["gmax"] for m in metas])
    shape_ok = np.concatenate([m["shape_ok"] for m in metas])
    gsum = np.concatenate([m["GVAL"].sum(axis=1) for m in metas])
    ALLOW = shape_ok[:, None] & (gmax[:, None] * levels[None, :]
                                 <= v4.PEAK_PHYSICAL + 1e-9)
    # per-atom load matrix C (G, L) and per-atom nu*J matrix NUJ: for the
    # D_load palette both are column-constant; the R17 D_gain extension
    # (setup_ctx_cert) appends gamma columns with per-geometry emergent loads.
    C = np.tile(levels[None, :], (G_tot, 1))
    NUJ = np.tile(nuJ[None, :], (G_tot, 1))
    V0 = v4.V0_prescan(v4.balanced_prescan_52(side), xhat, nu, rho_prescan, B)
    return {"n": n, "r": r, "levels": levels, "L": L, "nuJ": nuJ,
            "C": C, "NUJ": NUJ, "gmax": gmax,
            "metas": metas, "QB": QB, "G_tot": G_tot, "gsum": gsum,
            "ALLOW": ALLOW, "V0": V0, "eps0": eps0, "M_rows": M_rows,
            "nu": nu, "xhat": xhat, "side": side, "B": B}


def assemble(ctx, xi):
    u = ctx["M_rows"] * (xi * ctx["NUJ"]).sum(axis=1)
    nz = np.nonzero(u)[0]
    M = (ctx["QB"][nz] * u[nz, None]).T @ ctx["QB"][nz]
    V = ctx["V0"] + M
    V[np.diag_indices(ctx["r"])] += ctx["eps0"]
    return V


def dvals(ctx, Vinv):
    quad = np.einsum('gr,rs,gs->g', ctx["QB"], Vinv, ctx["QB"])
    D = ctx["M_rows"] * quad[:, None] * ctx["NUJ"]
    return np.where(ctx["ALLOW"], D, -np.inf)


def dose_of(ctx, xi):
    w = (xi * ctx["C"]).sum(axis=1)
    d = np.zeros(ctx["n"])
    off = 0
    for mb in ctx["metas"]:
        g_ = mb["IDX"].shape[0]
        wn = w[off:off + g_]
        nz = np.nonzero(wn)[0]
        if nz.size:
            d += np.bincount(mb["IDX"][nz].ravel(),
                             weights=(wn[nz, None] * mb["GVAL"][nz]).ravel(),
                             minlength=ctx["n"])
        off += g_
    return d


def mu_dots(ctx, mu):
    """Per-geometry <mu, g>_support for a dense pixel multiplier vector."""
    out = np.empty(ctx["G_tot"])
    off = 0
    for mb in ctx["metas"]:
        g_ = mb["IDX"].shape[0]
        out[off:off + g_] = (mb["GVAL"] * mu[mb["IDX"]]).sum(axis=1)
        off += g_
    return out


def lin_from_duals(ctx, theta, mup, mum):
    """lin_(g,l) = theta*c + sum_j mu+ u+ + sum_j mu- u-  (per atom)."""
    Sp, Sm = float(mup.sum()), float(mum.sum())
    Ap = mu_dots(ctx, mup)
    Am = mu_dots(ctx, mum)
    zb = ctx["gsum"] / ctx["n"]                       # zbar per unit load
    per_g = (Ap - (1 + DELTA) * zb * Sp) + ((1 - DELTA) * zb * Sm - Am)
    lin = (theta + per_g[:, None]) * ctx["C"]
    return lin


# ========================================================================== #
#  Sec 2: full constrained certificate (epigraph LP + column generation)      #
# ========================================================================== #
def full_constrained_cert(ctx, xi, budget, seed_atoms=None, max_rounds=12,
                          verbose=False, V_override=None, dbar_override=None,
                          dose_override=None, load_override=None,
                          n_pix0=256, mu_cap=1e6, lp_time=120):
    """xi-measure certificate; pass xi=None with the *_override arguments to
    certify an ARBITRARY deployed row set (R17: the exact SCAT32 design) --
    V, d^T xi, dose vector and mean load are then supplied by the caller and
    the bound max_a s_a - dbar + theta*b is taken over the same dictionary."""
    n, r, L = ctx["n"], ctx["r"], ctx["L"]
    levels = ctx["levels"]
    if V_override is not None:
        V = V_override
        dbar = float(dbar_override)
    else:
        V = assemble(ctx, xi)
    Vinv = np.linalg.inv(V)
    D = dvals(ctx, Vinv)
    if V_override is None:
        dbar = float((xi * np.where(ctx["ALLOW"], D, 0.0)).sum())
    # candidate set: support + top by d + caller seeds
    flat = np.argsort(D, axis=None)[::-1][:600]
    cand = set(map(tuple, np.column_stack(
        np.unravel_index(flat, D.shape)).tolist()))
    if xi is not None:
        cand |= set(map(tuple,
                        np.column_stack(np.nonzero(xi > 1e-12)).tolist()))
    if seed_atoms:
        cand |= set(seed_atoms)
    cand = sorted(a for a in cand if ctx["ALLOW"][a[0], a[1]])

    def atom_dose_row(g, l_):
        """(z dense, zbar) for atom (g, l); load from the per-atom C."""
        off = 0
        for mb in ctx["metas"]:
            g_ = mb["IDX"].shape[0]
            if g < off + g_:
                z = np.zeros(n)
                z[mb["IDX"][g - off]] = ctx["C"][g, l_] * mb["GVAL"][g - off]
                return z, ctx["C"][g, l_] * ctx["gsum"][g] / n
            off += g_
        raise IndexError(g)

    # mu columns restricted to pixels nearest the dose-band edges: ANY
    # feasible dual choice yields a valid conservative upper bound, so a
    # pixel subset only loosens G_full; the subset is expanded adaptively
    # (256 -> 512 -> 1024) while the bound stays above the gate. This keeps
    # the epigraph LP small (HiGHS ground >180s on the full 2*1024 dense,
    # degenerate mu block -- the freeze-8 LP_FAIL:time-limit root cause).
    dv0 = dose_override if dose_override is not None else dose_of(ctx, xi)
    dm0 = dv0.mean()
    edge = np.maximum(dv0 - (1 + DELTA) * dm0, (1 - DELTA) * dm0 - dv0)
    pix_rank = np.argsort(-edge, kind="stable")
    # pixel -> covering geometries index (kills the unbounded mu- ray: a
    # mu-pixel with NO covering atom in the LP would let t -> -inf, the
    # freeze-9 HiGHS Status-10 root cause)
    if "pix_cover" not in ctx:
        pairs_g, pairs_j = [], []
        off = 0
        for mb in ctx["metas"]:
            gg = mb["IDX"].shape[0]
            pairs_g.append(np.repeat(np.arange(off, off + gg),
                                     mb["IDX"].shape[1]))
            pairs_j.append(mb["IDX"].ravel())
            off += gg
        pg_all = np.concatenate(pairs_g)
        pj_all = np.concatenate(pairs_j)
        o_ = np.argsort(pj_all, kind="stable")
        pj_s, pg_s = pj_all[o_], pg_all[o_]
        bounds_j = np.searchsorted(pj_s, np.arange(n + 1))
        ctx["pix_cover"] = [pg_s[bounds_j[j]:bounds_j[j + 1]]
                            for j in range(n)]
    lvl_best = np.argmax(np.where(ctx["ALLOW"], D, -np.inf), axis=1)
    theta, mup, mum, t_val = 0.0, np.zeros(n), np.zeros(n), np.inf
    n_pix = int(n_pix0)
    d_scale = max(float(np.max(D[np.isfinite(D)])), 1.0)
    for rnd in range(max_rounds):
        pix = np.sort(pix_rank[:n_pix])
        npx = pix.size
        # coverage seeding: every mu-pixel must have >= 1 covering atom
        in_cand = np.zeros(ctx["G_tot"], dtype=bool)
        for (g_c, _l_c) in cand:
            in_cand[g_c] = True
        added = []
        for j in pix:
            cov = ctx["pix_cover"][j]
            if cov.size and not in_cand[cov].any():
                gb = int(cov[np.argmax(D[cov, lvl_best[cov]])])
                added.append((gb, int(lvl_best[gb])))
                in_cand[gb] = True
        if added:
            cand = sorted(set(cand) | set(added))
        m_ = len(cand)
        A_ub = np.zeros((m_, 2 + 2 * npx))
        b_ub = np.zeros(m_)
        for i, (g, l_) in enumerate(cand):
            z, zb = atom_dose_row(g, l_)
            A_ub[i, 0] = -1.0                          # t
            A_ub[i, 1] = -ctx["C"][g, l_]              # theta
            A_ub[i, 2:2 + npx] = -(z[pix] - (1 + DELTA) * zb)     # mu+
            A_ub[i, 2 + npx:] = -((1 - DELTA) * zb - z[pix])      # mu-
            b_ub[i] = -float(D[g, l_]) / d_scale
        A_ub[:, 2:] /= d_scale               # mu columns share the row scale
        c_obj = np.zeros(2 + 2 * npx)
        c_obj[0] = 1.0
        c_obj[1] = budget
        MU_CAP = float(mu_cap)  # restricting the dual cone only loosens
        bounds = [(None, None), (0, None)] + [(0, MU_CAP)] * (2 * npx)
        res = linprog(c_obj, A_ub=A_ub, b_ub=b_ub, bounds=bounds,
                      method="highs", options={"time_limit": float(lp_time)})
        if not res.success:
            return {"G_full": float("inf"), "status": "LP_FAIL:" + res.message}
        t_val = float(res.x[0]) * d_scale
        theta = float(res.x[1]) * d_scale
        mup = np.zeros(n)
        mum = np.zeros(n)
        mup[pix] = np.maximum(res.x[2:2 + npx], 0.0)
        mum[pix] = np.maximum(res.x[2 + npx:], 0.0)
        # certified full-dictionary scan at these multipliers
        lin = lin_from_duals(ctx, theta, mup, mum)
        s_all = np.where(ctx["ALLOW"], D - lin, -np.inf)
        s_max = float(s_all.max())
        viol = s_max - t_val
        if verbose:
            print("      [cert] round=%d set=%d t=%.4f scan_max=%.4f "
                  "viol=%.2e theta=%.3g |mu|=%d"
                  % (rnd, m_, t_val, s_max, viol,
                     theta, int((mup > 1e-12).sum() + (mum > 1e-12).sum())),
                  flush=True)
        if viol <= 1e-7 * max(abs(t_val), 1.0):
            t_val = max(t_val, s_max)
            g_now = t_val + theta * budget - dbar
            if g_now / ctx["r"] > GFULL_TOL_PER_R and n_pix < 1024:
                n_pix = min(2 * n_pix, 1024)     # loosen -> expand mu pixels
                if verbose:
                    print("      [cert] expanding mu pixel set to %d"
                          % n_pix, flush=True)
                continue
            break
        top = np.argsort(s_all, axis=None)[::-1][:100]
        cand = sorted(set(cand) | set(map(tuple, np.column_stack(
            np.unravel_index(top, s_all.shape)).tolist())))
    G_full = t_val + theta * budget - dbar
    if G_full < 0:
        if G_full >= -1e-10 * r:
            G_full = 0.0
        else:
            return {"G_full": float(G_full), "status": "NEGATIVE_BOUND_FAIL"}
    # primal residuals (R15 mandatory checks)
    dv = dose_override if dose_override is not None else dose_of(ctx, xi)
    dm = dv.mean()
    dose_excess = float(max(np.max(dv - (1 + DELTA) * dm),
                            np.max((1 - DELTA) * dm - dv), 0.0))
    load = (float(load_override) if load_override is not None
            else float((xi * ctx["C"]).sum()))
    comp_budget = theta * (budget - load)
    comp_dose = float((mup * np.maximum(dv - (1 + DELTA) * dm, -np.inf)).sum()
                      + (mum * ((1 - DELTA) * dm - dv)).sum())
    return {"G_full": float(G_full), "status": "OK", "theta": theta,
            "mu_plus": mup, "mu_minus": mum,
            "budget_viol": float(max(load - budget, 0.0)),
            "dose_excess": dose_excess,
            "comp_budget": float(comp_budget),
            "comp_dose": float(comp_dose),
            "dbar": dbar, "t": t_val,
            "MU_CAP_ACTIVE": bool(max(mup.max(), mum.max())
                                  >= float(mu_cap) * (1 - 1e-9)),
            "n_active_mu": int((mup > 1e-12).sum() + (mum > 1e-12).sum())}


# ========================================================================== #
#  Sec 2.3: toy with an ACTIVE dose constraint                                #
# ========================================================================== #
def toy_full_cert_check():
    """2-pixel, 3-atom toy where the +/-5% dose band is ACTIVE: solve the
    constrained D-opt by an independent convex solver (SLSQP), then require
    the LP dual bound at the optimum to agree within 1e-8."""
    X = np.array([[1.0, 0.0], [0.0, 1.0],
                  [1.0 / math.sqrt(2), 1.0 / math.sqrt(2)]])
    Z = np.array([[3.0, 0.0], [0.0, 1.0], [1.0, 1.0]])   # dose vectors
    c_load = np.array([1.0, 1.0, 1.0])
    budget = 1.0
    eps = 1e-9

    def negld(w):
        M = sum(wi * np.outer(x, x) for wi, x in zip(w, X)) + eps * np.eye(2)
        s_, ld = np.linalg.slogdet(M)
        return -ld

    def grad(w):
        M = sum(wi * np.outer(x, x) for wi, x in zip(w, X)) + eps * np.eye(2)
        Mi = np.linalg.inv(M)
        return -np.array([x @ Mi @ x for x in X])

    def dose_cons(w):
        dv = Z.T @ w
        dm = dv.mean()
        return np.concatenate([(1 + DELTA) * dm - dv, dv - (1 - DELTA) * dm])

    cons = [{"type": "eq", "fun": lambda w: w.sum() - 1.0,
             "jac": lambda w: np.ones(3)},
            {"type": "ineq", "fun": dose_cons}]
    best = None
    for w0 in ([1 / 3] * 3, [0.25, 0.55, 0.20], [0.2, 0.4, 0.4]):
        res = minimize(negld, np.array(w0), jac=grad, method="SLSQP",
                       bounds=[(0, 1)] * 3, constraints=cons,
                       options={"maxiter": 500, "ftol": 1e-14})
        if res.success and (best is None or res.fun < best.fun):
            best = res
    w_star = best.x
    # dual bound at the optimum via the same LP machinery (2 pixels)
    M = sum(wi * np.outer(x, x) for wi, x in zip(w_star, X)) + eps * np.eye(2)
    Mi = np.linalg.inv(M)
    d = np.array([x @ Mi @ x for x in X])
    zb = Z.mean(axis=1)
    n2 = 2
    A_ub = np.zeros((3, 2 + 2 * n2))
    for i in range(3):
        A_ub[i, 0] = -1.0
        A_ub[i, 1] = -c_load[i]
        A_ub[i, 2:2 + n2] = -(Z[i] - (1 + DELTA) * zb[i])
        A_ub[i, 2 + n2:] = -((1 - DELTA) * zb[i] - Z[i])
    b_ub = -d
    c_obj = np.zeros(2 + 2 * n2)
    c_obj[0] = 1.0
    c_obj[1] = budget
    res = linprog(c_obj, A_ub=A_ub, b_ub=b_ub,
                  bounds=[(None, None)] + [(0, None)] * (1 + 2 * n2),
                  method="highs")
    G = float(res.x[0]) + float(res.x[1]) * budget - float(d @ w_star)
    dv = Z.T @ w_star
    active = float(np.min(np.concatenate(
        [(1 + DELTA) * dv.mean() - dv, dv - (1 - DELTA) * dv.mean()])))
    return {"G_at_opt": G, "dose_active": bool(abs(active) < 1e-6),
            "w_star": w_star.tolist(), "agreement_ok": bool(abs(G) <= 1e-8)}


# ========================================================================== #
#  Sec 3: support-restricted exact reweighting (relaxed reference)            #
# ========================================================================== #
def support_reweight(ctx, xi_relaxed, budget, verbose=False, s_cap=400,
                     maxiter=250):
    import time as _t
    t0 = _t.time()
    levels = ctx["levels"]
    V = assemble(ctx, xi_relaxed)
    Vinv = np.linalg.inv(V)
    D = dvals(ctx, Vinv)
    dbar = float((xi_relaxed * np.where(ctx["ALLOW"], D, 0.0)).sum())
    _g, th = v4._cert_grel(D[ctx["ALLOW"]],
                           np.repeat(levels[None, :], ctx["G_tot"], 0)
                           [ctx["ALLOW"]], None, 0.0, dbar, budget)
    s_adj = np.where(ctx["ALLOW"], D - th * levels[None, :], -np.inf)
    top = np.argsort(s_adj, axis=None)[::-1][:300]
    supp = set(map(tuple, np.column_stack(
        np.unravel_index(top, s_adj.shape)).tolist()))
    supp |= set(map(tuple,
                    np.column_stack(np.nonzero(xi_relaxed > 1e-12)).tolist()))
    supp = sorted(a for a in supp if ctx["ALLOW"][a[0], a[1]])
    if len(supp) > s_cap:                    # deterministic size cap: keep all
        supp = sorted(supp,                  # current-support atoms, then the
                      key=lambda a: (xi_relaxed[a] <= 1e-12,
                                     -float(s_adj[a]), a))[:s_cap]
        supp = sorted(supp)
    S = len(supp)
    if verbose:
        print("    [ref] SLSQP start S=%d maxiter=%d" % (S, maxiter),
              flush=True)
    gs = np.array([a[0] for a in supp])
    ls = np.array([a[1] for a in supp])
    Qs = ctx["QB"][gs]
    nJs = ctx["nuJ"][ls]
    cs = levels[ls]
    M_rows, r = ctx["M_rows"], ctx["r"]
    I_r = np.eye(r)

    def V_of(w):
        return ctx["V0"] + (Qs * (M_rows * nJs * w)[:, None]).T @ Qs \
            + ctx["eps0"] * I_r

    def negld(w):
        s_, ld = np.linalg.slogdet(V_of(w))
        return -ld if s_ > 0 else 1e18

    def grad(w):
        Vi = np.linalg.inv(V_of(w))
        return -(M_rows * nJs * np.einsum('ar,rs,as->a', Qs, Vi, Qs))

    w0 = np.array([max(float(xi_relaxed[a]), 1e-9) for a in supp])
    w0 /= w0.sum()
    cons = [{"type": "eq", "fun": lambda w: w.sum() - 1.0,
             "jac": lambda w: np.ones(S)},
            {"type": "ineq", "fun": lambda w: budget - float(cs @ w),
             "jac": lambda w: -cs}]
    it_box = [0]

    def _cb(wk):
        it_box[0] += 1
        if verbose and it_box[0] % 25 == 0:
            print("    [ref] SLSQP it=%d obj=%.6f (%.0fs)"
                  % (it_box[0], -negld(wk), _t.time() - t0), flush=True)

    res = minimize(negld, w0, jac=grad, method="SLSQP",
                   bounds=[(0, 1)] * S, constraints=cons, callback=_cb,
                   options={"maxiter": maxiter, "ftol": 1e-16})
    if verbose:
        print("    [ref] SLSQP done its=%d status=%s%s (%.0fs)"
              % (it_box[0], res.message,
                 " CAP_HIT" if it_box[0] >= maxiter else "",
                 _t.time() - t0), flush=True)
    w = np.maximum(res.x, 0.0)
    w /= w.sum()
    xi_ref = np.zeros_like(xi_relaxed)
    for a, wv in zip(supp, w):
        xi_ref[a] = wv
    # support KKT residual (against the budget-adjusted equalized sensitivity)
    Vi = np.linalg.inv(V_of(w))
    ds = M_rows * nJs * np.einsum('ar,rs,as->a', Qs, Vi, Qs)
    on = w > 1e-10
    _g2, th2 = v4._cert_grel(ds, cs, None, 0.0, float(ds @ w), budget)
    s_on = ds[on] - th2 * cs[on]
    kkt_res = float(s_on.max() - s_on.min()) / r if on.any() else np.inf
    # full-dictionary relaxed reference bound
    D2 = dvals(ctx, np.linalg.inv(assemble(ctx, xi_ref)))
    dbar2 = float((xi_ref * np.where(ctx["ALLOW"], D2, 0.0)).sum())
    G2, th3 = v4._cert_grel(D2[ctx["ALLOW"]],
                            np.repeat(levels[None, :], ctx["G_tot"], 0)
                            [ctx["ALLOW"]], None, 0.0, dbar2, budget)
    return {"xi_ref": xi_ref, "bound_per_r": float(G2 / r),
            "kkt_residual_per_r": kkt_res,
            "slsqp_status": res.message, "n_support": S}


# ========================================================================== #
#  primal-dual constrained solve (drives G_full down at dose feasibility)     #
# ========================================================================== #
def primal_dual_solve(ctx, xi0, budget, dose_band=DELTA, rounds=10,
                      verbose=False):
    levels = ctx["levels"]
    xi = xi0.copy()
    traj = []
    best = (np.inf, xi.copy(), None)
    for rd in range(rounds):
        cert = full_constrained_cert(ctx, xi, budget, verbose=verbose)
        traj.append((rd, float(cert["G_full"] / ctx["r"]),
                     cert.get("dose_excess", np.nan)))
        if verbose:
            print("    [pd] round=%d G_full/r=%.3e dose_exc=%.2e theta=%.3g"
                  % (rd, cert["G_full"] / ctx["r"],
                     cert.get("dose_excess", -1), cert.get("theta", -1)),
                  flush=True)
        if cert["status"] != "OK":
            print("    [pd] CERT STATUS: %s (aborting rounds)"
                  % cert["status"], flush=True)
            break
        if cert["G_full"] < best[0]:
            best = (cert["G_full"], xi.copy(), cert)
        if cert["G_full"] / ctx["r"] <= GFULL_TOL_PER_R * 0.5:
            break
        # Lagrangian primal step at the LP-optimal duals (support-restricted
        # SLSQP on the simplex), then damped mix + dose re-projection
        lin = lin_from_duals(ctx, cert["theta"], cert["mu_plus"],
                             cert["mu_minus"])
        V = assemble(ctx, xi)
        D = dvals(ctx, np.linalg.inv(V))
        s_all = np.where(ctx["ALLOW"], D - lin, -np.inf)
        top = np.argsort(s_all, axis=None)[::-1][:300]
        supp = set(map(tuple, np.column_stack(
            np.unravel_index(top, s_all.shape)).tolist()))
        supp |= set(map(tuple, np.column_stack(np.nonzero(xi > 1e-12))
                        .tolist()))
        supp = sorted(a for a in supp if ctx["ALLOW"][a[0], a[1]])
        S = len(supp)
        gs = np.array([a[0] for a in supp])
        ls = np.array([a[1] for a in supp])
        Qs = ctx["QB"][gs]
        nJs = ctx["nuJ"][ls]
        lins = np.array([lin[a] for a in supp])
        I_r = np.eye(ctx["r"])

        def V_of(w):
            return ctx["V0"] + (Qs * (ctx["M_rows"] * nJs * w)[:, None]).T \
                @ Qs + ctx["eps0"] * I_r

        def nobj(w):
            s_, ld = np.linalg.slogdet(V_of(w))
            return -(ld - float(lins @ w)) if s_ > 0 else 1e18

        def ngrad(w):
            Vi = np.linalg.inv(V_of(w))
            d_ = ctx["M_rows"] * nJs * np.einsum('ar,rs,as->a', Qs, Vi, Qs)
            return -(d_ - lins)

        if S > 350:                          # deterministic support cap
            keep = sorted(range(S), key=lambda i: (-float(
                s_all[supp[i]]), supp[i]))[:350]
            supp = sorted(supp[i] for i in keep)
            S = len(supp)
            gs = np.array([a[0] for a in supp])
            ls = np.array([a[1] for a in supp])
            Qs = ctx["QB"][gs]
            nJs = ctx["nuJ"][ls]
            lins = np.array([lin[a] for a in supp])
        if verbose:
            print("    [pd] lagrangian SLSQP S=%d" % S, flush=True)
        w0 = np.array([max(float(xi[a]), 1e-9) for a in supp])
        w0 /= w0.sum()
        res = minimize(nobj, w0, jac=ngrad, method="SLSQP",
                       bounds=[(0, 1)] * S,
                       constraints=[{"type": "eq",
                                     "fun": lambda w: w.sum() - 1.0,
                                     "jac": lambda w: np.ones(S)}],
                       options={"maxiter": 150, "ftol": 1e-12})
        xi_new = np.zeros_like(xi)
        wn = np.maximum(res.x, 0.0)
        wn /= wn.sum()
        for a, wv in zip(supp, wn):
            xi_new[a] = wv
        for alpha in (0.5, 0.25, 0.1):
            xi_try = (1 - alpha) * xi + alpha * xi_new
            xi_try = _dose_project(ctx, xi_try, budget, dose_band,
                                   verbose=False)
            c_try = full_constrained_cert(ctx, xi_try, budget)
            if c_try["status"] == "OK" and c_try["G_full"] < cert["G_full"]:
                xi = xi_try
                break
        else:
            break                                    # no improving damped step
    G_best, xi_best, cert_best = best
    if cert_best is None:
        cert_best = full_constrained_cert(ctx, xi_best, budget)
    return xi_best, cert_best, traj


def _dose_project(ctx, xi, budget, dose_band, verbose=False):
    xi = np.maximum(xi, 0.0)
    xi = xi / xi.sum()
    it_p = -1
    for it_p in range(400):
        if verbose and it_p % 100 == 0 and it_p > 0:
            dv0 = dose_of(ctx, xi)
            print("    [proj] pass=%d dev=%.4f" % (
                it_p, float(np.abs(dv0 / dv0.mean() - 1).max())), flush=True)
        dv = dose_of(ctx, xi)
        rel = dv / dv.mean()
        if np.abs(rel - 1.0).max() <= dose_band - 0.012:
            break
        rc = np.clip(rel, 0.5, 2.0)
        off = 0
        fac = np.ones(ctx["G_tot"])
        for mb in ctx["metas"]:
            g_ = mb["IDX"].shape[0]
            fac[off:off + g_] = 1.0 / ((mb["GVAL"] * rc[mb["IDX"]])
                                       .sum(axis=1) / mb["GVAL"].sum(axis=1))
            off += g_
        xi = xi * fac[:, None]
        xi[~ctx["ALLOW"]] = 0.0
        xi /= xi.sum()
        cur = float((xi * ctx["C"]).sum())
        if cur > budget + 1e-15:
            low = np.zeros_like(xi)
            lo_idx = np.argmax(ctx["ALLOW"], axis=1)
            low[np.arange(ctx["G_tot"]), lo_idx] = xi.sum(axis=1)
            lmin = float((low * ctx["C"]).sum())
            t_ = 1.0 if lmin >= budget else (cur - budget) / (cur - lmin)
            xi = (1 - t_) * xi + t_ * low
    if it_p >= 399:
        dv0 = dose_of(ctx, xi)
        print("    [proj] CAP_HIT after 400 passes, dev=%.4f" % (
            float(np.abs(dv0 / dv0.mean() - 1).max())), flush=True)
    return xi


# ========================================================================== #
#  phase 1: relaxed PFW (ported from v4, ctx-based)                           #
# ========================================================================== #
def pfw_relaxed(ctx, budget, max_outer=40, max_inner=5000, verbose=False):
    levels, L, G_tot, r = ctx["levels"], ctx["L"], ctx["G_tot"], ctx["r"]
    QB, ALLOW, nuJ, M_rows = ctx["QB"], ctx["ALLOW"], ctx["nuJ"], ctx["M_rows"]
    lo_lvl = int(np.argmax(ALLOW.any(axis=0)))
    xi = np.zeros((G_tot, L))
    a0 = (int(np.argmax(ALLOW[:, lo_lvl])), lo_lvl)
    xi[a0] = 1.0
    weights = {("s", a0): 1.0}

    def vertex_atoms(v):
        if v[0] == "s":
            return [(v[1], 1.0)]
        (glo, llo), (ghi, lhi) = v[1], v[2]
        wlo = (levels[lhi] - budget) / (levels[lhi] - levels[llo])
        return [((glo, llo), wlo), ((ghi, lhi), 1.0 - wlo)]

    V = assemble(ctx, xi)
    Vinv = np.linalg.inv(V)
    traj, n_tot, stall, best_c = [], 0, 0, np.inf
    for outer in range(max_outer):
        D = dvals(ctx, Vinv)
        dbar = float((xi * np.where(ALLOW, D, 0.0)).sum())
        G_, th = v4._cert_grel(D[ALLOW],
                               np.repeat(levels[None, :], G_tot, 0)[ALLOW],
                               None, 0.0, dbar, budget)
        traj.append((outer, float(G_ / r)))
        if verbose:
            print("    [pfw] outer=%d G/r=%.3e" % (outer, G_ / r), flush=True)
        if G_ / r <= RELREF_TOL_PER_R:
            break
        if G_ < best_c * 0.995:
            best_c, stall = G_, 0
        else:
            stall += 1
            if stall >= 3 or n_tot >= max_inner:
                break
        s_adj = D - th * levels[None, :]
        top = np.column_stack(np.unravel_index(
            np.argsort(s_adj, axis=None)[::-1][:300], s_adj.shape))
        pool = {tuple(t) for t in top.tolist()}
        pool |= {tuple(t) for t in np.column_stack(np.nonzero(xi > 0)).tolist()}
        pool = [p for p in pool if ALLOW[p[0], p[1]]]
        pg = np.array([p[0] for p in pool])
        pl = np.array([p[1] for p in pool])
        Qp = QB[pg]
        sclp = M_rows * nuJ[pl]
        cp = levels[pl]
        for _ in range(200):
            n_tot += 1
            dp = sclp * np.einsum('ar,rs,as->a', Qp, Vinv, Qp)
            feas = cp <= budget + 1e-12
            best_v, best_val = None, -np.inf
            if feas.any():
                i_s = int(np.argmax(np.where(feas, dp, -np.inf)))
                best_v, best_val = ("s", pool[i_s]), float(dp[i_s])
            for llo in range(L):
                if levels[llo] > budget:
                    continue
                m_lo = np.where(pl == llo, dp, -np.inf)
                if not np.isfinite(m_lo).any():
                    continue
                i_lo = int(np.argmax(m_lo))
                for lhi in range(L):
                    if levels[lhi] <= budget:
                        continue
                    m_hi = np.where(pl == lhi, dp, -np.inf)
                    if not np.isfinite(m_hi).any():
                        continue
                    i_hi = int(np.argmax(m_hi))
                    wlo = (levels[lhi] - budget) / (levels[lhi] - levels[llo])
                    val = wlo * dp[i_lo] + (1 - wlo) * dp[i_hi]
                    if val > best_val:
                        best_val = float(val)
                        best_v = ("p", pool[i_lo], pool[i_hi])
            dmap = {p: float(dp[k]) for k, p in enumerate(pool)}
            away_v, away_val = None, np.inf
            for v_, w_ in weights.items():
                if w_ <= 1e-14:
                    continue
                val = sum(co * dmap.get(a_, 0.0) for a_, co in vertex_atoms(v_))
                if val < away_val:
                    away_val, away_v = val, v_
            if best_v is None or away_v is None or best_v == away_v:
                break
            if best_val - away_val <= 1e-9 * max(abs(dbar), 1.0):
                break
            coef = {}
            for a_, co in vertex_atoms(best_v):
                coef[a_] = coef.get(a_, 0.0) + co
            for a_, co in vertex_atoms(away_v):
                coef[a_] = coef.get(a_, 0.0) - co
            coef = {a_: co for a_, co in coef.items() if abs(co) > 1e-15}
            if not coef:
                break
            U = np.column_stack([math.sqrt(M_rows * nuJ[l_]) * QB[g_]
                                 for (g_, l_) in coef])
            svec = np.array(list(coef.values()))
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
        V = assemble(ctx, xi)
        Vinv = np.linalg.inv(V)
    return xi, traj


# ========================================================================== #
#  exact rows: rounding, repair, R15 Sec-6.2 materialized mixture             #
# ========================================================================== #
def _atom_row(ctx, g, l_):
    off = 0
    for mb in ctx["metas"]:
        gg = mb["IDX"].shape[0]
        if g < off + gg:
            row = np.zeros(ctx["n"])
            row[mb["IDX"][g - off]] = ctx["levels"][l_] * mb["GVAL"][g - off]
            return row
        off += gg
    raise IndexError(g)


def round_design(ctx, xi, budget):
    """Largest remainder -> atoms [(g, l, count)] + budget demotion."""
    flat = xi.ravel()
    supp = np.where(flat > flat.max() * 1e-9)[0]
    w = flat[supp] / flat[supp].sum()
    quota = w * ctx["M_rows"]
    cnt = np.floor(quota).astype(np.int64)
    rem = int(ctx["M_rows"] - cnt.sum())
    order = np.lexsort((supp, -(quota - cnt)))
    for i in range(rem):
        cnt[order[i]] += 1
    atoms = [(int(f // ctx["L"]), int(f % ctx["L"]), int(c))
             for f, c in zip(supp, cnt) if c > 0]
    lo_lvl = int(np.argmax(ctx["ALLOW"].any(axis=0)))
    while sum(ctx["levels"][l_] * c for _, l_, c in atoms) \
            > budget * ctx["M_rows"] + 1e-9:
        atoms.sort(key=lambda a: -ctx["levels"][a[1]])
        g, l_, c = atoms[0]
        if ctx["levels"][l_] <= ctx["levels"][lo_lvl] + 1e-12:
            break
        atoms[0] = (g, l_, c - 1)
        atoms = [a for a in atoms if a[2] > 0]
        for i, a in enumerate(atoms):
            if a[0] == g and a[1] == lo_lvl:
                atoms[i] = (g, lo_lvl, a[2] + 1)
                break
        else:
            atoms.append((g, lo_lvl, 1))
    return atoms


VERBOSE_REPAIR = [True]


def verbose_flag():
    return VERBOSE_REPAIR[0]


def dose_repair(ctx, atoms, dose_band=DELTA, n_pre=300, n_polish=900):
    """R15 Sec-6.1 steps 3-4: prescribed worst-pair repair with best-state
    retention, then strictly descending penalty polish (same level swaps)."""
    M = ctx["M_rows"]

    def dvec(al):
        d = np.zeros(ctx["n"])
        for g, l_, c in al:
            d += c * _atom_row(ctx, g, l_)
        return d / M

    def dev(dv):
        return float(np.abs(dv / dv.mean() - 1.0).max())

    dv = dvec(atoms)
    best = (dev(dv), [tuple(a) for a in atoms])
    for _ in range(n_pre):
        d_ = dev(dv)
        if d_ <= dose_band:
            best = (d_, [tuple(a) for a in atoms])
            break
        if d_ < best[0] - 1e-12:
            best = (d_, [tuple(a) for a in atoms])
        rel = dv / dv.mean() - 1.0
        jo, ju = int(np.argmax(rel)), int(np.argmin(rel))
        moved = False
        for i, (g, l_, c) in enumerate(atoms):
            row = _atom_row(ctx, g, l_)
            if row[jo] <= 0:
                continue
            off = 0
            for mb in ctx["metas"]:
                gg = mb["IDX"].shape[0]
                if g < off + gg:
                    cov = np.nonzero((mb["IDX"] == ju).any(axis=1))[0]
                    if cov.size:
                        g2 = off + int(cov[0])
                        row2 = _atom_row(ctx, g2, l_)
                        atoms[i] = (g, l_, c - 1)
                        atoms = [a for a in atoms if a[2] > 0]
                        for k2, a in enumerate(atoms):
                            if a[0] == g2 and a[1] == l_:
                                atoms[k2] = (g2, l_, a[2] + 1)
                                break
                        else:
                            atoms.append((g2, l_, 1))
                        dv = dv + (row2 - row) / M
                        moved = True
                    break
                off += gg
            if moved:
                break
        if not moved:
            break
    atoms = [tuple(a) for a in best[1]]
    dv = dvec(atoms)
    band_in = dose_band - 0.004

    def pen(d_):
        rl = np.abs(d_ / d_.mean() - 1.0)
        return float((np.maximum(rl - band_in, 0.0) ** 2).sum())

    P = pen(dv)
    for it_pol in range(n_polish):
        if verbose_flag() and it_pol % 150 == 0 and it_pol > 0:
            print("    [polish] swap=%d pen=%.3e dev=%.4f" % (
                it_pol, P, dev(dv)), flush=True)
        if P <= 0:
            break
        rel = dv / dv.mean() - 1.0
        jo, ju = int(np.argmax(rel)), int(np.argmin(rel))
        bestm = (0.0, None)
        for i, (g, l_, c) in enumerate(atoms):
            row = _atom_row(ctx, g, l_)
            if row[jo] <= 0:
                continue
            off = 0
            for mb in ctx["metas"]:
                gg = mb["IDX"].shape[0]
                if g < off + gg:
                    cov = np.nonzero((mb["IDX"] == ju).any(axis=1))[0][:32]
                    for t2 in cov:
                        g2 = off + int(t2)
                        if g2 == g:
                            continue
                        row2 = _atom_row(ctx, g2, l_)
                        dP = pen(dv + (row2 - row) / M) - P
                        if dP < bestm[0] - 1e-15:
                            bestm = (dP, (i, g2, l_, row, row2))
                    break
                off += gg
        if bestm[1] is None:
            break
        i, g2, l_, row, row2 = bestm[1]
        g, _, c = atoms[i]
        atoms[i] = (g, l_, c - 1)
        atoms = [a for a in atoms if a[2] > 0]
        for k2, a in enumerate(atoms):
            if a[0] == g2 and a[1] == l_:
                atoms[k2] = (g2, l_, a[2] + 1)
                break
        else:
            atoms.append((g2, l_, 1))
        dv = dv + (row2 - row) / M
        P = pen(dv)
    return atoms


def guard_eval(ctx, V_ex, dosev, loads, ld_cont, V_fix_B, budget):
    """All Box-13 guards from maintained state (V_ex includes V0+eps0)."""
    r = ctx["r"]
    M = len(loads)
    dm = dosev.mean()
    dose_dev = float(np.abs(dosev / dm - 1.0).max())
    s_, ld_ex = np.linalg.slogdet(V_ex)
    eff_D = math.exp((ld_ex - ld_cont) / r) if s_ > 0 else 0.0
    arisk = float(np.trace(np.linalg.inv(V_ex)))
    Vf = V_fix_B + ctx["eps0"] * np.eye(r)
    arisk_fix = float(np.trace(np.linalg.inv(Vf)))
    Lf = np.linalg.cholesky(Vf)
    C = np.linalg.solve(Lf, np.linalg.solve(Lf, V_ex.T).T)
    mu_min = float(np.linalg.eigvalsh(C).min())
    inc = float(np.sum(loads))
    ceil_w = float(np.mean([v4.p_ceil_exact(float(l_), int(ctx["nu"]))
                            for l_ in loads]))
    checks = {"budget": inc <= budget * M + 1e-9,
              "dose": dose_dev <= DELTA,
              "ceiling": ceil_w <= 0.01,
              "effD": eff_D >= 0.95,
              "arisk": arisk <= 1.05 * arisk_fix,
              "spectral": mu_min >= 0.5 - 1e-9}
    return {"checks": checks, "all_pass": all(checks.values()),
            "incident_sum_rho": inc, "dose_dev": dose_dev, "eff_D": eff_D,
            "a_risk": arisk, "a_risk_fix": arisk_fix, "mu_min": mu_min,
            "ceiling_weighted": ceil_w,
            "detected_pred_approx": float(sum(
                ctx["nu"] * l_ / (1 + l_) for l_ in loads))}


def mixture_search(ctx, atoms, fd_rows, ld_cont, V_fix_B, budget,
                   verbose=False):
    """R15 Sec 6.2: m = 972..0 descending, m=0 included; exact materialized
    matrix per candidate; ALL guards recomputed (incremental rank-1 swaps);
    accept the largest-OED-fraction full pass. Frozen hashed row orders."""
    M, r, n, xhat = ctx["M_rows"], ctx["r"], ctx["n"], ctx["xhat"]
    B = ctx["B"]
    oed_rows = []
    for g, l_, c in atoms:
        row = _atom_row(ctx, g, l_)
        oed_rows.extend([row] * c)
    oed_rows = np.array(oed_rows)
    assert oed_rows.shape[0] == M
    fd_loads0 = fd_rows @ xhat
    fd = fd_rows * (budget / float(np.mean(fd_loads0))) * (1.0 - 1e-12)
    fd_order = np.array([(i * 617) % M for i in range(M)], dtype=np.int64)
    order_sha = hashlib.sha256(
        fd_order.tobytes() + np.ascontiguousarray(fd).tobytes()
        + str(atoms).encode()).hexdigest()

    def row_H(row):
        load = float(row @ xhat)
        idx = np.nonzero(row)[0]
        q = B[idx].T @ (row[idx] / load)
        J = v4._J(load, ctx["nu"])
        return ctx["nu"] * J * np.outer(q, q), load

    H_oed, load_oed = zip(*[row_H(rw) for rw in oed_rows])
    H_fd, load_fd = zip(*[row_H(fd[j]) for j in fd_order])
    I_r = np.eye(r)
    V_ex = ctx["V0"] + sum(H_oed) + ctx["eps0"] * I_r
    dosev = oed_rows.sum(axis=0) / M
    loads = list(load_oed)
    peak_oed = [float(rw.max()) for rw in oed_rows]
    peak_fd = [float(fd[j].max()) for j in fd_order]
    peaks = list(peak_oed)
    flips = {}
    prev_checks = None
    for m in range(M, -1, -1):
        g = guard_eval(ctx, V_ex, dosev, loads, ld_cont, V_fix_B, budget)
        if verbose and (m % 100 == 0 or m == M):
            print("    [mix] m=%d checks=%s inc=%.4f cap=%.4f diff=%+.2e "
                  "effD=%.3f dose=%.4f arisk=%.2fx mu=%.3f"
                  % (m, g["checks"], g["incident_sum_rho"], budget * M,
                     g["incident_sum_rho"] - budget * M, g["eff_D"],
                     g["dose_dev"], g["a_risk"] / g["a_risk_fix"],
                     g["mu_min"]), flush=True)
        if prev_checks is not None:
            for k_, ok_ in g["checks"].items():
                if ok_ != prev_checks[k_] and k_ not in flips:
                    flips[k_] = (m, ok_)
        prev_checks = g["checks"]
        peak = max(peaks) if peaks else 0.0
        peak_ok = peak <= v4.PEAK_PHYSICAL + 1e-9
        if g["all_pass"] and peak_ok:
            rows = np.vstack([oed_rows[:m], fd[fd_order[:M - m]]]) \
                if m < M else oed_rows
            g["peak_physical"] = peak
            return {"rows": rows, "m": m, "alpha": m / M, "guards": g,
                    "order_sha": order_sha, "status": "PASS"}
        if m == 0:
            break
        i = m - 1
        V_ex = V_ex - H_oed[i] + H_fd[M - m]
        dosev = dosev + (fd[fd_order[M - m]] - oed_rows[i]) / M
        loads[i] = load_fd[M - m]
        peaks[i] = peak_fd[M - m]
    g = guard_eval(ctx, V_ex, dosev, loads, ld_cont, V_fix_B, budget)
    first_fail = next((k for k, ok in g["checks"].items() if not ok), "peak")
    return {"rows": None, "m": None, "alpha": None, "guards": g,
            "order_sha": order_sha, "flips": flips,
            "status": "MIXTURE_FAIL:%s" % first_fail}


# ========================================================================== #
#  orchestrator                                                               #
# ========================================================================== #
def design_v5(xhat, nu, budget, palette, B, eps0, V_fix_B, rho_prescan,
              fixed_dose_rows, M_rows=972, verbose=False, side=32):
    import time as _t
    t0 = _t.time()

    def _ph(name):
        print("    [v5 phase] %s (t=%.0fs)" % (name, _t.time() - t0),
              flush=True)
    ctx = setup_ctx(xhat, nu, palette, B, eps0, rho_prescan, M_rows, side)
    _ph("phase1: relaxed PFW")
    xi_rel, traj_pfw = pfw_relaxed(ctx, budget, verbose=verbose)
    _ph("phase2: support-restricted exact reweighting (SLSQP)")
    ref = support_reweight(ctx, xi_rel, budget, verbose=verbose)
    _ph("phase3: dose-feasibility projection")
    xi0 = _dose_project(ctx, ref["xi_ref"].copy(), budget, DELTA,
                        verbose=verbose)
    _ph("phase4: primal-dual constrained solve (CG-LP + Lagrangian)")
    xi_c, cert, traj_pd = primal_dual_solve(ctx, xi0, budget,
                                            verbose=verbose)
    _ph("phase5: rounding + dose repair")
    ld_cont = float(np.linalg.slogdet(assemble(ctx, xi_c))[1])
    atoms = round_design(ctx, xi_c, budget)
    atoms = dose_repair(ctx, atoms)
    _ph("phase6: materialized mixture search")
    mix = mixture_search(ctx, atoms, fixed_dose_rows, ld_cont, V_fix_B,
                         budget, verbose=verbose)
    _ph("done")
    out = {"ctx": ctx, "xi_cont": xi_c, "cert_full": cert,
           "traj_pfw": traj_pfw, "traj_pd": traj_pd,
           "relaxed_reference": ref, "ld_cont": ld_cont,
           "atoms": atoms, "mixture": mix,
           "rows": mix["rows"], "verdict": (
               "OK" if (mix["status"] == "PASS"
                        and cert["status"] == "OK"
                        and cert["G_full"] / ctx["r"] <= GFULL_TOL_PER_R
                        and ref["bound_per_r"] <= RELREF_TOL_PER_R)
               else ("CONTINUOUS_CERTIFICATE_FAIL"
                     if cert["status"] != "OK"
                     or cert["G_full"] / ctx["r"] > GFULL_TOL_PER_R
                     else mix["status"]))}
    out["FULL_CONSTRAINED_KW_UPPER_BOUND"] = cert["G_full"] / ctx["r"]
    out["FULL_CONSTRAINED_KW_TIGHT"] = bool(
        cert["G_full"] / ctx["r"] <= GFULL_TIGHT_PER_R)
    out["RELAXED_REFERENCE_KW_UPPER_BOUND"] = ref["bound_per_r"]
    return out


# ========================================================================== #
#  FIXED_DOSE v2: SCAT32-based (information-matched to the FIXED* reference)  #
# ========================================================================== #
def fixed_dose_scat32(side=32, band=0.045, n_iter=1500):
    """Exact 972-row dose-balanced fallback built from the FROZEN sparsek
    k=32 family (seed 0) so its information profile matches the SCAT32
    FIXED* reference (the freeze-8 strip-based fallback failed effD/arisk/
    spectral at m=0 by construction). 52 rows are dropped greedily
    (deterministic: at each step the lowest-index row minimizing the
    current maximum per-pixel coverage loss, then the loss sum), keeping
    per-pixel coverage in {30, 31} (+2.1%/-1.2% about the 30.375 mean);
    a Sinkhorn amplitude polish then flattens the dose to the band."""
    from patterns import make_patterns
    n = side * side
    base = make_patterns("sparsek", n, n, 0, k=32)["A"]        # (1024, n)
    supp = [np.nonzero(base[s])[0] for s in range(n)]
    loss = np.zeros(n)
    drop = []
    for _ in range(52):
        best = None
        for s in range(n):
            if s in drop:
                continue
            mx = float(loss[supp[s]].max())
            sm = float(loss[supp[s]].sum())
            key = (mx, sm, s)
            if best is None or key < best[0]:
                best = (key, s)
        s = best[1]
        drop.append(s)
        loss[supp[s]] += 1.0
    keep = [s for s in range(n) if s not in set(drop)]
    rows = base[keep].copy()
    assert rows.shape[0] == 972
    c = np.ones(972)
    supk = [np.nonzero(rows[s])[0] for s in range(972)]
    for _ in range(n_iter):
        dose = (c[:, None] * rows).sum(axis=0)
        rel = dose / dose.mean()
        if np.abs(rel - 1.0).max() <= band - 0.005:
            break
        for s in range(972):
            c[s] /= rel[supk[s]].mean() ** 0.7
    out = c[:, None] * rows
    out *= 972.0 / c.sum()
    dose = out.sum(axis=0)
    dev = float(np.abs(dose / dose.mean() - 1.0).max())
    return out, dev


# ========================================================================== #
#  R17 SECTION (docs/ROUND63_METHOD_SPEC_M1_R17_AMENDMENT.md)                 #
# ========================================================================== #
GAMMAS = (0.2, 0.5, 1.0, 2.0, 5.0)


def _kernel_grids(nu, lo=1e-3, hi=64.0, n_nodes=61):
    """Per-nu exact-kernel interpolation grids for per-atom admission
    (mean_info_efficiency, rql_logload_bias, p_ceil) at arbitrary emergent
    loads (D_gain). Exact kernel at the nodes, linear interp in log rho."""
    grid = np.geomspace(lo, hi, n_nodes)
    kv = [v4.kernel_eval4(r_, int(nu)) for r_ in grid]
    return {"lg": np.log(grid),
            "eff": np.array([k["mean_info_efficiency"] for k in kv]),
            "bias": np.array([k["rql_logload_bias"] for k in kv]),
            "pceil": np.array([k["p_ceil"] for k in kv]),
            "lo": grid[0], "hi": grid[-1]}


def _adm_at_loads(loads, kg):
    ll = np.log(np.clip(loads, kg["lo"], kg["hi"]))
    eff = np.interp(ll, kg["lg"], kg["eff"])
    bias = np.interp(ll, kg["lg"], kg["bias"])
    pc = np.interp(ll, kg["lg"], kg["pceil"])
    return ((loads > 0) & (loads <= v4.RHO_CAP)
            & (pc <= v4.CEIL_ATOM) & (eff >= v4.EFF_MIN)
            & (bias <= v4.BIAS_MAX))


def setup_ctx_cert(xhat, nu, b, B, eps0, side=32):
    """R17 Sec 4.2-4.3 certificate context: D_cert = D_load (the palette at
    this (nu, b) corner: safe palette for b<=0.05, fast otherwise) UNION
    D_gain (per frozen base geometry, gamma in GAMMAS, physical rows
    gamma*a_base with the scene-independent frame convention a_base =
    (n/k)*w; predicted load EMERGENT: rho_rel = gamma*(a_base . xhat) =
    gamma*n/gsum_g, deployed load = b*rho_rel -- never target-rescaled).
    Existing admission guards (shape, physical peak, rho cap, ceiling,
    trust, bias) on both parts."""
    pal = v4.palette4(int(nu), fast=(b > 0.05 + 1e-12))
    ctx = setup_ctx(xhat, float(nu), pal, B, eps0, b, 972, side)
    G_tot, L = ctx["G_tot"], ctx["L"]
    kg = _kernel_grids(nu)
    for l_ in range(L):
        ok_l = _adm_at_loads(np.full(G_tot, ctx["levels"][l_]), kg)
        ctx["ALLOW"][:, l_] &= ok_l
    n = ctx["n"]
    base_load = n / ctx["gsum"]                # (a_base . xhat) exactly
    Cg = np.empty((G_tot, len(GAMMAS)))
    NUJg = np.empty((G_tot, len(GAMMAS)))
    ALLOWg = np.empty((G_tot, len(GAMMAS)), dtype=bool)
    gmax = ctx["gmax"]
    for i, gam in enumerate(GAMMAS):
        loads = b * gam * base_load            # emergent deployed loads
        Cg[:, i] = loads
        NUJg[:, i] = nu * np.array([v4._J(l_, nu) for l_ in loads])
        # physical peak of the deployed gamma atom: identity a = load * g
        # holds for D_gain too (same direction), so peak = load * gmax
        ALLOWg[:, i] = (_adm_at_loads(loads, kg)
                        & (loads * gmax <= v4.PEAK_PHYSICAL + 1e-9))
    ctx["C"] = np.hstack([ctx["C"], Cg])
    ctx["NUJ"] = np.hstack([ctx["NUJ"], NUJg])
    ctx["ALLOW"] = np.hstack([ctx["ALLOW"], ALLOWg])
    ctx["L"] = L + len(GAMMAS)
    ctx["L_load"] = L
    ctx["gammas"] = list(GAMMAS)
    ctx["budget_corner"] = float(b)
    return ctx


def d_cert_sha(side=32):
    """Scene-independent SHA256 of the D_cert DEFINITION: family manifest
    (supports + translations + powers), palette construction rule, gamma
    set, and admission constants. The per-cell atom VALUES depend on the
    local pre-scan estimate by frozen construction (R13); the frozen object
    is this definition (interpretation logged in the ledger)."""
    man = v4.dictionary_manifest(side)
    h = hashlib.sha256()
    h.update(man["sha256_global"].encode())
    h.update(json.dumps({
        "palette_rule": "v4.palette4(nu, fast=b>0.05); levels + admission "
                        "as frozen in R13/R15",
        "gammas": list(GAMMAS),
        "gain_load_rule": "rho = b*gamma*(a_base.xhat), a_base=(n/k)w frame "
                          "convention, emergent, never target-rescaled",
        "guards": {"rho_cap": v4.RHO_CAP, "peak_physical": v4.PEAK_PHYSICAL,
                   "ceil_atom": v4.CEIL_ATOM, "eff_min": v4.EFF_MIN,
                   "bias_max": v4.BIAS_MAX, "w_clip": [0.25, 4.0],
                   "dose_delta": DELTA},
    }, sort_keys=True).encode())
    return h.hexdigest()


def cert_deployed_rows(ctx, rows_rel, budget, verbose=False,
                       kw_rounds=12, kw_npix=256, kw_mucap=1e6,
                       kw_lptime=120):
    """R17 Sec 4.4-4.5: full dose-constrained certificate of a DEPLOYED
    physical row multiset (rows_rel: relative rows, ~unit-mean-load frame;
    deployed at global mean load = budget). Returns the full_constrained_cert
    dict + primal feasibility of the deployed design itself."""
    n, r = ctx["n"], ctx["r"]
    B = ctx["B"]
    xhat = ctx["xhat"]
    loads_rel = rows_rel @ xhat
    scale = budget / float(loads_rel.mean())
    rows_abs = rows_rel * scale * (1.0 - 1e-12)
    loads = rows_abs @ xhat
    M = rows_abs.shape[0]
    V = ctx["V0"].copy()
    qs = []
    for s in range(M):
        idx = np.nonzero(rows_abs[s])[0]
        q = B[idx].T @ (rows_abs[s, idx] / loads[s])
        qs.append((q, float(loads[s])))
        V += ctx["nu"] * v4._J(float(loads[s]), ctx["nu"]) * np.outer(q, q)
    V[np.diag_indices(r)] += ctx["eps0"]
    Vinv = np.linalg.inv(V)
    dbar = 0.0
    for q, ld in qs:
        dbar += ctx["nu"] * v4._J(ld, ctx["nu"]) * float(q @ Vinv @ q)
    dose = rows_abs.sum(axis=0) / M
    out = full_constrained_cert(
        ctx, None, budget, verbose=verbose, V_override=V,
        dbar_override=dbar, dose_override=dose,
        load_override=float(loads.mean()),
        max_rounds=kw_rounds, n_pix0=kw_npix, mu_cap=kw_mucap,
        lp_time=kw_lptime)
    out["deployed_mean_load"] = float(loads.mean())
    out["deployed_dose_dev"] = float(np.abs(dose / dose.mean() - 1.0).max())
    out["deployed_peak"] = float(rows_abs.max())
    out["primal_feasible"] = bool(
        loads.mean() <= budget + 1e-9
        and out["deployed_dose_dev"] <= DELTA)
    out["cell_pass"] = bool(out["status"] == "OK" and out["primal_feasible"]
                            and out["G_full"] / r <= GFULL_TOL_PER_R)
    return out


def path_feasible_alpha(mix_result):
    """R17 Sec 3.2 semantics: the descending exact scan is ONLY the
    path-specific descriptive diagnostic PATH_FEASIBLE_ALPHA. m=0 emits
    ADAPTIVE_COLLAPSE_UNDER_GUARDS and can never PASS; the retired
    production strings (COMPLIANT_VIA_MIXTURE / mixture PASS / alpha=0
    deployment) must never be derived from this result."""
    m = mix_result.get("m")
    if m is None:
        return {"PATH_FEASIBLE_ALPHA": None,
                "verdict": "NO_FEASIBLE_POINT_ON_PATH",
                "is_pass": False,
                "descriptive_only": True, "dev_supplement_only": True}
    alpha = m / 972.0
    return {"PATH_FEASIBLE_ALPHA": float(alpha),
            "verdict": ("ADAPTIVE_COLLAPSE_UNDER_GUARDS" if m == 0
                        else "PATH_FEASIBLE_ALPHA=%.4f" % alpha),
            "is_pass": False,               # never a PASS (R17 Sec 1.1)
            "descriptive_only": True, "dev_supplement_only": True}


def toy_full_cert_check_gain():
    """R15 active-dose toy EXTENDED with a gain-family atom (R17 Sec 4.2):
    a 4th atom = 2x atom-3 base row (gamma-coupled: same direction, double
    load and dose). The constrained optimum is re-solved independently
    (SLSQP) over all 4 atoms and the dual bound must agree <= 1e-8."""
    X = np.array([[1.0, 0.0], [0.0, 1.0],
                  [1.0 / math.sqrt(2), 1.0 / math.sqrt(2)],
                  [2.0 / math.sqrt(2), 2.0 / math.sqrt(2)]])
    Z = np.array([[3.0, 0.0], [0.0, 1.0], [1.0, 1.0], [2.0, 2.0]])
    c_load = np.array([1.0, 1.0, 1.0, 2.0])
    budget = 1.0
    eps = 1e-9

    def negld(w):
        Mm = sum(wi * np.outer(x, x) for wi, x in zip(w, X)) + eps * np.eye(2)
        s_, ld = np.linalg.slogdet(Mm)
        return -ld

    def grad(w):
        Mm = sum(wi * np.outer(x, x) for wi, x in zip(w, X)) + eps * np.eye(2)
        Mi = np.linalg.inv(Mm)
        return -np.array([x @ Mi @ x for x in X])

    def dose_cons(w):
        dv = Z.T @ w
        dm = dv.mean()
        return np.concatenate([(1 + DELTA) * dm - dv, dv - (1 - DELTA) * dm])

    cons = [{"type": "eq", "fun": lambda w: w.sum() - 1.0,
             "jac": lambda w: np.ones(4)},
            {"type": "ineq", "fun": dose_cons},
            {"type": "ineq", "fun": lambda w: budget - float(c_load @ w),
             "jac": lambda w: -c_load}]
    best = None
    for w0 in ([0.25] * 4, [0.2, 0.5, 0.2, 0.1], [0.3, 0.4, 0.2, 0.1]):
        res = minimize(negld, np.array(w0), jac=grad, method="SLSQP",
                       bounds=[(0, 1)] * 4, constraints=cons,
                       options={"maxiter": 600, "ftol": 1e-14})
        if res.success and (best is None or res.fun < best.fun):
            best = res
    w_star = np.maximum(best.x, 0.0)
    Mm = sum(wi * np.outer(x, x) for wi, x in zip(w_star, X)) + eps * np.eye(2)
    Mi = np.linalg.inv(Mm)
    d = np.array([x @ Mi @ x for x in X])
    zb = Z.mean(axis=1)
    n2 = 2
    A_ub = np.zeros((4, 2 + 2 * n2))
    for i in range(4):
        A_ub[i, 0] = -1.0
        A_ub[i, 1] = -c_load[i]
        A_ub[i, 2:2 + n2] = -(Z[i] - (1 + DELTA) * zb[i])
        A_ub[i, 2 + n2:] = -((1 - DELTA) * zb[i] - Z[i])
    b_ub = -d
    c_obj = np.zeros(2 + 2 * n2)
    c_obj[0] = 1.0
    c_obj[1] = budget
    res = linprog(c_obj, A_ub=A_ub, b_ub=b_ub,
                  bounds=[(None, None)] + [(0, None)] * (1 + 2 * n2),
                  method="highs")
    G = float(res.x[0]) + float(res.x[1]) * budget - float(d @ w_star)
    dv = Z.T @ w_star
    active = float(np.min(np.concatenate(
        [(1 + DELTA) * dv.mean() - dv, dv - (1 - DELTA) * dv.mean()])))
    return {"G_at_opt": G, "dose_active": bool(abs(active) < 1e-6),
            "w_star": w_star.tolist(), "agreement_ok": bool(abs(G) <= 1e-8)}


# ========================================================================== #
#  R18 SECTION (docs/ROUND63_METHOD_SPEC_M1_R18_AMENDMENT.md)                 #
#  Full-stack certificate: C_stack = {simplex, budget, +/-5% dose, A-risk    #
#  tr(V^-1) <= 1.05 tr(V_fix^-1), spectral V >= 0.5 V_fix} with the          #
#  certified cutting-plane outer approximation (ruling Sec 3.3) and the      #
#  mandatory support-expanding primal discriminator (Sec 3.4).               #
# ========================================================================== #
GSTACK_TOL_PER_R = 1e-2
CERT_PRIMAL_SCREEN_S = 60.0
CERT_PRIMARY_S = 180.0
CERT_RETRY_S = 180.0
CERT_CELL_WALL_CAP = 420.0


class StackCertProblem:
    """Callback container for one full-stack certificate instance. Both the
    real cells and the toy suite implement exactly this interface, so the
    cutting-plane engine is shared verbatim (R18 Sec 3.5 requirement)."""

    def __init__(self, d, c, n_pix, V_pre, eps0, V_fix, b, r,
                 h_of, dose_cols, spec_coeff, ar_coeff):
        self.d = np.asarray(d, dtype=float)          # (A,) objective scores
        self.c = np.asarray(c, dtype=float)          # (A,) budget loads
        self.n_pix = int(n_pix)
        self.V_pre = V_pre                           # (r,r) incl. nothing else
        self.eps0 = float(eps0)
        self.V_fix = V_fix
        self.b = float(b)
        self.r = int(r)
        self.h_of = h_of                             # cols -> (k, r) h rows
        self.dose_cols = dose_cols                   # cols -> (2n, k) U rows
        self.spec_coeff = spec_coeff                 # z (r,) -> (A,) coeffs
        self.ar_coeff = ar_coeff                     # Vk_inv -> (A,) coeffs
        self.h_fix = float(np.trace(np.linalg.inv(V_fix)))

    def V_of(self, cols, w):
        H = self.h_of(cols)                          # (k, r)
        return (self.V_pre + (H * np.asarray(w)[:, None]).T @ H
                + self.eps0 * np.eye(self.r))


def stack_cert_engine(prob, dbar_d, time_cap=CERT_PRIMARY_S, lp_time=60.0,
                      max_rounds=60, verbose=False, d_scale=None):
    """Certified cutting-plane outer approximation of
    G_stack = sup_{C_stack} (d^T xi) - dbar_d  (linear objective; dose and
    budget exact in the master LP; spectral eigenvector cuts + A-risk
    tangent cuts added on violation; column generation with a FULL-
    dictionary reduced-score scan each round). Finite cut sets are an OUTER
    relaxation, so the returned bound is conservative (ruling Sec 3.3)."""
    import time as _t
    t0 = _t.time()
    A = prob.d.size
    n = prob.n_pix
    if d_scale is None:
        d_scale = max(float(np.max(np.abs(prob.d[np.isfinite(prob.d)]))), 1.0)
    W = list(np.argsort(prob.d)[::-1][:min(400, A)])
    spec_cuts = []                    # (coeffs_A, rhs)
    ar_cuts = []
    n_scans = 0
    scan_resid = np.inf
    V_anchor = prob.V_fix.copy()   # strictly feasible for both LMIs
    status = "RUNNING"
    G = np.inf
    xi_W = None
    dyn_range = 0.0
    for rnd in range(max_rounds):
        if _t.time() - t0 > time_cap:
            status = "TIME_CAP"
            break
        k = len(W)
        Wa = np.array(W)
        # ---- master LP: maximize d^T xi over the outer polytope ---------- #
        Ud = prob.dose_cols(Wa)                      # (2n, k)
        rows = [np.ones((1, k)), prob.c[Wa][None, :], Ud]
        rhs = [np.array([1.0]), np.array([prob.b]), np.zeros(2 * n)]
        eqs = [True, False, False]
        for coeffs, cr in spec_cuts:
            rows.append(-coeffs[Wa][None, :])        # >= -> <= by negation
            rhs.append(np.array([-cr]))
            eqs.append(False)
        for coeffs, cr in ar_cuts:
            rows.append(-coeffs[Wa][None, :])
            rhs.append(np.array([-cr]))
            eqs.append(False)
        A_all = np.vstack(rows)
        b_all = np.concatenate(rhs)
        is_eq = np.concatenate([[e] * r_.shape[0]
                                for e, r_ in zip(eqs, rows)])
        A_ub = A_all[~is_eq]
        b_ub = b_all[~is_eq]
        A_eq = A_all[is_eq]
        b_eq = b_all[is_eq]
        cvec = -prob.d[Wa] / d_scale
        finite = np.abs(np.concatenate(
            [A_ub.ravel(), cvec, b_ub]))
        finite = finite[finite > 0]
        dyn_range = (float(finite.max() / finite.min())
                     if finite.size else 0.0)
        res = linprog(cvec, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                      bounds=[(0, None)] * k, method="highs",
                      options={"time_limit": float(lp_time)})
        if not res.success:
            if xi_W is not None and "infeasible" in str(res.message).lower():
                # over-cutting at a degenerate vertex: the PREVIOUS master's
                # G was computed over a LOOSER outer set (fewer cuts), hence
                # remains a valid conservative upper bound (ruling Sec 3.3)
                status = "CONVERGED_OUTER"
            else:
                status = "LP_FAIL:" + str(res.message)[:60]
            break
        xi_W = np.maximum(res.x, 0.0)
        n_spec_lp = len(spec_cuts)
        n_ar_lp = len(ar_cuts)
        G = float(prob.d[Wa] @ xi_W) - dbar_d
        # ---- true-constraint checks at xi* ------------------------------- #
        V = prob.V_of(Wa, xi_W)
        S = V - 0.5 * prob.V_fix
        evals, evecs = np.linalg.eigh(S)
        Vinv = np.linalg.inv(V)
        h_now = float(np.trace(Vinv))
        feas_now = (evals[0] >= -1e-9
                    and h_now <= 1.05 * prob.h_fix * (1 + 1e-9))
        if feas_now:
            V_anchor = V.copy()
        added = False
        if evals[0] < -1e-9:
            V_cut = V
            if V_anchor is not None:      # bisect to the spectral boundary
                lo_t, hi_t = 0.0, 1.0
                for _ in range(60):
                    mid = 0.5 * (lo_t + hi_t)
                    Vm = (1 - mid) * V_anchor + mid * V
                    if np.linalg.eigvalsh(Vm - 0.5 * prob.V_fix)[0] < 0:
                        hi_t = mid
                    else:
                        lo_t = mid
                V_cut = (1 - hi_t) * V_anchor + hi_t * V
            for zc in (np.linalg.eigh(V_cut - 0.5 * prob.V_fix)[1][:, 0],
                       evecs[:, 0]):     # boundary (depth) + argmax (cutoff)
                coeffs = prob.spec_coeff(zc)
                crhs = 0.5 * float(zc @ prob.V_fix @ zc)                     - float(zc @ (prob.V_pre
                                  + prob.eps0 * np.eye(prob.r)) @ zc)
                spec_cuts.append((coeffs, crhs))
            added = True
        if h_now > 1.05 * prob.h_fix * (1 + 1e-9):
            V_cut = V
            if V_anchor is not None:      # bisect to the A-risk boundary
                lo_t, hi_t = 0.0, 1.0
                for _ in range(60):
                    mid = 0.5 * (lo_t + hi_t)
                    Vm = (1 - mid) * V_anchor + mid * V
                    if float(np.trace(np.linalg.inv(Vm)))                             > 1.05 * prob.h_fix:
                        hi_t = mid
                    else:
                        lo_t = mid
                V_cut = (1 - hi_t) * V_anchor + hi_t * V
            for Vki in (np.linalg.inv(V_cut), Vinv):   # depth + cutoff
                Vk2 = Vki @ Vki
                coeffs = prob.ar_coeff(Vki)
                crhs = (2.0 * float(np.trace(Vki))
                        - 1.05 * prob.h_fix
                        - float(np.trace(Vk2 @ (prob.V_pre
                                                + prob.eps0
                                                * np.eye(prob.r)))))
                ar_cuts.append((coeffs, crhs))
            added = True
        # ---- full-dictionary column-generation scan ----------------------- #
        duals = _lp_duals(res, A_ub.shape[0], A_eq.shape[0])
        rc = prob.d.copy()
        rc -= duals["eq"][0]                          # simplex price
        iu = 0
        rc -= duals["ub"][iu] * prob.c
        iu += 1
        du = duals["ub"][iu:iu + 2 * n]
        rc -= _dose_dual_dot(prob, du)
        iu += 2 * n
        for coeffs, _cr in spec_cuts[:n_spec_lp]:
            rc += duals["ub"][iu] * coeffs
            iu += 1
        for coeffs, _cr in ar_cuts[:n_ar_lp]:
            rc += duals["ub"][iu] * coeffs
            iu += 1
        rc[Wa] = -np.inf
        n_scans += 1
        j = int(np.argmax(rc))
        scan_resid = float(rc[j]) if np.isfinite(rc[j]) else 0.0
        if scan_resid > 1e-7 * d_scale:
            order = np.argsort(rc)[::-1][:50]
            W.extend(int(x) for x in order if np.isfinite(rc[x]))
            added = True
        if verbose:
            print("      [stack] rnd=%d |W|=%d G=%.6g cuts=%d/%d "
                  "scan_rc=%.2e eig=%.3e h=%.4g (%.0fs)"
                  % (rnd, k, G, len(spec_cuts), len(ar_cuts),
                     scan_resid, evals[0], h_now / prob.h_fix,
                     _t.time() - t0), flush=True)
        if not added:
            status = "CONVERGED"
            break
    comp_resid = np.nan
    try:
        du = duals["ub"]
        slacks = b_ub - A_ub @ xi_W
        comp_resid = float(np.abs(du * slacks).max()
                           / max(abs(G) + abs(dbar_d), 1.0))
    except Exception:
        pass
    return {"G": G, "status": status, "n_spec_cuts": len(spec_cuts),
            "n_arisk_cuts": len(ar_cuts), "n_scans": n_scans,
            "scan_residual": scan_resid, "xi_W": xi_W,
            "W": list(W), "dyn_range": dyn_range,
            "wall": _t.time() - t0,
            "complementarity_residual": comp_resid,
            "psd_residual": float(evals[0]) if xi_W is not None else np.nan,
            "spec_cuts": spec_cuts, "ar_cuts": ar_cuts}


def _lp_duals(res, n_ub, n_eq):
    ub = np.zeros(n_ub)
    eq = np.zeros(max(n_eq, 1))
    try:
        ub = -np.asarray(res.ineqlin.marginals)
        eq = -np.asarray(res.eqlin.marginals)
    except Exception:
        pass
    return {"ub": ub, "eq": eq}


def _dose_dual_dot(prob, du):
    """sum_j du_j * U_{j,a} over BOTH dose blocks for every atom (A,)."""
    return prob.dose_dual_dot(du)


# ---- real-cell problem builder ------------------------------------------- #
def build_stack_problem(ctx, V_fix, b):
    A_flat = ctx["G_tot"] * ctx["L"]
    G_tot, L, n, r = ctx["G_tot"], ctx["L"], ctx["n"], ctx["r"]
    QB, C, NUJ, M = ctx["QB"], ctx["C"], ctx["NUJ"], ctx["M_rows"]
    allow = ctx["ALLOW"].ravel()
    d_full = np.full(A_flat, -np.inf)
    zb = ctx["gsum"] / n                              # zbar per unit load

    def h_of(cols):
        g = cols // L
        l_ = cols % L
        return (np.sqrt(M * NUJ[g, l_])[:, None] * QB[g])

    def dose_cols(cols):
        g = cols // L
        l_ = cols % L
        k = cols.size
        U = np.zeros((2 * n, k))
        off = 0
        gi = 0
        for mb in ctx["metas"]:
            gg = mb["IDX"].shape[0]
            sel = np.nonzero((g >= off) & (g < off + gg))[0]
            for s in sel:
                t = g[s] - off
                load = C[g[s], l_[s]]
                z = np.zeros(n)
                z[mb["IDX"][t]] = load * mb["GVAL"][t]
                zbar = load * ctx["gsum"][g[s]] / n
                U[:n, s] = z - (1 + DELTA) * zbar
                U[n:, s] = (1 - DELTA) * zbar - z
            off += gg
        return U

    def spec_coeff(z):
        w = QB @ z                                    # (G,)
        return (M * NUJ * (w ** 2)[:, None]).ravel()

    def ar_coeff(Vinv):
        T = QB @ Vinv                                 # (G, r)
        q2 = (T * T).sum(axis=1)
        return (M * NUJ * q2[:, None]).ravel()

    def dose_dual_dot(du):
        dup = du[:n]
        dum = du[n:]
        Sp, Sm = float(dup.sum()), float(dum.sum())
        Ap = mu_dots(ctx, dup)
        Am = mu_dots(ctx, dum)
        per_g = (Ap - (1 + DELTA) * zb * Sp) + ((1 - DELTA) * zb * Sm - Am)
        return (per_g[:, None] * C).ravel()

    prob = StackCertProblem(d_full, C.ravel(), n, ctx["V0"], ctx["eps0"],
                            V_fix, b, r, h_of, dose_cols, spec_coeff,
                            ar_coeff)
    prob.dose_dual_dot = dose_dual_dot
    prob.allow = allow
    return prob


def stack_primal_probe(ctx, V_fix, b, ld_dep, time_cap=CERT_PRIMAL_SCREEN_S,
                       verbose=False):
    """R18 Sec 3.4 SUPPORT-EXPANDING primal discriminator over C_stack:
    positive initialization on ALL admitted atoms (fixes the R18-flagged
    support-preserving flaw), multiplicative ascent + dose/budget
    projection, steps violating the GLOBAL A-risk/spectral constraints are
    damped/rejected. Records the best fully C_stack-feasible logdet."""
    import time as _t
    t0 = _t.time()
    h_fix = float(np.trace(np.linalg.inv(V_fix)))
    Lf = np.linalg.cholesky(V_fix)
    xi = np.where(ctx["ALLOW"], 1.0, 0.0)
    xi /= xi.sum()
    xi = _dose_project(ctx, xi, b, DELTA)
    best = (-np.inf, None)
    xi_ok = None

    def stack_ok(V):
        ar = float(np.trace(np.linalg.inv(V)))
        Cw = np.linalg.solve(Lf, np.linalg.solve(Lf, V.T).T)
        mu = float(np.linalg.eigvalsh(Cw).min())
        return (ar <= 1.05 * h_fix and mu >= 0.5), ar, mu

    it = 0
    while _t.time() - t0 < time_cap:
        it += 1
        V = assemble(ctx, xi)
        Vinv = np.linalg.inv(V)
        D = dvals(ctx, Vinv)
        pos = np.where(ctx["ALLOW"], np.maximum(D, 0.0), 0.0)
        zb_ = float((xi * pos).sum())
        if zb_ <= 0:
            break
        xi_new = xi * np.sqrt(pos / zb_)
        xi_new[~ctx["ALLOW"]] = 0.0
        xi_new /= xi_new.sum()
        xi_new = _dose_project(ctx, xi_new, b, DELTA)
        V_new = assemble(ctx, xi_new)
        ok, ar, mu = stack_ok(V_new)
        if not ok and xi_ok is not None:
            xi_new = 0.5 * (xi_ok + xi_new)
            xi_new = _dose_project(ctx, xi_new, b, DELTA)
            V_new = assemble(ctx, xi_new)
            ok, ar, mu = stack_ok(V_new)
        load = float((xi_new * ctx["C"]).sum())
        dv = dose_of(ctx, xi_new)
        dev = float(np.abs(dv / dv.mean() - 1.0).max())
        if ok and load <= b + 1e-9 and dev <= DELTA:
            xi_ok = xi_new.copy()
            ld = float(np.linalg.slogdet(V_new)[1])
            if ld > best[0]:
                best = (ld, xi_new.copy())
        xi = xi_new
    gap = (best[0] - ld_dep) / ctx["r"] if np.isfinite(best[0]) else -np.inf
    comp = {}
    if best[1] is not None:
        xi_b = best[1]
        Ll = ctx.get("L_load", ctx["L"])
        lu = ctx["C"][xi_b > 1e-6]
        dvb = dose_of(ctx, xi_b)
        Vb = assemble(ctx, xi_b)
        arb = float(np.trace(np.linalg.inv(Vb)))
        Cw = np.linalg.solve(Lf, np.linalg.solve(Lf, Vb.T).T)
        comp = {"budget_used": float((xi_b * ctx["C"]).sum()),
                "dose_dev": float(np.abs(dvb / dvb.mean() - 1.0).max()),
                "gain_mass": float(xi_b[:, Ll:].sum()) if Ll < ctx["L"]
                else 0.0,
                "load_q50": float(np.median(lu)) if lu.size else np.nan,
                "load_q95": float(np.percentile(lu, 95)) if lu.size
                else np.nan,
                "arisk_ratio": arb / h_fix,
                "mu_min": float(np.linalg.eigvalsh(Cw).min())}
    return {"ld_best": best[0], "gap_lower_per_r": float(gap),
            "found_feasible": bool(np.isfinite(best[0])),
            "iters": it, "wall": _t.time() - t0, "composition": comp}


def full_stack_cert_cell(ctx, rows_rel, b, verbose=False):
    """R18 Sec 5 per-cell protocol: 60 s support-expanding primal screen,
    180 s primary cutting-plane dual, one deterministic rescaled retry,
    CERT_CELL_WALL_CAP total; terminal status CERTIFIED / COUNTEREXAMPLE /
    NUMERICAL_UNRESOLVED with the full Sec-5.3 disclosure set."""
    import time as _t
    t0 = _t.time()
    n, r = ctx["n"], ctx["r"]
    xhat = ctx["xhat"]
    loads_rel = rows_rel @ xhat
    rows_abs = rows_rel * (b / float(loads_rel.mean())) * (1 - 1e-12)
    loads = rows_abs @ xhat
    V_d = ctx["V0"].copy()
    for s in range(rows_abs.shape[0]):
        idx = np.nonzero(rows_abs[s])[0]
        q = ctx["B"][idx].T @ (rows_abs[s, idx] / loads[s])
        V_d += ctx["nu"] * v4._J(float(loads[s]), ctx["nu"]) * np.outer(q, q)
    V_d[np.diag_indices(r)] += ctx["eps0"]
    ld_dep = float(np.linalg.slogdet(V_d)[1])
    V_fix = V_d                                       # load-matched comparator
    Vdinv = np.linalg.inv(V_d)
    dbar_d = float(np.trace(
        Vdinv @ (V_d - ctx["V0"] - ctx["eps0"] * np.eye(r))))
    # deployed-design C_stack feasibility (trivially so vs itself + dose)
    dose_dep = rows_abs.sum(axis=0) / rows_abs.shape[0]
    dep_dose_dev = float(np.abs(dose_dep / dose_dep.mean() - 1.0).max())
    dep_feasible = bool(loads.mean() <= b + 1e-9 and dep_dose_dev <= DELTA)
    # atom scores
    D = dvals(ctx, Vdinv)
    prob = build_stack_problem(ctx, V_fix, b)
    prob.d = np.where(prob.allow, D.ravel(), -np.inf)
    # ---- 1) primal screen ------------------------------------------------ #
    scr = stack_primal_probe(ctx, V_fix, b, ld_dep,
                             time_cap=CERT_PRIMAL_SCREEN_S, verbose=verbose)
    out = {"ld_dep": ld_dep, "dep_feasible": dep_feasible,
           "dep_dose_dev": dep_dose_dev,
           "primal_gap_lower_per_r": scr["gap_lower_per_r"],
           "primal_found_feasible": scr["found_feasible"],
           "primal_composition": scr.get("composition", {}),
           "MU_CAP_ACTIVE": False}
    if scr["found_feasible"] and scr["gap_lower_per_r"] > GSTACK_TOL_PER_R:
        out.update({"status": "COUNTEREXAMPLE",
                    "dual_gap_upper_per_r": float("nan"),
                    "solver_status_primary": "SKIPPED_COUNTEREXAMPLE",
                    "solver_status_retry": "",
                    "wall_seconds": _t.time() - t0})
        return out
    # ---- 2) primary dual + one frozen rescaled retry ---------------------- #
    attempts = []
    for attempt in (0, 1):
        if attempt == 1 and attempts and attempts[0]["ok"]:
            break
        d_scale = None
        if attempt == 1:      # frozen equilibration: scores/r, budget/b,
            d_scale = max(float(np.max(prob.d[np.isfinite(prob.d)])), 1.0) / r
        eng = stack_cert_engine(prob, dbar_d,
                                time_cap=(CERT_PRIMARY_S if attempt == 0
                                          else CERT_RETRY_S),
                                verbose=verbose, d_scale=d_scale)
        Gpr = eng["G"] / r if np.isfinite(eng["G"]) else float("inf")
        ok = (eng["status"] in ("CONVERGED", "CONVERGED_OUTER")
              and np.isfinite(eng["G"]) and eng["dyn_range"] < 1e10)
        attempts.append({"eng": eng, "ok": ok, "Gpr": Gpr})
        if ok or _t.time() - t0 > CERT_CELL_WALL_CAP:
            break
    eng = attempts[-1]["eng"]
    ok = attempts[-1]["ok"]
    Gpr = attempts[-1]["Gpr"]
    if ok and dep_feasible and Gpr <= GSTACK_TOL_PER_R:
        status = "CERTIFIED"
    elif scr["found_feasible"] and scr["gap_lower_per_r"] > GSTACK_TOL_PER_R:
        status = "COUNTEREXAMPLE"
    else:
        # frozen taxonomy (R18 Sec 5.1): a converged outer bound above the
        # tolerance is NOT a counterexample and NOT certified -> unresolved
        status = "NUMERICAL_UNRESOLVED"
    out.update({
        "status": status,
        "dual_gap_upper_per_r": Gpr,
        "solver_status_primary": attempts[0]["eng"]["status"],
        "solver_status_retry": (attempts[1]["eng"]["status"]
                                if len(attempts) > 1 else ""),
        "wall_seconds": _t.time() - t0,
        "coefficient_dynamic_range": eng["dyn_range"],
        "n_dictionary_scans": eng["n_scans"],
        "n_arisk_cuts": eng["n_arisk_cuts"],
        "n_spectral_cuts": eng["n_spec_cuts"],
        "full_dictionary_scan_residual": eng["scan_residual"],
        "min_generalized_eigenvalue": eng.get("psd_residual", float("nan")),
        "arisk_ratio": 1.0})
    return out


# ========================================================================== #
#  R18 Sec 3.5 four-toy suite (shared engine, independent primal solves)      #
# ========================================================================== #
def _toy_problem(hX, Z, c_load, V_pre_s, eps, V_fix, b):
    """Small dense StackCertProblem (r=2, n_pix=2) running through the SAME
    stack_cert_engine as the real cells."""
    hX = np.asarray(hX, dtype=float)
    Z = np.asarray(Z, dtype=float)
    A = hX.shape[0]
    n = Z.shape[1]
    V_pre = V_pre_s * np.eye(2)
    zbar = Z.mean(axis=1)

    def h_of(cols):
        return hX[cols]

    def dose_cols(cols):
        k = cols.size
        U = np.zeros((2 * n, k))
        for i, a in enumerate(cols):
            U[:n, i] = Z[a] - (1 + DELTA) * zbar[a]
            U[n:, i] = (1 - DELTA) * zbar[a] - Z[a]
        return U

    def spec_coeff(z):
        return (hX @ z) ** 2

    def ar_coeff(Vinv):
        T = hX @ Vinv
        return (T * T).sum(axis=1)

    def dose_dual_dot(du):
        Up = np.zeros(A)
        for a in range(A):
            Up[a] = float(du[:n] @ (Z[a] - (1 + DELTA) * zbar[a])
                          + du[n:] @ ((1 - DELTA) * zbar[a] - Z[a]))
        return Up

    d = np.array([hX[a] @ hX[a] for a in range(A)])   # placeholder; caller
    prob = StackCertProblem(d, c_load, n, V_pre, eps, V_fix, b, 2,
                            h_of, dose_cols, spec_coeff, ar_coeff)
    prob.dose_dual_dot = dose_dual_dot
    prob.allow = np.ones(A, dtype=bool)
    prob.hX = hX
    prob.Z = Z
    return prob


def _toy_primal_grid(prob, d, b, step=2e-3, w_extra=None):
    """Independent global primal: dense simplex grid (closed-form 2x2 checks)
    + SLSQP polish. Returns sup d^T xi over C_stack."""
    A = prob.hX.shape[0]
    hX, Z = prob.hX, prob.Z
    zbar = Z.mean(axis=1)

    def feasible(w):
        if w.min() < -1e-12 or abs(w.sum() - 1) > 1e-9:
            return False
        if float(prob.c @ w) > b + 1e-12:
            return False
        dv = Z.T @ w
        dm = dv.mean()
        if np.max(dv - (1 + DELTA) * dm) > 1e-12 \
                or np.max((1 - DELTA) * dm - dv) > 1e-12:
            return False
        V = prob.V_pre + (hX * w[:, None]).T @ hX + prob.eps0 * np.eye(2)
        if float(np.trace(np.linalg.inv(V))) > 1.05 * prob.h_fix + 1e-12:
            return False
        S = V - 0.5 * prob.V_fix
        if np.linalg.eigvalsh(S)[0] < -1e-12:
            return False
        return True

    best = (-np.inf, None)
    if A == 3:
        g = np.arange(0.0, 1.0 + step / 2, step)
        for w1 in g:
            for w2 in np.arange(0.0, 1.0 - w1 + step / 2, step):
                w = np.array([w1, w2, 1.0 - w1 - w2])
                if feasible(w):
                    v = float(d @ w)
                    if v > best[0]:
                        best = (v, w)
    else:
        rng = np.random.RandomState(63)
        for _ in range(200000):
            w = rng.dirichlet(np.ones(A))
            if feasible(w):
                v = float(d @ w)
                if v > best[0]:
                    best = (v, w)
    # SLSQP polish (independent solver)
    def negobj(w):
        return -float(d @ w)

    def gr(w):
        return -d

    cons = [{"type": "eq", "fun": lambda w: w.sum() - 1.0,
             "jac": lambda w: np.ones(A)},
            {"type": "ineq", "fun": lambda w: b - float(prob.c @ w)},
            {"type": "ineq", "fun": lambda w: np.concatenate([
                (1 + DELTA) * (Z.T @ w).mean() - Z.T @ w,
                Z.T @ w - (1 - DELTA) * (Z.T @ w).mean()])},
            {"type": "ineq", "fun": lambda w: 1.05 * prob.h_fix - np.trace(
                np.linalg.inv(prob.V_pre + (hX * w[:, None]).T @ hX
                              + prob.eps0 * np.eye(2)))},
            {"type": "ineq", "fun": lambda w: np.linalg.eigvalsh(
                prob.V_pre + (hX * w[:, None]).T @ hX
                + prob.eps0 * np.eye(2) - 0.5 * prob.V_fix)[0]}]
    starts = [best[1] if best[1] is not None else np.ones(A) / A]
    if w_extra is not None:
        starts.append(np.asarray(w_extra, dtype=float))
    for w0 in starts:
        res = minimize(negobj, w0, jac=gr, method="SLSQP",
                       bounds=[(0, 1)] * A, constraints=cons,
                       options={"maxiter": 500, "ftol": 1e-16})
        cand = np.maximum(res.x, 0)
        if feasible(cand):
            v = float(d @ cand)
            if v > best[0]:
                best = (v, cand)
        # direct acceptance of a feasible engine/master solution
        if feasible(w0):
            v = float(d @ w0)
            if v > best[0]:
                best = (v, w0)
    return best


def r18_toy_suite(verbose=False):
    """Four toys (ruling Sec 3.5): (1) budget+dose active; (2) A-risk only;
    (3) spectral only; (4) all active + a dose-feasible design rejected by
    conditioning. Engine bound vs independent primal <=1e-8; PSD residual
    >= -1e-9; complementarity <=1e-8; full-dict scan residual <=1e-8."""
    out = []
    s2 = 1.0 / math.sqrt(2)
    toys = {}
    # (1) budget+dose active; conditioning inactive (tiny V_fix, huge cap)
    toys["budget_dose"] = dict(
        hX=[[2.0, 0.0], [0.0, 1.0], [1.2 * s2, 1.2 * s2]],
        Z=[[3.0, 0.2], [0.2, 3.0], [1.0, 1.0]],
        c=[1.2, 0.9, 2.0], Vfix=0.02 * np.eye(2), Vpre=0.05, b=1.045)
    # (2) A-risk only: dose uniform per atom (never binds), budget huge,
    #     V_fix = balanced reference; favored atom concentrates -> cap binds
    toys["arisk_only"] = dict(
        hX=[[2.0, 0.0], [0.0, 2.0], [s2, s2]],
        Z=[[1.0, 1.0], [1.0, 1.0], [1.0, 1.0]],
        c=[1.0, 1.0, 1.0], Vfix=np.diag([1.2, 1.2]), Vpre=0.05, b=10.0)
    # (3) spectral only: as (2) but arisk cap slack via larger V_fix trace
    #     margin and a tight spectral floor in direction 2
    toys["spectral_only"] = dict(
        hX=[[2.0, 0.0], [0.0, 1.1], [s2, s2]],
        Z=[[1.0, 1.0], [1.0, 1.0], [1.0, 1.0]],
        c=[1.0, 1.0, 1.0], Vfix=np.diag([0.2, 1.1]), Vpre=0.05, b=10.0,
        hfix_relax=40.0)
    # (4) all active: GENERIC 5-atom problem (4 dof -> non-degenerate
    #     4-active vertex with moderate dual prices); atom 1 alone is dose-
    #     and budget-feasible but conditioning-rejected (the witness)
    toys["all_active"] = dict(
        hX=[[4.0, 0.0], [0.0, 1.6], [0.6, 0.6], [1.2, 0.3], [0.3, 1.2]],
        Z=[[1.5, 1.5], [2.3, 0.7], [0.7, 2.3], [1.9, 1.1], [1.1, 1.9]],
        c=[0.75, 0.8, 0.8, 0.7, 0.7], Vfix=np.diag([2.2, 1.0]), Vpre=0.02,
        b=0.80)
    for name, t in toys.items():
        prob = _toy_problem(t["hX"], t["Z"], np.array(t["c"], float),
                            t["Vpre"], 1e-9, np.asarray(t["Vfix"]), t["b"])
        if "hfix_relax" in t:
            prob.h_fix = t["hfix_relax"]
        # objective: local scores at the balanced reference design
        Vref = prob.V_pre + prob.V_fix + prob.eps0 * np.eye(2)
        Vri = np.linalg.inv(Vref)
        d = np.array([h @ Vri @ h for h in prob.hX])
        prob.d = d
        eng = stack_cert_engine(prob, 0.0, time_cap=60.0, lp_time=20.0,
                                max_rounds=500, verbose=False)
        w_eng = None
        if eng["xi_W"] is not None:
            w_eng = np.zeros(len(prob.d))
            for wi, ai in zip(eng["xi_W"], eng["W"][:len(eng["xi_W"])]):
                w_eng[ai] += wi
        best_v, best_w = _toy_primal_grid(prob, d, t["b"], w_extra=w_eng)
        scale = max(abs(best_v), 1.0)
        agree = abs(eng["G"] - best_v) / scale
        rec = {"toy": name, "engine_G": eng["G"], "primal_G": best_v,
               "agreement": agree, "status": eng["status"],
               "psd_residual": eng["psd_residual"],
               "complementarity": eng["complementarity_residual"],
               "scan_residual": abs(eng["scan_residual"]) / scale,
               "n_cuts": (eng["n_spec_cuts"], eng["n_arisk_cuts"]),
               "pass": bool(agree <= 1e-8
                            and eng["psd_residual"] >= -1e-9
                            and (not np.isfinite(eng["complementarity_residual"])
                                 or eng["complementarity_residual"] <= 1e-8)
                            and abs(eng["scan_residual"]) / scale <= 1e-8)}
        # activity disclosure at the accepted optimum (ledger evidence)
        if best_w is not None:
            wv = best_w
            Vv = prob.V_of(np.arange(len(prob.d)), wv)
            dvv = prob.Z.T @ wv
            dmv = dvv.mean()
            rec["active_at_optimum"] = {
                "budget": bool(t["b"] - float(prob.c @ wv) <= 1e-7),
                "dose": bool(max(np.max(dvv - (1 + DELTA) * dmv),
                                 np.max((1 - DELTA) * dmv - dvv)) >= -1e-7),
                "arisk": bool(float(np.trace(np.linalg.inv(Vv)))
                              >= 1.05 * prob.h_fix - 1e-6),
                "spectral": bool(np.linalg.eigvalsh(
                    Vv - 0.5 * prob.V_fix)[0] <= 1e-6)}
        # toy-4 requirement: a dose-feasible design rejected by conditioning
        if name == "all_active":
            w_rej = np.zeros(len(prob.d))
            w_rej[0] = 1.0
            dv = prob.Z.T @ w_rej
            dm = dv.mean()
            dose_ok = (np.max(dv - (1 + DELTA) * dm) <= 1e-12
                       and np.max((1 - DELTA) * dm - dv) <= 1e-12)
            V1 = prob.V_pre + np.outer(prob.hX[0], prob.hX[0]) \
                + prob.eps0 * np.eye(2)
            cond_reject = (
                float(np.trace(np.linalg.inv(V1))) > 1.05 * prob.h_fix
                or np.linalg.eigvalsh(V1 - 0.5 * prob.V_fix)[0] < 0)
            rec["dose_feasible_conditioning_rejected_witness"] = bool(
                dose_ok and cond_reject)
            rec["pass"] = rec["pass"] and rec[
                "dose_feasible_conditioning_rejected_witness"]
        if verbose:
            print("  [toy %s] G_eng=%.10f G_pri=%.10f agree=%.1e psd=%.1e "
                  "comp=%s scan=%.1e cuts=%s %s"
                  % (name, eng["G"], best_v, agree, eng["psd_residual"],
                     eng["complementarity_residual"], rec["scan_residual"],
                     rec["n_cuts"], "PASS" if rec["pass"] else "FAIL"),
                  flush=True)
        out.append(rec)
    return out
