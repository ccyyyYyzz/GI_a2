"""ROUND63 lambda_TV selection + model-adequacy audit — F1-FROZEN rule.

Source of record: docs/ROUND63_GPT_ROUND4_DIGEST.md (GPT round-4 ruling,
"outer AUDIT split + coherent refit bootstrap"). This file supersedes the
round-3 six-step draft rule; the changes and their reasons:

  * OUTER 80/20 DEV/AUDIT split per cell (logical groups; hadpair pairs are
    atomic). AUDIT never touches lam_max / eta selection / initializers.
  * eta* per arm from a grouped K=5 fold cross-fit INSIDE DEV, one-SE rule on
    per-fold MEAN squared calibrated renewal residuals. GOF is fully OUT of
    the acceptance set: adequacy can never change the estimator.
  * Model adequacy = MODEL_FAIL_PREDICTIVE, tested ONCE per cell by RQL (not
    per baseline arm): a coherent plugin scene x_plugin (RQL fit on DEV at
    eta_min, 25 iters) predicts the untouched AUDIT; B=39 refit-per-replicate
    parametric bootstrap replays the same pipeline coherently (generator =
    A @ x_plugin + dark, never concatenated cross-fitted lambdas); exact
    Monte-Carlo p-value, MODEL_FAIL_PREDICTIVE <=> p <= 0.025 (level 1/40).
    Semantics: "the calibrated renewal detector model plus the frozen RQL-TV
    reconstruction procedure cannot explain independent audit data" — NOT a
    pure likelihood-class test (at M < n, detector misfit and null-space
    under-recovery are not strictly separable; stated in the paper).
  * Fixed-lambda-hat bootstrap survives only as the LOWER-tail leakage
    diagnostic (LEAKAGE_SUSPECT, p_low <= 0.01). mean-r / load-correlation are
    independent WARNINGS, never part of the MODEL_FAIL union.
  * MODEL_FAIL never alters eta*, never removes a cell, never consults PSNR.
  * All RNG keys are integers (the round-4 audit caught int(round(T)) == 0 for
    physical T, which had different-nu cells sharing bootstrap streams).

Frozen constants (F1): AUDIT_FRAC=0.2, MIN_AUDIT_GROUPS=128, K=5,
ETA_GRID as below, TV_NULL_REL=1e-3, lam_max expansion cap 40 / 26 log
bisections, N_SEL=60 path iters, N_AUDIT=25 audit iters, B_AUDIT=39
(early-stop after 19 all-below), P_FAIL=0.025, B_LEAK=199, P_LEAK=0.01.

RNG stream layout (all-integer keys, disjoint by the trailing tag):
  outer split      rng_for(*cell_key, seed, 63, 4, 0)
  inner DEV folds  rng_for(*cell_key, seed, 63, 4, 1)
  audit bootstrap  rng_for(*cell_key, seed, 63, 4, 2, tau_ns)
  leak bootstrap   rng_for(*cell_key, seed, 63, 4, 3, tau_ns)
cell_key = (kind_id, rho_milli, nu_int, M, side) — integers built by the
caller (campaign.run_cell / probes); tau_ns = int(round(tau * 1e9)).

API (campaign.py contract):
  split_dev_audit(M, meta, cell_key, seed) -> dict split
  select_eta_arm(arm_name, A, b, ctx, cell_key, seed, split) -> sel dict
  audit_cell(A, b, ctx, cell_key, seed, split, sel_rql) -> audit dict
  select_eta_and_fit(arm_name, A, b, ctx, cell_key=None, seed=0,
                     run_audit=None) -> (x_final, info)
    run_audit defaults to (arm_name == "RQL"). info keys kept from the old
    rule where they still exist (eta_star, lam_tv, lam_max_arm, eta_grid,
    d_bar_curve, SE_min, ...) plus the audit block on RQL.
float64 throughout.
"""
import os
import sys
from dataclasses import replace

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))

from gi_core.utils import rng_for
from physics import qmle_mean_var, simulate_counts
from solvers import (tv_fista, tv_value, init_gi_flux, run_arm,
                     _arm_factory, _qmle_irls_factory)

# ---------------------------------------------------------------- F1 constants
ETA_GRID = (1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0)
K_FOLDS = 5
AUDIT_FRAC = 0.20
MIN_AUDIT_GROUPS = 128
TV_NULL_REL = 1e-3
N_SEL = 60          # selection-path / lam_max fits (outer iterations)
N_AUDIT = 25        # plugin + per-replicate audit fits
B_AUDIT = 39        # refit bootstrap replicates (exact 1/40 MC level)
B_EARLY = 19        # early-stop block: any exceedance in first 19 => PASS
P_FAIL = 0.025
B_LEAK = 199
P_LEAK = 0.01
LAM_FLOOR_REL = 1e-6


def _default_cell_key(ctx, M):
    """Fallback integer cell key when the caller supplies none (probes/smoke).
    Production (campaign.run_cell) always passes the real grid coordinates."""
    nu_int = int(round(ctx.T / ctx.det.tau))
    return (0, 0, nu_int, int(M), int(ctx.side))


def _tau_ns(ctx):
    return max(1, int(round(ctx.det.tau * 1e9)))


# ---------------------------------------------------------------- step 1: split
def _logical_groups(M, meta):
    """Logical measurement groups: one physical pattern row (bern50/gam4) or one
    complementary hadpair pair (atomic — a fold/split never separates it)."""
    pair_indices = meta.get("pair_indices") if isinstance(meta, dict) else None
    if pair_indices is not None and len(pair_indices) > 0:
        pairs = np.asarray(pair_indices, dtype=np.int64)
        if pairs.ndim != 2 or pairs.shape[1] != 2:
            raise ValueError("meta['pair_indices'] must be (P,2)")
        groups = [pairs[g] for g in range(pairs.shape[0])]
        covered = np.zeros(M, dtype=bool)
        covered[pairs.ravel()] = True
        groups += [np.array([i]) for i in np.nonzero(~covered)[0]]
        return groups
    return [np.array([i]) for i in range(M)]


def split_dev_audit(M, meta, cell_key, seed):
    """Frozen outer split: first 80% of hash-permuted logical groups -> DEV,
    rest -> AUDIT. Shared by ALL arms of the cell. AUDIT takes no part in
    lam_max, eta selection, or initializer construction."""
    if isinstance(meta, dict) and (meta.get("continuous")
                                   or meta.get("fold_mode") == "block"):
        raise NotImplementedError(
            "continuous acquisition has no independent AUDIT (unbounded "
            "afterpulse tails); per the F1 ruling such cells inherit eta* from "
            "the corresponding active-start cell with GOF_NA_DEPENDENT.")
    groups = _logical_groups(M, meta)
    G = len(groups)
    rng = rng_for(*cell_key, seed, 63, 4, 0)
    perm = rng.permutation(G)
    n_dev = int(round((1.0 - AUDIT_FRAC) * G))
    n_dev = min(max(n_dev, 1), G - 1) if G >= 2 else G
    dev_mask = np.zeros(M, dtype=bool)
    audit_mask = np.zeros(M, dtype=bool)
    for gi in perm[:n_dev]:
        dev_mask[groups[gi]] = True
    for gi in perm[n_dev:]:
        audit_mask[groups[gi]] = True
    return {"dev_mask": dev_mask, "audit_mask": audit_mask,
            "n_groups": G, "n_dev_groups": int(n_dev),
            "n_audit_groups": int(G - n_dev),
            "underpowered": bool((G - n_dev) < MIN_AUDIT_GROUPS),
            "dev_groups": [groups[gi] for gi in perm[:n_dev]]}


def _dev_folds(split, cell_key, seed, K=K_FOLDS):
    """Grouped K folds INSIDE DEV (frame-index -> fold id; -1 outside DEV)."""
    dev_groups = split["dev_groups"]
    P = len(dev_groups)
    K = max(2, min(K, P))
    rng = rng_for(*cell_key, seed, 63, 4, 1)
    perm = rng.permutation(P)
    M = split["dev_mask"].size
    folds = np.full(M, -1, dtype=np.int64)
    for k, blk in enumerate(np.array_split(perm, K)):
        for gi in blk:
            folds[dev_groups[gi]] = k
    return folds, K


# ------------------------------------------------- lam_max (DEV-only by caller)
def lam_max_arm(f_grad, x0, side, n_iter=N_SEL, tv_null_rel=TV_NULL_REL,
                max_expand=40, n_bisect=26):
    """Smallest lam_TV whose solution is TV-null (TV < tv_null_rel * TV(x0)).
    Anchors the dimensionless eta path. F1 NOTE: the caller must hand this
    DEV-only data — a full-data lam_max leaks AUDIT through the penalty scale.
    Frozen: anchor ||grad f(x0)||, geometric expand (cap 40), 26 log bisections,
    all fits at the selection budget."""
    x0 = np.asarray(x0, dtype=np.float64).ravel()
    tv0 = tv_value(x0, side)
    _, g0 = f_grad(x0)
    lam0 = float(np.linalg.norm(g0)) + 1e-30
    if tv0 <= 0:
        return lam0
    thresh = tv_null_rel * tv0

    def is_null(lam):
        x, _ = tv_fista(f_grad, x0, float(lam), n_iter=n_iter, side=side)
        return tv_value(x, side) < thresh

    hi = lam0
    if is_null(hi):
        lo = hi
        for _ in range(max_expand):
            lo *= 0.5
            if not is_null(lo):
                break
        else:
            return hi
    else:
        lo = hi
        for _ in range(max_expand):
            hi *= 2.0
            if is_null(hi):
                break
        else:
            return hi
    for _ in range(n_bisect):
        mid = float(np.sqrt(lo * hi))
        if is_null(mid):
            hi = mid
        else:
            lo = mid
    return hi


# ------------------------------------------------------------- shared plumbing
def _make_factory(arm_name, ctx, x_ref):
    if arm_name == "QMLE":
        return _qmle_irls_factory(ctx, x_ref)
    return _arm_factory(arm_name, ctx)


def _std_residuals(b, lam_hat, ctx):
    """Standardized calibrated non-paralyzable renewal residuals."""
    mu, v = qmle_mean_var(lam_hat, ctx.T, ctx.det.tau, ctx.sigma_b)
    return (np.asarray(b, dtype=np.float64) - mu) / np.sqrt(v)


def _corr(r, rho):
    r = np.asarray(r)
    rho = np.asarray(rho)
    if r.size < 2 or r.std() <= 0 or rho.std() <= 0:
        return float("nan")
    return float(np.corrcoef(r, rho)[0, 1])


def _predict_lam(A_sub, x_hat, ctx):
    lam_floor = LAM_FLOOR_REL * ctx.Phi
    return np.maximum(ctx.Phi * (A_sub @ x_hat) + ctx.det.dark, lam_floor)


# ------------------------------------------------------- step 2: eta* per arm
def select_eta_arm(arm_name, A, b, ctx, cell_key, seed, split,
                   eta_grid=ETA_GRID):
    """DEV-only eta selection for one arm. Grouped K=5 cross-fit inside DEV,
    per-fold MEAN squared calibrated renewal residual, one-SE rule
    (SE = fold-dispersion heuristic, not a confidence interval). GOF plays no
    part here. Returns the sel dict consumed by audit_cell / the final refit."""
    A = np.asarray(A, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    dev = split["dev_mask"]
    A_dev, b_dev = A[dev], b[dev]
    folds_full, K = _dev_folds(split, cell_key, seed)
    folds = folds_full[dev]                    # fold ids aligned to DEV frames

    x0_dev = ctx.x0
    if x0_dev is None:
        x0_dev = init_gi_flux(A_dev, b_dev, ctx.Phi, ctx.det.dark, ctx.T,
                              ctx.det.tau)
    fg_dev = _make_factory(arm_name, ctx, x0_dev)(A_dev, b_dev)
    lam_max = lam_max_arm(fg_dev, x0_dev, ctx.side, n_iter=N_SEL)

    etas = list(eta_grid)
    d_k = np.full((len(etas), K), np.nan)      # per-fold mean squared residual
    for i, eta in enumerate(etas):
        lam_tv = float(eta) * lam_max
        for k in range(K):
            val = folds == k
            if not val.any():
                d_k[i, k] = np.nan
                continue
            tr = ~val
            A_tr, b_tr = A_dev[tr], b_dev[tr]
            x0k = init_gi_flux(A_tr, b_tr, ctx.Phi, ctx.det.dark, ctx.T,
                               ctx.det.tau)
            fg = _make_factory(arm_name, ctx, x0k)(A_tr, b_tr)
            x_hat, _ = tv_fista(fg, x0k, lam_tv, n_iter=N_SEL, side=ctx.side)
            r = _std_residuals(b_dev[val], _predict_lam(A_dev[val], x_hat, ctx),
                               ctx)
            d_k[i, k] = float(np.mean(r ** 2))

    d_bar = np.nanmean(d_k, axis=1)
    i_min = int(np.argmin(d_bar))              # first occurrence = smallest eta
    dk_min = d_k[i_min][np.isfinite(d_k[i_min])]
    SE_min = (float(np.std(dk_min, ddof=1) / np.sqrt(dk_min.size))
              if dk_min.size >= 2 else 0.0)
    thresh = d_bar[i_min] + SE_min
    admissible = [i for i in range(len(etas)) if d_bar[i] <= thresh]
    i_star = max(admissible)                   # max regularization one-SE pick
    return {
        "arm": arm_name, "eta_grid": [float(e) for e in etas],
        "eta_star": float(etas[i_star]), "eta_min": float(etas[i_min]),
        "i_star": i_star, "i_min": i_min,
        "lam_max_arm": float(lam_max),
        "lam_tv": float(etas[i_star]) * float(lam_max),
        "lam_plugin": float(etas[i_min]) * float(lam_max),
        "d_bar_curve": [float(x) for x in d_bar],
        "SE_min": SE_min, "one_se_threshold": float(thresh),
        "K": int(K), "x0_dev": x0_dev,
    }


# ---------------------------------------------------- steps 3-7: cell adequacy
def audit_cell(A, b, ctx, cell_key, seed, split, sel_rql):
    """MODEL_FAIL_PREDICTIVE audit — once per cell, carried by RQL.
    Coherent plugin scene at eta_min (never eta*: one-SE's deliberate smoothing
    bias must not trigger the flag) -> predict untouched AUDIT -> B=39 coherent
    refit bootstrap -> exact MC p-value. Flag-only: nothing downstream changes."""
    A = np.asarray(A, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    dev, aud = split["dev_mask"], split["audit_mask"]
    if split["underpowered"]:
        return {"GOF_STATUS": "GOF_UNDERPOWERED", "MODEL_FAIL_PREDICTIVE": None,
                "n_audit_groups": split["n_audit_groups"]}

    A_dev, b_dev, A_aud, b_aud = A[dev], b[dev], A[aud], b[aud]
    lam_plugin = sel_rql["lam_plugin"]
    x0_dev = sel_rql["x0_dev"]
    fg_dev = _arm_factory("RQL", ctx)(A_dev, b_dev)
    x_plugin, _ = tv_fista(fg_dev, x0_dev, lam_plugin, n_iter=N_AUDIT,
                           side=ctx.side)

    lam_aud = _predict_lam(A_aud, x_plugin, ctx)
    r_obs = _std_residuals(b_aud, lam_aud, ctx)
    D_obs = float(np.sum(r_obs ** 2))
    mr_obs = float(np.mean(r_obs))
    c_obs = _corr(r_obs, lam_aud * ctx.det.tau)

    # coherent generator scene over ALL patterns (dark added inside simulate)
    u_scene = np.maximum(A @ x_plugin, 0.0)
    tau_ns = _tau_ns(ctx)

    # steps 4-5: refit-per-replicate bootstrap, early stop, exact MC p-value
    rng_boot = rng_for(*cell_key, seed, 63, 4, 2, tau_ns)
    D_star, mr_star, c_star = [], [], []
    n_exceed = 0
    for bi in range(B_AUDIT):
        b_star, _ = simulate_counts(u_scene, ctx.Phi, ctx.T, ctx.det, rng_boot,
                                    sigma_b=ctx.sigma_b)
        b_star = np.asarray(b_star, dtype=np.float64)
        x0s = init_gi_flux(A_dev, b_star[dev], ctx.Phi, ctx.det.dark, ctx.T,
                           ctx.det.tau)
        fgs = _arm_factory("RQL", ctx)(A_dev, b_star[dev])
        x_s, _ = tv_fista(fgs, x0s, lam_plugin, n_iter=N_AUDIT, side=ctx.side)
        lam_s = _predict_lam(A_aud, x_s, ctx)
        r_s = _std_residuals(b_star[aud], lam_s, ctx)
        Ds = float(np.sum(r_s ** 2))
        D_star.append(Ds)
        mr_star.append(float(np.mean(r_s)))
        c_star.append(_corr(r_s, lam_s * ctx.det.tau))
        if Ds >= D_obs:
            n_exceed += 1
        # early stop: p = (1+n_exceed)/(B+1) <= 0.025 requires n_exceed == 0
        # over all 39, so ANY exceedance decides PASS immediately; and after the
        # first 19 all-below we must finish the block (binary-equivalent).
        if n_exceed > 0 and (bi + 1) >= 1:
            break
    B_used = len(D_star)
    model_fail = (n_exceed == 0)               # <=> p = 1/40 <= 0.025
    p_value = (1 + n_exceed) / (B_AUDIT + 1.0) if B_used == B_AUDIT \
        else (1 + n_exceed) / (B_used + 1.0)   # early-stopped: report bound
    # NOTE: on early stop the decision is exact; the numeric p is then only an
    # upper-bound style summary (>= 2/40), recorded as such.

    # step 6: fixed-lam_hat lower-tail leakage diagnostic (no refits)
    rng_leak = rng_for(*cell_key, seed, 63, 4, 3, tau_ns)
    mu_a, v_a = qmle_mean_var(lam_aud, ctx.T, ctx.det.tau, ctx.sigma_b)
    inv_sd = 1.0 / np.sqrt(v_a)
    u_aud = np.maximum((lam_aud - ctx.det.dark) / ctx.Phi, 0.0)
    n_low = 0
    for _ in range(B_LEAK):
        b0, _ = simulate_counts(u_aud, ctx.Phi, ctx.T, ctx.det, rng_leak,
                                sigma_b=ctx.sigma_b)
        D0 = float(np.sum(((np.asarray(b0, dtype=np.float64) - mu_a)
                           * inv_sd) ** 2))
        if D0 <= D_obs:
            n_low += 1
    p_low = (1 + n_low) / (B_LEAK + 1.0)

    # step 7: residual-structure warnings (never in the MODEL_FAIL union)
    mean_warn = bool(abs(mr_obs) > max(abs(m) for m in mr_star))
    c_finite = [abs(c) for c in c_star if np.isfinite(c)]
    corr_warn = (bool(abs(c_obs) > max(c_finite))
                 if (np.isfinite(c_obs) and c_finite) else None)

    return {
        "GOF_STATUS": "GOF_OK", "MODEL_FAIL_PREDICTIVE": bool(model_fail),
        "p_value": float(p_value), "B_used": int(B_used),
        "D_obs": D_obs, "D_star_max": float(max(D_star)),
        "D_star_mean": float(np.mean(D_star)),
        "LEAKAGE_SUSPECT": bool(p_low <= P_LEAK), "p_low": float(p_low),
        "MEAN_RESIDUAL_WARN": mean_warn, "LOAD_CORR_WARN": corr_warn,
        "mean_r_obs": mr_obs, "corr_obs": c_obs,
        "n_audit_frames": int(aud.sum()),
        "n_audit_groups": split["n_audit_groups"],
        "lam_plugin": float(lam_plugin),
    }


# --------------------------------------------------------------- orchestrator
def select_eta_and_fit(arm_name, A, b, ctx, cell_key=None, seed=0,
                       run_audit=None, eta_grid=ETA_GRID, split=None):
    """Full F1 rule for one arm: DEV/AUDIT split -> DEV eta* selection ->
    (RQL only) adequacy audit -> final refit on ALL frames at eta*.

    MODEL_FAIL_PREDICTIVE is an adequacy FLAG: it never changes eta*, never
    suppresses the reconstruction (the old E=empty -> eta_min fallback is
    deliberately gone). Returns (x_final, info)."""
    A = np.asarray(A, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    M = A.shape[0]
    if cell_key is None:
        cell_key = _default_cell_key(ctx, M)
    cell_key = tuple(int(v) for v in cell_key)
    if run_audit is None:
        run_audit = (arm_name == "RQL")
    if getattr(ctx.det, "start_mode", "active") == "continuous":
        raise NotImplementedError(
            "continuous acquisition: inherit eta* from the corresponding "
            "active-start cell (GOF_NA_DEPENDENT) at the campaign layer.")

    if split is None:
        split = split_dev_audit(M, ctx.meta, cell_key, seed)
    sel = select_eta_arm(arm_name, A, b, ctx, cell_key, seed, split,
                         eta_grid=eta_grid)

    audit = None
    if run_audit:
        sel_rql = sel if arm_name == "RQL" else \
            select_eta_arm("RQL", A, b, ctx, cell_key, seed, split,
                           eta_grid=eta_grid)
        audit = audit_cell(A, b, ctx, cell_key, seed, split, sel_rql)

    ctx_final = replace(ctx, lam_tv=sel["lam_tv"], n_iter=int(ctx.n_iter),
                        x0=None)
    x_final, fit_info = run_arm(arm_name, A, b, ctx_final)

    info = {
        "arm": arm_name, "rule": "f1_audit_split_coherent_bootstrap",
        "cell_key": list(cell_key), "seed": int(seed),
        "eta_star": sel["eta_star"], "eta_min": sel["eta_min"],
        "lam_tv": sel["lam_tv"], "lam_max_arm": sel["lam_max_arm"],
        "eta_grid": sel["eta_grid"], "d_bar_curve": sel["d_bar_curve"],
        "SE_min": sel["SE_min"], "one_se_threshold": sel["one_se_threshold"],
        "K": sel["K"], "split": {k: split[k] for k in
                                 ("n_groups", "n_dev_groups", "n_audit_groups",
                                  "underpowered")},
        "audit": audit,
        "MODEL_FAIL": (audit or {}).get("MODEL_FAIL_PREDICTIVE"),
        "fit_info": fit_info,
    }
    return x_final, info


# ------------------------------------------------------------------ smoke test
def _smoke():
    import time
    from physics import Detector
    from solvers import ArmContext

    side, n, M = 32, 32 * 32, 1500
    tau, nu, rho = 1.0, 500.0, 0.6
    T = nu * tau
    img = np.full((side, side), 0.03)
    yy, xx = np.mgrid[0:side, 0:side]
    img[(yy - 9) ** 2 + (xx - 9) ** 2 <= 25] = 1.0
    img[20:28, 4:12] = 0.55
    img[6:9, 16:30] = 1.0
    x_true = img.ravel() / img.sum()

    rng_p = rng_for(0, 63, 1, 1)
    A = 2.0 * (rng_p.random((M, n)) < 0.5)
    u = A @ x_true
    det = Detector(tau=tau, dark=0.0, start_mode="active")
    Phi = rho / (tau * float(u.mean()))
    b, N = simulate_counts(u, Phi, T, det, rng_for(0, 63, 3, 1))
    print("[select_eta F1 smoke] side=%d M=%d rho_emp=%.4f mean_counts=%.1f"
          % (side, M, float(np.mean(Phi * u) * tau), float(N.mean())),
          flush=True)

    ctx = ArmContext(Phi=Phi, det=det, T=T, side=side, n_iter=150,
                     select_iter=N_SEL, pattern_kind="bern50",
                     meta={"kind": "bern50"})
    cell_key = (1, 600, 500, M, side)
    ok = True
    for arm in ("RQL", "POISSON-LIN"):
        t0 = time.time()
        x, info = select_eta_and_fit(arm, A, b, ctx, cell_key=cell_key, seed=0)
        dt = time.time() - t0
        finite = bool(np.all(np.isfinite(x)))
        ok = ok and finite and np.isfinite(info["lam_tv"])
        au = info["audit"]
        print("  %-11s eta*=%.3g lam_tv=%.3g audit=%s (%.1fs)"
              % (arm, info["eta_star"], info["lam_tv"],
                 (None if au is None else
                  {k: au.get(k) for k in ("GOF_STATUS",
                                          "MODEL_FAIL_PREDICTIVE", "p_value",
                                          "B_used", "LEAKAGE_SUSPECT",
                                          "MEAN_RESIDUAL_WARN")}), dt),
              flush=True)
    print("[select_eta F1 smoke] %s" % ("PASS" if ok else "FAIL"), flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(_smoke())
