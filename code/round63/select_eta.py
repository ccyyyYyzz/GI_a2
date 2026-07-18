"""ROUND63 lambda_TV rule + descriptive measurement audit — F1-FROZEN
(`analytic_score_concentration`, GPT round-6 ruling; audit layer = round-5).

Rule history (each superseded layer probed to death before replacement —
docs/ROUND63_GPT_ROUND{4,5,6}_*.md):
  round-3 six-step discrepancy      -> probe: 10/10 null false alarms (gate)
  round-4 coherent-bootstrap gate   -> probe: 7/20 size + 0/10 power
  round-5 descriptive audit (kept)  -> no binary gate, flag-free diagnostics
  round-6 selection replacement     -> S1: one-SE discrepancy collapses to
    eta*=1 (TV-null) everywhere; mechanism = bucket-SPI multiplex
    disadvantage (per-frame pattern contrast sd(u) = sqrt(sum x^2) ~ 1.6%
    << shot noise, so held-out renewal fit is FLAT in eta at M <= n).
    Truth-free per-image search is structurally powerless here — replaced by
    an ANALYTIC noise-scaled TV weight with one development-calibrated
    concentration threshold.

THE FROZEN RULE (round-6 digest §1, all constants final):
  lambda_TV,a = c_i * sigma_{g,a} * sqrt(2 ln n)
  sigma_{g,a} = Phi * sqrt(kappa_A * v_{s,a} / M)
    kappa_A   = max_j (1/M) sum_i a_ij^2          (actual pattern matrix)
    v_{s,a}   = Var_{N ~ exact NP renewal(lam_bar, T, tau_hat)}[s_a(N)]
                (physics.score_variance — exact PMF enumeration; CLT moments
                are FORBIDDEN, they break at nu = 5, 10)
    lam_bar   = Phi + dark_hat                    (E[u] = 1 analytic)
  c_i two-bin: 0.50 if C_hat_i <= C0 else 0.25   (no continuous maps)
  C_hat = clip[ n (S_N^2 - V0)_+ / ((Phi mu'_0)^2 omega_A), 1, 64 ]
    S_N^2 from DEV raw counts; mu0, V0, mu'_0 exact-PMF moments at lam_bar;
    omega_A = (1/(M_dev n)) sum_{i in DEV} sum_j (a_ij - abar_j)^2
    S_N^2 <= V0 -> C_hat = 1 (smooth bin; no negative extrapolation)
  C0 is THE one development-calibrated constant (Pass A, round-6 §1.6/§5):
  candidates {0,1.5,2,3,4,6,8,12,16,inf}, endpoint-oracle regret objective
  J(C0) = Q0.90 over dev images, tie -> larger C0 within 0.02 dB. Frozen in
  C0_FROZEN.json next to this file after Pass A; before that callers must
  pass C0 explicitly. ALL arms share the same c_i.
  Paper wording: "analytically noise-scaled TV with one development-
  calibrated concentration threshold" — never "training-free".

The DEV/AUDIT 80/20 split SURVIVES: DEV -> the concentration statistic;
AUDIT -> the round-5 descriptive audit (unchanged semantics: continuous
ranks, no adequacy gate, nothing may affect reconstruction or inference).
The old production path (lam_max bisection, K-fold cross-fit, one-SE) is
deleted; git history keeps it.

RNG stream layout (integer keys): outer split rng_for(*cell_key, seed, 63,
4, 0); audit bootstrap (..., 63, 4, 2, tau_ns); leak bootstrap (..., 63, 4,
3, tau_ns). cell_key = (kind_id, rho_milli, nu_int, M, side).

API (campaign.py contract):
  split_dev_audit(M, meta, cell_key, seed) -> split dict
  concentration_stat(A, b, ctx, split) -> dict (C_hat + logged ingredients)
  analytic_lambda(arm_name, A, ctx, C_hat, C0) -> dict (lam_tv + ingredients)
  audit_cell(A, b, ctx, cell_key, seed, split, lam_plugin) -> audit dict
  select_eta_and_fit(arm_name, A, b, ctx, cell_key=None, seed=0,
                     run_audit=None, split=None, C0=None) -> (x, info)
"""
import json
import os
import sys
from dataclasses import replace

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))

from gi_core.utils import rng_for
from physics import (qmle_mean_var, simulate_counts, exact_moments,
                     score_variance)
from solvers import tv_fista, init_gi_flux, run_arm, _arm_factory

# ---------------------------------------------------------------- F1 constants
AUDIT_FRAC = 0.20
MIN_AUDIT_GROUPS = 128
C_BIN_LOW = 0.50            # c for C_hat <= C0  (smooth / diffuse scenes)
C_BIN_HIGH = 0.25           # c for C_hat >  C0  (concentrated scenes)
C_HAT_CLIP = (1.0, 64.0)
C0_CANDIDATES = (0.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0, float("inf"))
C0_FILE = os.path.join(HERE, "C0_FROZEN.json")
N_AUDIT = 25                # plugin + per-replicate audit fits (round-5)
B_DIAG = 39                 # ALWAYS all 39 (round-5: no early stop)
B_LEAK = 199
P_LEAK = 0.01
LAM_FLOOR_REL = 1e-6


def _default_cell_key(ctx, M):
    nu_int = int(round(ctx.T / ctx.det.tau))
    return (0, 0, nu_int, int(M), int(ctx.side))


def _tau_ns(ctx):
    return max(1, int(round(ctx.det.tau * 1e9)))


def frozen_C0():
    """The Pass-A calibrated concentration threshold; None before Pass A."""
    if os.path.exists(C0_FILE):
        with open(C0_FILE) as f:
            return float(json.load(f)["C0"])
    return None


# ---------------------------------------------------------------- outer split
def _logical_groups(M, meta):
    """One physical pattern row (bern50/gam4) or one complementary hadpair
    pair (atomic) per logical group."""
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
    """Frozen 80/20 logical-group split (hash-permuted). DEV feeds only the
    concentration statistic; AUDIT only the descriptive audit; neither runs
    any reconstruction hyperparameter search (round-6 §1.1)."""
    if isinstance(meta, dict) and (meta.get("continuous")
                                   or meta.get("fold_mode") == "block"):
        raise NotImplementedError(
            "continuous acquisition has no independent AUDIT; such cells "
            "inherit lambda from the corresponding active-start cell with "
            "AUDIT_STATUS=NA_DEPENDENT (campaign layer).")
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
            "underpowered": bool((G - n_dev) < MIN_AUDIT_GROUPS)}


# ------------------------------------------- the concentration statistic (DEV)
def concentration_stat(A, b, ctx, split):
    """C_hat = clip[ n (S_N^2 - V0)_+ / ((Phi mu'_0)^2 omega_A), 1, 64 ]
    from DEV raw counts + exact-PMF moments at lam_bar = Phi + dark (round-6
    §1.4). Deployment-legal: no truth, no reconstruction."""
    A = np.asarray(A, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    dev = split["dev_mask"]
    b_dev = b[dev]
    A_dev = A[dev]
    M_dev = int(dev.sum())
    n = A.shape[1]
    S_N2 = float(np.var(b_dev, ddof=1)) if M_dev >= 2 else 0.0
    lam_bar = ctx.Phi + ctx.det.dark
    mu0, V0, mup0 = exact_moments(lam_bar, ctx.T, ctx.det.tau)
    abar = A_dev.mean(axis=0)
    omega_A = float(np.mean((A_dev - abar) ** 2))
    denom = (ctx.Phi * mup0) ** 2 * omega_A
    if denom <= 0 or S_N2 <= V0:
        C_hat = 1.0
    else:
        C_hat = float(np.clip(n * (S_N2 - V0) / denom,
                              C_HAT_CLIP[0], C_HAT_CLIP[1]))
    return {"C_hat": C_hat, "S_N2": S_N2, "V0_exact": V0, "mu0_exact": mu0,
            "muprime0_exact": mup0, "omega_A": omega_A, "M_dev": M_dev}


# --------------------------------------------------- the analytic lambda rule
def analytic_lambda(arm_name, A, ctx, C_hat, C0):
    """lambda_TV,a = c(C_hat; C0) * Phi * sqrt(kappa_A v_{s,a} / M) *
    sqrt(2 ln n) (round-6 §1.2/§1.5). Same c for every arm."""
    A = np.asarray(A, dtype=np.float64)
    M, n = A.shape
    kappa_A = float(np.max(np.mean(A ** 2, axis=0)))
    lam_bar = ctx.Phi + ctx.det.dark
    v_s = score_variance(arm_name, lam_bar, ctx.T, ctx.det.tau)
    sigma_g = ctx.Phi * np.sqrt(kappa_A * v_s / M)
    c_used = C_BIN_LOW if C_hat <= C0 else C_BIN_HIGH
    return {"c_used": c_used, "sigma_grad_arm": float(sigma_g),
            "kappa_A": kappa_A, "v_s_arm": float(v_s),
            "lambda_tv_arm": float(c_used * sigma_g * np.sqrt(2.0 * np.log(n)))}


# --------------------------------------- descriptive measurement audit (cell)
def _std_residuals(b, lam_hat, ctx):
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


def audit_cell(A, b, ctx, cell_key, seed, split, lam_plugin):
    """Round-5 DESCRIPTIVE audit, unchanged semantics (continuous ranks, no
    adequacy gate, nothing downstream may change). One round-6 adaptation:
    the plugin smoothing level is now the production analytic lambda for RQL
    (the old eta_min * lam_max anchor no longer exists); recorded in the
    output for provenance."""
    A = np.asarray(A, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    dev, aud = split["dev_mask"], split["audit_mask"]
    if split["underpowered"]:
        return {"AUDIT_STATUS": "UNDERPOWERED",
                "n_audit_groups": split["n_audit_groups"]}

    A_dev, b_dev, A_aud, b_aud = A[dev], b[dev], A[aud], b[aud]
    x0_dev = init_gi_flux(A_dev, b_dev, ctx.Phi, ctx.det.dark, ctx.T,
                          ctx.det.tau)
    fg_dev = _arm_factory("RQL", ctx)(A_dev, b_dev)
    x_plugin, _ = tv_fista(fg_dev, x0_dev, lam_plugin, n_iter=N_AUDIT,
                           side=ctx.side)

    lam_aud = _predict_lam(A_aud, x_plugin, ctx)
    r_obs = _std_residuals(b_aud, lam_aud, ctx)
    D_obs = float(np.sum(r_obs ** 2))
    mr_obs = float(np.mean(r_obs))
    c_obs = _corr(r_obs, lam_aud * ctx.det.tau)

    u_scene = np.maximum(A @ x_plugin, 0.0)
    tau_ns = _tau_ns(ctx)

    rng_boot = rng_for(*cell_key, seed, 63, 4, 2, tau_ns)
    D_star = np.empty(B_DIAG)
    mr_star = np.empty(B_DIAG)
    c_star = np.empty(B_DIAG)
    for bi in range(B_DIAG):
        b_star, _ = simulate_counts(u_scene, ctx.Phi, ctx.T, ctx.det, rng_boot,
                                    sigma_b=ctx.sigma_b)
        b_star = np.asarray(b_star, dtype=np.float64)
        x0s = init_gi_flux(A_dev, b_star[dev], ctx.Phi, ctx.det.dark, ctx.T,
                           ctx.det.tau)
        fgs = _arm_factory("RQL", ctx)(A_dev, b_star[dev])
        x_s, _ = tv_fista(fgs, x0s, lam_plugin, n_iter=N_AUDIT, side=ctx.side)
        lam_s = _predict_lam(A_aud, x_s, ctx)
        r_s = _std_residuals(b_star[aud], lam_s, ctx)
        D_star[bi] = float(np.sum(r_s ** 2))
        mr_star[bi] = float(np.mean(r_s))
        c_star[bi] = _corr(r_s, lam_s * ctx.det.tau)

    q_D = (1 + int(np.sum(D_star >= D_obs))) / (B_DIAG + 1.0)
    q_mean = (1 + int(np.sum(np.abs(mr_star) >= abs(mr_obs)))) / (B_DIAG + 1.0)
    c_fin = np.abs(c_star[np.isfinite(c_star)])
    q_corr = ((1 + int(np.sum(c_fin >= abs(c_obs)))) / (B_DIAG + 1.0)
              if (np.isfinite(c_obs) and c_fin.size) else None)

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
    q_low = (1 + n_low) / (B_LEAK + 1.0)

    min_rank = 1.0 / (B_DIAG + 1.0)
    return {
        "AUDIT_STATUS": "OK",
        "D_obs": D_obs, "D_star_mean": float(np.mean(D_star)),
        "D_star_sd": float(np.std(D_star, ddof=1)),
        "D_ratio": float(D_obs / np.mean(D_star)),
        "plugin_upper_rank": float(q_D),           # NOT a calibrated p-value
        "mean_r_obs": mr_obs, "q_mean": float(q_mean),
        "MEAN_RESIDUAL_WARN": bool(abs(q_mean - min_rank) < 1e-12),
        "corr_obs": c_obs,
        "q_corr": (float(q_corr) if q_corr is not None else None),
        "LOAD_CORR_WARN": (bool(abs(q_corr - min_rank) < 1e-12)
                           if q_corr is not None else None),
        "q_low": float(q_low), "LEAKAGE_SUSPECT": bool(q_low <= P_LEAK),
        "B_diag": int(B_DIAG),
        "n_audit_frames": int(aud.sum()),
        "n_audit_groups": split["n_audit_groups"],
        "lam_plugin": float(lam_plugin),
    }


# --------------------------------------------------------------- orchestrator
def select_eta_and_fit(arm_name, A, b, ctx, cell_key=None, seed=0,
                       run_audit=None, split=None, C0=None, c_force=None):
    """Analytic lambda + final fit for one arm; RQL additionally carries the
    once-per-cell descriptive audit. c_force overrides the two-bin c (used
    ONLY by the Pass-A calibration to compute the .25/.50 endpoints)."""
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
            "continuous acquisition: inherit lambda from the corresponding "
            "active-start cell (AUDIT_STATUS=NA_DEPENDENT) at the campaign "
            "layer.")
    if split is None:
        split = split_dev_audit(M, ctx.meta, cell_key, seed)

    conc = concentration_stat(A, b, ctx, split)
    if c_force is None:
        if C0 is None:
            C0 = frozen_C0()
        if C0 is None:
            raise RuntimeError(
                "C0 not frozen yet (Pass A pending) and no explicit C0/"
                "c_force given — refusing to guess the concentration "
                "threshold.")
        sel = analytic_lambda(arm_name, A, ctx, conc["C_hat"], C0)
    else:
        sel = analytic_lambda(arm_name, A, ctx, conc["C_hat"],
                              float("inf") if c_force == C_BIN_LOW else 0.0)
        sel["c_used"] = float(c_force)
        sel["lambda_tv_arm"] = float(
            c_force * sel["sigma_grad_arm"]
            * np.sqrt(2.0 * np.log(A.shape[1])))

    audit = None
    if run_audit:
        if arm_name == "RQL":
            lam_plugin = sel["lambda_tv_arm"]
        else:
            rql = analytic_lambda("RQL", A, ctx, conc["C_hat"],
                                  C0 if C0 is not None else float("inf"))
            lam_plugin = rql["lambda_tv_arm"]
        audit = audit_cell(A, b, ctx, cell_key, seed, split, lam_plugin)

    ctx_final = replace(ctx, lam_tv=sel["lambda_tv_arm"],
                        n_iter=int(ctx.n_iter), x0=None)
    x_final, fit_info = run_arm(arm_name, A, b, ctx_final)

    info = {
        "arm": arm_name, "rule": "analytic_score_concentration",
        "cell_key": list(cell_key), "seed": int(seed),
        "lam_tv": sel["lambda_tv_arm"],
        "C0": (float(C0) if C0 is not None else None),
        **conc, **sel,
        "split": {k: split[k] for k in
                  ("n_groups", "n_dev_groups", "n_audit_groups",
                   "underpowered")},
        "audit": audit,
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
    print("[select F1v2 smoke] side=%d M=%d rho_emp=%.4f mean_counts=%.1f"
          % (side, M, float(np.mean(Phi * u) * tau), float(N.mean())),
          flush=True)

    ctx = ArmContext(Phi=Phi, det=det, T=T, side=side, n_iter=200,
                     pattern_kind="bern50", meta={"kind": "bern50"})
    cell_key = (1, 600, 500, M, side)
    ok = True
    for arm in ("RQL", "POISSON-LIN"):
        t0 = time.time()
        x, info = select_eta_and_fit(arm, A, b, ctx, cell_key=cell_key,
                                     seed=0, C0=4.0)
        dt = time.time() - t0
        finite = bool(np.all(np.isfinite(x)))
        ok = ok and finite and np.isfinite(info["lam_tv"])
        au = info["audit"]
        print("  %-11s C_hat=%.2f c=%.2f sigma_g=%.3g lam_tv=%.3g audit=%s"
              " (%.1fs)"
              % (arm, info["C_hat"], info["c_used"], info["sigma_grad_arm"],
                 info["lam_tv"],
                 (None if au is None else
                  {k: au.get(k) for k in ("AUDIT_STATUS", "D_ratio",
                                          "plugin_upper_rank",
                                          "LEAKAGE_SUSPECT")}), dt),
              flush=True)
    print("[select F1v2 smoke] %s" % ("PASS" if ok else "FAIL"), flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(_smoke())
