"""ROUND63 lam_TV selection — cross-fitted common renewal discrepancy rule
(GPT round-3 digest §4, six steps, followed EXACTLY).

This REPLACES the deprecated own-NLL rule (`solvers.select_lam_tv`): instead of
scoring each arm on its OWN data-term NLL (which lets an over-confident model
prefer its own optimistic fit), every arm is scored on ONE COMMON, physically
grounded, deployment-legal evaluator — the standardized renewal residual
    r_i = (N_i - mu_NP(lam_hat_i)) / sqrt(V_NP(lam_hat_i)),
where mu_NP / V_NP are the non-paralyzable renewal moments
(`physics.qmle_mean_var`) evaluated with the CALIBRATED detector params carried
by `ArmContext.det` (tau_est, dark_est) — NOT the truth. The reconstruction
lam_hat is CROSS-FITTED (K folds; each frame is scored by a model that never saw
it), so the discrepancy cannot be gamed by in-sample over-fitting.

The six steps (digest §4):
  1. make_folds        : K=5 frozen hash folds; complementary Hadamard pairs land
                         in the SAME fold; a documented NotImplemented guard for
                         continuous-acquisition (state-carrying) frames.
  2. lam_max_arm       : per-arm TV scale via a frozen bisection on the smallest
                         lam_TV that drives the solution TV-null; the eta path is
                         eta = lam_TV / lam_max in ETA_GRID (dimensionless).
  3. cross-fit         : for each eta, fit the arm on K-1 folds with its OWN data
                         fidelity, lam_TV = eta * lam_max, predict the held fold.
  4. common evaluator  : held-out r_i, D(eta) = sum r_i^2, mean(r), corr(r, rho).
  5. GOF + one-SE      : exact renewal parametric bootstrap at eta_min gives the
                         [2.5%, 97.5%] acceptance band for D and mean(r);
                         E = {eta: D(eta) <= D(eta_min) + SE_min AND GOF pass};
                         eta_star = max(E); E empty -> eta_min with MODEL_FAIL
                         (NEVER look at PSNR to rescue).
  6. refit             : refit on ALL frames at the frozen eta_star.

The theory backing is the Morozov discrepancy principle; SURE/GSURE are appendix
sensitivity only. The truth-oracle lam is a DIAGNOSTIC (ORACLE-NOT-DEPLOYABLE)
and is never used here.

API:  select_eta_and_fit(arm_name, A, b, ctx, K=5, seed=0) -> (x_final, info).
float64 throughout; selection fits use ctx.select_iter, the final refit ctx.n_iter.
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

# dimensionless regularization path eta = lam_TV / lam_max_arm (digest §4 step 2)
ETA_GRID = (1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0)

# RNG sub-streams (disjoint from patterns=1, images=2, buckets=3): folds=4,
# selection bootstrap=5 (SEED0 system, spec §5, utils.rng_for).
FOLD_SUBSTREAM = 4
BOOT_SUBSTREAM = 5

# frozen bisection rule for lam_max_arm: a solution is "TV-null" when its
# isotropic TV has collapsed below this fraction of the initializer's TV.
TV_NULL_REL = 1e-3


# ----------------------------------------------------------------------
# step 1: deterministic hash-based K folds (complementary pairs co-located)
# ----------------------------------------------------------------------
def make_folds(M, K=5, seed=0, meta=None):
    """Deterministic fold assignment for M frames.

    A single frozen permutation of the frame (or pair) indices is drawn from
    rng_for(seed, 63, FOLD_SUBSTREAM, M, K) — the SEED0 hash system — and split
    into K near-equal contiguous blocks; block k -> fold k. Because the whole
    assignment is a pure function of (seed, M, K) it is reproducible and audit-
    stable (no wall-clock, no set iteration order).

    Complementary-pair constraint (hadpair): if `meta` carries "pair_indices"
    (the REAL key emitted by patterns.make_patterns — a (P,2) list of the two
    physical row indices of each signed measurement), the PAIRS are permuted and
    folded, and BOTH physical rows of a pair inherit their pair's fold. This
    keeps a signed measurement whole: a fold never splits a (+,-) pair, so the
    held-out discrepancy is computed on complete measurements.

    Continuous-acquisition guard: for the 'continuous' start mode the detector
    carries dead-time / afterpulse state ACROSS frames, so frames are NOT
    exchangeable and per-frame (or per-pair) hold-out is invalid — a correct fold
    must instead cut CONTIGUOUS time blocks with a discarded guard frame at each
    seam. That path is intentionally left as a documented NotImplemented guard;
    the S2 main arms are per-frame independent (active start), so it is never hit
    in production. Trigger it by passing meta={"continuous": True} (or a
    "fold_mode" == "block").

    Returns an int array `folds` of length M with values in 0..K-1.
    """
    M = int(M)
    K = int(K)
    if meta is not None and (meta.get("continuous")
                             or meta.get("fold_mode") == "block"):
        raise NotImplementedError(
            "continuous-acquisition contiguous-block folding is not implemented "
            "(frames carry detector state across the seam; per-frame hold-out is "
            "invalid). S2 main arms use the 'active' start mode and are per-frame "
            "independent.")
    if M <= 0:
        return np.zeros(0, dtype=np.int64)
    K = max(1, min(K, M))
    rng = rng_for(seed, 63, FOLD_SUBSTREAM, M, K)
    folds = np.empty(M, dtype=np.int64)

    pair_indices = meta.get("pair_indices") if isinstance(meta, dict) else None
    if pair_indices is not None and len(pair_indices) > 0:
        pairs = np.asarray(pair_indices, dtype=np.int64)
        if pairs.ndim != 2 or pairs.shape[1] != 2:
            raise ValueError("meta['pair_indices'] must be (P,2)")
        P = pairs.shape[0]
        Kp = max(1, min(K, P))
        perm = rng.permutation(P)
        fold_of_pair = np.empty(P, dtype=np.int64)
        for k, blk in enumerate(np.array_split(perm, Kp)):
            fold_of_pair[blk] = k
        # both physical rows of a pair share the pair's fold
        folds[:] = -1
        folds[pairs[:, 0]] = fold_of_pair
        folds[pairs[:, 1]] = fold_of_pair
        if (folds < 0).any():
            # any physical rows not covered by a pair (should not happen for a
            # well-formed hadpair meta) get folded individually as a fallback
            loose = np.nonzero(folds < 0)[0]
            for k, blk in enumerate(np.array_split(rng.permutation(loose.size),
                                                   Kp)):
                folds[loose[blk]] = k
        return folds

    perm = rng.permutation(M)
    for k, blk in enumerate(np.array_split(perm, K)):
        folds[blk] = k
    return folds


# ----------------------------------------------------------------------
# step 2: per-arm TV scale via frozen bisection (smallest TV-null lam_TV)
# ----------------------------------------------------------------------
def lam_max_arm(f_grad, x0, side, n_iter=120, tv_null_rel=TV_NULL_REL,
                max_expand=40, n_bisect=26):
    """Smallest lam_TV for which `tv_fista(f_grad, x0, lam_TV)` returns a TV-null
    solution, i.e. isotropic TV(x) < tv_null_rel * TV(x0) (a near-constant image,
    the fully-regularized endpoint). This anchors the dimensionless eta path.

    Frozen rule (audit-stable):
      * anchor lam0 = ||grad f_data(x0)||_2 (the arm's own data-misfit scale);
      * EXPAND: double lam upward (<= max_expand times) until the solution is
        TV-null -> an upper bracket `hi`; if never null, return the last lam;
      * if lam0 was already TV-null, halve (<= max_expand times) to find a
        non-null `lo`;
      * BISECT in log-space (geometric midpoint) n_bisect times, keeping `hi` the
        smallest known TV-null lam; return `hi`.
    Log-space bisection matches the multi-decade lam scale. All fits use the
    reduced `n_iter` (selection budget) — lam_max only needs to be approximate.
    """
    x0 = np.asarray(x0, dtype=np.float64).ravel()
    tv0 = tv_value(x0, side)
    _, g0 = f_grad(x0)
    lam0 = float(np.linalg.norm(g0)) + 1e-30
    if tv0 <= 0:
        # flat initializer: no TV to collapse; fall back to the data-misfit scale
        return lam0

    thresh = tv_null_rel * tv0

    def is_null(lam):
        x, _ = tv_fista(f_grad, x0, float(lam), n_iter=n_iter, side=side)
        return tv_value(x, side) < thresh

    hi = lam0
    if is_null(hi):
        # already null: descend to find a non-null lower bracket
        lo = hi
        for _ in range(max_expand):
            lo *= 0.5
            if not is_null(lo):
                break
        else:
            return hi  # null all the way down — degenerate; hi is a valid bound
    else:
        # ascend to find a null upper bracket
        lo = hi
        for _ in range(max_expand):
            hi *= 2.0
            if is_null(hi):
                break
        else:
            return hi  # never null within the expansion cap; return the cap

    for _ in range(n_bisect):
        mid = float(np.sqrt(lo * hi))          # geometric midpoint (log-space)
        if is_null(mid):
            hi = mid
        else:
            lo = mid
    return hi


# ----------------------------------------------------------------------
# arm data-fidelity factory (own fidelity, per digest §4 step 3)
# ----------------------------------------------------------------------
def _make_factory(arm_name, ctx, x_ref):
    """f_grad_factory(A_sub, b_sub) -> f_grad for `arm_name` using the arm's OWN
    data fidelity (reuses the solvers.py machinery). QMLE uses Wedderburn
    quasi-score weights frozen at lam(x_ref) (truth-free); all other iterative
    arms (RQL / POISSON-LIN / SAT-POISSON / PRECORRECT) use `solvers._arm_factory`.
    """
    if arm_name == "QMLE":
        return _qmle_irls_factory(ctx, x_ref)
    return _arm_factory(arm_name, ctx)


# ----------------------------------------------------------------------
# common physical evaluator (digest §4 step 4)
# ----------------------------------------------------------------------
def _std_residuals(b, lam_hat, ctx):
    """Standardized non-paralyzable renewal residuals r_i, with the CALIBRATED
    detector moments (physics.qmle_mean_var at ctx.det.tau / ctx.sigma_b)."""
    mu, v = qmle_mean_var(lam_hat, ctx.T, ctx.det.tau, ctx.sigma_b)
    return (np.asarray(b, dtype=np.float64) - mu) / np.sqrt(v)


def _corr(r, rho):
    r = np.asarray(r); rho = np.asarray(rho)
    if r.std() <= 0 or rho.std() <= 0:
        return float("nan")
    return float(np.corrcoef(r, rho)[0, 1])


def _cross_fit(arm_name, A, b, ctx, folds, K, lam_max, eta, n_iter):
    """Cross-fitted held-out lam_hat over K folds at lam_TV = eta * lam_max.

    Each fold's model is fit on the OTHER K-1 folds (own data fidelity, own
    truth-free per-fold GI-flux init) and predicts its held-out frames, so every
    lam_hat_i comes from a model that never saw frame i. Returns (lam_hat, D,
    mean_r, corr, Dk, nk).
    """
    A = np.asarray(A, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    M = A.shape[0]
    lam_tv = float(eta) * float(lam_max)
    lam_floor = 1e-6 * ctx.Phi
    lam_hat = np.empty(M, dtype=np.float64)
    for k in range(K):
        val = folds == k
        if not val.any():
            continue
        tr = ~val
        A_tr, b_tr, A_val = A[tr], b[tr], A[val]
        x0 = init_gi_flux(A_tr, b_tr, ctx.Phi, ctx.det.dark, ctx.T, ctx.det.tau)
        fg = _make_factory(arm_name, ctx, x0)(A_tr, b_tr)
        x_hat, _ = tv_fista(fg, x0, lam_tv, n_iter=n_iter, side=ctx.side)
        lam_hat[val] = np.maximum(ctx.Phi * (A_val @ x_hat) + ctx.det.dark,
                                  lam_floor)
    r = _std_residuals(b, lam_hat, ctx)
    rho = lam_hat * ctx.det.tau
    Dk = np.array([float(np.sum(r[folds == k] ** 2)) for k in range(K)])
    nk = np.array([int(np.sum(folds == k)) for k in range(K)])
    return (lam_hat, float(np.sum(r ** 2)), float(np.mean(r)),
            _corr(r, rho), Dk, nk)


# ----------------------------------------------------------------------
# step 5: exact renewal parametric bootstrap acceptance band
# ----------------------------------------------------------------------
def _bootstrap_band(lam_hat, ctx, seed, B=200):
    """Exact renewal parametric bootstrap at the fitted lam_hat. Simulate B exact
    renewal count vectors from lam_hat with the CALIBRATED detector
    (physics.simulate_counts) and score each with the SAME CLT residual as the
    data. Returns ([D_lo,D_hi], [mr_lo,mr_hi], D_boot_mean).

    Because the residual uses the CLT moments while the counts are drawn from the
    exact renewal law, the band captures the finite-count CLT-vs-exact gap — the
    honest null for the discrepancy statistic.
    """
    mu, v = qmle_mean_var(lam_hat, ctx.T, ctx.det.tau, ctx.sigma_b)
    inv_sd = 1.0 / np.sqrt(v)
    u_boot = np.maximum((lam_hat - ctx.det.dark) / ctx.Phi, 0.0)
    rng = rng_for(seed, 63, BOOT_SUBSTREAM, int(round(ctx.T)), lam_hat.size)
    D_boot = np.empty(B, dtype=np.float64)
    mr_boot = np.empty(B, dtype=np.float64)
    for j in range(B):
        b_boot, _ = simulate_counts(u_boot, ctx.Phi, ctx.T, ctx.det, rng,
                                    sigma_b=ctx.sigma_b)
        rb = (b_boot - mu) * inv_sd
        D_boot[j] = float(np.sum(rb ** 2))
        mr_boot[j] = float(np.mean(rb))
    D_band = [float(np.quantile(D_boot, 0.025)),
              float(np.quantile(D_boot, 0.975))]
    mr_band = [float(np.quantile(mr_boot, 0.025)),
               float(np.quantile(mr_boot, 0.975))]
    return D_band, mr_band, float(np.mean(D_boot))


def _bootstrap_band_refit(arm_name, A, ctx, folds, K, lam_max, eta_min,
                          lam_hat_min, seed, B=30, n_iter=60):
    """Refit-per-replicate parametric bootstrap at eta_min (candidate-A GOF,
    round-4 adjudication).

    The fixed-lam_hat band above tests 'lam_hat == lam' (counting noise only);
    but the cross-fitted D additionally carries the held-out estimation error
    sum_i (mu(lam_i)-mu(lam_hat_i))^2 / V_i — Theta(M)-scale on underdetermined
    cells — so the fixed band under-covers and MODEL_FAILs the CORRECT model
    (observed: D=1881 vs band ~1604 at M=1500 on a same-model smoke). Here each
    replicate re-simulates exact-renewal counts from lam_hat(eta_min) AND reruns
    the SAME K-fold cross-fit at eta_min, so the null distribution of D*
    includes estimation variance + smoothing bias by construction. B is small
    (each replicate costs K fits): the gate widens the raw quantiles with a
    moment band (mean +/- 2.5 sd) to de-noise the small-B tail; both recorded.
    Caveat (probe-checked): the generator is the SMOOTHED lam_hat, so replicate
    reconstruction error can understate the true-lam error — coverage is
    verified end-to-end by probe_gof_coverage.py, not assumed.
    """
    u_boot = np.maximum((lam_hat_min - ctx.det.dark) / ctx.Phi, 0.0)
    rng = rng_for(seed, 63, BOOT_SUBSTREAM, 7001, lam_hat_min.size)
    D_b = np.empty(B, dtype=np.float64)
    mr_b = np.empty(B, dtype=np.float64)
    for j in range(B):
        b_j, _ = simulate_counts(u_boot, ctx.Phi, ctx.T, ctx.det, rng,
                                 sigma_b=ctx.sigma_b)
        _, Dj, mrj, _, _, _ = _cross_fit(arm_name, A, b_j, ctx, folds, K,
                                         lam_max, eta_min, n_iter)
        D_b[j] = Dj
        mr_b[j] = mrj

    def _band(v):
        q = [float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975))]
        m, s = float(np.mean(v)), float(np.std(v, ddof=1))
        return {"q": q, "mean": m, "sd": s, "B": int(B),
                "gate": [min(q[0], m - 2.5 * s), max(q[1], m + 2.5 * s)]}

    return _band(D_b), _band(mr_b)


# ----------------------------------------------------------------------
# driver: the full six-step rule
# ----------------------------------------------------------------------
def select_eta_and_fit(arm_name, A, b, ctx, K=5, seed=0, B=200,
                       eta_grid=ETA_GRID, gof_mode="refit", B_refit=30):
    """Choose lam_TV by the cross-fitted common renewal discrepancy rule and
    refit on all frames at the frozen choice (digest §4).

    gof_mode (round-4 adjudication; final rule frozen at F1):
      "fixed" — literal digest §4: gate band from the fixed-lam_hat bootstrap.
                KNOWN DEFECT: band covers counting noise only -> correct-model
                MODEL_FAIL misfire on underdetermined cells (see
                _bootstrap_band_refit docstring + probe_gof_coverage.py).
      "refit" — candidate A: gate band from the refit-per-replicate bootstrap
                at eta_min (B_refit replicates x K fits). The fixed band is
                still computed and recorded; D_obs BELOW its lower edge raises
                the 'overfit_flag' diagnostic (leakage), never MODEL_FAIL.
      "off"   — candidate C: no absolute band; E from the one-SE rule alone
                (selection unaffected); bands recorded as diagnostics only.

    Parameters
    ----------
    arm_name : one iterative arm ('RQL' | 'POISSON-LIN' | 'SAT-POISSON' |
               'QMLE' | 'PRECORRECT').
    A, b     : (M,n) illumination matrix and (M,) buckets (raw counts; b == N at
               sigma_b = 0). For hadpair these are the 2*M physical rows.
    ctx      : solvers.ArmContext with the CALIBRATED detector ctx.det (tau_est,
               dark_est), ctx.Phi, ctx.T, ctx.side, ctx.sigma_b, ctx.select_iter
               (selection budget), ctx.n_iter (final refit budget), and
               ctx.meta (pattern meta; its 'pair_indices' co-locates hadpair
               folds).
    K, seed  : fold count and frozen fold/bootstrap seed.

    Returns
    -------
    (x_final, info) where info carries eta_star, lam_tv, MODEL_FAIL, the full D
    curve, the GOF band + details, the fold spec and the acceptance set E.
    """
    A = np.asarray(A, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    M = A.shape[0]
    side = ctx.side
    if getattr(ctx.det, "start_mode", "active") == "continuous":
        raise NotImplementedError(
            "select_eta_and_fit: continuous acquisition needs contiguous-block "
            "folds (see make_folds guard); not implemented for S2.")

    K = max(1, min(int(K), M))
    folds = make_folds(M, K=K, seed=seed, meta=ctx.meta)
    n_sel = int(ctx.select_iter)

    # step 2: per-arm lam_max on the FULL data (frozen eta scale), full-data init
    x0_full = ctx.x0
    if x0_full is None:
        x0_full = init_gi_flux(A, b, ctx.Phi, ctx.det.dark, ctx.T, ctx.det.tau)
    fg_full = _make_factory(arm_name, ctx, x0_full)(A, b)
    lam_max = lam_max_arm(fg_full, x0_full, side, n_iter=n_sel)

    # step 3-4: cross-fitted discrepancy along the eta path
    etas = list(eta_grid)
    D = np.empty(len(etas)); mean_r = np.empty(len(etas))
    corr = np.empty(len(etas)); lam_hats = []
    Dk_all = []; nk_all = []
    for i, eta in enumerate(etas):
        lh, Di, mri, ci, Dk, nk = _cross_fit(arm_name, A, b, ctx, folds, K,
                                             lam_max, eta, n_sel)
        D[i] = Di; mean_r[i] = mri; corr[i] = ci
        lam_hats.append(lh); Dk_all.append(Dk); nk_all.append(nk)

    # step 5a: one-SE band at the discrepancy minimizer
    i_min = int(np.argmin(D))
    Dk_min = Dk_all[i_min]; nk_min = np.maximum(nk_all[i_min], 1)
    m_k = Dk_min / nk_min                      # per-fold mean squared residual
    if K >= 2:
        SE_min = float(M * np.std(m_k, ddof=1) / np.sqrt(K))
    else:
        SE_min = 0.0

    # step 5b: exact renewal parametric bootstrap GOF band (at eta_min).
    # The fixed-lam_hat band is ALWAYS computed (cheap, B sims, no refits): in
    # "fixed" mode it is the gate; otherwise it is the overfit/leakage floor.
    D_band, mr_band, D_boot_mean = _bootstrap_band(lam_hats[i_min], ctx, seed, B)

    refit_band = None
    if gof_mode == "refit":
        refit_band = _bootstrap_band_refit(
            arm_name, A, ctx, folds, K, lam_max, etas[i_min], lam_hats[i_min],
            seed, B=B_refit, n_iter=n_sel)
        gate_D, gate_mr = refit_band[0]["gate"], refit_band[1]["gate"]
    elif gof_mode == "fixed":
        gate_D, gate_mr = D_band, mr_band
    elif gof_mode == "off":
        gate_D = gate_mr = None
    else:
        raise ValueError("unknown gof_mode %r" % (gof_mode,))

    def gof_pass(i):
        if gate_D is None:
            return True
        return (gate_D[0] <= D[i] <= gate_D[1]
                and gate_mr[0] <= mean_r[i] <= gate_mr[1])

    gof = [bool(gof_pass(i)) for i in range(len(etas))]

    # step 5c: one-SE + max-regularization acceptance set
    thresh = D[i_min] + SE_min
    E = [i for i in range(len(etas)) if D[i] <= thresh and gof[i]]
    if E:
        i_star = max(E)                        # largest eta = most regularized
        model_fail = False
    else:
        i_star = i_min                         # fall back to the D-argmin
        model_fail = True
    eta_star = float(etas[i_star])
    lam_tv_star = eta_star * lam_max

    # step 6: refit on ALL frames at frozen eta_star (production machinery)
    ctx_final = replace(ctx, lam_tv=lam_tv_star, n_iter=int(ctx.n_iter),
                        x0=x0_full)
    x_final, fit_info = run_arm(arm_name, A, b, ctx_final)

    info = {
        "arm": arm_name, "rule": "cross_fitted_common_renewal_discrepancy",
        "eta_star": eta_star, "lam_tv": float(lam_tv_star),
        "lam_max_arm": float(lam_max), "MODEL_FAIL": bool(model_fail),
        "eta_min": float(etas[i_min]), "i_star": i_star, "i_min": i_min,
        "K": K, "B": int(B), "seed": int(seed),
        "eta_grid": [float(e) for e in etas],
        "D_curve": [float(x) for x in D],
        "mean_r_curve": [float(x) for x in mean_r],
        "corr_r_rho_curve": [float(x) for x in corr],
        "SE_min": SE_min, "one_se_threshold": float(thresh),
        "gof_pass": gof, "accept_set_E": [int(i) for i in E],
        "gof_mode": gof_mode,
        # D at the CHOSEN eta below the counting-noise floor => held-out
        # residuals are too good to be true (leakage/overfit) — diagnostic only
        "overfit_flag": bool(D[i_star] < D_band[0]),
        "gof": {"D_band": D_band, "mean_r_band": mr_band,
                "D_boot_mean": D_boot_mean,
                "refit": refit_band,
                "D_at_eta_min": float(D[i_min]),
                "mean_r_at_eta_min": float(mean_r[i_min])},
        "fold_spec": {"K": K, "seed": seed, "substream": FOLD_SUBSTREAM,
                      "fold_sizes": [int(np.sum(folds == k)) for k in range(K)],
                      "hadpair_pairs_colocated": bool(
                          isinstance(ctx.meta, dict)
                          and ctx.meta.get("pair_indices") is not None)},
        "fit_info": fit_info,
    }
    return x_final, info


# ----------------------------------------------------------------------
# smoke: 32x32, M=1500, rho=0.6, arms RQL and POISSON-LIN
# ----------------------------------------------------------------------
def _smoke():
    import time
    from physics import Detector
    from solvers import ArmContext

    side, n, M = 32, 32 * 32, 1500
    tau, nu, rho = 1.0, 500.0, 0.6
    T = nu * tau
    # high-contrast phantom (dead-time compression varies frame-to-frame)
    img = np.full((side, side), 0.03)
    yy, xx = np.mgrid[0:side, 0:side]
    img[(yy - 9) ** 2 + (xx - 9) ** 2 <= 25] = 1.0
    img[20:28, 4:12] = 0.55
    img[6:9, 16:30] = 1.0
    x_true = (img.ravel() / img.sum())

    rng_p = rng_for(0, 63, 1, 1)
    A = 2.0 * (rng_p.random((M, n)) < 0.5)
    u = A @ x_true
    det = Detector(tau=tau, dark=0.0, start_mode="active")
    Phi = rho / (tau * float(u.mean()))
    b, N = simulate_counts(u, Phi, T, det, rng_for(0, 63, 3, 1))
    print("[select_eta smoke] side=%d M=%d rho_emp=%.4f mean_counts=%.1f"
          % (side, M, float(np.mean(Phi * u) * tau), float(N.mean())), flush=True)

    ctx = ArmContext(Phi=Phi, det=det, T=T, side=side, n_iter=150,
                     select_iter=40, pattern_kind="bern50",
                     meta={"kind": "bern50"})
    ok = True
    for arm in ("RQL", "POISSON-LIN"):
        t0 = time.time()
        x, info = select_eta_and_fit(arm, A, b, ctx, K=5, seed=0)
        dt = time.time() - t0
        xp = np.maximum(x, 0.0)
        s = xp.sum()
        xs = xp * (x_true.sum() / s) if s > 0 else xp
        mse = float(np.mean((xs - x_true) ** 2))
        psnr = np.inf if mse <= 0 else 10.0 * np.log10(float(x_true.max()) ** 2 / mse)
        finite = bool(np.all(np.isfinite(x)))
        ok = ok and finite and np.isfinite(info["lam_tv"])
        print("  %-11s eta_star=%.3g lam_tv=%.3g lam_max=%.3g MODEL_FAIL=%s "
              "|E|=%d D_min=%.1f D_band=[%.1f,%.1f] PSNR=%.2f finite=%s (%.1fs)"
              % (arm, info["eta_star"], info["lam_tv"], info["lam_max_arm"],
                 info["MODEL_FAIL"], len(info["accept_set_E"]),
                 info["gof"]["D_at_eta_min"], info["gof"]["D_band"][0],
                 info["gof"]["D_band"][1], psnr, finite, dt), flush=True)
    print("[select_eta smoke] %s" % ("PASS" if ok else "FAIL"), flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(_smoke())
