#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
DLGI CONFIRMATORY CAMPAIGN -- the Neyman machinery (R34 sec.1, exact).

eta = (log t_c, log CV).  For every candidate eta_0 on the frozen grid and photon-
SNR stratum, the frozen full-pipeline PROFILE statistic is

    W(Y; eta_0) = 2 { ell_joint(eta_hat, x_hat) - ell_joint(eta_0, x_hat_eta0) }.

The exact-likelihood approximation is the FROZEN log-Whittle/profile construction
from the feasibility probe (dl_bar4_final.whittle_nll_grid + _periodogram, imported
VERBATIM).  ell_joint(eta, x_hat) = -whittle_nll(eta) on the time-ordered log-gain
residual series z(x_hat) with the white floor sv2 profiled; x_hat_eta0 ~= x_hat (the
joint reconstruction) so W over the (log t_c, log CV) grid is the JOINT (2-DOF)
profile LR surface

    W_grid(Y) = 2 ( nll_grid - min nll_grid ).

Calibration replaces the asymptotic chi^2_2 threshold by a SIMULATION-CALIBRATED
critical value: every calibration replicate REGENERATES the gain path + Poisson
noise and RERUNS the entire scene+medium pipeline (joint_dual_ledger) -- holding
the reconstructed scene fixed is forbidden (R34 sec.1).  Then

  1. c_0.95(eta_0) = conservative Monte-Carlo order statistic of W(Y_rep; eta_0=true)
     over the calibration replicates at each declared (t_c,CV,SNR) knot;
  2. a FROZEN monotone upper-envelope interpolation in (log t_c, log CV, log SNR)
     extends c_0.95 to every grid cell (enclosing-box max, hold-flat extrapolation);
  3. Neyman inversion  C_0.95(Y) = { eta_0 : W(Y;eta_0) <= c_0.95(eta_0) };
  4. projections of the joint region give the t_c and CV intervals.

These are simulation-calibrated, model-conditional intervals under the validated
scalar-gain/OU model -- not distribution-free guarantees (R34 wording).

This module provides the machinery + an end-to-end SMOKE unit-test (__main__) on
ONE cell (tc64_cv40, nominal SNR) where the probe's profile likelihood already
worked.  Read-only on the frozen probe; writes only the smoke JSON.  CPU.
"""
import os
import sys
import json
import math
import time
import platform
from datetime import datetime, timezone

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dlgi_common as K
import dl_common as C                                   # frozen estimator (via K's path)
# frozen Whittle construction, imported VERBATIM from the probe:
from dl_bar4_final import (LOG_TC, LOG_CV, LOG_SV2, whittle_nll_grid, _periodogram)

CAMPAIGN = K.CAMPAIGN
LEVEL = 0.95


def _mem_safe_nproc():
    """Cap worker count so 32 scipy-importing processes do not exhaust Windows
    commit charge (the '页面文件太小' DLL-load failure).  ~0.4 GB/worker."""
    ncpu = os.cpu_count() or 4
    try:
        import psutil
        avail = psutil.virtual_memory().available / 1e9
        return max(2, min(ncpu, int((avail - 1.0) / 0.4)))
    except Exception:
        return min(ncpu, 10)

# fine-grid parameter axes (exp of the frozen log grids)
TC_AXIS = np.exp(LOG_TC)                                 # (58,)
CV_AXIS = np.exp(LOG_CV)                                 # (46,)


# ======================================================= W-statistic surface
def residual_series(Y, x_hat, slot):
    """Time-ordered demeaned log-gain residual series z and validity mask, from a
    record and its joint scene estimate (reuses the frozen _z_residuals)."""
    zc, R, valid, mu = C._z_residuals(Y, x_hat)
    order = np.argsort(slot)
    return zc[order], valid[order]


def W_surface(zc, valid):
    """Return (W_grid, nll_grid) on (LOG_TC x LOG_CV); sv2 profiled.  W_grid is the
    joint profile LR statistic 2(nll - nll_min)."""
    Ik, wk = _periodogram(zc, valid)
    nll = whittle_nll_grid(Ik, wk)                       # (58,46), frozen
    return 2.0 * (nll - nll.min()), nll


def grid_index(tc, cv):
    """Nearest (i,j) fine-grid index to a parameter point (for W-at-true)."""
    i = int(np.argmin(np.abs(np.log(TC_AXIS) - np.log(tc))))
    j = int(np.argmin(np.abs(np.log(CV_AXIS) - np.log(cv))))
    return i, j


def W_at_true(W_grid, tc, cv):
    i, j = grid_index(tc, cv)
    return float(W_grid[i, j])


# ======================================================= conservative crit
def wilson_ci(k, n, z=1.96):
    """Wilson 95% interval for a binomial proportion (coverage CI)."""
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def conservative_critical(W_samples, level=LEVEL):
    """Conservative Monte-Carlo order-statistic critical value: the k-th order
    statistic with k = ceil((n+1) * level), clamped to n.  (n+1) makes the
    estimated quantile conservative -- the region errs toward over-coverage.)"""
    w = np.asarray(W_samples, float)
    w = w[np.isfinite(w)]
    n = len(w)
    if n < 20:
        return float("nan"), n
    k = min(int(math.ceil((n + 1) * level)), n)
    return float(np.sort(w)[k - 1]), n


# ======================================================= critical-value surface
class CritSurface:
    """FROZEN monotone upper-envelope interpolation of the MC critical values in
    (u,v,w) = (log t_c, log CV, log SNR).  Knots are the calibration (t_c,CV,SNR)
    cells.  Rule (frozen before confirmatory): clamp the query into the calibration
    knot box (hold-flat extrapolation), locate the enclosing box on each axis, and
    return the MAX critical value over that box's corner knots (an upper envelope,
    hence conservative / coverage >= nominal between knots)."""

    def __init__(self, tc_knots, cv_knots, ph_knots, cvals):
        self.u = np.log(np.asarray(tc_knots, float))
        self.v = np.log(np.asarray(cv_knots, float))
        self.w = np.log(np.asarray(ph_knots, float))
        self.c = np.asarray(cvals, float)                # shape (len u, len v, len w)
        assert self.c.shape == (len(self.u), len(self.v), len(self.w))

    @staticmethod
    def _box(axis, q):
        """indices [lo,hi] of the knot cell enclosing clamped q (hold-flat outside)."""
        q = min(max(q, axis[0]), axis[-1])
        hi = int(np.searchsorted(axis, q, side="left"))
        hi = min(max(hi, 1), len(axis) - 1)
        return hi - 1, hi

    def eval(self, tc, cv, ph):
        iu = self._box(self.u, math.log(tc))
        iv = self._box(self.v, math.log(cv))
        iw = self._box(self.w, math.log(ph))
        block = self.c[iu[0]:iu[1] + 1, iv[0]:iv[1] + 1, iw[0]:iw[1] + 1]
        return float(np.nanmax(block))

    def grid_at(self, ph):
        """c_0.95 evaluated on the full (LOG_TC x LOG_CV) fine grid at one SNR."""
        out = np.empty((len(TC_AXIS), len(CV_AXIS)))
        for i, tc in enumerate(TC_AXIS):
            for j, cv in enumerate(CV_AXIS):
                out[i, j] = self.eval(tc, cv, ph)
        return out

    def to_dict(self):
        return dict(rule=("enclosing-box upper envelope in (log t_c, log CV, log "
                          "SNR); hold-flat extrapolation outside the calibration box"),
                    tc_knots=np.exp(self.u).tolist(), cv_knots=np.exp(self.v).tolist(),
                    ph_knots=np.exp(self.w).tolist(), cvals=self.c.tolist())


def flat_surface(cval, tc_knots=K.TC_PRIMARY, cv_knots=K.CV_PRIMARY,
                 ph_knots=K.PHOTON_PRIMARY):
    """Degenerate single-value surface (used only to exercise the inversion code
    path in the single-cell smoke, where the full 27-knot table does not yet
    exist)."""
    c = np.full((len(tc_knots), len(cv_knots), len(ph_knots)), cval, float)
    return CritSurface(tc_knots, cv_knots, ph_knots, c)


# ======================================================= Neyman inversion
def neyman_region(W_grid, c_grid):
    """C_0.95 = {eta_0 : W(eta_0) <= c_0.95(eta_0)} on the fine grid, projected to
    t_c and CV intervals.  Reports boundedness (region does not touch a grid edge)
    and connectedness (included indices form a contiguous run per projected axis)."""
    inside = W_grid <= c_grid
    if not inside.any():
        return dict(empty=True, tc_lo=float("nan"), tc_hi=float("nan"),
                    cv_lo=float("nan"), cv_hi=float("nan"),
                    bounded=False, connected=False, n_cells=0)
    ii = np.where(inside.any(axis=1))[0]                 # tc indices present
    jj = np.where(inside.any(axis=0))[0]                 # cv indices present
    tc_lo, tc_hi = float(TC_AXIS[ii.min()]), float(TC_AXIS[ii.max()])
    cv_lo, cv_hi = float(CV_AXIS[jj.min()]), float(CV_AXIS[jj.max()])
    connected = bool(np.array_equal(ii, np.arange(ii.min(), ii.max() + 1)) and
                     np.array_equal(jj, np.arange(jj.min(), jj.max() + 1)))
    bounded = bool(ii.min() > 0 and ii.max() < len(TC_AXIS) - 1 and
                   jj.min() > 0 and jj.max() < len(CV_AXIS) - 1)
    return dict(empty=False, tc_lo=tc_lo, tc_hi=tc_hi, cv_lo=cv_lo, cv_hi=cv_hi,
                bounded=bounded, connected=connected, n_cells=int(inside.sum()))


# ======================================================= parallel worker
def _worker(task):
    """One full-pipeline replicate -> W_grid + point estimates (R34: regenerate gain
    + Poisson + rerun joint_dual_ledger; the reconstructed scene is NOT held fixed)."""
    bank, ckey, scene_id, x, rep, tc, sig, ph = task
    grng, prng, sched = K.record_rngs(bank, ckey, scene_id, rep)
    Y, a_time, slot = K.simulate_record_scaled(x, grng, prng, tc, sig, ph,
                                               schedule="paired", sched_seed=sched)
    r = C.joint_dual_ledger(Y, slot=slot, n_outer=2)     # full pipeline, scene re-fit
    zc, valid = residual_series(Y, r["x_hat"], slot)
    W_grid, _ = W_surface(zc, valid)
    return dict(W_grid=W_grid.astype(np.float32),
                tc_hat=float(r["med"]["tc_hat"]), cv_hat=float(r["med"]["cv_hat"]))


def run_records(bank, ckey, scenes, reps, tc, sig, ph, nproc=None):
    """Full-pipeline replicates for one cell across `scenes` (dict id->x), `reps`
    seeds each.  Returns list of worker dicts.  Multiprocess."""
    tasks = [(bank, ckey, sid, scenes[sid], rep, tc, sig, ph)
             for sid in scenes for rep in reps]
    if nproc is None:
        nproc = _mem_safe_nproc()
    from concurrent.futures import ProcessPoolExecutor
    out = []
    with ProcessPoolExecutor(max_workers=nproc) as ex:
        for res in ex.map(_worker, tasks, chunksize=4):
            out.append(res)
    return out


# ======================================================= SMOKE (unit test)
def smoke(n_cal=180, n_test=180, tc=64.0, cv=0.40, ph=1.0, nproc=None):
    """End-to-end machinery test on ONE cell.  Held-out split: n_cal calibration
    replicates set c_0.95 at the true knot; n_test DISJOINT replicates measure the
    calibrated coverage.  Verifies the pointwise JOINT coverage lands near nominal
    (0.95) and exercises the inversion + projection code path."""
    t0 = time.time()
    sig = C.sigma_l_of_cv(cv)
    ckey = K.cell_key(tc, cv, ph)
    cal_scenes, _ = K.load_bank("calibration")
    i_true, j_true = grid_index(tc, cv)
    # DISCIPLINE: the smoke keeps the confirmatory bank SEALED.  Both splits draw
    # from the CALIBRATION bank with DISJOINT replicate-index ranges (independent
    # SeedSequence streams) -- a genuine held-out coverage check on independent
    # records, with no confirmatory scene opened before the calibration table exists.
    TEST_OFFSET = 1_000_000

    def alloc(n, nsc):
        base = n // nsc
        return [base + (1 if k < n % nsc else 0) for k in range(nsc)]

    cal_ids = list(cal_scenes)
    cal_reps = {sid: list(range(a)) for sid, a in zip(cal_ids, alloc(n_cal, len(cal_ids)))}
    conf_reps = {sid: list(range(TEST_OFFSET, TEST_OFFSET + a))
                 for sid, a in zip(cal_ids, alloc(n_test, len(cal_ids)))}
    conf_scenes = cal_scenes

    print(f"[smoke] cell {ckey}: {n_cal} calibration + {n_test} held-out test "
          f"(disjoint seeds, calibration bank; confirmatory stays sealed) "
          f"full-pipeline replicates (~1.26 s each)", flush=True)

    # -- calibration: W(Y;eta_0=true) null distribution -> conservative c_0.95 --
    cal = _run_split("calibration", ckey, cal_scenes, cal_reps, tc, sig, ph, nproc)
    Wcal = np.array([W_at_true(d["W_grid"].astype(np.float64), tc, cv) for d in cal])
    c95, n_used = conservative_critical(Wcal, LEVEL)
    chi2_2 = 5.991                                       # asymptotic 2-DOF threshold
    print(f"[smoke] c_0.95(true) = {c95:.3f} from {n_used} reps "
          f"(asymptotic chi2_2 = {chi2_2:.3f}; MC/asymptotic = {c95/chi2_2:.2f})",
          flush=True)

    # -- held-out coverage: pointwise joint + projected intervals (flat surface) --
    test = _run_split("calibration", ckey, conf_scenes, conf_reps, tc, sig, ph, nproc)
    Wtest = np.array([W_at_true(d["W_grid"].astype(np.float64), tc, cv) for d in test])
    n_hit = int(np.sum(Wtest <= c95))
    joint_cov = float(n_hit / len(Wtest))
    cov_lo, cov_hi = wilson_ci(n_hit, len(Wtest))
    asym_cov = float(np.mean(Wtest <= chi2_2))

    surf = flat_surface(c95)
    c_grid = surf.grid_at(ph)
    tc_hit = cv_hit = bounded = connected = 0
    widths_tc = []; widths_cv = []
    for d in test:
        reg = neyman_region(d["W_grid"].astype(np.float64), c_grid)
        if reg["empty"]:
            continue
        tc_hit += int(reg["tc_lo"] <= tc <= reg["tc_hi"])
        cv_hit += int(reg["cv_lo"] <= cv <= reg["cv_hi"])
        bounded += int(reg["bounded"]); connected += int(reg["connected"])
        if np.isfinite(reg["tc_hi"]) and reg["tc_lo"] > 0:
            widths_tc.append(math.log(reg["tc_hi"]) - math.log(reg["tc_lo"]))
        if np.isfinite(reg["cv_hi"]) and reg["cv_lo"] > 0:
            widths_cv.append(math.log(reg["cv_hi"]) - math.log(reg["cv_lo"]))
    m = len(test)
    result = dict(
        cell=ckey, tc=tc, cv=cv, ph=ph, n_cal=n_used, n_test=m,
        c_0_95_true=c95, asymptotic_chi2_2=chi2_2, mc_over_asymptotic=c95 / chi2_2,
        pointwise_joint_coverage=joint_cov,
        pointwise_joint_coverage_wilson95=[cov_lo, cov_hi],
        asymptotic_joint_coverage=asym_cov,
        projected_tc_interval_coverage_flatsurf=tc_hit / m,
        projected_cv_interval_coverage_flatsurf=cv_hit / m,
        frac_bounded=bounded / m, frac_connected=connected / m,
        median_logwidth_tc=float(np.median(widths_tc)) if widths_tc else float("nan"),
        median_logwidth_cv=float(np.median(widths_cv)) if widths_cv else float("nan"),
        pf_tc_hat_median=float(np.median([d["tc_hat"] for d in test])),
        pf_cv_hat_median=float(np.median([d["cv_hat"] for d in test])),
        runtime_s=time.time() - t0,
        note=("pointwise_joint_coverage is the real calibration check (held-out: "
              "independent replicate seeds within the calibration bank; confirmatory "
              "bank stays SEALED); projected-interval coverage uses a FLAT single-knot "
              "surface only to exercise the inversion path -- the true projected "
              "coverage needs the full 27-knot calibration table (separate GO)."))
    # PASS if the coverage Wilson CI overlaps the nominal band [0.92,0.98] AND the
    # MC calibration is far better than the asymptotic threshold on this cell.
    overlaps_band = (cov_hi >= 0.92) and (cov_lo <= 0.98)
    beats_asym = joint_cov - asym_cov > 0.10
    verdict = "PASS" if (overlaps_band and beats_asym) else "CHECK"
    result["verdict"] = verdict
    result["verdict_basis"] = ("Wilson95 overlaps [0.92,0.98] AND MC-calibrated "
                               "coverage exceeds asymptotic-chi2 coverage by >0.10")
    meta = dict(utc=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ"),
                python=platform.python_version(), numpy=np.__version__,
                whittle_source="dl_bar4_final (frozen probe, verbatim import)",
                grids=dict(n_log_tc=len(LOG_TC), n_log_cv=len(LOG_CV),
                           n_log_sv2=len(LOG_SV2)))
    blob = dict(meta=meta, smoke=result)
    os.makedirs(os.path.join(CAMPAIGN, "smoke"), exist_ok=True)
    fn = os.path.join(CAMPAIGN, "smoke", "smoke_tc64_cv40.json")
    json.dump(blob, open(fn, "w"), indent=2)
    print(f"[smoke] pointwise JOINT coverage (MC-calibrated)   = {joint_cov:.3f}  "
          f"[Wilson95 {cov_lo:.3f}-{cov_hi:.3f}]  (target ~0.95)", flush=True)
    print(f"[smoke] pointwise JOINT coverage (asymptotic chi2) = {asym_cov:.3f}",
          flush=True)
    print(f"[smoke] projected t_c interval coverage (flat surf)= "
          f"{tc_hit/m:.3f}  cv = {cv_hit/m:.3f}", flush=True)
    print(f"[smoke] bounded={bounded/m:.2f} connected={connected/m:.2f}  "
          f"median log-width t_c={result['median_logwidth_tc']:.3f} "
          f"cv={result['median_logwidth_cv']:.3f}", flush=True)
    print(f"[smoke] wrote {fn}  ({result['runtime_s']:.1f}s)", flush=True)
    print(f"[smoke] VERDICT: {verdict}  (Wilson95 overlaps [0.92,0.98]: {overlaps_band}; "
          f"MC beats asymptotic by {joint_cov-asym_cov:+.2f})", flush=True)
    return result


def _run_split(bank, ckey, scenes, reps, tc, sig, ph, nproc):
    tasks = [(bank, ckey, sid, scenes[sid], rep, tc, sig, ph)
             for sid in scenes for rep in reps.get(sid, [])]
    if nproc is None:
        nproc = _mem_safe_nproc()
    from concurrent.futures import ProcessPoolExecutor
    out = []
    with ProcessPoolExecutor(max_workers=nproc) as ex:
        for res in ex.map(_worker, tasks, chunksize=4):
            out.append(res)
    return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--ncal", type=int, default=180)
    ap.add_argument("--ntest", type=int, default=180)
    ap.add_argument("--nproc", type=int, default=None)
    a = ap.parse_args()
    smoke(n_cal=a.ncal, n_test=a.ntest, nproc=a.nproc)
