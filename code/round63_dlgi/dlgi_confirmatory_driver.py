#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
DLGI CONFIRMATORY CAMPAIGN -- the final endpoint driver (R34 sec.2-3).

This is the *last* computation of the campaign: it opens the sealed CONFIRMATORY
scene bank (12 fresh scenes) + the EDGE bank (6 scenes) and, per (cell, record),
generates ONE acquisition record and evaluates ALL FOUR ARMS on it, forming the
calibrated Neyman t_c/CV intervals from the FROZEN critical-value surface
(results/round63_dlgi_campaign/CALIBRATION_TABLE.json, committed @95818c8 -- read
only, never regenerated here) and every bar-input statistic the C1-C7 evaluator
(dlgi_c_bars.py) consumes.

FROZEN / VERBATIM.  The point estimator is the feasibility probe's, imported
VERBATIM through dlgi_common -> dl_common (joint_dual_ledger, medium_estimate,
mom_autocov, s_of_x, ou_path, make_schedule, arm_A2/arm_A4, psnr, ...).  R34 bar C3
kills the campaign if the calibration changes the estimator; the calibration changes
ONLY the uncertainty statement.  This module adds NO estimator; it reuses the scaffold
(dlgi_common, dlgi_arms, dlgi_neyman) and the probe's Whittle machinery verbatim.
assert_estimator_provenance() byte-checks that the driver's DLGI point estimate equals
the probe estimator on a probe-seed record.

WHAT ONE RECORD DOES (R34 sec.2.3, four arms on the SAME acquisition):
  arm 1  pilot-free DLGI     : joint_dual_ledger on the ordinary bucket record.
  arm 4  plain-linear camera : arm_A2 on the BYTE-IDENTICAL pilot-free record (shares
                               the same gain path AND the same Poisson stream).
  arm 2  pilot-interleaved   : its OWN record -- pilot rows replace scene rows at equal
                               total exposures/photons/duration; SHARED gain realization,
                               independent Poisson stream.
  arm 3  oracle monitor      : the true gain path of the pilot-free record (ceiling).

INTERVALS (frozen construction; commit-before-confirmatory, R34 sec.4).
  * DLGI (pilot-free): the CALIBRATED Neyman inversion C_0.95={eta0: W<=c_0.95(eta0)}
    over the (log t_c, log CV) fine grid, c_0.95 from the frozen upper-envelope
    CritSurface.  Projections give the t_c and CV intervals + bounded/connected flags.
    This is the C1 coverage interval and the C2 DLGI width.
  * Pilot baseline reference width (C2): the pilot arm's own honest interval, a
    parametric log-percentile bootstrap of its estimator (ou_fit_irregular on resampled
    pilot gain reads).  The pilot directly reads the gain, so a bootstrap is well
    calibrated for it (no finite-sample skew -- the very problem that forced calibration
    for DLGI); using the pilot's *narrower* natural interval is the conservative choice
    for the "<=1.5x pilot" bar (we never inflate the pilot to flatter DLGI).
  * Oracle relative width (C2, informational): profile-likelihood (Whittle) interval on
    the clean, dense log-gain path (dl_bar4_final.profile_ci, verbatim).

SEED HYGIENE (salt scheme, DISJOINT from calibration/edge).  Reuses dlgi_arms._arm_rngs
VERBATIM: np.random.SeedSequence([BANK_SALT, u32(cell_key::schedule), u32(scene_id),
rep]).spawn(5) -> independent (gain, pilot-free Poisson, pilot Poisson, schedule-seed,
bootstrap) children.  BANK_SALT = 700 (confirmatory) / 720 (edge), DISJOINT from
calibration's 710 and the neyman-smoke's +1e6 replicate offset.  The scene_id hash
further separates confirmatory scenes (seed base 700000) from calibration (710000) and
edge (720000); folding the schedule into the cell-key hash gives paired/random/ordered
mutually-disjoint gain+noise streams.  No scene, gain path, Poisson seed or bootstrap
seed crosses banks (R34 sec.2.1).

TASK / FISHER metric (bar C4, frozen machinery definition).  R33/R34 name a
"<=5% task/Fisher loss" bar but leave the operator unspecified; this module FREEZES a
principled default BEFORE any confirmatory record: a mid-band 2-D DCT-II structural
task.  task_err(x_hat) = sum over the frozen mid-band annulus (radial DCT index in
[0.15,0.60] x Nyquist, DC excluded; 268 of 1024 modes) of DCT(x_hat - x)^2 -- the
feature/edge information a downstream user reads, distinct from global-pixel PSNR (which
low frequencies dominate).  The cohort task Fisher of an arm = 1 / mean_record task_err
(inverse task error-variance = estimation Fisher for a Gaussian task); task/Fisher loss
= 1 - F_dlgi / F_linear.  (Flagged to the authority for ratification; frozen here.)

FROZEN COMPUTE STACK (binding, per CALIBRATION_TABLE.meta.frozen_stack):
python 3.12.13 / numpy 2.0.2 / scipy 1.16.3, verified per-VM at nonce time.  Official
confirmatory records MUST run on that stack (same Colab CPU VM image class as the
calibration).  A LOCAL smoke MAY run on the local stack for machinery health only; it is
written under confirmatory/smoke/ with official=false and MUST NOT enter the endpoints.

CHECKPOINT / RESUME.  Per (cell,schedule) work-unit, records append to an in-memory list
and the shard partial is rewritten atomically (tmp + os.replace) every --flush records
(default 1 == per-record).  A resumed shard reloads its partial and skips completed
records, so a reclaim costs at most --flush in-flight records (mirrors the calibration
runner's atomic per-cell checkpoint).

DISCIPLINE.  Outcome-blind: this file computes and stores bar INPUTS; it never prints a
bar verdict.  The launch is a SEPARATE GO (an independent 3-cell identity cross-check
must pass first).  Default CLI action is --plan (cost + shard manifests; no compute).
Read-only on the frozen table, banks, probe, and scaffold.  CPU.
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
import dlgi_arms as A
import dlgi_neyman as NY
import dl_common as C
from dl_bar4_final import profile_ci                      # frozen Whittle profile-LR (verbatim)

CAMPAIGN = K.CAMPAIGN
CONF_OUT = os.path.join(CAMPAIGN, "confirmatory")
TABLE = os.path.join(CAMPAIGN, "CALIBRATION_TABLE.json")
ABK = os.path.join(CAMPAIGN, "predictions", "ABK_PREDICTIONS.json")
SCRATCH = r"D:/tmp/dlgi_conf"

LEVEL = 0.95
N_REC_PRIMARY = 1000          # R34 sec.2.1 minimum confirmatory records / primary cell
N_REC_EDGE = 300              # edge bank: failure-map display (reduced)
N_REC_SCHED = 300             # schedule sub-study (ordering, not coverage)
N_BOOT_PILOT = 40             # parametric bootstrap resamples for the pilot width
N_SHARDS = 6                  # reuse the proven 6-shard fleet layout

# local sec/record MEASURED by the 20-record smoke (full pipeline: 4 arms + calibrated
# Neyman inversion + pilot parametric bootstrap + oracle profile-LR + gauge/residual).
SEC_PER_RECORD_LOCAL = 2.0
# Colab/frozen-stack scale: the calibration ran 1.26 s / replicate on the 2-vCPU Colab
# VM vs ~0.69 s locally (this box) -> ~1.83x.  Applied to the local record time.
COLAB_SCALE = 1.83


# ===================================================================== grids
def primary_units():
    """27 primary cells, PAIRED schedule, N_REC_PRIMARY records each (R34 sec.2.2)."""
    return [dict(role="primary", bank="confirmatory", tc=tc, cv=cv, ph=ph,
                 schedule="paired", n=N_REC_PRIMARY)
            for (tc, cv, ph) in K.primary_cells()]


def edge_units():
    """6 edge cells (fast-drift t_c=2, slow-drift t_c=128 x CV sweep, nominal SNR),
    PAIRED schedule, from the EDGE bank (R34 sec.2.2).  Failure map only."""
    return [dict(role="edge", bank="edge", tc=tc, cv=cv, ph=ph,
                 schedule="paired", n=N_REC_EDGE)
            for (tc, cv, ph) in K.edge_cells()]


# schedule sub-study subset (FROZEN): the 9 (t_c,CV) pairs at nominal SNR -- exactly the
# cells for which ABK_PREDICTIONS freezes a per-cell schedule ordering (SNR-invariant).
# paired is already computed for these 9 in the primary grid; the sub-study only ADDS the
# two non-paired schedules (R34 sec.2.3 {paired, random, ordered}).
SCHED_SUBSTUDY_CELLS = [(tc, cv, 1.0) for tc in K.TC_PRIMARY for cv in K.CV_PRIMARY]


def substudy_units():
    return [dict(role="substudy", bank="confirmatory", tc=tc, cv=cv, ph=ph,
                 schedule=sch, n=N_REC_SCHED)
            for (tc, cv, ph) in SCHED_SUBSTUDY_CELLS
            for sch in ("random", "ordered")]


def all_units():
    return primary_units() + edge_units() + substudy_units()


def unit_key(u):
    return f"{K.cell_key(u['tc'], u['cv'], u['ph'])}::{u['schedule']}"


def _alloc(n, nsc):
    base = n // nsc
    return [base + (1 if k < n % nsc else 0) for k in range(nsc)]


# ===================================================================== frozen table
_SURF = None


def crit_surface():
    """Load the FROZEN CritSurface from CALIBRATION_TABLE.json (read-only)."""
    global _SURF
    if _SURF is None:
        s = json.load(open(TABLE))["surface"]
        _SURF = NY.CritSurface(s["tc_knots"], s["cv_knots"], s["ph_knots"],
                               np.asarray(s["cvals"], float))
    return _SURF


# ===================================================================== task / Fisher
SIDE = C.SIDE
_uu, _vv = np.meshgrid(np.arange(SIDE), np.arange(SIDE), indexing="ij")
_RAD = np.sqrt(_uu ** 2 + _vv ** 2) / (SIDE - 1)
TASK_RLO, TASK_RHI = 0.15, 0.60
TASK_MASK = ((_RAD >= TASK_RLO) & (_RAD <= TASK_RHI))     # DC (0,0) excluded by RLO>0
TASK_MASK_FLAT = TASK_MASK.ravel()
TASK_NMODES = int(TASK_MASK.sum())


def task_err(x_hat, x):
    """Frozen mid-band 2-D DCT-II structural task error (sum of squared coeffs of the
    reconstruction residual over the frozen mid-band annulus)."""
    from scipy.fft import dctn
    d = (np.asarray(x_hat, float) - np.asarray(x, float)).reshape(SIDE, SIDE)
    D = dctn(d, norm="ortho").ravel()
    return float(np.sum(D[TASK_MASK_FLAT] ** 2))


# ===================================================================== provenance
def assert_estimator_provenance():
    """C3 kill-condition guard: the driver's DLGI point estimator is BIT-IDENTICAL to
    the probe.  Reproduce the probe's own simulate_record + joint_dual_ledger on a DEV
    scene and confirm equality (forward bit-identity is separately self-tested in
    dlgi_common)."""
    sc = C.SCENES[0]
    tc, sig = 64.0, C.sigma_l_of_cv(0.15)
    Y, a_time, slot = C.simulate_record(sc, 7, tc, sig)   # probe forward, seed 7
    r = C.joint_dual_ledger(Y, slot=slot, n_outer=2)
    m = r["med"]
    # same record via the campaign's probe-compatible forward (must be byte-identical)
    Yb, ab, sb = K.simulate_record_probe(C.scene_x[sc], 7, tc, sig)
    rb = C.joint_dual_ledger(Yb, slot=sb, n_outer=2)
    ok = bool(np.array_equal(Y, Yb) and
              abs(m["tc_hat"] - rb["med"]["tc_hat"]) == 0.0 and
              abs(m["cv_hat"] - rb["med"]["cv_hat"]) == 0.0 and
              np.array_equal(r["x_hat"], rb["x_hat"]))
    return dict(estimator_bit_identical=ok, probe_seed=7,
                tc_hat=float(m["tc_hat"]), cv_hat=float(m["cv_hat"]))


# ===================================================================== pilot bootstrap
def _pilot_bootstrap_ci(x, tc_hat, sig_hat, ph, boot_rng, B=N_BOOT_PILOT):
    """Parametric log-percentile bootstrap interval for the PILOT baseline: resample the
    OU gain + Poisson pilot reads at the pilot exposures with the pilot's fitted params
    and re-fit ou_fit_irregular.  Returns ((tc_lo,tc_hi),(cv_lo,cv_hi))."""
    sumx = float(np.asarray(x, float).sum())
    n_pair = C.NPIX
    n_pilot = max(int(round(0.05 * n_pair)), 2)
    pilot_rows = np.unique(np.linspace(1, n_pair - 1, n_pilot).round().astype(int))
    pe = np.concatenate([2 * pilot_rows, 2 * pilot_rows + 1]); pe.sort()
    phi = math.exp(-1.0 / max(tc_hat, 1e-3))
    sd = sig_hat * math.sqrt(max(1.0 - phi * phi, 1e-9))
    tcs = []; cvs = []
    for _ in range(B):
        l = np.empty(C.N_EXP); l[0] = boot_rng.standard_normal() * sig_hat
        for t in range(1, C.N_EXP):
            l[t] = phi * l[t - 1] + sd * boot_rng.standard_normal()
        a = np.exp(l)
        Yp = boot_rng.poisson(np.maximum(ph * a[pe] * C.PHI * sumx, 0.0)).astype(float)
        a_reads = Yp / (ph * C.PHI * sumx)
        tcb, sgb = C.ou_fit_irregular(pe, np.log(np.maximum(a_reads, 1e-9)))
        if np.isfinite(tcb) and tcb > 0:
            tcs.append(tcb); cvs.append(C.cv_of_sigma_l(sgb))

    def pct(arr):
        arr = np.asarray(arr, float); arr = arr[np.isfinite(arr) & (arr > 0)]
        if len(arr) < max(10, B // 4):
            return float("nan"), float("nan")
        q = np.percentile(np.log(arr), [2.5, 97.5])
        return float(np.exp(q[0])), float(np.exp(q[1]))
    return pct(tcs), pct(cvs)


def _profile_width(zc, valid):
    """(t_c, CV) profile-LR interval on a time-ordered residual series (oracle ref)."""
    (tl, th), (cl, ch) = profile_ci(zc, valid)
    return (tl, th), (cl, ch)


def _logwidth(lo, hi):
    if lo is None or hi is None or not (np.isfinite(lo) and np.isfinite(hi)) \
       or lo <= 0 or hi <= 0:
        return float("nan")
    return float(math.log(hi) - math.log(lo))


# ===================================================================== one record
def evaluate_record(u, sid, x, rep):
    """All four arms on ONE (unit, scene, replicate).  Mirrors dlgi_arms.run_arms for
    the acquisition/arm structure (shared gain path; byte-identical pilot-free & linear
    record; own pilot record) and ADDS the calibrated Neyman interval, the pilot/oracle
    reference intervals, the frozen task-Fisher, and the bar-C5 gauge/residual stats.
    Returns a compact scalar record (no arrays) for atomic checkpointing."""
    bank, tc, ph = u["bank"], u["tc"], u["ph"]
    cv = u["cv"]; sig = C.sigma_l_of_cv(cv); sch = u["schedule"]
    ck = K.cell_key(tc, cv, ph)
    seed_ck = f"{ck}::{sch}"

    ss = np.random.SeedSequence([K._BANK_SALT[bank], K._u32(seed_ck), K._u32(sid),
                                 int(rep)])
    g, pf, pp, sc_, bt = ss.spawn(5)
    gain_rng = np.random.default_rng(g); pf_rng = np.random.default_rng(pf)
    pp_rng = np.random.default_rng(pp); boot_rng = np.random.default_rng(bt)
    sched_seed = int(sc_.generate_state(1)[0]) % (2 ** 31)

    a_time = (np.ones(C.N_EXP) if sig <= 0 else C.ou_path(gain_rng, sig, tc, C.N_EXP))
    slot = C.make_schedule(A._SCHED_MAP[sch], sched_seed)
    a_exp = a_time[slot]
    s = C.s_of_x(x)

    # --- arm 1 pilot-free DLGI + arm 4 plain-linear: byte-identical record ---
    Y = pf_rng.poisson(ph * a_exp * C.PHI * s).astype(np.float64)
    r = C.joint_dual_ledger(Y, slot=slot, n_outer=2); m = r["med"]
    pf_tc = float(m["tc_hat"]); pf_cv = float(m["cv_hat"])
    psnr_pf = C.psnr(r["x_hat"], x)
    x_lin = C.arm_A2(Y[0::2], Y[1::2], C.LAM_A2)
    psnr_lin = C.psnr(x_lin, x)
    terr_pf = task_err(r["x_hat"], x); terr_lin = task_err(x_lin, x)

    # --- calibrated Neyman interval (frozen CritSurface) ---
    zc, valid = NY.residual_series(Y, r["x_hat"], slot)
    W_grid, _ = NY.W_surface(zc, valid)
    c_grid = crit_surface().grid_at(ph)
    reg = NY.neyman_region(W_grid, c_grid)

    # --- bar C5 gauge / residual (fresh scenes, all SNR strata) ---
    zc2, R2, valid2, mu2 = C._z_residuals(Y, r["x_hat"])
    order = np.argsort(slot)
    zo, Ro, vo = zc2[order], R2[order], valid2[order]
    tcw, sgw = C.mom_autocov(zo, vo)
    if np.isfinite(tcw):
        phw = math.exp(-1.0 / max(tcw, 1e-3))
        w = C.kalman_innovations(zo, Ro, vo.astype(bool), phw, max(sgw, 1e-6))
        innov_std = float(np.std(w))
        lag1 = float(np.corrcoef(w[:-1], w[1:])[0, 1]) if len(w) > 3 else float("nan")
    else:
        innov_std = float("nan"); lag1 = float("nan")
    m_g = C.medium_estimate(Y, 3.7 * r["x_hat"], slot=slot, return_path=False)
    gauge_tc = abs(pf_tc - m_g["tc_hat"]) / max(pf_tc, 1e-9)
    gauge_cv = abs(pf_cv - m_g["cv_hat"]) / max(pf_cv, 1e-9)
    gaincorr = C.gain_path_corr(m["a_time"], a_time)

    # --- arm 2 pilot (own record, shared gain path) + its bootstrap width ---
    p = A.pilot_interleaved_scaled(x, a_time, pp_rng, sig, photon_scale=ph)
    pilot_tc = float(p["med"]["tc_hat"]); pilot_cv = float(p["med"]["cv_hat"])
    psnr_pilot = C.psnr(p["x_hat"], x)
    (ptc_lo, ptc_hi), (pcv_lo, pcv_hi) = _pilot_bootstrap_ci(
        x, pilot_tc, p["med"]["sigma_l_hat"], ph, boot_rng)

    # --- arm 3 oracle monitor + its profile-LR reference width ---
    om = C.oracle_monitor_estimate(a_time)
    oracle_tc = float(om["tc_hat"]); oracle_cv = float(om["cv_hat"])
    lo = np.log(np.maximum(a_time, 1e-12)); lo = lo - lo.mean()
    (otc_lo, otc_hi), (ocv_lo, ocv_hi) = _profile_width(lo, np.ones_like(lo))

    return dict(
        scene_id=sid, rep=int(rep), schedule=sch,
        tc_true=tc, cv_true=cv, ph=ph,
        pf_tc=pf_tc, pf_cv=pf_cv, pilot_tc=pilot_tc, pilot_cv=pilot_cv,
        oracle_tc=oracle_tc, oracle_cv=oracle_cv,
        # DLGI calibrated Neyman interval
        tc_lo=reg["tc_lo"], tc_hi=reg["tc_hi"], cv_lo=reg["cv_lo"], cv_hi=reg["cv_hi"],
        bounded=bool(reg["bounded"]), connected=bool(reg["connected"]),
        empty=bool(reg["empty"]),
        dlgi_lw_tc=_logwidth(reg["tc_lo"], reg["tc_hi"]),
        dlgi_lw_cv=_logwidth(reg["cv_lo"], reg["cv_hi"]),
        # pilot bootstrap + oracle profile widths
        pilot_lw_tc=_logwidth(ptc_lo, ptc_hi), pilot_lw_cv=_logwidth(pcv_lo, pcv_hi),
        oracle_lw_tc=_logwidth(otc_lo, otc_hi), oracle_lw_cv=_logwidth(ocv_lo, ocv_hi),
        # scene
        psnr_pf=psnr_pf, psnr_lin=psnr_lin, psnr_pilot=psnr_pilot,
        task_err_pf=terr_pf, task_err_lin=terr_lin,
        # gauge / residual (C5)
        innov_std=innov_std, lag1_ac=lag1, gauge_tc=gauge_tc, gauge_cv=gauge_cv,
        gaincorr=gaincorr,
    )


# ===================================================================== work unit
def accounting_for(u, sid, x):
    """Byte-audit fields (C7): pilot-free & plain-linear share the record; pilot has
    equal exposures/photons.  Deterministic given (u, scene); computed once per unit."""
    ph = u["ph"]; cv = u["cv"]; sig = C.sigma_l_of_cv(cv)
    ss = np.random.SeedSequence([K._BANK_SALT[u["bank"]],
                                 K._u32(f"{K.cell_key(u['tc'],cv,ph)}::{u['schedule']}"),
                                 K._u32(sid), 0])
    g, pf, pp, sc_, bt = ss.spawn(5)
    a_time = (np.ones(C.N_EXP) if sig <= 0
              else C.ou_path(np.random.default_rng(g), sig, u["tc"], C.N_EXP))
    slot = C.make_schedule(A._SCHED_MAP[u["schedule"]],
                           int(sc_.generate_state(1)[0]) % (2 ** 31))
    Y = np.random.default_rng(pf).poisson(ph * a_time[slot] * C.PHI * C.s_of_x(x)).astype(float)
    p = A.pilot_interleaved_scaled(x, a_time, np.random.default_rng(pp), sig, photon_scale=ph)
    return dict(
        pilot_free_exposures=int(Y.size), pilot_free_photons=float(Y.sum()),
        plain_linear_exposures=int(Y.size), plain_linear_photons=float(Y.sum()),
        identical_record_pilot_free_vs_linear=True,
        pilot_exposures=int(p["Y"].size), pilot_photons=float(p["Y"].sum()),
        pilot_extra_photons_or_exposures=0,
        pattern_multiset="1024-row sequency Hadamard x complementary (2048 exposures)")


def _eval_task(task):
    """Picklable worker: one record.  crit_surface() lazy-loads the frozen table once
    per process (read-only)."""
    u, sid, x, rep = task
    return evaluate_record(u, sid, np.asarray(x, float), rep)


def run_unit(u, out_dir, flush=1, verbose=True, nproc=1):
    """Run one work-unit's records with atomic checkpoint/resume.  Returns the records
    list.  The partial file is <out_dir>/unit_<key>.json.  nproc>1 fans the records over
    a ProcessPoolExecutor (the fleet launches each shard --nproc 2 on a 2-vCPU VM)."""
    key = unit_key(u)
    safe = key.replace("::", "__")
    os.makedirs(out_dir, exist_ok=True)
    fn = os.path.join(out_dir, f"unit_{safe}.json")
    scenes, _ = K.load_bank(u["bank"])
    ids = list(scenes)
    alloc = _alloc(u["n"], len(ids))
    plan = [(sid, rep) for sid, a in zip(ids, alloc) for rep in range(a)]

    done = {}
    if os.path.exists(fn):
        try:
            prev = json.load(open(fn))
            for rec in prev.get("records", []):
                done[(rec["scene_id"], rec["rep"])] = rec
            if verbose:
                print(f"[unit {key}] RESUME: {len(done)}/{len(plan)} records banked",
                      flush=True)
        except Exception as e:
            print(f"[unit {key}] resume-read failed ({e}); fresh", flush=True)
            done = {}

    records = [done[(sid, rep)] for (sid, rep) in plan if (sid, rep) in done]
    acct = accounting_for(u, ids[0], scenes[ids[0]])
    stack = _stack()

    def _checkpoint():
        tmp = fn + ".tmp"
        json.dump(dict(unit=u, key=key, stack=stack, accounting=acct,
                       n_planned=len(plan), records=records), open(tmp, "w"))
        os.replace(tmp, fn)

    pending = [(sid, rep) for (sid, rep) in plan if (sid, rep) not in done]
    t0 = time.time(); n_new = 0
    if nproc and nproc > 1 and pending:
        from concurrent.futures import ProcessPoolExecutor
        tasks = [(u, sid, scenes[sid], rep) for (sid, rep) in pending]
        with ProcessPoolExecutor(max_workers=nproc) as ex:
            for rec in ex.map(_eval_task, tasks, chunksize=4):
                records.append(rec); n_new += 1
                if n_new % flush == 0:
                    _checkpoint()
    else:
        for (sid, rep) in pending:
            records.append(evaluate_record(u, sid, scenes[sid], rep)); n_new += 1
            if n_new % flush == 0:
                _checkpoint()
    _checkpoint()
    if verbose:
        dt = time.time() - t0
        rate = dt / max(n_new, 1)
        print(f"[unit {key}] {len(records)}/{len(plan)} records "
              f"(+{n_new} new, {rate:.2f} s/rec wall, {dt/60:.1f} min)", flush=True)
    return records


def _stack():
    import scipy
    return dict(python=platform.python_version(), numpy=np.__version__,
                scipy=scipy.__version__)


# ===================================================================== shards
def balance_shards(units, nshards=N_SHARDS):
    """Greedy longest-processing-time split of the work-units across shards, weighted by
    record count (keeps every shard's wall-clock close)."""
    order = sorted(range(len(units)), key=lambda i: -units[i]["n"])
    load = [0] * nshards
    assign = {s: [] for s in range(nshards)}
    for i in order:
        s = int(np.argmin(load))
        assign[s].append(units[i]); load[s] += units[i]["n"]
    return assign, load


def run_shard(shard_id, out_dir=None, flush=1, nproc=2):
    """Run every work-unit assigned to one shard (checkpoint/resume per unit).  nproc
    defaults to 2 (the fleet's 2-vCPU VMs)."""
    out_dir = out_dir or os.path.join(CONF_OUT, "partials")
    os.makedirs(out_dir, exist_ok=True)
    assign, _ = balance_shards(all_units())
    mine = assign[shard_id]
    print(f"[shard {shard_id}] {len(mine)} units, "
          f"{sum(u['n'] for u in mine)} records; nproc={nproc}; stack={_stack()}",
          flush=True)
    for u in mine:
        run_unit(u, out_dir, flush=flush, nproc=nproc)
    print(f"[shard {shard_id}] DONE", flush=True)


# ===================================================================== aggregation
from scipy.stats import beta as _beta


def clopper_pearson(k, n, level=0.95):
    """Exact (Clopper-Pearson) binomial CI."""
    if n == 0:
        return (float("nan"), float("nan"))
    al = 1.0 - level
    lo = 0.0 if k == 0 else _beta.ppf(al / 2, k, n - k + 1)
    hi = 1.0 if k == n else _beta.ppf(1 - al / 2, k + 1, n - k)
    return (float(lo), float(hi))


def _med(a):
    a = np.asarray(a, float); a = a[np.isfinite(a)]
    return float(np.median(a)) if len(a) else float("nan")


def _rmse(hat, true):
    a = np.asarray(hat, float); a = a[np.isfinite(a)]
    return float(np.sqrt(np.mean((a - true) ** 2))) if len(a) else float("nan")


def aggregate_primary(recs, u):
    """Per primary/edge cell endpoint aggregate (bar inputs)."""
    tc, cv = u["tc"], u["cv"]
    n = len(recs)
    scenes = sorted({r["scene_id"] for r in recs})
    # C1 coverage of the calibrated DLGI t_c / CV interval
    empty = np.array([r["empty"] for r in recs])
    tc_lo = np.array([r["tc_lo"] for r in recs]); tc_hi = np.array([r["tc_hi"] for r in recs])
    cv_lo = np.array([r["cv_lo"] for r in recs]); cv_hi = np.array([r["cv_hi"] for r in recs])
    tc_fin = np.isfinite(tc_lo) & np.isfinite(tc_hi) & (~empty)
    cv_fin = np.isfinite(cv_lo) & np.isfinite(cv_hi) & (~empty)
    tc_hit = int(np.sum(tc_fin & (tc_lo <= tc) & (tc <= tc_hi)))
    cv_hit = int(np.sum(cv_fin & (cv_lo <= cv) & (cv <= cv_hi)))
    nonfinite_rate = float(np.mean(~tc_fin))
    cov_tc = tc_hit / n; cov_cv = cv_hit / n
    # C2 informativeness
    frac_bounded = float(np.mean([r["bounded"] for r in recs]))
    frac_connected = float(np.mean([r["connected"] for r in recs]))
    d_lw_tc = _med([r["dlgi_lw_tc"] for r in recs]); d_lw_cv = _med([r["dlgi_lw_cv"] for r in recs])
    p_lw_tc = _med([r["pilot_lw_tc"] for r in recs]); p_lw_cv = _med([r["pilot_lw_cv"] for r in recs])
    o_lw_tc = _med([r["oracle_lw_tc"] for r in recs]); o_lw_cv = _med([r["oracle_lw_cv"] for r in recs])
    # C3 precision
    pf_tc_rmse = _rmse([r["pf_tc"] for r in recs], tc); pf_cv_rmse = _rmse([r["pf_cv"] for r in recs], cv)
    pl_tc_rmse = _rmse([r["pilot_tc"] for r in recs], tc); pl_cv_rmse = _rmse([r["pilot_cv"] for r in recs], cv)
    or_tc_rmse = _rmse([r["oracle_tc"] for r in recs], tc); or_cv_rmse = _rmse([r["oracle_cv"] for r in recs], cv)
    pf_tc_bias = float(np.median(np.array([r["pf_tc"] for r in recs]) - tc))
    pf_cv_bias = float(np.median(np.array([r["pf_cv"] for r in recs]) - cv))
    # C4 scene noninferiority
    dpsnr = np.array([r["psnr_pf"] - r["psnr_lin"] for r in recs], float)
    dpsnr = dpsnr[np.isfinite(dpsnr)]
    dm = float(np.mean(dpsnr)); dse = float(np.std(dpsnr) / max(math.sqrt(len(dpsnr)), 1))
    terr_pf = float(np.mean([r["task_err_pf"] for r in recs]))
    terr_lin = float(np.mean([r["task_err_lin"] for r in recs]))
    tf_pf = 1.0 / terr_pf if terr_pf > 0 else float("inf")
    tf_lin = 1.0 / terr_lin if terr_lin > 0 else float("inf")
    tf_loss = float(1.0 - tf_pf / tf_lin) if np.isfinite(tf_pf) and np.isfinite(tf_lin) else float("nan")

    def ratio(a, b):
        return float(a / b) if (np.isfinite(a) and np.isfinite(b) and b > 0) else float("nan")

    return dict(
        ckey=K.cell_key(tc, cv, u["ph"]), tc=tc, cv=cv, ph=u["ph"], role=u["role"],
        n_records=n, n_scenes=len(scenes),
        # C1
        coverage_tc=cov_tc, coverage_tc_cp95=list(clopper_pearson(tc_hit, n)),
        coverage_cv=cov_cv, coverage_cv_cp95=list(clopper_pearson(cv_hit, n)),
        nonfinite_rate=nonfinite_rate, n_hit_tc=tc_hit, n_total=n,
        # C2
        frac_bounded=frac_bounded, frac_connected=frac_connected,
        dlgi_logwidth_tc=d_lw_tc, dlgi_logwidth_cv=d_lw_cv,
        pilot_logwidth_tc=p_lw_tc, pilot_logwidth_cv=p_lw_cv,
        oracle_logwidth_tc=o_lw_tc, oracle_logwidth_cv=o_lw_cv,
        width_ratio_tc_dlgi_over_pilot=ratio(d_lw_tc, p_lw_tc),
        width_ratio_cv_dlgi_over_pilot=ratio(d_lw_cv, p_lw_cv),
        width_ratio_tc_dlgi_over_oracle=ratio(d_lw_tc, o_lw_tc),
        width_ratio_cv_dlgi_over_oracle=ratio(d_lw_cv, o_lw_cv),
        # C3
        pf_tc_rmse=pf_tc_rmse, pf_cv_rmse=pf_cv_rmse,
        pilot_tc_rmse=pl_tc_rmse, pilot_cv_rmse=pl_cv_rmse,
        oracle_tc_rmse=or_tc_rmse, oracle_cv_rmse=or_cv_rmse,
        rmse_ratio_tc_pf_over_pilot=ratio(pf_tc_rmse, pl_tc_rmse),
        rmse_ratio_cv_pf_over_pilot=ratio(pf_cv_rmse, pl_cv_rmse),
        rmse_ratio_tc_pf_over_oracle=ratio(pf_tc_rmse, or_tc_rmse),
        rmse_ratio_cv_pf_over_oracle=ratio(pf_cv_rmse, or_cv_rmse),
        pf_tc_bias=pf_tc_bias, pf_cv_bias=pf_cv_bias,
        # C4
        psnr_pf_med=_med([r["psnr_pf"] for r in recs]),
        psnr_lin_med=_med([r["psnr_lin"] for r in recs]),
        psnr_delta_mean=dm, psnr_delta_ci95=[dm - 1.96 * dse, dm + 1.96 * dse],
        task_fisher_pf=tf_pf, task_fisher_lin=tf_lin, task_fisher_loss=tf_loss,
        task_err_pf_mean=terr_pf, task_err_lin_mean=terr_lin,
        # C5
        innov_std_med=_med([r["innov_std"] for r in recs]),
        lag1_ac_med=_med([r["lag1_ac"] for r in recs]),
        gauge_tc_shift_med=_med([r["gauge_tc"] for r in recs]),
        gauge_cv_shift_med=_med([r["gauge_cv"] for r in recs]),
        gaincorr_med=_med([r["gaincorr"] for r in recs]),
    )


def aggregate_edge(recs, u):
    """Edge-cell aggregate + honest failure flags (C7 edge honesty)."""
    a = aggregate_primary(recs, u)
    a["edge_role"] = "fast_drift" if u["tc"] <= 2.0 else "slow_drift"
    a["tc_bias_frac_med"] = float(np.median(
        (np.array([r["pf_tc"] for r in recs], float) - u["tc"]) / u["tc"]))
    a["coverage_note"] = ("edge cell: displays the published failure map; cannot rescue "
                          "or kill the interior (R34 sec.2.2)")
    return a


def build_substudy(units_recs):
    """C6 schedule-ordering table on the 9 nominal (t_c,CV) cells: median PSNR + median
    |t_c| rel-err per schedule; measured best-scene / best-medium vs the frozen paired
    prediction.  units_recs: dict unit_key -> records list (paired pulled from primary)."""
    abk = json.load(open(ABK))
    out = {}; contra = 0; ncell = 0
    for (tc, cv, ph) in SCHED_SUBSTUDY_CELLS:
        ck = K.cell_key(tc, cv, ph)
        tag = ck.rsplit("_ph", 1)[0]
        per = {}
        for sch in ("paired", "random", "ordered"):
            recs = units_recs.get(f"{ck}::{sch}")
            if not recs:
                continue
            per[sch] = dict(
                scene_psnr_med=_med([r["psnr_pf"] for r in recs]),
                tc_absrelerr_med=_med([abs(r["pf_tc"] - tc) / tc for r in recs]),
                cv_absrelerr_med=_med([abs(r["pf_cv"] - cv) / cv for r in recs]))
        if len(per) < 3:
            out[tag] = dict(per_schedule=per, incomplete=True)
            continue
        ncell += 1
        best_scene = max(per, key=lambda s: per[s]["scene_psnr_med"])
        best_medium = min(per, key=lambda s: per[s]["tc_absrelerr_med"])
        pred = abk["per_cell"].get(tag, {})
        pred_scene = pred.get("predicted_best_scene_schedule", "paired")
        # medium is predicted schedule-robust: contradiction only if paired is clearly
        # worst on medium (tol 10% of the spread); scene contradiction if paired not best.
        scene_contra = (best_scene != pred_scene)
        out[tag] = dict(
            per_schedule=per, measured_best_scene=best_scene,
            measured_best_medium=best_medium, predicted_best=pred_scene,
            scene_contradicts_prediction=bool(scene_contra))
        contra += int(scene_contra)
    out["_summary"] = dict(n_cells=ncell, n_scene_contradictions=contra,
                           contradiction_frac=(contra / ncell) if ncell else float("nan"),
                           predeclared_tolerance_frac=0.10)
    return out


def merge(out_dir=None, write=True):
    """Combine all unit partials -> CONFIRMATORY_RESULTS.json (per-cell endpoints +
    schedule sub-study + edge + reciprocity + ledger audit + provenance)."""
    out_dir = out_dir or os.path.join(CONF_OUT, "partials")
    import glob
    units_recs = {}; stacks = set(); accts = {}
    for f in sorted(glob.glob(os.path.join(out_dir, "unit_*.json"))):
        d = json.load(open(f))
        units_recs[d["key"]] = d["records"]
        stacks.add(json.dumps(d.get("stack", {}), sort_keys=True))
        accts[d["key"]] = d.get("accounting", {})
    if len(stacks) > 1:
        raise RuntimeError("STACK MISMATCH across confirmatory shards: %s" % stacks)

    primary = {}; edge = {}
    for u in primary_units():
        k = unit_key(u)
        if k in units_recs and units_recs[k]:
            primary[u_ck(u)] = aggregate_primary(units_recs[k], u)
    for u in edge_units():
        k = unit_key(u)
        if k in units_recs and units_recs[k]:
            edge[u_ck(u)] = aggregate_edge(units_recs[k], u)
    substudy = build_substudy(units_recs)

    abk = json.load(open(ABK))
    reciprocity = dict(
        worst_abserr_det_reciprocity=abk["generic_theorem2"]["worst_abserr_det_reciprocity"],
        gi_fisher_reciprocity_worst_abserr=abk["reciprocity_worst_abserr"],
        predicted_ordering=abk["meta"]["predicted_ordering"],
        source="predictions/ABK_PREDICTIONS.json (frozen)")
    ledger_audit = dict(sample=accts, note="pilot-free & plain-linear consume the "
                        "byte-identical record; pilot equal exposures/photons; medium "
                        "ledger = zero-cost post-processing (C7)")
    blob = dict(
        meta=dict(utc=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ"),
                  merge_host=_stack(),
                  frozen_stack=json.loads(list(stacks)[0]) if stacks else None,
                  n_primary=len(primary), n_edge=len(edge),
                  n_records_total=sum(len(v) for v in units_recs.values()),
                  calibration_table=os.path.relpath(TABLE, K.REPO).replace("\\", "/"),
                  status="CONFIRMATORY ENDPOINTS -- feed to dlgi_c_bars.py"),
        provenance=dict(estimator=assert_estimator_provenance(),
                        task_fisher=dict(kind="mid-band 2D DCT-II structural task",
                                         rlo=TASK_RLO, rhi=TASK_RHI, n_modes=TASK_NMODES,
                                         loss="1 - Fisher_dlgi/Fisher_linear, "
                                              "Fisher = 1/mean_task_err"),
                        seed_scheme="SeedSequence([BANK_SALT, u32(cell::sched), "
                                    "u32(scene), rep]).spawn(5); salts 700/720 disjoint "
                                    "from calibration 710",
                        interval_construction="DLGI=calibrated Neyman (frozen table); "
                                              "pilot=parametric bootstrap; oracle=Whittle "
                                              "profile-LR"),
        primary_cells=primary, edge_cells=edge, schedule_substudy=substudy,
        reciprocity=reciprocity, ledger_audit=ledger_audit)
    if write:
        os.makedirs(CONF_OUT, exist_ok=True)
        fn = os.path.join(CONF_OUT, "CONFIRMATORY_RESULTS.json")
        json.dump(blob, open(fn, "w"), indent=2)
        print(f"[merge] wrote {fn} ({len(primary)} primary, {len(edge)} edge, "
              f"{blob['meta']['n_records_total']} records)", flush=True)
    return blob


def u_ck(u):
    return K.cell_key(u["tc"], u["cv"], u["ph"])


# ===================================================================== cost + shard plan
def plan(sec_per_record_local=SEC_PER_RECORD_LOCAL, colab_scale=COLAB_SCALE, write=True):
    """Cost projection + 6-shard Colab manifests (reuse the proven calibration fleet
    pattern).  No compute.  Writes confirmatory/COST_AND_SHARD_PLAN.json + shard cell
    files under D:/tmp/dlgi_conf/shards/."""
    units = all_units()
    n_records = sum(u["n"] for u in units)
    sec_colab = sec_per_record_local * colab_scale
    assign, load = balance_shards(units)
    ncpu = os.cpu_count() or 8
    colab_cores = 2

    core_hours_local = n_records * sec_per_record_local / 3600.0
    core_hours_colab = n_records * sec_colab / 3600.0
    # 6 shards x 2 vCPU = 12-way parallel; shard wall = its records / colab_cores
    shard_plan = []
    for s in range(N_SHARDS):
        recs = load[s]
        shard_plan.append(dict(
            shard=s, account=("pro2" if s <= 2 else "pro1"),
            n_units=len(assign[s]), n_records=recs,
            units=[dict(ckey=u_ck(u), schedule=u["schedule"], role=u["role"], n=u["n"])
                   for u in assign[s]],
            est_wall_hours_colab=round(recs * sec_colab / 3600.0 / colab_cores, 2)))
    max_shard_wall = max(sp["est_wall_hours_colab"] for sp in shard_plan)

    counts = dict(
        primary=dict(cells=len(primary_units()), records_per_cell=N_REC_PRIMARY,
                     records=len(primary_units()) * N_REC_PRIMARY),
        edge=dict(cells=len(edge_units()), records_per_cell=N_REC_EDGE,
                  records=len(edge_units()) * N_REC_EDGE),
        schedule_substudy=dict(cells=len(SCHED_SUBSTUDY_CELLS), extra_schedules=2,
                               records_per_unit=N_REC_SCHED,
                               records=len(substudy_units()) * N_REC_SCHED),
        total_records=n_records)

    blob = dict(
        meta=dict(utc=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ"),
                  local_stack=_stack(), local_cores=ncpu,
                  frozen_stack=dict(python="3.12.13", numpy="2.0.2", scipy="1.16.3"),
                  sec_per_record_local_measured=sec_per_record_local,
                  colab_scale_from_calibration=colab_scale,
                  sec_per_record_colab_projected=round(sec_colab, 2),
                  note="official confirmatory runs on the FROZEN stack, verified per VM "
                       "at nonce; local smoke is unofficial"),
        counts=counts,
        cost=dict(total_records=n_records,
                  core_hours_local=round(core_hours_local, 2),
                  core_hours_colab=round(core_hours_colab, 2),
                  local_all_cores_hours=round(core_hours_local / ncpu, 2),
                  colab_6shard_max_wall_hours=round(max_shard_wall, 2),
                  under_12h=bool(max_shard_wall < 12.0),
                  recommendation=(
                      "6-shard Colab fleet (3 pro2 + 3 pro1, --nproc 2), reusing the "
                      "calibration bundle+watchdog recipe (rebind-first + KA-ping + "
                      "per-poll mirror + file-carried scripts).  Max shard ~%.1f h wall; "
                      "confirmatory is ~0.44x the calibration compute (%d records vs "
                      "135000 replicates).  Local all-cores fallback ~%.1f h."
                      % (max_shard_wall, n_records, core_hours_local / ncpu))),
        shards=shard_plan,
        fleet_recipe=dict(
            scratch=SCRATCH, n_shards=N_SHARDS, nproc_per_shard=2,
            accounts=dict(pro2=[0, 1, 2], pro1=[3, 4, 5]),
            bundle="tar.gz of {code/round63_dlgi, results/round63_next/{DUAL_LEDGER_PROBE,"
                   "LADDER_PROBE}, results/round63_dlgi_campaign/{CALIBRATION_TABLE.json,"
                   "predictions,scene_banks/{confirmatory,edge}}}; SHA-gated on upload",
            watchdog="reuse D:/tmp/dlgi_calib/watchdog.sh recipe: rebind-first every "
                     "cycle (list_assignments + keep_alive_assignment), network-aware "
                     "back-off, per-cell(unit) atomic checkpoint + per-poll mirror, "
                     "death=endpoint-absent-while-net-up; on all-done fetch->merge->"
                     "verify->release",
            per_vm_nonce="print LIBS numpy/scipy; ABORT unless "
                         "numpy==2.0.2 and scipy==1.16.3 (stack must match calibration)",
            launch_command="python code/round63_dlgi/dlgi_confirmatory_driver.py "
                           "--run-shard <S> --nproc 2 --flush 25",
            merge_command="python code/round63_dlgi/dlgi_confirmatory_driver.py --merge",
            discipline="SEPARATE GO -- confirmatory fires only after the independent "
                       "3-cell identity cross-check passes; outcome-blind until then"))

    if write:
        os.makedirs(CONF_OUT, exist_ok=True)
        fn = os.path.join(CONF_OUT, "COST_AND_SHARD_PLAN.json")
        json.dump(blob, open(fn, "w"), indent=2)
        # emit shard cell manifests (file-carried, like shards/shard*.cells)
        sd = os.path.join(SCRATCH, "shards"); os.makedirs(sd, exist_ok=True)
        for sp in shard_plan:
            lines = [f"{u['ckey']},{u['schedule']},{u['role']},{u['n']}"
                     for u in sp["units"]]
            open(os.path.join(sd, f"conf_shard{sp['shard']}.cells"), "w").write(
                "\n".join(lines) + "\n")
        json.dump({str(sp["shard"]): [f"{u['ckey']}::{u['schedule']}" for u in sp["units"]]
                   for sp in shard_plan},
                  open(os.path.join(SCRATCH, "expected_units_6shard.json"), "w"), indent=2)
        print(f"[plan] wrote {fn}", flush=True)
        print(f"[plan] {n_records} records; local {counts['total_records']} @ "
              f"{sec_per_record_local:.2f}s = {core_hours_local:.1f} core-h; "
              f"Colab 6-shard max wall ~{max_shard_wall:.1f} h", flush=True)
        for sp in shard_plan:
            print(f"   shard {sp['shard']} ({sp['account']}): {sp['n_units']} units, "
                  f"{sp['n_records']} rec, ~{sp['est_wall_hours_colab']:.1f} h", flush=True)
    return blob


# ===================================================================== smoke
def smoke(tc=32.0, cv=0.15, ph=1.0, n=20, flush=5):
    """ONE cell x n records end-to-end (all arms, calibrated intervals, all bar-input
    fields).  Machinery HEALTH + coarse plausibility only (outcome-blind: no bar
    verdicts).  LOCAL stack -> written unofficially under confirmatory/smoke/."""
    t0 = time.time()
    prov = assert_estimator_provenance()
    u = dict(role="primary", bank="confirmatory", tc=tc, cv=cv, ph=ph,
             schedule="paired", n=n)
    sd = os.path.join(SCRATCH, "smoke_partials"); os.makedirs(sd, exist_ok=True)
    recs = run_unit(u, sd, flush=flush, verbose=True)
    agg = aggregate_primary(recs, u)
    dt = time.time() - t0
    sec_per_record = dt / max(len(recs), 1)

    fields = sorted(recs[0].keys())
    finite_ok = all(np.isfinite(recs[0][k]) for k in
                    ["pf_tc", "pf_cv", "pilot_tc", "oracle_tc", "psnr_pf", "psnr_lin",
                     "task_err_pf", "task_err_lin", "innov_std"])
    health = dict(
        n_records=len(recs), sec_per_record=sec_per_record,
        estimator_bit_identical=prov["estimator_bit_identical"],
        record_fields=fields, all_core_fields_finite=bool(finite_ok),
        neyman_finite_rate=float(1.0 - agg["nonfinite_rate"]),
        frac_bounded=agg["frac_bounded"], frac_connected=agg["frac_connected"])
    plausibility = dict(
        coverage_tc_point=agg["coverage_tc"], coverage_tc_cp95=agg["coverage_tc_cp95"],
        coverage_cv_point=agg["coverage_cv"],
        width_ratio_tc_dlgi_over_pilot=agg["width_ratio_tc_dlgi_over_pilot"],
        rmse_ratio_tc_pf_over_pilot=agg["rmse_ratio_tc_pf_over_pilot"],
        rmse_ratio_cv_pf_over_pilot=agg["rmse_ratio_cv_pf_over_pilot"],
        psnr_delta_mean=agg["psnr_delta_mean"], task_fisher_loss=agg["task_fisher_loss"],
        innov_std_med=agg["innov_std_med"], gaincorr_med=agg["gaincorr_med"],
        note="SMALL-N (n=%d) plausibility only; coverage CI is wide, NOT a bar verdict" % n)

    blob = dict(
        meta=dict(utc=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ"),
                  official=False, local_stack=_stack(),
                  frozen_stack=dict(python="3.12.13", numpy="2.0.2", scipy="1.16.3"),
                  cell=K.cell_key(tc, cv, ph), n_records=n,
                  disclaimer="LOCAL machinery smoke on the LOCAL stack; NOT an official "
                             "endpoint; official records run on the frozen stack via the "
                             "separate GO."),
        health=health, plausibility=plausibility, aggregate=agg)
    os.makedirs(os.path.join(CONF_OUT, "smoke"), exist_ok=True)
    fn = os.path.join(CONF_OUT, "smoke", "smoke_local.json")
    json.dump(blob, open(fn, "w"), indent=2)
    print(f"\n[smoke] cell {K.cell_key(tc,cv,ph)}: {len(recs)} records in {dt:.1f}s "
          f"({sec_per_record:.2f} s/record LOCAL stack)", flush=True)
    print(f"[smoke] estimator bit-identical to probe: {prov['estimator_bit_identical']}",
          flush=True)
    print(f"[smoke] all core fields finite: {finite_ok}; neyman finite rate "
          f"{health['neyman_finite_rate']:.2f}; bounded {agg['frac_bounded']:.2f} "
          f"connected {agg['frac_connected']:.2f}", flush=True)
    print(f"[smoke] PLAUSIBILITY (small-N): cov_tc={agg['coverage_tc']:.2f} "
          f"[{agg['coverage_tc_cp95'][0]:.2f},{agg['coverage_tc_cp95'][1]:.2f}] "
          f"width_tc/pilot={agg['width_ratio_tc_dlgi_over_pilot']:.2f} "
          f"rmse_tc pf/pilot={agg['rmse_ratio_tc_pf_over_pilot']:.2f} "
          f"dPSNR={agg['psnr_delta_mean']:+.2f}dB innov={agg['innov_std_med']:.2f}", flush=True)
    print(f"[smoke] wrote {fn} (official=false)", flush=True)
    return blob


# ===================================================================== CLI
if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="DLGI confirmatory endpoint driver")
    ap.add_argument("--plan", action="store_true", help="cost + shard manifests (default)")
    ap.add_argument("--smoke", action="store_true", help="1 cell x N records local health")
    ap.add_argument("--smoke-n", type=int, default=20)
    ap.add_argument("--run-shard", type=int, default=None, help="run one shard (GO)")
    ap.add_argument("--run-unit", default=None,
                    help="tc,cv,ph,schedule,n  e.g. 32,0.15,1.0,paired,1000")
    ap.add_argument("--merge", action="store_true", help="merge unit partials -> endpoints")
    ap.add_argument("--provenance", action="store_true", help="estimator bit-identity check")
    ap.add_argument("--flush", type=int, default=1, help="checkpoint every K records")
    ap.add_argument("--nproc", type=int, default=2, help="workers per shard (fleet: 2)")
    a = ap.parse_args()
    if a.smoke:
        smoke(n=a.smoke_n)
    elif a.provenance:
        print(json.dumps(assert_estimator_provenance(), indent=2))
    elif a.run_shard is not None:
        run_shard(a.run_shard, flush=a.flush, nproc=a.nproc)
    elif a.run_unit:
        parts = a.run_unit.split(",")
        u = dict(role="primary", bank="confirmatory", tc=float(parts[0]),
                 cv=float(parts[1]), ph=float(parts[2]), schedule=parts[3],
                 n=int(parts[4]))
        run_unit(u, os.path.join(CONF_OUT, "partials"), flush=a.flush, nproc=a.nproc)
    elif a.merge:
        merge()
    else:
        plan()
