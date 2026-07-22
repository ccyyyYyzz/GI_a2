"""ROUND63 DEV bridge — six-arm cell harness (T-D, phase 1).

NORMATIVE specs
  docs/ROUND63_BRIDGE_BUILD_PLAN.md          (T-D entry, interface contract)
  docs/ROUND63_GPT_ROUND24_RULING_RAW.md  §3 (cohort / physics / arms / gates)
  docs/ROUND63_GPT_ROUND25_RULING_RAW.md §11-12 (Gate C replacement; A_j)

Six arms per R24 §3.3 (one cell = scene x arm x nu x c x seed):
  1 SCAT32-060   fixed L0 rows, global load 0.60
  2 RIDGE-SCAT32 L0 rows, global multiplier -> rho_R(nu) (m1_runner ridge logic)
  3 TRUE-X FW    R18 support-expanding primal on the TRUE scene (oracle; nondeployable)
  4 XHAT FW      same solver on the 52-pattern pre-scan estimate (diagnostic)
  5 RLMI         the full run_rlmi pipeline (rlmi.run_rlmi)
  6 ORACLE-LIB   8 banks materialized pure; best by true PSNR_rad (DEV-only label)

Physics per R24 §3.2:
  M = 52 pre-scan + 972 main = 1024 (pre-scan CHARGED to every arm and INCLUDED
  in reconstruction); tau = 50 ns; nu in {200,2000}; c in {0,0.05}; 5 paired
  seeds; adaptive-arm incident budget capped by the RIDGE arm's per-cell budget;
  identical dwell and pattern counts everywhere.  The c=0.05 cells use the
  jittered-hold forward path (iid lognormal holds, mean tau, CV = c) reproduced
  from the study-2 jitter machinery (jitter_score_diag.simulate); the c=0 path is
  deterministic-hold NP-renewal counting, verified equal to physics.simulate_counts.

Reconstruction: the frozen RQL renewal-quasi-likelihood + isotropic-TV pipeline
  (solvers.run_arm("RQL", ...), truth-free lam_TV via solvers.select_lam_tv — the
  campaign.run_cell RQL path's deployment-legal equivalent).  PSNR_rad primary.

-------------------------------------------------------------------------------
INTERPRETATION REGISTER (load-bearing conventions; quote + choice)
-------------------------------------------------------------------------------
SCALE.  The bridge scenes are raw [0,1] images (sum ~ 244).  The campaign truth
  convention is "unit-mean patterns x SUM-NORMALIZED truth, E[u]=1 analytic"
  (campaign.py docstring).  We sum-normalize every scene to sum == 1 and compute
  PSNR_rad on that normalized truth, exactly as campaign.run_cell does.

LOAD UNITS.  The T-B bank rows are load-profile rows: <row, xhat> == operating
  load rho (atom_load in the bank npz).  We measure the whole harness in load
  units: Phi == 1/tau, so lambda_i = <row_i, x>/tau and rho_i = <row_i, x>.  The
  52 balanced pre-scan rows P (v4.balanced_prescan_52, unit-mean) are scaled to
  PRESCAN_LOAD*P so <P_load_i, x> is also an operating load; this reproduces
  m1.prescan_estimate (rho=0.60, Phi=0.60/tau) exactly (lambda identical).

ESTIMATOR.  The deployed RQL estimator is the c=0 renewal quasi-likelihood; it
  does NOT know the hidden hold CV c.  The c>0 truth/estimator mismatch is the
  physical fact under test — never fed to the estimator.

ARM BUDGET.  The RIDGE-SCAT32 arm defines the per-cell incident budget
  B_inc = sum over its 972 main rows of row-sum (scene-independent light cost).
  The library/adaptive arms that MATERIALIZE (TRUE-X FW, XHAT FW, RLMI, and each
  ORACLE-LIB bank) round under incident_budget = B_inc (rlmi.materialize).

FW SOLVE GRANULARITY (R24 §3.2 / T-D brief; phase-1.5 coordinator ruling 3).
  TRUE-X FW uses only the true scene and is deterministic given (scene, nu) —
  and c-independent (the surrogate uses the c=0 v4._J kernel) — so it is solved
  ONCE per (scene, nu) and cached (reused across both c and all 5 seeds; only the
  final acquisition is per-(c, seed)).  XHAT FW is a DISCLOSED DIAGNOSTIC ECONOMY
  (phase-1.5 ruling 3): it feeds only decomposition gap 2 (XHAT vs TRUE-X, the
  pre-scan/localization loss) and gates nothing, so it is solved at SEED 0 ONLY,
  once per (scene, nu, c) = 64 solves, and its seed-0 design is reused for the
  acquisition at seeds 1-4.  RLMI conditions on the per-seed pre-scan estimate and
  is solved per (scene, nu, c, seed).  ORACLE-LIB bank materializations are
  scene/xhat-free designs but acquired per seed.

FW ITERATIONS (phase-1.5 ruling 3).  primal_probe fw_iters=40 (the R18 primal's
  best feasible log-det is flat from it~20 on the smoke cell; 40 keeps margin).

RIDGE at c=0.05 (phase-1.5 ruling 4).  The RIDGE-SCAT32 arm targets the frozen
  c=0 production ridge rho_R(nu) (22.25 at nu=2000); the c=0.05 information
  optimum is ~5.7, so RIDGE runs "too hot" at c=0.05 and its PSNR is EXPECTED to
  be low there.  This is not a bug: the paired composite comparator
  Q_base = max(SCAT32-060, RIDGE-SCAT32) (R24 §3.3) already absorbs it (SCAT32-060
  dominates at c=0.05).  Do NOT "fix" the ridge target to the c=0.05 optimum —
  that would break the frozen M1 calibration and the composite-baseline contract.

MATERIALIZER v2 (phase-1.5 ruling 1).  rlmi.materialize is the R25 §8/§9 two-stage
  materializer: integer largest-remainder + deterministic dose/incident exchange,
  THEN a bounded per-row source-power repair gamma_r in [0.94, 1.06] (R25 §8
  per-row source settings) that pulls per-pixel dose into the +-5% band; guards +
  realized-regret recomputed WITH the trimmed powers.  Genuine (non-fallback)
  designs carry the gamma vector; only solver/materialization/guard failure
  returns the flagged L0 fallback.  R27 §2: pure library banks use their frozen
  exact-972 admission WITNESS verbatim (rlmi._materialize_from_witness).

GATE A -> LIBRARY_REACHABILITY_PASS (R28 §1).  The M2 science go/no-go is the
  finite-library image oracle: per stress scene, k*_j = min arg max_k five-seed-
  mean PSNR over the 8 GENUINE ORACLE-LIB bank materializations, dQ^A_j =
  ORACLE-LIB - max(SCAT32-060, RIDGE).  It NEVER reads an FW arm.  Gate B's 60%
  capture denominator is this ORACLE-LIB gain (R28 §2).

FW ARMS -- k_eff>=32 restriction + DOSE-RELAXED diagnostic (R28 §3; supersedes
  the R27 dose-compliant witness path, which smoke-4 showed destroys the oracle,
  9.5 dB vs 22.8 baseline, because dose-uniformization removes the FW direction).
  (1) the online primal_probe dictionary is restricted to k_eff>=32 super-atoms
      (_restrict_ctx_keff32; +-5% at 972 rows is infeasible for k=16, R27 §1.1);
  (2) rlmi.materialize_fw_relaxed keeps exactly 972 rows, the incident-budget cap,
      row-wise peak/load/admission/nonneg guards, k>=32 dict, dwell and RQL, but
      removes ONLY the +-5% dose band -- the FW information direction is preserved
      at its frozen amplitudes (NOT dose-uniformized).  Both arms are labeled
      NONDEPLOYABLE_DOSE_RELAXED_FW_DIAGNOSTIC and record the full realized dose
      profile / max deviation / power distribution / budget / peak / loads.  NO
      Gate A-C statistic reads them (asserted by bridge_gates.gates_read_no_fw +
      a unit test).  They feed only the descriptive four-quantity decomposition
      (R28 §4): q3 = XHAT-FW - TRUE-X-FW, q4 = TRUE-X-FW - ORACLE-LIB.

L7 PRUNE (phase-1.5 ruling 2b).  load_banks() prunes L7 to 99% cumulative weight
  mass ONLY when the bank carries no exact-972 admission witness (R27 §2 witnessed
  banks are consumed verbatim).

BLINDNESS.  smoke-1..3 reconstruct ONLY bridge_control_0; smoke-4 additionally
  reconstructs bridge_twopop_0 for the four FW cells (coordinator-authorized,
  materialization validation only).  The other 14 stress scenes stay blind until
  the checkpoint-gated grid.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from dataclasses import dataclass

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
ROOT = os.path.dirname(os.path.dirname(HERE))

import rlmi
import physics
import solvers
import m1_runner as m1
import oed_design_v4 as v4
import oed_design_v5 as v5
import dev_gap_probe as dgp

# --------------------------------------------------------------------------- #
#  frozen constants
# --------------------------------------------------------------------------- #
TAU = 50e-9
SIDE = 32
N_PIX = SIDE * SIDE
M_PRE = 52
M_MAIN = 972
M_TOTAL = M_PRE + M_MAIN
PRESCAN_LOAD = 0.60
PEAK_CAP = 1536.0
DOSE_BAND = 0.05
NU_LIST = (200.0, 2000.0)
C_LIST = (0.0, 0.05)
SEEDS = (0, 1, 2, 3, 4)
ARMS = ("SCAT32-060", "RIDGE-SCAT32", "TRUE-X-FW", "XHAT-FW", "RLMI", "ORACLE-LIB")

SCENES_DIR = os.path.join(ROOT, "data", "r63_bridge_scenes")
LIB_DIR = os.path.join(ROOT, "results", "round63_bridge", "library")
OUT_DIR = os.path.join(ROOT, "results", "round63_bridge")
CELLS_DIR = os.path.join(OUT_DIR, "cells")
CACHE_DIR = os.path.join(OUT_DIR, "cache")

# base RNG streams (disjoint): pre-scan (1), main acquisition (2)
_STREAM_PRESCAN = 1
_STREAM_MAIN = 2


# --------------------------------------------------------------------------- #
#  1.  jittered-hold forward path (study-2 lognormal-hold renewal counting)
# --------------------------------------------------------------------------- #
def _m_max_for(nu, rho_max):
    """Column budget above the largest attainable count (jitter_score_diag._m_max)."""
    mean_ct = nu * rho_max / (1.0 + rho_max)
    return int(mean_ct + 12.0 * np.sqrt(max(nu, 1.0)) + 50.0)


def simulate_counts_jitter(loads, nu, c, rng, chunk=4000):
    """Per-frame detector counts under iid lognormal dead-time holds.

    loads : (M,) operating loads rho_i = <row_i, x> (>= 0).
    Model (jitter_score_diag.simulate, per-frame lambda): active start; live
    waits W ~ Exp(lam_i); hold B ~ tau*LogNormal(mean tau, CV = c) for c>0, else
    deterministic tau; the m-th detection is counted iff
        (sum of m waits) + (sum of m-1 holds) <= T,   T = nu*tau.
    c=0 reduces to exact active-start non-paralyzable renewal counting (verified
    equal to physics._counts_nonpar_vec on a smoke cell).  Returns int counts N.
    """
    loads = np.asarray(loads, dtype=np.float64)
    M = loads.shape[0]
    T = nu * TAU
    lam = np.maximum(loads, 0.0) / TAU
    rho_max = float(loads.max()) if M else 0.0
    m_max = _m_max_for(nu, rho_max)
    N = np.empty(M, dtype=np.int64)
    if c > 0:
        sig2 = np.log1p(c * c)
        mu = -0.5 * sig2
    pos = 0
    while pos < M:
        nb = min(chunk, M - pos)
        A = np.cumsum(rng.exponential(1.0, size=(nb, m_max)), axis=1)
        if c > 0:
            B = TAU * np.exp(mu + np.sqrt(sig2) * rng.standard_normal((nb, m_max)))
        else:
            B = np.full((nb, m_max), TAU)
        D = np.concatenate([np.zeros((nb, 1)), np.cumsum(B, axis=1)[:, :-1]], axis=1)
        slack = T - D                                  # lam-independent
        lam_c = lam[pos:pos + nb][:, None]
        n = (A <= lam_c * slack).sum(axis=1)
        if M and int(n.max()) >= m_max:
            raise RuntimeError("count saturated m_max=%d (rho_max=%.3f nu=%g); "
                               "widen _m_max_for" % (m_max, rho_max, nu))
        N[pos:pos + nb] = n
        pos += nb
    return N


def fano_verification(nu=2000.0, c=0.05, loads=(0.60, 10.0), n_mc=40000, seed=777):
    """Verify the c>0 forward path against the long-window renewal prediction.

    Long-window (renewal) Fano of the count over a window with iid intervals
    W+B (E[W]=1/lam, Var W = 1/lam^2; E[B]=tau, Var B = (c*tau)^2) is
        Fano = Var(interval)/Mean(interval)^2 = (1 + c^2 rho^2)/(1 + rho)^2,
    and E[N] -> nu*rho/(1+rho).  Reports empirical vs predicted at each load, and
    the c=0 vs physics.simulate_counts cross-check at the first load.
    """
    out = {"nu": nu, "c": c, "loads": list(loads), "n_mc": n_mc, "per_load": []}
    for rho in loads:
        rng = np.random.default_rng(seed + int(1000 * rho))
        Nj = simulate_counts_jitter(np.full(n_mc, rho), nu, c, rng)
        mean_pred = nu * rho / (1.0 + rho)
        fano_pred = (1.0 + c * c * rho * rho) / (1.0 + rho) ** 2
        out["per_load"].append({
            "rho": rho,
            "mean_emp": float(Nj.mean()), "mean_pred": float(mean_pred),
            "fano_emp": float(Nj.var() / Nj.mean()), "fano_pred": float(fano_pred),
            "fano_pred_c0": float(1.0 / (1.0 + rho) ** 2),
        })
    # c=0 cross-check against physics.simulate_counts at the first load
    rho0 = loads[0]
    Nj0 = simulate_counts_jitter(np.full(n_mc, rho0), nu, 0.0,
                                 np.random.default_rng(seed + 1))
    det = physics.Detector(tau=TAU)
    _, Np = physics.simulate_counts(np.full(n_mc, rho0), 1.0 / TAU, nu * TAU, det,
                                    np.random.default_rng(seed + 2))
    out["c0_crosscheck"] = {
        "rho": rho0, "jitter_mean": float(Nj0.mean()), "physics_mean": float(Np.mean()),
        "jitter_fano": float(Nj0.var() / Nj0.mean()),
        "physics_fano": float(Np.var() / Np.mean())}
    return out


# --------------------------------------------------------------------------- #
#  2.  scenes + banks
# --------------------------------------------------------------------------- #
def scene_manifest():
    with open(os.path.join(SCENES_DIR, "BRIDGE_SCENES.json")) as f:
        return json.load(f)


def scene_groups():
    """{scene_id: group}; stress = contour|twopop|microtex (12), control = aligned (4)."""
    man = scene_manifest()
    return {s["scene_id"]: s["group"] for s in man["scenes"]}


def load_scene(scene_id):
    """Sum-normalized truth x (n,), sum == 1 (campaign convention)."""
    d = np.load(os.path.join(SCENES_DIR, "%s.npz" % scene_id))
    x = np.asarray(d["x"], dtype=np.float64).ravel()
    s = x.sum()
    return x / s if s > 0 else np.full(N_PIX, 1.0 / N_PIX)


_BANKS = None
L7_PRUNE_MASS = 0.99                     # phase-1.5 ruling 2b


def _prune_measure(bk, mass):
    """Prune a bank measure to `mass` cumulative weight; renormalize weights.
    Guards (peak/admission per row) are carried per-atom, so they re-verify by
    construction after the renorm (linear budget/dose guards inherit)."""
    o = np.argsort(bk.weights)[::-1]
    cum = np.cumsum(bk.weights[o])
    keep = np.sort(o[:int(np.searchsorted(cum, mass)) + 1])
    w = bk.weights[keep] / bk.weights[keep].sum()
    return rlmi.Bank(name=bk.name, rows=bk.rows[keep], weights=w,
                     load=bk.load[keep], incident=bk.incident[keep],
                     peak=bk.peak[keep], admission=bk.admission[keep],
                     is_knob=bk.is_knob, meta=dict(bk.meta, pruned_to_mass=mass,
                                                   pruned_atoms=int(keep.size)))


def load_banks():
    """The 8 mixable measures L0..L7 (L0 = knob).  Cached module-level.

    L7 is pruned to L7_PRUNE_MASS cumulative weight mass (phase-1.5 ruling 2b)
    ONLY when it does NOT carry an exact-972 admission witness.  The R27 §2 bank
    rebuild regenerates LIBRARY_MANIFEST.json with per-bank admission witnesses
    (an exact 972-row realization proving ±5% feasibility); when a bank carries
    that witness, pruning is skipped (the witness is defined on the full measure)
    so smoke-3 consumes the witnessed measures verbatim."""
    global _BANKS
    if _BANKS is None:
        b = [rlmi.load_bank(os.path.join(LIB_DIR, "L%d.npz" % k),
                            name="L%d" % k, is_knob=(k == 0)) for k in range(8)]
        if not _has_admission_witness(b[7]) and b[7].rows.shape[0] > 5000:
            b[7] = _prune_measure(b[7], L7_PRUNE_MASS)
        _BANKS = b
    return _BANKS


def _has_admission_witness(bank):
    """True if the bank carries an R27 §2 exact-972 admission witness
    (schema R27_S2_V1: meta['witness_mult'] loaded from npz admission_witness_mult,
    or meta['admission']['witness_multiplicities'])."""
    m = bank.meta or {}
    if m.get("witness_mult") is not None:
        return True
    adm = m.get("admission")
    return bool(isinstance(adm, dict) and adm.get("witness_multiplicities"))


def reset_caches():
    """Drop the module bank cache + the rlmi scenario/oracle cache so a rerun
    (smoke-3) picks up a freshly regenerated library.  Call after the R27 §2
    T-B rebuild lands a new LIBRARY_MANIFEST.json."""
    global _BANKS
    _BANKS = None
    try:
        rlmi._SM_CACHE.clear()
    except Exception:
        pass


def prescan_matrix_load():
    """52 balanced pre-scan rows scaled to load PRESCAN_LOAD (load units)."""
    return PRESCAN_LOAD * v4.balanced_prescan_52(SIDE)


# --------------------------------------------------------------------------- #
#  3.  pre-scan estimate + reconstruction + metric
# --------------------------------------------------------------------------- #
def _rng(seed, stream, *tags):
    """Deterministic per-(seed, stream, tags) generator (SHA-seeded)."""
    key = "|".join([str(seed), str(stream)] + [str(t) for t in tags])
    h = hashlib.sha256(key.encode()).digest()
    return np.random.default_rng(int.from_bytes(h[:8], "little"))


def run_prescan(x, nu, c, seed, scene_id=""):
    """Acquire the 52 balanced rows (charged) and form xhat (GI + 4x4 block
    smooth + floor + sum-normalize; reproduces m1.prescan_estimate).  Returns
    (xhat, P_load (52,n), N_pre (52,))."""
    P_load = prescan_matrix_load()
    loads = P_load @ x
    rng = _rng(seed, _STREAM_PRESCAN, scene_id, int(nu), int(round(c * 1000)))
    N_pre = simulate_counts_jitter(loads, nu, c, rng)
    det_est = physics.Detector(tau=TAU, dark=0.0)
    ctx = solvers.ArmContext(Phi=1.0 / TAU, det=det_est, T=nu * TAU, side=SIDE,
                             sigma_b=0.0, n_iter=200, select_iter=60,
                             pattern_kind="bridge_prescan",
                             meta={"n_physical_rows": M_PRE})
    xh, _ = solvers.run_arm("GI", P_load, N_pre.astype(np.float64), ctx)
    xh = np.maximum(np.asarray(xh, dtype=np.float64).ravel(), 0.0)
    blk = xh.reshape(SIDE // 4, 4, SIDE // 4, 4).mean(axis=(1, 3))
    xhat = np.repeat(np.repeat(blk, 4, axis=0), 4, axis=1).ravel()
    s = xhat.sum()
    xhat = xhat / s if s > 0 else np.full(N_PIX, 1.0 / N_PIX)
    xhat = np.maximum(xhat, 0.05 / N_PIX)
    return xhat / xhat.sum(), P_load, N_pre


def reconstruct(A_full, N_full, nu, lam_tv=None):
    """Frozen RQL + isotropic-TV reconstruction (truth-free lam_TV).  A_full in
    load units, Phi=1/tau.  Returns (xhat, lam_tv, runtime_s)."""
    det_est = physics.Detector(tau=TAU, dark=0.0)
    ctx = solvers.ArmContext(Phi=1.0 / TAU, det=det_est, T=nu * TAU, side=SIDE,
                             sigma_b=0.0, n_iter=200, select_iter=60,
                             pattern_kind="bridge_main",
                             meta={"n_physical_rows": A_full.shape[0]},
                             lam_tv=lam_tv)
    t0 = time.perf_counter()
    xh, info = solvers.run_arm("RQL", A_full, N_full.astype(np.float64), ctx)
    return (np.asarray(xh, dtype=np.float64).ravel(),
            float(info.get("lam_tv", np.nan)), time.perf_counter() - t0)


def psnr_rad(xh, x):
    """Radiometric PSNR, NO rescaling (campaign spec D2 §4 PRIMARY metric)."""
    xp = np.maximum(xh, 0.0)
    mse = float(np.mean((xp - x) ** 2))
    if mse <= 0:
        return float("inf")
    return 10.0 * np.log10(float(x.max()) ** 2 / mse)


def _load_stats(loads):
    q = np.percentile(loads, [5, 50, 95])
    return {"q5": float(q[0]), "q50": float(q[1]), "q95": float(q[2]),
            "max": float(loads.max()), "mean": float(loads.mean())}


def acquire_and_reconstruct(x, main_rows, P_load, N_pre, nu, c, seed, scene_id,
                            tag, lam_tv=None):
    """Charge pre-scan + main rows, simulate main counts (jittered path), and
    reconstruct on [P_load; main_rows] with [N_pre; N_main].  Returns a result
    dict (PSNR_rad, PSNR, achieved loads, incident, detected, wall)."""
    A_full = np.vstack([P_load, main_rows])
    loads_main = main_rows @ x
    rng = _rng(seed, _STREAM_MAIN, scene_id, tag, int(nu), int(round(c * 1000)))
    N_main = simulate_counts_jitter(loads_main, nu, c, rng)
    N_full = np.concatenate([N_pre, N_main])
    xh, lam, rt = reconstruct(A_full, N_full, nu, lam_tv=lam_tv)
    incident = float(main_rows.sum()) + float(P_load.sum())   # total light cost
    detected = int(N_full.sum())
    return {
        "PSNR_rad": psnr_rad(xh, x),
        "PSNR_flux": _psnr_flux(xh, x),
        "loads_main": _load_stats(loads_main),
        "incident_main": float(main_rows.sum()),
        "incident_total": incident,
        "detected_total": detected,
        "detected_main": int(N_main.sum()),
        "lam_tv": lam, "recon_wall_s": rt,
        "n_rows": int(A_full.shape[0]),
    }, xh


def _psnr_flux(xh, x):
    """Flux-matched PSNR (rescale recon to truth flux) — secondary disclosure."""
    xp = np.maximum(xh, 0.0)
    s = xp.sum()
    if s <= 0:
        return float("nan")
    xr = xp * (x.sum() / s)
    mse = float(np.mean((xr - x) ** 2))
    return 10.0 * np.log10(float(x.max()) ** 2 / mse) if mse > 0 else float("inf")


# --------------------------------------------------------------------------- #
#  4.  arm designs (each returns 972 main rows in load units + metadata)
# --------------------------------------------------------------------------- #
def _l0_rows():
    return load_banks()[0].rows                       # (972, n) at ref load 0.60


def arm_scat32_060(xhat):
    """L0 rows scaled by one global multiplier so predicted mean load == 0.60."""
    rows = _l0_rows()
    m = 0.60 / max(float((rows @ xhat).mean()), 1e-12)
    return rows * m, {"global_multiplier": float(m), "target_load": 0.60}


def _kernel_grid_ok(m_, base_load, kg):
    loads = m_ * base_load
    ll = np.log(np.clip(loads, kg["lo"], kg["hi"]))
    pc = np.interp(ll, kg["lg"], kg["pceil"])
    return (pc.mean() <= v4.CEIL_TARGET and pc.max() <= v4.CEIL_ATOM
            and np.interp(ll, kg["lg"], kg["eff"]).min() >= v4.EFF_MIN
            and np.interp(ll, kg["lg"], kg["bias"]).max() <= v4.BIAS_MAX)


def arm_ridge_scat32(xhat, nu):
    """L0 rows with a global multiplier targeting rho_R(nu), then the frozen
    global safety clip (kernel-grid guards).  Reproduces m1_runner.
    ridge_scat32_calibration's per-dwell multiplier + bisection safety clip
    (documented reimplementation; m1's version caches by m1 image-id and asserts
    an m1/m1_dev imageset, so it is not directly callable on bridge scenes)."""
    rows = _l0_rows()
    base_load = rows @ xhat
    rho_R = float(v4.ridge_target4(int(nu))["rho_R_production"])
    kg = v5._kernel_grids(nu)
    m0 = rho_R / max(float(base_load.mean()), 1e-12)
    clipped = False
    m_ = m0
    if not _kernel_grid_ok(m_, base_load, kg):
        clipped = True
        lo, hi = 0.0, m_
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if _kernel_grid_ok(mid, base_load, kg):
                lo = mid
            else:
                hi = mid
        m_ = lo
    achieved = float(m_ * base_load.mean())
    main = rows * m_
    info = {"global_multiplier": float(m_), "rho_R": rho_R,
            "achieved_load": achieved, "clip_applied": bool(clipped),
            "guard_clipped": bool(achieved < 0.90 * rho_R),
            "incident_budget": float(main.sum())}
    return main, info


# ---- FW oracle machinery (dev_gap_probe primal on the certificate context) -- #
def _floor_design(x):
    """Floor a design image away from exact zero (0.05/N, the m1.prescan_estimate
    floor) and sum-normalize.  The v5 dictionary divides by the local scene load
    (build_blocks: g = w/(w.x)), so a geometry whose support is entirely dark
    yields inf/nan directions.  xhat is already floored in run_prescan; the
    TRUE-X oracle scene (raw truth) must be floored the same way for the DESIGN
    context only — the acquisition still uses the true (possibly-zero) scene."""
    xf = np.maximum(np.asarray(x, dtype=np.float64), 0.05 / N_PIX)
    return xf / xf.sum()


FW_KEFF_MIN = 32                         # R27 (TRUE-X/XHAT hole ruling) point 1


def _restrict_ctx_keff32(ctx, k_min=FW_KEFF_MIN):
    """R27 (FW-hole ruling §1): restrict the online primal_probe dictionary to
    the SAME k_eff>=32 physical super-atom class T-B used for L7.  The +-5% dose
    band at 972 rows is arithmetically infeasible for k=16 atoms (R27 §1.1) for
    ANY arm; the restriction is conservative (it can only WEAKEN the oracle, so
    Gate A cannot be inflated by it).  Implemented by zeroing ALLOW for every
    geometry whose physical support has fewer than k_min pixels, so primal_probe
    never selects a sub-32 atom."""
    k_per_geom = np.concatenate([np.full(mb["IDX"].shape[0], mb["IDX"].shape[1],
                                         dtype=np.int64) for mb in ctx["metas"]])
    ctx["ALLOW"] = ctx["ALLOW"] & (k_per_geom >= k_min)[:, None]
    ctx["k_per_geom"] = k_per_geom
    return ctx


def _fw_ctx(x_design, nu, budget=0.60):
    """Build the R18 certificate context on a design image (x_true for TRUE-X,
    xhat for XHAT), RESTRICTED to k_eff>=32 super-atoms (R27 FW-hole ruling §1).
    info_matrix_full(deployed SCAT32 + pre-scan) -> fixed-star subspace B ->
    setup_ctx_cert (D_load U D_gain palette under the full safety cone) ->
    k_eff>=32 restriction.  The design scene is floored (see _floor_design)."""
    xd = _floor_design(x_design)
    rows_dep, _ = m1.deployed_scat32()
    P = v4.balanced_prescan_52(SIDE)
    V_full = v4.info_matrix_full(rows_dep, xd, int(nu), budget, P=P)
    B, eps0, _ = v4.subspace_from_fixedstar(V_full)
    ctx = v5.setup_ctx_cert(xd, float(nu), budget, B, eps0, SIDE)
    return _restrict_ctx_keff32(ctx)


def _ctx_atom_rows(ctx, gl_pairs):
    """Physical rows of certificate atoms (g_global, l): row = C[g,l]*scatter(
    GVAL[g], IDX[g])  (GVAL is the load-normalized direction, <GVAL,xhat>=1, so
    the row's operating load at xhat equals C[g,l]).  Same identity used by
    v5.dose_of for D_load and setup_ctx_cert for D_gain."""
    metas = ctx["metas"]
    n = ctx["n"]
    C = ctx["C"]
    offs, acc = [], 0
    for mb in metas:
        offs.append(acc)
        acc += mb["IDX"].shape[0]
    offs = np.asarray(offs)
    rows = np.zeros((len(gl_pairs), n))
    for i, (g, l) in enumerate(gl_pairs):
        bi = int(np.searchsorted(offs, g, side="right") - 1)
        mb = metas[bi]
        gl = g - offs[bi]
        rows[i, mb["IDX"][gl]] = C[g, l] * mb["GVAL"][gl]
    return rows


def fw_design_bank(ctx, xi, name, x_design):
    """Convert the primal_probe design measure xi (G_tot x L) into a single
    mixable Bank of support-atom physical rows (weights = xi), ready for the
    rlmi.materialize single-measure case (R25 SS9)."""
    gl = list(map(tuple, np.column_stack(np.nonzero(xi > 1e-9)).tolist()))
    rows = _ctx_atom_rows(ctx, gl)
    w = np.array([xi[g, l] for (g, l) in gl], dtype=np.float64)
    w = w / w.sum()
    load = rows @ x_design
    return rlmi.Bank(name=name, rows=rows, weights=w, load=load,
                     incident=rows.sum(axis=1), peak=rows.max(axis=1),
                     admission=np.ones(rows.shape[0], dtype=bool), is_knob=False)


def fw_design_directions(ctx, xi):
    """FW design as UNIT geometry directions for the T-B witness construction
    (R27 FW-hole ruling §2).  Aggregates the design measure over load LEVELS per
    geometry (weight_g = sum_l xi[g,l]) and returns the unit-amplitude direction
    patterns GVAL[g] scattered onto IDX[g] (NOT the load-scaled C[g,l] rows) --
    the fixed_dose construction sets the per-atom nominal power via the Sinkhorn,
    so the atoms must enter at unit amplitude (base pattern = n/k on support),
    exactly as make_patterns('sparsek', k=32) does in fixed_dose_scat32.  Returns
    (rows (G_sel, n), weights (G_sel,))."""
    metas = ctx["metas"]
    n = ctx["n"]
    wg = xi.sum(axis=1)                                 # weight per geometry
    gsel = np.where(wg > 1e-12)[0]
    offs, acc = [], 0
    for mb in metas:
        offs.append(acc)
        acc += mb["IDX"].shape[0]
    offs = np.asarray(offs)
    rows = np.zeros((gsel.size, n))
    for i, g in enumerate(gsel):
        bi = int(np.searchsorted(offs, g, side="right") - 1)
        mb = metas[bi]
        gl = int(g - offs[bi])
        rows[i, mb["IDX"][gl]] = mb["GVAL"][gl]        # unit direction (g.xhat=1)
    w = wg[gsel]
    return rows, w / w.sum()


FW_LABEL = "NONDEPLOYABLE_DOSE_RELAXED_FW_DIAGNOSTIC"   # R28 §3


def arm_fw(x_design, xhat_mat, nu, budget, B_inc, name, iters=80):
    """R28 §3 dose-RELAXED FW diagnostic (NONDEPLOYABLE).  Run the R18 support-
    expanding primal on the k_eff>=32-restricted dictionary (R27 §1), then
    materialize the per-cell design at its FROZEN FW amplitudes via
    rlmi.materialize_fw_relaxed: exact 972 rows + incident-budget cap + peak/load/
    admission/nonneg guards, but NO +-5% dose band (the FW information direction is
    preserved, not dose-uniformized).  Records the full realized dose profile /
    power distribution / peak / budget / loads.  NO Gate A-C statistic reads this
    arm.  On a peak/incident/count guard violation -> L0 fallback."""
    t0 = time.perf_counter()
    ctx = _fw_ctx(x_design, nu, budget)
    ld, xi, it = dgp.primal_probe(ctx, budget, iters=iters)
    t_solve = time.perf_counter() - t0
    if xi is None:
        rows, sc = arm_scat32_060(xhat_mat)
        return rows, {"fw_solve_wall_s": t_solve, "fw_logdet": None,
                      "fw_support": 0, "fallback": "L0", "fw_label": FW_LABEL, **sc}
    bank = fw_design_bank(ctx, xi, name, xhat_mat)     # frozen FW amplitudes
    t1 = time.perf_counter()
    mat = rlmi.materialize_fw_relaxed(bank.rows, bank.weights, M_rows=M_MAIN,
                                      peak_cap=PEAK_CAP, incident_budget=B_inc)
    t_mat = time.perf_counter() - t1
    if mat is None:
        rows, sc = arm_scat32_060(xhat_mat)
        return rows, {"fw_solve_wall_s": t_solve, "fw_materialize_wall_s": t_mat,
                      "fw_logdet": float(ld), "fw_support": int((xi > 1e-9).sum()),
                      "materialize_flag": "GUARD_FAIL", "fallback": "L0",
                      "fw_label": FW_LABEL, **sc}
    return mat.rows, {"fw_solve_wall_s": t_solve, "fw_materialize_wall_s": t_mat,
                      "fw_logdet": float(ld), "fw_iters": int(it),
                      "fw_support": int((xi > 1e-9).sum()),
                      "materialize_flag": "OK", "fw_label": FW_LABEL,
                      "materialize_guards": mat.guards}


def arm_rlmi(xhat, nu, c, B_inc):
    """The full R25 Step-3 pipeline (rlmi.run_rlmi) over the 8 banks.  Returns
    (main_rows, disclosures)."""
    banks = load_banks()
    rows_dep, _ = m1.deployed_scat32()
    P = v4.balanced_prescan_52(SIDE)
    V_full = v4.info_matrix_full(rows_dep, xhat, int(nu), 0.60, P=P)
    B, eps0, _ = v4.subspace_from_fixedstar(V_full)
    H0_full = rlmi.build_H0_full(xhat, SIDE, float(nu), float(c),
                                 lambda_TV=1.0, eps=1e-6)
    cfg = rlmi.RLMIConfig(n_scenarios=16, nu=float(nu), c=float(c), side=SIDE,
                          M_rows=M_MAIN, dose_band=DOSE_BAND, peak_cap=PEAK_CAP,
                          incident_budget=float(B_inc))
    disc, mat, sm, scen, mm = rlmi.run_rlmi(xhat, banks, B, H0_full, cfg)
    return mat.rows, disc


def arm_oracle_lib(x, xhat, P_load, N_pre, nu, c, seed, scene_id, B_inc):
    """Each of the 8 banks materialized pure (w = e_k) under B_inc, acquired and
    reconstructed; best by TRUE PSNR_rad (DEV-only label).  Returns (best_result,
    per_bank list, best_rows)."""
    banks = load_banks()
    per_bank = []
    best = None
    best_rows = None
    for k, bk in enumerate(banks):
        w = np.zeros(len(banks))
        w[k] = 1.0
        mat = rlmi.materialize(banks, w, xhat, M_rows=M_MAIN, dose_band=DOSE_BAND,
                               incident_budget=float(B_inc), peak_cap=PEAK_CAP,
                               side=SIDE)
        res, _ = acquire_and_reconstruct(
            x, mat.rows, P_load, N_pre, nu, c, seed, scene_id,
            tag="oraclelib_L%d" % k)
        rec = {"bank": bk.name, "PSNR_rad": res["PSNR_rad"],
               "materialize_flag": mat.flag,
               "incident_main": res["incident_main"]}
        per_bank.append(rec)
        if best is None or res["PSNR_rad"] > best["PSNR_rad"]:
            best = {**res, "best_bank": bk.name}
            best_rows = mat.rows
    return best, per_bank, best_rows


# --------------------------------------------------------------------------- #
#  5.  cell runner (one (scene, nu, c, seed) group through requested arms)
# --------------------------------------------------------------------------- #
def _truex_cache_path(scene_id, nu, c):
    # TRUE-X FW design is c-INDEPENDENT (setup_ctx_cert / primal_probe /
    # info_matrix_full take no c; the c=0 v4._J kernel is used for the surrogate)
    # AND seed-independent (uses x_true) -> one solve per (scene, nu) for the
    # whole grid; reused across both c and all 5 seeds.
    return os.path.join(CACHE_DIR, "truex_%s_nu%d.npz" % (scene_id, int(nu)))


def _xhatfw_cache_path(scene_id, nu, c):
    # XHAT FW = disclosed diagnostic economy (phase-1.5 ruling 3): solved at
    # SEED 0 only, once per (scene, nu, c) = 64 solves; the seed-0 design is
    # reused for the acquisition at seeds 1-4 (it feeds only decomposition gap 2
    # and gates nothing).
    return os.path.join(CACHE_DIR, "xhatfw_%s_nu%d_c%03d.npz"
                        % (scene_id, int(nu), int(round(c * 1000))))


def run_group(scene_id, nu, c, seed, arms=ARMS, fw_iters=40, save=True,
              verbose=False):
    """Run one (scene, nu, c, seed) group.  Shares the charged pre-scan across
    all arms; caches the TRUE-X FW design per (scene, nu, c).  Returns a dict
    {arm: result}.  Writes per-(scene,arm,nu,c,seed) JSON when save=True."""
    x = load_scene(scene_id)
    t_group = time.perf_counter()
    xhat, P_load, N_pre = run_prescan(x, nu, c, seed, scene_id)

    # RIDGE first (defines the incident budget for the adaptive/library arms)
    ridge_main, ridge_info = arm_ridge_scat32(xhat, nu)
    B_inc = ridge_info["incident_budget"]

    results = {}

    def _emit(arm, res, extra=None):
        rec = {"scene": scene_id, "group": scene_groups().get(scene_id),
               "arm": arm, "nu": nu, "c": c, "seed": seed, **res}
        if extra:
            rec["arm_info"] = extra
        results[arm] = rec
        if save:
            os.makedirs(CELLS_DIR, exist_ok=True)
            fn = os.path.join(CELLS_DIR, "%s_%s_nu%d_c%03d_s%d.json"
                              % (scene_id, arm, int(nu), int(round(c * 1000)), seed))
            with open(fn, "w") as f:
                json.dump(_jsonable(rec), f, indent=1)
        if verbose:
            print("   [%s] PSNR_rad=%.3f" % (arm, res.get("PSNR_rad", float("nan"))),
                  flush=True)

    if "SCAT32-060" in arms:
        main, sc = arm_scat32_060(xhat)
        res, _ = acquire_and_reconstruct(x, main, P_load, N_pre, nu, c, seed,
                                         scene_id, "scat060")
        _emit("SCAT32-060", res, sc)

    if "RIDGE-SCAT32" in arms:
        res, _ = acquire_and_reconstruct(x, ridge_main, P_load, N_pre, nu, c, seed,
                                         scene_id, "ridge")
        _emit("RIDGE-SCAT32", res, ridge_info)

    if "TRUE-X-FW" in arms:
        cpath = _truex_cache_path(scene_id, nu, c)
        if os.path.exists(cpath):
            dd = np.load(cpath, allow_pickle=True)
            main = dd["rows"]
            fw_info = json.loads(str(dd["info"]))
            fw_info["cache_hit"] = True
        else:
            main, fw_info = arm_fw(x, x, nu, 0.60, B_inc, "TRUEXFW", iters=fw_iters)
            fw_info["cache_hit"] = False
            if save:
                os.makedirs(CACHE_DIR, exist_ok=True)
                np.savez_compressed(cpath, rows=main,
                                    info=json.dumps(_jsonable(fw_info)))
        res, _ = acquire_and_reconstruct(x, main, P_load, N_pre, nu, c, seed,
                                         scene_id, "truexfw")
        _emit("TRUE-X-FW", res, fw_info)

    if "XHAT-FW" in arms:
        cpath = _xhatfw_cache_path(scene_id, nu, c)
        if os.path.exists(cpath):
            dd = np.load(cpath, allow_pickle=True)
            main = dd["rows"]
            fw_info = json.loads(str(dd["info"]))
            fw_info["cache_hit"] = True
        else:
            # seed-0 pre-scan estimate defines the (disclosed) diagnostic design
            xhat0 = xhat if seed == 0 else run_prescan(x, nu, c, 0, scene_id)[0]
            main, fw_info = arm_fw(xhat0, xhat0, nu, 0.60, B_inc, "XHATFW",
                                   iters=fw_iters)
            fw_info["cache_hit"] = False
            fw_info["design_seed"] = 0
            if save:
                os.makedirs(CACHE_DIR, exist_ok=True)
                np.savez_compressed(cpath, rows=main,
                                    info=json.dumps(_jsonable(fw_info)))
        res, _ = acquire_and_reconstruct(x, main, P_load, N_pre, nu, c, seed,
                                         scene_id, "xhatfw")
        _emit("XHAT-FW", res, fw_info)

    if "RLMI" in arms:
        main, disc = arm_rlmi(xhat, nu, c, B_inc)
        res, _ = acquire_and_reconstruct(x, main, P_load, N_pre, nu, c, seed,
                                         scene_id, "rlmi")
        res["rlmi_disclosures"] = disc
        _emit("RLMI", res, {"A_realized": disc.get("A_realized"),
                            "route_latency_s": disc.get("alloc_time_s"),
                            "flag": disc.get("flag"),
                            "w_realized": disc.get("w_realized")})

    if "ORACLE-LIB" in arms:
        best, per_bank, _ = arm_oracle_lib(x, xhat, P_load, N_pre, nu, c, seed,
                                           scene_id, B_inc)
        best["per_bank"] = per_bank
        _emit("ORACLE-LIB", best, {"best_bank": best.get("best_bank")})

    if verbose:
        print("  group (%s,nu=%g,c=%g,s=%d) done %.1fs"
              % (scene_id, nu, c, seed, time.perf_counter() - t_group), flush=True)
    return results


def _jsonable(o):
    if isinstance(o, dict):
        return {k: _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return o


# --------------------------------------------------------------------------- #
#  6.  manifest generator (shardable; same conventions as the M1 colab infra)
# --------------------------------------------------------------------------- #
def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def build_grid_cells():
    """The full DEV-bridge grid as a flat cell list (unit = one (scene,nu,c,seed)
    group; each group emits the six arm outputs).  320 groups."""
    scenes = [s["scene_id"] for s in scene_manifest()["scenes"]]
    cells = []
    for sc in scenes:
        for nu in NU_LIST:
            for c in C_LIST:
                for seed in SEEDS:
                    cells.append({
                        "cell_id": "%s_nu%d_c%03d_s%d"
                                   % (sc, int(nu), int(round(c * 1000)), seed),
                        "scene": sc, "nu": nu, "c": c, "seed": seed,
                        "arms": list(ARMS)})
    return cells


def frozen_inputs():
    """sha256-verified inputs (scenes + banks), repo-relative paths."""
    fi = []
    for s in scene_manifest()["scenes"]:
        p = os.path.join(SCENES_DIR, "%s.npz" % s["scene_id"])
        fi.append({"path": os.path.relpath(p, ROOT).replace("\\", "/"),
                   "sha256": _sha256(p)})
    for k in range(8):
        p = os.path.join(LIB_DIR, "L%d.npz" % k)
        fi.append({"path": os.path.relpath(p, ROOT).replace("\\", "/"),
                   "sha256": _sha256(p)})
    # code provenance (R27 smoke-3/4: harness + allocator file hashes so a code
    # change invalidates the frozen manifest)
    for fn in ("bridge_harness.py", "rlmi.py", "bridge_gates.py"):
        p = os.path.join(HERE, fn)
        if os.path.exists(p):
            fi.append({"path": os.path.relpath(p, ROOT).replace("\\", "/"),
                       "sha256": _sha256(p), "role": "code"})
    return fi


# Route layouts (phase-1.75 sequencing note).  Primary assumes the R16
# diagnostic has released its 2 pro2 slots by smoke-3 time (all 3 pro2 + 2 pro1
# free); the 3-route fallback assumes R16 still holds 2 pro2 slots (only pro2_a
# free) — so the grid runs on pro2_a + pro1_a + pro1_b.
ROUTE_LAYOUTS = {
    "primary_5route": ["pro2_a", "pro2_b", "pro2_c", "pro1_a", "pro1_b"],
    "fallback_3route": ["pro2_a", "pro1_a", "pro1_b"],
}


def _emit_layout(name, routes, cells, fi, keys, out_dir):
    ldir = os.path.join(out_dir, name)
    os.makedirs(ldir, exist_ok=True)
    R = len(routes)
    buckets = {r: [] for r in routes}
    for i, key in enumerate(keys):                     # round-robin (scene, nu) keys
        buckets[routes[i % R]].append(key)
    shards = []
    for r in routes:
        keyset = set(buckets[r])
        scells = [c for c in cells if (c["scene"], c["nu"]) in keyset]
        shard = {"shard_id": "bridge_%s_%s" % (name, r), "stage": "DEV_BRIDGE",
                 "route": r, "layout": name, "frozen_inputs": fi, "cells": scells,
                 "output_dir": "results/round63_bridge/cells",
                 "output_csv": "results/round63_bridge/bridge_%s_%s.csv" % (name, r)}
        p = os.path.join(ldir, "bridge_%s.json" % r)
        with open(p, "w") as f:
            json.dump(shard, f, indent=1)
        shards.append({"route": r, "n_keys": len(buckets[r]),
                       "n_groups": len(scells), "manifest": p})
    return {"layout": name, "routes": routes, "shards": shards}


def generate_manifests(out_dir=None, layouts=None):
    """Write shard manifests (schema mirrors shard_runner.py) partitioning the
    320 groups across Colab routes, for BOTH the primary 5-route layout and the
    3-route R16-contention fallback (phase-1.75 sequencing note).  Groups are
    assigned so that all 5 seeds and both c of a (scene, nu) stay in one shard
    (TRUE-X FW design reuse across seeds; no cross-shard cache dependence)."""
    out_dir = out_dir or os.path.join(OUT_DIR, "manifests")
    os.makedirs(out_dir, exist_ok=True)
    layouts = layouts or ROUTE_LAYOUTS
    cells = build_grid_cells()
    fi = frozen_inputs()
    scenes = [s["scene_id"] for s in scene_manifest()["scenes"]]
    keys = [(sc, nu) for sc in scenes for nu in NU_LIST]        # 32 keys
    emitted = {name: _emit_layout(name, routes, cells, fi, keys, out_dir)
               for name, routes in layouts.items()}
    idx = {"n_groups_total": len(cells), "n_keys": len(keys),
           "layouts": emitted,
           "note": ("primary_5route assumes R16 released its 2 pro2 slots; "
                    "fallback_3route assumes R16 still holds pro2_b/pro2_c "
                    "(~1-2h) so only pro2_a + pro1_a + pro1_b are free.")}
    with open(os.path.join(out_dir, "INDEX.json"), "w") as f:
        json.dump(idx, f, indent=1)
    return idx


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--fano", action="store_true")
    ap.add_argument("--manifests", action="store_true")
    ap.add_argument("--smoke", type=str, default=None,
                    help="scene_id to smoke through all six arms")
    ap.add_argument("--nu", type=float, default=2000.0)
    ap.add_argument("--c", type=float, default=0.05)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    if args.fano:
        print(json.dumps(fano_verification(), indent=1))
    if args.manifests:
        print(json.dumps(generate_manifests(), indent=1))
    if args.smoke:
        res = run_group(args.smoke, args.nu, args.c, args.seed, verbose=True)
        for arm, r in res.items():
            print("%-12s PSNR_rad=%.3f" % (arm, r["PSNR_rad"]))
