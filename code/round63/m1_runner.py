"""M1 method-campaign runner (docs/ROUND63_METHOD_SPEC_M1.md; R10+R11 frozen).

Modeled on study2_runner.py. Six arms, each reconstructed by the frozen RQL
production path through campaign.run_cell on the M1 imageset:

  SCAT16      pattern 'sparsek' k=16          (Study-2 bridge; no gate)
  SCAT32      pattern 'sparsek' k=32          (ladder rung; no gate)
  LBLOB16     pattern 'lblob16'               (strongest fixed baseline)
  RIDGE-FIXED pattern 'lblob16', per-dwell global multiplier from
              oed_design_v3.ridge_fixed_design on the cell's own pre-scan
              estimate (fast rho_bar = achieved mean load rho_ach(nu))
  OED-EQLOAD  v3 optimizer, ALL atoms at the arm rho (kernel ablation)
  OED-DT      v3 optimizer, full R11 palette + incident budget + dose
              machinery (the ONLY gate-carrying arm)

DESIGN GRANULARITY (frozen semantics)
-------------------------------------
ONE design per (image, seed, arm), computed at the TERMINAL dwell nu=2000
palette and that arm's operating rho; the SAME pattern matrix is reused
across the whole nu sweep (the Q90 time-to-quality endpoint compares fixed
patterns swept over dwell). RIDGE-FIXED keeps fixed supports and varies only
its GLOBAL multiplier per dwell (R11 §2: a ridge-tracking policy, not a
servo); that multiplier enters the cell as its fast rho_bar. Designs are
cached to results/round63_m1/designs/<image>_<seed>_<arm>.npz with a sha256
and the full R11 disclosure block (rho_R, clip reasons, budget slack, dose
residual, alpha-mixture, certificate gap).

PRE-SCAN accounting (R10 Q7 §4)
-------------------------------
Every arm's measurement of every cell = the frozen 52-row balanced multiscale
pre-scan (v3's prescan_52_supports: 32 even-parity 4x4 blocks + 16 8x8
blocks + 4 quadrants, each row unit-mean-load normalized (n/k_row) *
indicator) at the cell's own (rho, nu), STACKED with the arm's 972 main
patterns: A = vstack(prescan, main), 1024 physical exposures total; the
pre-scan counts are folded into the reconstruction data for ALL arms (they
are simply the first 52 rows of the cell's A). Only OED/RIDGE arms use the
pre-scan RECONSTRUCTION (GI on the 52 rows, clipped nonneg, 4x4-block
smoothed, sum-normalized) to parameterize their designs; that design-time
estimate is simulated once per (image, seed) at the terminal (0.60, 2000)
operating point on the dedicated rng_for(seed, 63, 9, ...) stream.

Safe reference for EVERY arm = that arm's own patterns at rho = 0.05
(spec §6). Fast operating point: 0.60 for the budget arms; rho_R(nu)
achieved (post safety clip) for RIDGE-FIXED.

Ambiguity ledger (A6..A12; resolved toward the spec text)
---------------------------------------------------------
A6  Design-time pre-scan estimate: the spec's "the cell's own 52-pattern
    pre-scan estimate" conflicts with the frozen design granularity (one
    design per (image, seed, arm) at terminal dwell) if read per-cell.
    Resolved: ONE estimate per (image, seed), simulated at (0.60, 2000) on
    the design stream (63, 9); every cell still FOLDS its own per-cell
    pre-scan counts into reconstruction (data path unaffected).
A7  Fixed arms' 972 main rows: the frozen square families have 1024 rows;
    main = the FIRST 972 rows in the family's frozen order (deterministic;
    the dropped 52 are the accounting cost of the pre-scan slots).
A8  See campaign.py M1 appendix: pattern kinds 'm1pat:<ARM>:<image>' are
    resolved by campaign at run time so the frozen Colab shard infra
    (shard_runner -> campaign.run_cell) runs M1 manifests unchanged.
A9  OED pattern matrices are handed to run_cell normalized to unit MEAN
    predicted load (A_main / mean_i(a_i . xhat)), so the campaign's
    rho_bar semantics (Phi = rho/tau at E[u] = 1) hold; the R11 relative
    per-row load profile is preserved. At fast rho_bar = 0.60 the realized
    incident budget sits exactly AT the cap 972*0.60 (<=, allowed) even
    when the continuous design left slack.
A10 CNR column: campaign only routes frozen ROIs for detail32*; M1 rows
    carry blank cnr (the analyzer uses PSNR_rad; ROI wiring would touch
    non-additive campaign code).
A11 OED-EQLOAD palette: v3.design_v3 hardwires load_palette; the EQLOAD
    single-level palette is injected by a scoped monkey-patch
    (try/finally) — runtime-only, no file modified.
A12 R11 disclosure columns are injected into the LOCAL per-arm CSVs by this
    runner; Colab shard CSVs carry the core run_cell rows and the
    disclosures are joined at merge time from the sha-frozen design caches
    (the caches are declared frozen_inputs of their shards).

Usage (repo root, py311):
  python code/round63/m1_runner.py --designs [--arms OED-DT,OED-EQLOAD,RIDGE-FIXED]
  python code/round63/m1_runner.py --arm OED-DT       # local resume-safe sweep
  python code/round63/m1_runner.py --smoke            # item-6 end-to-end smoke
"""
import argparse
import csv
import hashlib
import json
import os
import sys
import time
import traceback
import zlib

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
ROOT = os.path.dirname(os.path.dirname(HERE))

import oed_design as v1                     # noqa: E402
import oed_design_v2 as v2                  # noqa: E402
import oed_design_v3 as v3                  # noqa: E402
from gi_core.utils import rng_for           # noqa: E402
from patterns import make_patterns          # noqa: E402
from physics import Detector, simulate_counts  # noqa: E402
from solvers import ArmContext, run_arm     # noqa: E402
from select_eta import frozen_C0            # noqa: E402
import pilot_s1                             # noqa: E402

# ---- frozen M1 geometry (spec §2) ----------------------------------------- #
SIDE = 32
N_PIX = SIDE * SIDE
M_TOTAL = 1024
M_PRE = 52
M_MAIN = 972
TAU = 50e-9
NU_FULL = [5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0, 1000.0, 2000.0]
RHO_SAFE = 0.05
RHO_FAST = 0.60
SEEDS5 = [0, 1, 2, 3, 4]
NU_TERMINAL = 2000.0

ARMS_ALL = ["SCAT16", "SCAT32", "LBLOB16", "MATCH1", "RIDGE-FIXED",
            "OED-EQLOAD", "OED-DT"]
DESIGN_ARMS = ["MATCH1", "RIDGE-FIXED", "OED-EQLOAD", "OED-DT"]

OUT_DIR = os.path.join(ROOT, "results", "round63_m1")
DESIGN_DIR = os.path.join(OUT_DIR, "designs")
RIDGE_JSON = os.path.join(DESIGN_DIR, "ridge_targets.json")

# ---- CSV schema ----------------------------------------------------------- #
MECH_COLS = ["cnr", "C_u", "Gamma", "S_det", "S_inc", "k_occupancy"]
PROV_COLS = ["m1_arm", "stage", "cell_id"]
DISC_COLS = ["rho_R_production", "ridge_clip_reason", "requested_mean_load",
             "achieved_mean_load", "RIDGE_GUARD_CLIPPED", "min_trust",
             "mean_J_exact", "incident_sum_rho", "budget_slack",
             "budget_active", "detected_pred_exact", "dose_resid",
             "baseline_dose_dev", "alpha_mixture", "cert_gap",
             "cert_gap_v3_loop", "theta_degenerate", "mean_load_pred",
             "rho_5", "rho_50", "rho_95", "rho_max", "design_sha"]
M1_FIELDNAMES = list(pilot_s1.FIELDNAMES) + MECH_COLS + PROV_COLS + DISC_COLS
RESUME_COLS = ["pattern", "rho_bar", "nu", "M", "seed", "image", "arm"]


# =========================================================================== #
# frozen 52-row pre-scan matrix (unit-mean-load rows)                          #
# =========================================================================== #
def prescan_matrix(side=SIDE):
    """(52 x n) pre-scan rows. R13 Section 2 REPLACEMENT: the balanced,
    scene-independent construction (32 paired-4x4-block rows + 16 8x8 rows +
    4 quadrant rows) whose row average is exactly the all-ones vector and
    whose per-pixel cumulative dose is exactly equal. Supersedes the v2/v3
    unbalanced multiscale set for every M1 arm."""
    import oed_design_v4 as v4
    return v4.balanced_prescan_52(side)


# =========================================================================== #
# design-time pre-scan estimate  (A6: one per (image, seed))                   #
# =========================================================================== #
def _img_tag(image):
    return zlib.adler32(image.encode("utf-8")) % 100000


def prescan_estimate(x_true, image, seed, rho=RHO_FAST, nu=NU_TERMINAL):
    """GI on the simulated 52-row pre-scan, clipped nonneg, 4x4-block
    smoothed, sum-normalized (spec item 3). Design stream (63, 9, ...)."""
    P = prescan_matrix()
    T = nu * TAU
    Phi = rho / TAU
    det = Detector(tau=TAU, dark=0.0)
    rng = rng_for(int(seed), 63, 9, _img_tag(image))
    u = P @ x_true
    b, _N = simulate_counts(u, Phi, T, det, rng, sigma_b=0.0)
    ctx = ArmContext(Phi=Phi, det=det, T=T, side=SIDE, sigma_b=0.0,
                     n_iter=200, select_iter=60, pattern_kind="m1prescan",
                     meta={"kind": "m1prescan", "n_physical_rows": M_PRE})
    xh, _info = run_arm("GI", P, b, ctx)
    xh = np.maximum(np.asarray(xh, dtype=float).ravel(), 0.0)
    blk = xh.reshape(SIDE // 4, 4, SIDE // 4, 4).mean(axis=(1, 3))
    xhat = np.repeat(np.repeat(blk, 4, axis=0), 4, axis=1).ravel()
    s = xhat.sum()
    xhat = xhat / s if s > 0 else np.full(N_PIX, 1.0 / N_PIX)
    # A13: R10's chain rule is g = a / max(a.xhat, floor); the v3 machinery
    # divides by support loads, so a real (clipped) GI estimate with zero
    # blocks would produce zero-load atoms -> singular design. Floor every
    # pixel at 5% of the uniform level, then renormalize. Design-time only;
    # the reconstruction data path never sees this estimate.
    xhat = np.maximum(xhat, 0.05 / N_PIX)
    return xhat / xhat.sum()


# =========================================================================== #
# ridge-target cache (global per nu)                                           #
# =========================================================================== #
def ridge_targets(nu_list=NU_FULL):
    """rho_R(nu) records for the dwell grid, cached to a JSON side file."""
    os.makedirs(DESIGN_DIR, exist_ok=True)
    cache = {}
    if os.path.exists(RIDGE_JSON):
        with open(RIDGE_JSON) as f:
            cache = json.load(f)
    changed = False
    for nu in nu_list:
        key = "%g" % nu
        if key not in cache:
            t0 = time.time()
            rr = v3.ridge_target(int(nu))
            rr["compute_s"] = round(time.time() - t0, 1)
            cache[key] = rr
            changed = True
            print("[m1 ridge] nu=%g rho_R=%.4f clip=%s (%.1fs)"
                  % (nu, rr.get("rho_R_production", float("nan")),
                     rr.get("ridge_clip_reason"), rr["compute_s"]), flush=True)
    if changed:
        with open(RIDGE_JSON, "w") as f:
            json.dump(cache, f, indent=1, sort_keys=True)
    return cache


# =========================================================================== #
# design caches                                                                #
# =========================================================================== #
def _design_path(image, seed, arm):
    return os.path.join(DESIGN_DIR, "%s_%d_%s.npz" % (image, int(seed), arm))


def _atoms_arrays(atoms):
    names = np.array([a[0] for a in atoms])
    t = np.array([a[1] for a in atoms], dtype=np.int64)
    p = np.array([a[2] for a in atoms], dtype=np.float64)
    rho = np.array([a[3] for a in atoms], dtype=np.float64)
    cnt = np.array([a[4] for a in atoms], dtype=np.int64)
    return names, t, p, rho, cnt


def _sha_arrays(*arrays):
    h = hashlib.sha256()
    for a in arrays:
        h.update(np.ascontiguousarray(a).tobytes())
    return h.hexdigest()


def rebuild_A_main(names, t, p, rho, cnt, xhat, side=SIDE):
    """Deterministically rebuild the design's 972 main rows from its atom
    list + xhat (mirrors v2._build_geometries' G4 pipeline exactly):
    w = renorm(clip((xh_on/mean)^p, 1/4, 4)); g = w/(w.xh_on); row = rho*g."""
    n = side * side
    shapes = v1._default_shapes(side)
    rows = np.zeros((int(cnt.sum()), n))
    r = 0
    for nm, t_, p_, rho_, c_ in zip(names, t, p, rho, cnt):
        sidx = v1._shape_support(shapes[str(nm)], side)[int(t_)]
        xh_on = xhat[sidx]
        w0 = (xh_on / xh_on.mean()) ** float(p_)
        w1 = np.clip(w0, 0.25, 4.0)
        w = w1 / w1.mean()
        g = w / (w * xh_on).sum()
        row = np.zeros(n)
        row[sidx] = float(rho_) * g
        rows[r:r + int(c_)] = row
        r += int(c_)
    return rows


def posthoc_cert(xhat, atoms, levels, budget=RHO_FAST, side=SIDE):
    """A14: post-hoc constrained-polytope FW certificate of the DEPLOYED
    rounded design (this is the ledger's cert_gap; the v3 in-loop Lagrangian
    gap is kept as provenance only, because its theta comes from a damped
    level-assignment bisection that degenerates to a 1e18 sentinel on real
    scenes -> gap = inf, which certifies nothing).

    gap = (LMO - <xi, d>) / <xi, d> with d(a) = nu*J(rho_a) * g_a^T V^-1 g_a,
    V = V0 + sum_rows H(row) at the deployed counts, and the LMO taken over
    the R11-feasible polytope {xi in simplex, <rho, xi> <= budget}: the best
    single atom with rho <= budget or the best two-level boundary mixture.
    By concavity of log det, gap upper-bounds the relative directional
    improvement available inside the frozen dictionary + budget (dose
    constraints are NOT priced -- same caveat as v3, stated in the ledger).

    Also returns the scatter16@0.6-equivalent BASELINE dose deviation, the
    evidence field for the A15 finding (on real scenes the +/-5% per-pixel
    dose band is violated by the load-normalized BASELINE itself: equal
    detector load makes physical dose anti-correlated with scene brightness,
    so G5-as-frozen is scene-dependent-infeasible, flagged for R13)."""
    n = side * side
    shapes = v1._default_shapes(side)
    levels = [float(l_) for l_ in levels]
    rho_arr = np.asarray(levels)
    IDX, GVAL, names_g, t_g, p_g, _A0, _g3, _g4 = v2._build_geometries(
        xhat, shapes, (0, 1), side, (0.25, 4.0), np.inf, rho_arr)
    sumg = GVAL.sum(axis=1)
    geo_ok = GVAL.max(axis=1) * (16.0 / sumg) <= 4.0 + 1e-9   # A2 shape guard
    nuJ = np.array([NU_TERMINAL * v3.kernel_eval(l_, int(NU_TERMINAL))["J_exact"]
                    for l_ in levels])
    # V = V0 + deployed information (one exposure per row)
    V0, _ = v2.build_V0(xhat, NU_TERMINAL, RHO_FAST, v2.get_J_source(), side)
    V = V0.copy()
    lvl_of = {round(l_, 12): i for i, l_ in enumerate(levels)}
    deployed = []                                   # (sidx, g_on, level, count)
    for (nm, t_, p_, rho_, c_) in atoms:
        sidx = v1._shape_support(shapes[str(nm)], side)[int(t_)]
        xh_on = xhat[sidx]
        w0 = (xh_on / xh_on.mean()) ** float(p_)
        w = np.clip(w0, 0.25, 4.0)
        w = w / w.mean()
        g_on = w / (w * xh_on).sum()
        li = lvl_of[round(float(rho_), 12)]
        V[np.ix_(sidx, sidx)] += c_ * nuJ[li] * np.outer(g_on, g_on)
        deployed.append((sidx, g_on, li, int(c_)))
    eps = 1e-9 * np.trace(V) / n
    V[np.diag_indices(n)] += eps
    Vinv = np.linalg.inv(V)
    sub = Vinv[IDX[:, :, None], IDX[:, None, :]]
    dgeo = np.einsum('ai,aij,aj->a', GVAL, sub, GVAL)
    D = np.where(geo_ok[:, None], dgeo[:, None] * nuJ[None, :], -np.inf)
    M_dep = sum(c_ for _, _, _, c_ in deployed)
    dbar = sum((c_ / M_dep) * nuJ[li] *
               float(g_on @ Vinv[np.ix_(sidx, sidx)] @ g_on)
               for (sidx, g_on, li, c_) in deployed)
    # LMO over the budget polytope
    best = -np.inf
    dmax_l = D.max(axis=0)                          # best atom per level
    for i, l_ in enumerate(levels):
        if l_ <= budget + 1e-12:
            best = max(best, dmax_l[i])
    for i, lo in enumerate(levels):
        for j, hi in enumerate(levels):
            if lo <= budget < hi:
                wlo = (hi - budget) / (hi - lo)
                best = max(best, wlo * dmax_l[i] + (1 - wlo) * dmax_l[j])
    gap = (best - dbar) / dbar if dbar > 0 else float("inf")
    # baseline dose deviation (scatter16 p=0, level nearest 0.6, uniform)
    sel = (np.asarray(names_g) == "scatter16") & (p_g == 0.0) & geo_ok
    l6 = int(np.argmin(np.abs(rho_arr - 0.6)))
    w_geo = np.where(sel, rho_arr[l6] / max(sel.sum(), 1), 0.0)
    dose_fix = np.bincount(IDX.ravel(), weights=(w_geo[:, None] * GVAL).ravel(),
                           minlength=n)
    base_dev = float(np.abs(dose_fix / dose_fix.mean() - 1.0).max())
    return {"gap": float(gap), "dbar": float(dbar), "lmo": float(best),
            "baseline_dose_dev": base_dev}


def _oed_disclosure(out, ridge_rec, xhat):
    g5 = out["guards"]["G5"]
    bud = out["guards"]["budget"]
    loads = np.concatenate([[a[3]] * a[4] for a in out["atoms"]])
    v3_gap = out["cert"]["final_gap_adjusted"]
    cert = posthoc_cert(xhat, out["atoms"], out["cert"]["levels"])
    return {
        "rho_R_production": ridge_rec["rho_R_production"],
        "ridge_clip_reason": ridge_rec["ridge_clip_reason"],
        "incident_sum_rho": bud["incident_sum_rho"],
        "budget_slack": bud["slack"],
        "budget_active": bool(bud["slack"]
                              <= 1e-6 * max(bud["cap_sum_rho"], 1.0)),
        "detected_pred_exact": bud["detected_pred_exact"],
        "dose_resid": g5["max_rel_dev_rounded"],
        "baseline_dose_dev": cert["baseline_dose_dev"],
        "alpha_mixture": (g5["alpha_mixture"] if g5["alpha_mixture"]
                          is not None else ""),
        "cert_gap": cert["gap"],
        "cert_gap_v3_loop": (v3_gap if np.isfinite(v3_gap) else "inf"),
        "theta_degenerate": bool(not np.isfinite(v3_gap)),
        "mean_load_pred": float(loads.mean()),
        "rho_5": float(np.percentile(loads, 5)),
        "rho_50": float(np.percentile(loads, 50)),
        "rho_95": float(np.percentile(loads, 95)),
        "rho_max": float(loads.max()),
        "mean_J_exact": float(np.mean(
            [v3.kernel_eval(r_, int(NU_TERMINAL))["J_exact"]
             for r_ in sorted(set(loads.tolist()))])),
        "levels": out["cert"]["levels"],
        "G5_pass_direct": g5["pass"],
        "cert_note": ("cert_gap = post-hoc constrained-polytope FW gap of "
                      "the DEPLOYED rounded design (A14); dose constraints "
                      "not priced. cert_gap_v3_loop = v3 in-loop Lagrangian "
                      "gap (inf when the damped level-assignment theta "
                      "degenerated; provenance only)."),
    }


def build_design(image, seed, arm, x_true, force=False):
    """Build (or load) the cached design for one (image, seed, arm)."""
    path = _design_path(image, seed, arm)
    if os.path.exists(path) and not force:
        return dict(np.load(path, allow_pickle=False)), path
    os.makedirs(DESIGN_DIR, exist_ok=True)
    t0 = time.time()
    xhat = prescan_estimate(x_true, image, seed)
    ridge_rec = ridge_targets([NU_TERMINAL])["%g" % NU_TERMINAL]

    if arm == "MATCH1":
        arrs = dict(xhat=xhat)
        disc = {"kind": "MATCH1 xhat cache (matched-intensity weights)",
                "design_wall_s": round(time.time() - t0, 1)}
    elif arm == "RIDGE-FIXED":
        nu_grid, rho_fast, mult, req, clipped, bis = [], [], [], [], [], []
        detail = {}
        for nu in NU_FULL:
            rf = v3.ridge_fixed_design(xhat, int(nu), side=SIDE)
            nu_grid.append(float(nu))
            rho_fast.append(rf["achieved_mean_load"])
            mult.append(rf["global_multiplier"])
            req.append(rf["requested_mean_load"])
            clipped.append(bool(rf["RIDGE_GUARD_CLIPPED"]))
            bis.append(bool(rf["bisection_clip_applied"]))
            detail["%g" % nu] = {
                "ridge": {k_: rf["ridge"][k_] for k_ in
                          ("rho_ridge_exact_unconstrained", "rho_R_production",
                           "ridge_clip_reason", "J_exact_at_target",
                           "J_q_at_target", "p_ceiling_scalar")},
                "mean_p_ceil": rf["mean_p_ceil"], "max_p_ceil": rf["max_p_ceil"],
                "min_trust": rf["min_trust"], "mean_J_exact": rf["mean_J_exact"],
                "rho_quantiles": rf["rho_quantiles"]}
        arrs = dict(xhat=xhat, nu_grid=np.array(nu_grid),
                    rho_fast=np.array(rho_fast), multiplier=np.array(mult),
                    requested=np.array(req),
                    guard_clipped=np.array(clipped, dtype=np.int64),
                    bisection_applied=np.array(bis, dtype=np.int64))
        disc = {"per_nu": detail, "design_wall_s": round(time.time() - t0, 1)}
    else:
        if arm == "OED-DT":
            out = v3.design_v3(xhat, nu=NU_TERMINAL, M_rows=M_MAIN,
                               budget_mean=RHO_FAST)
        elif arm == "OED-EQLOAD":
            orig = v3.load_palette

            def _single(nu, ridge=None):                      # A11 scoped patch
                pal = orig(nu, ridge)
                recs = [r for r in pal["records"]
                        if abs(r["rho"] - RHO_FAST) < 1e-9] or pal["records"]
                return {"nu": pal["nu"], "levels": [RHO_FAST],
                        "records": recs, "ridge": pal["ridge"]}

            v3.load_palette = _single
            try:
                out = v3.design_v3(xhat, nu=NU_TERMINAL, M_rows=M_MAIN,
                                   budget_mean=RHO_FAST)
            finally:
                v3.load_palette = orig
        else:
            raise ValueError("unknown design arm %r" % (arm,))
        names, t, p, rho, cnt = _atoms_arrays(out["atoms"])
        assert int(cnt.sum()) == M_MAIN
        loads = np.concatenate([[r_] * int(c_) for r_, c_ in zip(rho, cnt)])
        arrs = dict(xhat=xhat, atom_shape=names.astype("U16"), atom_t=t,
                    atom_p=p, atom_rho=rho, atom_count=cnt,
                    mean_load_pred=np.array(float(loads.mean())))
        disc = _oed_disclosure(out, ridge_rec, xhat)
        disc["design_wall_s"] = round(time.time() - t0, 1)

    sha = _sha_arrays(*[arrs[k] for k in sorted(arrs)])
    arrs["disclosure_json"] = np.array(json.dumps(disc, sort_keys=True))
    arrs["sha256"] = np.array(sha)
    np.savez_compressed(path, **arrs)
    print("[m1 design] %s seed=%d %s -> %s (%.1fs, sha %s..)"
          % (image, seed, arm, os.path.basename(path),
             time.time() - t0, sha[:12]), flush=True)
    return dict(np.load(path, allow_pickle=False)), path


def design_disclosure(image, seed, arm, nu=None):
    """Disclosure column dict for CSV injection ('' for fixed arms)."""
    blank = {c: "" for c in DISC_COLS}
    if arm not in DESIGN_ARMS:
        return blank
    path = _design_path(image, seed, arm)
    if not os.path.exists(path):
        return blank
    d = dict(np.load(path, allow_pickle=False))
    disc = json.loads(str(d["disclosure_json"]))
    out = dict(blank)
    out["design_sha"] = str(d["sha256"])[:16]
    if arm == "RIDGE-FIXED":
        key = "%g" % (nu if nu is not None else NU_TERMINAL)
        per = disc["per_nu"].get(key, {})
        rid = per.get("ridge", {})
        i = list(d["nu_grid"]).index(float(key)) if float(key) in d["nu_grid"] \
            else None
        out.update({
            "rho_R_production": rid.get("rho_R_production", ""),
            "ridge_clip_reason": rid.get("ridge_clip_reason", ""),
            "requested_mean_load": (float(d["requested"][i])
                                    if i is not None else ""),
            "achieved_mean_load": (float(d["rho_fast"][i])
                                   if i is not None else ""),
            "RIDGE_GUARD_CLIPPED": (int(d["guard_clipped"][i])
                                    if i is not None else ""),
            "min_trust": per.get("min_trust", ""),
            "mean_J_exact": per.get("mean_J_exact", ""),
            "rho_5": per.get("rho_quantiles", {}).get("rho_5", ""),
            "rho_50": per.get("rho_quantiles", {}).get("rho_50", ""),
            "rho_95": per.get("rho_quantiles", {}).get("rho_95", ""),
            "rho_max": per.get("rho_quantiles", {}).get("rho_max", "")})
    else:
        for c in DISC_COLS:
            if c in disc:
                out[c] = disc[c]
    return out


# =========================================================================== #
# pattern loader (called from campaign._patterns via the M1 appendix)          #
# =========================================================================== #
def load_m1_pattern(kind, M, n, seed):
    """'m1pat:<ARM>:<image>' -> {A, exposures_per_row, meta}: the 52-row
    pre-scan stacked over the arm's 972 main rows (A7/A9 conventions)."""
    _tag, arm, image = kind.split(":", 2)
    if M != M_TOTAL or n != N_PIX:
        raise ValueError("M1 cells are frozen at M=%d, n=%d (got %d, %d)"
                         % (M_TOTAL, N_PIX, M, n))
    P = prescan_matrix()
    meta = {"kind": kind, "M_signed": M, "n": n, "seed": int(seed),
            "m1_arm": arm, "image": image, "prescan_rows": M_PRE,
            "main_rows": M_MAIN, "nonneg": True,
            "n_physical_rows": M_TOTAL, "total_exposures": M_TOTAL,
            "pixel_mean_target": 1.0}
    if arm in ("SCAT16", "SCAT32"):
        k = 16 if arm == "SCAT16" else 32
        base = make_patterns("sparsek", N_PIX, N_PIX, seed, k=k)
        A_main = base["A"][:M_MAIN]
        meta["main_construction"] = base["meta"]["construction"]
        meta["k"] = k
    elif arm in ("LBLOB16", "RIDGE-FIXED"):
        base = make_patterns("lblob16", N_PIX, N_PIX, seed)
        A_main = base["A"][:M_MAIN]
        meta["main_construction"] = base["meta"]["construction"]
        meta["k"] = 16
    elif arm == "MATCH1":
        # R10 mandatory matched-intensity no-gate arm: LBLOB16 first-972
        # multiset, p=1 weights from the cached pre-scan estimate,
        # G4-clipped [1/4,4] + support-renormalized; no servo.
        path = _design_path(image, seed, arm)
        if not os.path.exists(path):
            raise FileNotFoundError(
                "MATCH1 xhat cache missing: %s (run m1_runner.py --designs)"
                % path)
        d = dict(np.load(path, allow_pickle=False))
        xhat_m = np.asarray(d["xhat"], dtype=float)
        base = make_patterns("lblob16", N_PIX, N_PIX, seed)["A"][:M_MAIN]
        A_main = np.zeros_like(base)
        for s_ in range(M_MAIN):
            idx = np.nonzero(base[s_])[0]
            w0 = xhat_m[idx] / xhat_m[idx].mean()
            w = np.clip(w0, 0.25, 4.0)
            w = w / w.mean()
            A_main[s_, idx] = (N_PIX / 16.0) * w
        meta["main_construction"] = ("MATCH1: lblob16 first-972 x G4-clipped "
                                     "p=1 matched weights from cached xhat")
        meta["k"] = 16
    elif arm in ("OED-DT", "OED-EQLOAD"):
        path = _design_path(image, seed, arm)
        if not os.path.exists(path):
            raise FileNotFoundError(
                "M1 design cache missing: %s (run m1_runner.py --designs)"
                % path)
        d = dict(np.load(path, allow_pickle=False))
        A_rows = rebuild_A_main(d["atom_shape"], d["atom_t"], d["atom_p"],
                                d["atom_rho"], d["atom_count"], d["xhat"])
        A_main = A_rows / float(d["mean_load_pred"])          # A9 unit mean
        meta["main_construction"] = ("v3 design cache %s (sha %s..), rows "
                                     "rebuilt from atoms, mean-load "
                                     "normalized" % (os.path.basename(path),
                                                     str(d["sha256"])[:12]))
        meta["design_sha256"] = str(d["sha256"])
    else:
        raise ValueError("unknown M1 arm in pattern kind %r" % (kind,))
    A = np.vstack([P, A_main])
    assert A.shape == (M_TOTAL, N_PIX)
    meta["exposures_per_row"] = 1
    return {"A": A, "exposures_per_row": 1, "meta": meta}


# =========================================================================== #
# cell grid                                                                    #
# =========================================================================== #
def m1_images():
    import m1_scenes
    return [name for name, _f, _s in m1_scenes._m1_conf_table()]


def fast_rho_for(arm, nu, image, seed):
    if arm != "RIDGE-FIXED":
        return RHO_FAST
    path = _design_path(image, seed, arm)
    if not os.path.exists(path):
        return None                     # design cache required
    d = np.load(path, allow_pickle=False)
    grid = list(d["nu_grid"])
    return float(d["rho_fast"][grid.index(float(nu))])


def m1_cell(arm, image, seed, rho, nu, C0=None, imageset="m1"):
    kind = "m1pat:%s:%s" % (arm, image)
    cid = "M1_%s_%s_s%d_r%.4g_nu%g" % (arm, image, seed, rho, nu)
    return dict(side=SIDE, pattern=kind, rho_bar=float(rho), nu=float(nu),
                M=M_TOTAL, seed=int(seed), arms=["RQL"], images=[image],
                imageset=imageset, tau=TAU, sigma_b=0.0, fista_iters=200,
                select_iter=60, select_rule="discrepancy", use_lpips=False,
                audit=False, C0=C0, cell_id=cid, stage="M1_%s" % arm)


def cells_for_arm(arm, C0, images=None, seeds=None):
    cells = []
    for image in (images or m1_images()):
        for seed in (seeds or SEEDS5):
            for nu in NU_FULL:
                cells.append(m1_cell(arm, image, seed, RHO_SAFE, nu, C0))
            for nu in NU_FULL:
                rf = fast_rho_for(arm, nu, image, seed)
                if rf is None:
                    raise SystemExit(
                        "[m1] RIDGE-FIXED design cache missing for %s seed %d "
                        "(run --designs first)" % (image, seed))
                cells.append(m1_cell(arm, image, seed, rf, nu, C0))
    return cells


# =========================================================================== #
# local resume-safe sweep (modeled on study2_runner.run_mode)                  #
# =========================================================================== #
def _resume_key_row(r):
    return tuple(str(r[c]) for c in RESUME_COLS)


def _expected_row_keys(cell):
    return [(cell["pattern"], str(float(cell["rho_bar"])),
             str(float(cell["nu"])), str(int(cell["M"])),
             str(int(cell["seed"])), img, arm)
            for img in cell["images"] for arm in cell["arms"]]


def run_arm_grid(arm, out_dir=OUT_DIR, images=None, seeds=None):
    from campaign import run_cell
    os.makedirs(out_dir, exist_ok=True)
    C0 = frozen_C0()
    cells = cells_for_arm(arm, C0, images=images, seeds=seeds)
    path = os.path.join(out_dir, "M1_%s_rows.csv" % arm.replace("-", "_"))
    errlog = os.path.join(out_dir, "M1_%s_errors.log" % arm.replace("-", "_"))
    done = set()
    header_present = os.path.exists(path) and os.path.getsize(path) > 0
    if header_present:
        with open(path, newline="") as f:
            rd = csv.DictReader(f)
            if rd.fieldnames is not None and list(rd.fieldnames) != M1_FIELDNAMES:
                raise SystemExit("[m1:%s] incompatible header in %s" % (arm, path))
            for r in rd:
                done.add(_resume_key_row(r))
    t0 = time.time()
    with open(path, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=M1_FIELDNAMES, extrasaction="ignore")
        if not header_present:
            w.writeheader()
            f.flush()
        for ci, cell in enumerate(cells):
            if all(k_ in done for k_ in _expected_row_keys(cell)):
                continue
            tc = time.time()
            try:
                rows = run_cell(cell)
            except Exception:
                with open(errlog, "a") as ef:
                    ef.write("[%s] %s\n%s\n" % (
                        time.strftime("%Y-%m-%d %H:%M:%S"), cell["cell_id"],
                        traceback.format_exc()))
                print("[m1:%s] ERROR %d/%d %s (logged)"
                      % (arm, ci + 1, len(cells), cell["cell_id"]), flush=True)
                continue
            disc = design_disclosure(cell["images"][0], cell["seed"], arm,
                                     nu=cell["nu"])
            n_new = 0
            for r in rows:
                r.update(disc)
                r["m1_arm"] = arm
                r["stage"] = cell["stage"]
                r["cell_id"] = cell["cell_id"]
                k_ = _resume_key_row(r)
                if k_ in done:
                    continue
                w.writerow(r)
                done.add(k_)
                n_new += 1
            f.flush()
            os.fsync(f.fileno())
            print("[m1:%s] cell %d/%d %s -> +%d rows (%.1fs, total %.0fs)"
                  % (arm, ci + 1, len(cells), cell["cell_id"], n_new,
                     time.time() - tc, time.time() - t0), flush=True)
    print("[m1:%s] SWEEP DONE wall=%.0fs -> %s" % (arm, time.time() - t0, path),
          flush=True)
    return path


# =========================================================================== #
# design-cache builder                                                         #
# =========================================================================== #
def build_all_designs(arms=None, images=None, seeds=None, force=False):
    import campaign
    arms = arms or DESIGN_ARMS
    images = images or m1_images()
    seeds = seeds if seeds is not None else SEEDS5
    imgs = campaign._images(SIDE, "all", imageset="m1")
    ridge_targets(NU_FULL)
    for image in images:
        for seed in seeds:
            for arm in arms:
                build_design(image, seed, arm, imgs[image], force=force)


# =========================================================================== #
# smoke (deliverable item 6)                                                   #
# =========================================================================== #
def smoke():
    import campaign
    from campaign import run_cell
    t0 = time.time()
    C0 = frozen_C0()
    imgs = campaign._images(SIDE, "all", imageset="m1")
    image = m1_images()[0]
    seed = 0
    print("[m1 smoke] image=%s seed=%d C0=%s" % (image, seed, C0), flush=True)

    print("\n[m1 smoke] === designs (OED-DT + RIDGE-FIXED) ===", flush=True)
    for arm in ("OED-DT", "RIDGE-FIXED"):
        d, path = build_design(image, seed, arm, imgs[image])
        disc = json.loads(str(d["disclosure_json"]))
        print("  %s cache: %s" % (arm, os.path.basename(path)))
        if arm == "OED-DT":
            for k_ in ("rho_R_production", "ridge_clip_reason",
                       "incident_sum_rho", "budget_slack", "budget_active",
                       "detected_pred_exact", "dose_resid",
                       "baseline_dose_dev", "alpha_mixture", "cert_gap",
                       "cert_gap_v3_loop", "theta_degenerate",
                       "mean_load_pred", "rho_5", "rho_50",
                       "rho_95", "rho_max", "mean_J_exact", "levels",
                       "G5_pass_direct", "design_wall_s"):
                print("    %-22s = %s" % (k_, disc.get(k_)))
        else:
            for nu in (20.0, 2000.0):
                per = disc["per_nu"]["%g" % nu]
                i = list(d["nu_grid"]).index(nu)
                print("    nu=%-5g rho_R=%.4f clip=%s requested=%.4f "
                      "achieved=%.4f mult=%.4f GUARD_CLIPPED=%d min_trust=%.4f"
                      % (nu, per["ridge"]["rho_R_production"],
                         per["ridge"]["ridge_clip_reason"],
                         float(d["requested"][i]), float(d["rho_fast"][i]),
                         float(d["multiplier"][i]),
                         int(d["guard_clipped"][i]), per["min_trust"]))

    print("\n[m1 smoke] === cells through campaign.run_cell (imageset m1) ===",
          flush=True)
    cells = []
    for nu in (20.0, 2000.0):
        cells.append(("OED-DT safe", m1_cell("OED-DT", image, seed,
                                             RHO_SAFE, nu, C0)))
        cells.append(("OED-DT fast", m1_cell("OED-DT", image, seed,
                                             RHO_FAST, nu, C0)))
    rf2000 = fast_rho_for("RIDGE-FIXED", 2000.0, image, seed)
    cells.append(("RIDGE-FIXED fast", m1_cell("RIDGE-FIXED", image, seed,
                                              rf2000, 2000.0, C0)))
    cells.append(("SCAT16 fast", m1_cell("SCAT16", image, seed,
                                         RHO_FAST, 2000.0, C0)))
    ok = True
    for label, cell in cells:
        tc = time.time()
        rows = run_cell(cell)
        arm = cell["stage"].replace("M1_", "")
        disc = design_disclosure(image, seed, arm, nu=cell["nu"])
        r = rows[0]
        print("  %-17s nu=%-5g rho=%.4g  PSNR_rad=%s  PSNR=%s  "
              "mean_counts=%s  (%.1fs)"
              % (label, cell["nu"], cell["rho_bar"], r["PSNR_rad"], r["PSNR"],
                 r["mean_counts"], time.time() - tc), flush=True)
        keys = [k_ for k_ in DISC_COLS if disc.get(k_) not in ("", None)]
        if keys:
            print("      disclosure: " + ", ".join(
                "%s=%s" % (k_, disc[k_]) for k_ in keys))
        ok = ok and (r["PSNR_rad"] != "")
    print("\n[m1 smoke] wall=%.1fs  RESULT: %s"
          % (time.time() - t0, "SMOKE PASS" if ok else "*** SMOKE FAIL ***"),
          flush=True)
    return 0 if ok else 1


# =========================================================================== #
def main(argv=None):
    ap = argparse.ArgumentParser(description="M1 method-campaign runner.")
    ap.add_argument("--designs", action="store_true")
    ap.add_argument("--arm", type=str, default=None,
                    help="run one arm's full local grid (%s)" % "|".join(ARMS_ALL))
    ap.add_argument("--arms", type=str, default=None,
                    help="comma list for --designs (default all design arms)")
    ap.add_argument("--images", type=str, default=None)
    ap.add_argument("--seeds", type=str, default=None)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    a = ap.parse_args(sys.argv[1:] if argv is None else argv)
    images = a.images.split(",") if a.images else None
    seeds = [int(s) for s in a.seeds.split(",")] if a.seeds else None
    if a.smoke:
        return smoke()
    if a.designs:
        build_all_designs(arms=(a.arms.split(",") if a.arms else None),
                          images=images, seeds=seeds, force=a.force)
        return 0
    if a.arm:
        if a.arm not in ARMS_ALL:
            raise SystemExit("unknown arm %r" % a.arm)
        run_arm_grid(a.arm, images=images, seeds=seeds)
        return 0
    ap.print_help()
    return 2


if __name__ == "__main__":
    sys.exit(main())
