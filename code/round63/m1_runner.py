"""M1 campaign runner — R17 re-architecture (OPERATIVE).

Authority: docs/ROUND63_GPT_ROUND17_RULING_RAW.md (issue #9) via
docs/ROUND63_METHOD_SPEC_M1_R17_AMENDMENT.md. Supersedes the pre-R17
runner (the OED-DT / OED-EQLOAD / MATCH1 / RIDGE-FIXED arms and the
retired adaptive speed gate are removed from every production path; the frozen
rule m=0 => ADAPTIVE_COLLAPSE_UNDER_GUARDS != PASS lives in
oed_design_v5.path_feasible_alpha).

Architecture (amendment §B):
  Deployed spatial design  = the balanced exact 972-row SCAT32 multiset
  (oed_design_v5.fixed_dose_scat32; scene-independent; SHA-frozen),
  identical in all operating-point comparisons; plus the common balanced
  52-row pre-scan charged to every arm (52 + 972 = 1024 rows per cell).

  Endpoint-carrying modes:  SCAT32-SAFE (global load 0.05),
  SCAT32-060 (0.60), RIDGE-SCAT32 (one global source multiplier per dwell,
  calibrated at run time from the common pre-scan estimate so the predicted
  mean main-pattern load hits the exact production ridge rho_R(nu), then
  only the frozen global safety clip; no per-pattern servo).
  Context arms (no gate, descriptive): SCAT16, LBLOB16 at load 0.60.

Endpoints (amendment §C): RIDGE_OPERATING_PASS (primary, paired terminal-
dwell), RIDGE_SPEED_PASS (nine-dwell Q90 secondary), DOSE_SAFE_CERT_PASS
(480-cell expanded-class certificate) — emitted by m1_analyze_r17.

Certificate cells (amendment §C.3) run through run_cert_cell — shardable
via the existing manifest machinery (campaign.run_cell dispatches cells
carrying m1_cert=True to run_cert_cell; one CSV row per cell).

HARD RULE (R17 §5): no confirmatory (633000+) image may be loaded and no
confirmatory pre-scan count generated before the m1-freeze launch. Guarded
here: imageset "m1" requires the M1_FREEZE_LAUNCHED environment gate (set
by the post-freeze Colab launcher); DEV ("m1_dev") and synthetic paths are
always allowed.

Interpretations logged for the refreeze audit (I1..I5):
  I1 context arms run at load 0.60 across all nine dwells (descriptive
     occupancy-ladder columns; the amendment fixes only the load).
  I2 ridge calibration uses ONE common pre-scan estimate per (image, seed)
     at the (0.60, nu=2000) anchor for all dwells (design-granularity
     precedent; §2.1 "calibrated from the common pre-scan estimate").
  I3 certificate cells use their own (nu, b)-specific pre-scan realization
     (stream rng_for(seed, 63, 9, img, nu, b)) — §4.5 "every actual
     pre-scan realization" enumerates exactly 24x5x2x2 = 480 realizations.
  I4 the declared r=200 task subspace is the top-200 eigenspace of the
     DEPLOYED SCAT32 information matrix (pre-scan V0 + 972 rows) at the
     cell's own budget corner — the campaign's reference design under R17;
     eps0 = 1e-9 tr(B^T V B)/r as frozen.
  I5 "PSNR_rad mean over the five frozen measurement seeds" pairs seeds
     via identical pattern matrices and identical noise streams keyed by
     the shared pattern kind prefix (only the global multiplier differs).
"""
import argparse
import csv
import hashlib
import json
import math
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

import oed_design_v4 as v4                  # noqa: E402
import oed_design_v5 as v5                  # noqa: E402
from gi_core.utils import rng_for           # noqa: E402
from patterns import make_patterns          # noqa: E402
from physics import Detector, simulate_counts  # noqa: E402
from solvers import ArmContext, run_arm     # noqa: E402
from select_eta import frozen_C0            # noqa: E402
import pilot_s1                             # noqa: E402

# ---- frozen geometry ------------------------------------------------------ #
SIDE = 32
N_PIX = SIDE * SIDE
M_TOTAL, M_PRE, M_MAIN = 1024, 52, 972
TAU = 50e-9
NU_FULL = [5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0, 1000.0, 2000.0]
NU_TERMINAL = 2000.0
SEEDS5 = [0, 1, 2, 3, 4]

# ---- R17 arm roster ------------------------------------------------------- #
MODE_ARMS = ["SCAT32-SAFE", "SCAT32-060", "RIDGE-SCAT32"]
CONTEXT_ARMS = ["SCAT16", "LBLOB16"]
ARMS_ALL = MODE_ARMS + CONTEXT_ARMS
MODE_LOADS = {"SCAT32-SAFE": 0.05, "SCAT32-060": 0.60,
              "SCAT16": 0.60, "LBLOB16": 0.60}
RETIRED_ARMS = {"OED-DT", "OED-EQLOAD", "MATCH1", "RIDGE-FIXED"}
CERT_NUS = [200.0, 2000.0]
CERT_BUDGETS = [0.05, 0.60]

OUT_DIR = os.path.join(ROOT, "results", "round63_m1")
DESIGN_DIR = os.path.join(OUT_DIR, "designs")
RIDGE_JSON = os.path.join(DESIGN_DIR, "ridge_targets_r17.json")
SCAT32_NPZ = os.path.join(DESIGN_DIR, "scat32_deployed.npz")

MECH_COLS = ["cnr", "C_u", "Gamma", "S_det", "S_inc", "k_occupancy"]
PROV_COLS = ["m1_arm", "stage", "cell_id", "descriptive_context"]
DISC_COLS = ["rho_R_production", "ridge_clip_reason", "requested_mean_load",
             "achieved_mean_load", "RIDGE_GUARD_CLIPPED",
             "incident_dose_ratio", "detected_count_ratio",
             "load_q5", "load_q50", "load_q95", "load_max",
             "physical_peak", "ceiling_fraction_pred", "deployed_sha"]
M1_FIELDNAMES = list(pilot_s1.FIELDNAMES) + MECH_COLS + PROV_COLS + DISC_COLS
RESUME_COLS = ["pattern", "rho_bar", "nu", "M", "seed", "image", "arm"]

CERT_FIELDNAMES = ["image", "seed", "nu", "b", "G_full_per_r", "cell_pass",
                   "status", "primal_feasible", "theta", "n_active_mu",
                   "MU_CAP_ACTIVE", "budget_viol", "dose_excess",
                   "comp_budget", "comp_dose", "deployed_mean_load",
                   "deployed_dose_dev", "deployed_peak", "d_cert_sha",
                   "wall_s", "cell_id", "stage"]


def _confirmatory_guard(imageset):
    """R17 §5 hard rule: confirmatory scenes only after the freeze launch."""
    if imageset == "m1" and not os.environ.get("M1_FREEZE_LAUNCHED"):
        raise RuntimeError(
            "R17 launch policy: confirmatory imageset 'm1' is locked until "
            "the m1-freeze tag launch sets M1_FREEZE_LAUNCHED (use 'm1_dev' "
            "for development).")


# ---- balanced pre-scan + estimate ----------------------------------------- #
def prescan_matrix(side=SIDE):
    """The R13 Sec-2 balanced 52-row pre-scan (exact identities)."""
    return v4.balanced_prescan_52(side)


def _img_tag(image):
    return zlib.adler32(image.encode("utf-8")) % 100000


def prescan_estimate(x_true, image, seed, rho, nu, per_cell=False):
    """GI on the simulated balanced 52-row pre-scan, clipped nonneg,
    4x4-block smoothed, floored (A13) and sum-normalized. Stream: (63, 9,
    img) for the common per-(image, seed) estimate; (63, 9, img, nu, b) for
    per-cell certificate realizations (I3)."""
    P = prescan_matrix()
    T = nu * TAU
    Phi = rho / TAU
    det = Detector(tau=TAU, dark=0.0)
    if per_cell:
        rng = rng_for(int(seed), 63, 9, _img_tag(image), int(nu),
                      int(round(rho * 1000)))
    else:
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
    xhat = np.maximum(xhat, 0.05 / N_PIX)
    return xhat / xhat.sum()


# ---- deployed SCAT32 multiset (frozen + hashed) --------------------------- #
_SCAT32_CACHE = {}


def deployed_scat32():
    """The balanced exact 972-row SCAT32 multiset (R17 §2.1) + SHA256."""
    if "rows" in _SCAT32_CACHE:
        return _SCAT32_CACHE["rows"], _SCAT32_CACHE["sha"]
    if os.path.exists(SCAT32_NPZ):
        d = np.load(SCAT32_NPZ)
        rows, sha = d["rows"], str(d["sha256"])
    else:
        rows, dev = v5.fixed_dose_scat32(SIDE)
        sha = hashlib.sha256(np.ascontiguousarray(rows).tobytes()).hexdigest()
        os.makedirs(DESIGN_DIR, exist_ok=True)
        np.savez_compressed(SCAT32_NPZ, rows=rows, sha256=np.array(sha),
                            dose_dev=np.array(dev))
    _SCAT32_CACHE["rows"] = rows
    _SCAT32_CACHE["sha"] = sha
    return rows, sha


# ---- exact ridge targets (A1-amended kernel) ------------------------------ #
def ridge_targets(nu_list=NU_FULL):
    os.makedirs(DESIGN_DIR, exist_ok=True)
    cache = {}
    if os.path.exists(RIDGE_JSON):
        with open(RIDGE_JSON) as f:
            cache = json.load(f)
    changed = False
    for nu in nu_list:
        key = "%g" % nu
        if key not in cache:
            rr = v4.ridge_target4(int(nu))
            cache[key] = rr
            changed = True
            print("[m1 ridge] nu=%g rho_R=%.4f clip=%s"
                  % (nu, rr.get("rho_R_production", float("nan")),
                     rr.get("ridge_clip_reason")), flush=True)
    if changed:
        with open(RIDGE_JSON, "w") as f:
            json.dump(cache, f, indent=1, sort_keys=True)
    return cache


# ---- RIDGE-SCAT32 runtime calibration ------------------------------------- #
def _ridge_cache_path(image, seed):
    return os.path.join(DESIGN_DIR, "ridge_scat32_%s_%d.npz" % (image, seed))


def ridge_scat32_calibration(image, seed, imageset="m1_dev", x_true=None):
    """Per (image, seed): the common pre-scan estimate (I2) + per-dwell
    global multiplier so predicted mean main load = rho_R(nu), then ONLY
    the frozen global safety clip (kernel-grid guards). Cached npz."""
    path = _ridge_cache_path(image, seed)
    if os.path.exists(path):
        return dict(np.load(path, allow_pickle=False))
    _confirmatory_guard(imageset)
    if x_true is None:
        import campaign
        x_true = campaign._images(SIDE, "all", imageset=imageset)[image]
    xhat = prescan_estimate(x_true, image, seed, MODE_LOADS["SCAT32-060"],
                            NU_TERMINAL)
    rows, sha = deployed_scat32()
    base_load = rows @ xhat                       # relative predicted loads
    rt = ridge_targets(NU_FULL)
    nu_grid, mult, req, ach, clip, flags = [], [], [], [], [], []
    for nu in NU_FULL:
        kg = v5._kernel_grids(nu)
        rho_R = rt["%g" % nu]["rho_R_production"]
        m0 = rho_R / float(base_load.mean())

        def ok(m_):
            loads = m_ * base_load
            ll = np.log(np.clip(loads, kg["lo"], kg["hi"]))
            pc = np.interp(ll, kg["lg"], kg["pceil"])
            return (pc.mean() <= v4.CEIL_TARGET and pc.max() <= v4.CEIL_ATOM
                    and np.interp(ll, kg["lg"], kg["eff"]).min() >= v4.EFF_MIN
                    and np.interp(ll, kg["lg"], kg["bias"]).max()
                    <= v4.BIAS_MAX)

        clipped = False
        m_ = m0
        if not ok(m_):
            clipped = True
            lo, hi = 0.0, m_
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                if ok(mid):
                    lo = mid
                else:
                    hi = mid
            m_ = lo
        a_mean = float(m_ * base_load.mean())
        nu_grid.append(float(nu))
        mult.append(float(m_))
        req.append(float(rho_R))
        ach.append(a_mean)
        clip.append(bool(clipped))
        flags.append(bool(a_mean < 0.90 * rho_R))
    arrs = dict(xhat=xhat, nu_grid=np.array(nu_grid),
                multiplier=np.array(mult), requested=np.array(req),
                achieved=np.array(ach),
                clip_applied=np.array(clip, dtype=np.int64),
                guard_clipped=np.array(flags, dtype=np.int64),
                deployed_sha=np.array(sha))
    os.makedirs(DESIGN_DIR, exist_ok=True)
    np.savez_compressed(path, **arrs)
    return dict(np.load(path, allow_pickle=False))


def ridge_scat32_rho(image, seed, nu, imageset="m1_dev", x_true=None):
    cal = ridge_scat32_calibration(image, seed, imageset, x_true)
    i = list(cal["nu_grid"]).index(float(nu))
    return float(cal["achieved"][i]), cal


# ---- pattern loader (via the campaign M1 appendix) ------------------------ #
def load_m1_pattern(kind, M, n, seed):
    """'m1pat:<ARM>:<image>' -> 52-row balanced pre-scan + the arm's 972
    main rows. R17: the three SCAT32 modes share ONE deployed multiset."""
    _tag, arm, image = kind.split(":", 2)
    if arm in RETIRED_ARMS:
        raise ValueError("arm %r RETIRED_BY_R17 (issue #9); no production "
                         "path may build its patterns" % arm)
    if M != M_TOTAL or n != N_PIX:
        raise ValueError("M1 cells are frozen at M=%d, n=%d" % (M_TOTAL, N_PIX))
    P = prescan_matrix()
    meta = {"kind": kind, "M_signed": M, "n": n, "seed": int(seed),
            "m1_arm": arm, "image": image, "prescan_rows": M_PRE,
            "main_rows": M_MAIN, "nonneg": True,
            "n_physical_rows": M_TOTAL, "total_exposures": M_TOTAL,
            "pixel_mean_target": 1.0}
    if arm in ("SCAT32-SAFE", "SCAT32-060", "RIDGE-SCAT32"):
        rows, sha = deployed_scat32()
        A_main = rows
        meta["main_construction"] = ("deployed balanced exact 972-row "
                                     "SCAT32 multiset sha %s.." % sha[:12])
        meta["deployed_sha256"] = sha
        meta["k"] = 32
    elif arm == "SCAT16":
        A_main = make_patterns("sparsek", N_PIX, N_PIX, seed, k=16)["A"][:972]
        meta["k"] = 16
    elif arm == "LBLOB16":
        A_main = make_patterns("lblob16", N_PIX, N_PIX, seed)["A"][:972]
        meta["k"] = 16
    else:
        raise ValueError("unknown M1 arm %r" % arm)
    A = np.vstack([P, A_main])
    assert A.shape == (M_TOTAL, N_PIX)
    meta["exposures_per_row"] = 1
    return {"A": A, "exposures_per_row": 1, "meta": meta}


# ---- cells ---------------------------------------------------------------- #
def m1_images():
    import m1_scenes
    return [name for name, _f, _s in m1_scenes._m1_conf_table()]


def m1_cell(arm, image, seed, rho, nu, C0=None, imageset="m1"):
    kind = "m1pat:%s:%s" % (arm, image)
    cid = "M1_%s_%s_s%d_nu%g" % (arm, image, seed, nu)
    cell = dict(side=SIDE, pattern=kind, rho_bar=float(rho), nu=float(nu),
                M=M_TOTAL, seed=int(seed), arms=["RQL"], images=[image],
                imageset=imageset, tau=TAU, sigma_b=0.0, fista_iters=200,
                select_iter=60, select_rule="discrepancy", use_lpips=False,
                audit=False, C0=C0, cell_id=cid, stage="M1_%s" % arm)
    if arm == "RIDGE-SCAT32":
        cell["m1_ridge_dynamic"] = True       # rho_bar resolved at run time
    return cell


def cert_cell(image, seed, nu, b, imageset="m1"):
    return {"m1_cert": True, "image": image, "seed": int(seed),
            "nu": float(nu), "b": float(b), "imageset": imageset,
            "images": [image],
            "cell_id": "M1CERT_%s_s%d_nu%g_b%g" % (image, seed, nu, b),
            "stage": "M1_CERT"}


def cells_for_arm(arm, C0, images=None, seeds=None, imageset="m1"):
    cells = []
    for image in (images or m1_images()):
        for seed in (seeds or SEEDS5):
            for nu in NU_FULL:
                rho = MODE_LOADS.get(arm, -1.0)   # ridge: sentinel, dynamic
                cells.append(m1_cell(arm, image, seed, rho, nu, C0,
                                     imageset=imageset))
    return cells


def cert_cells_confirmatory():
    return [cert_cell(img, s, nu, b)
            for img in m1_images() for s in SEEDS5
            for nu in CERT_NUS for b in CERT_BUDGETS]


# ---- certificate cell runner (dispatched from campaign.run_cell) ---------- #
def run_cert_cell(cell):
    """One R17 §C.3 certificate cell -> [one CSV row dict]."""
    t0 = time.time()
    imageset = cell.get("imageset", "m1")
    _confirmatory_guard(imageset)
    import campaign
    image, seed = cell["image"], int(cell["seed"])
    nu, b = float(cell["nu"]), float(cell["b"])
    x_true = campaign._images(SIDE, "all", imageset=imageset)[image]
    xhat = prescan_estimate(x_true, image, seed, b, nu, per_cell=True)  # I3
    rows, sha_dep = deployed_scat32()
    # I4: declared r=200 subspace from the DEPLOYED design's info matrix
    V_full = v4.info_matrix_full(rows, xhat, int(nu), b, P=prescan_matrix())
    B, eps0, _tr = v4.subspace_from_fixedstar(V_full)
    ctx = v5.setup_ctx_cert(xhat, nu, b, B, eps0, SIDE)
    out = v5.cert_deployed_rows(ctx, rows, b)
    row = {"image": image, "seed": seed, "nu": nu, "b": b,
           "G_full_per_r": (out["G_full"] / ctx["r"]
                            if np.isfinite(out.get("G_full", np.inf))
                            else ""),
           "cell_pass": bool(out.get("cell_pass", False)),
           "status": out.get("status", ""),
           "primal_feasible": out.get("primal_feasible", ""),
           "theta": out.get("theta", ""),
           "n_active_mu": out.get("n_active_mu", ""),
           "MU_CAP_ACTIVE": out.get("MU_CAP_ACTIVE", ""),
           "budget_viol": out.get("budget_viol", ""),
           "dose_excess": out.get("dose_excess", ""),
           "comp_budget": out.get("comp_budget", ""),
           "comp_dose": out.get("comp_dose", ""),
           "deployed_mean_load": out.get("deployed_mean_load", ""),
           "deployed_dose_dev": out.get("deployed_dose_dev", ""),
           "deployed_peak": out.get("deployed_peak", ""),
           "d_cert_sha": v5.d_cert_sha()[:16],
           "wall_s": round(time.time() - t0, 1),
           "cell_id": cell.get("cell_id", ""),
           "stage": cell.get("stage", "M1_CERT")}
    return [row]


# ---- R17 analyzer --------------------------------------------------------- #
def _family_of(image):
    parts = image.split("_")
    return parts[1] if len(parts) >= 3 else image


def _nested_boot_lb(per_image, families, B=10000, q=2.5, seed_tag=13):
    """10k nested family-stratified bootstrap LB of the median (families
    resampled with replacement; images within family resampled)."""
    fam_map = {}
    for img, v_ in per_image.items():
        fam_map.setdefault(families[img], []).append(v_)
    fams = sorted(fam_map)
    rng = rng_for(0, 63, seed_tag)
    meds = np.empty(B)
    for b_ in range(B):
        fsel = rng.integers(0, len(fams), size=len(fams))
        vals = []
        for fi in fsel:
            arr = fam_map[fams[fi]]
            isel = rng.integers(0, len(arr), size=len(arr))
            vals.extend(arr[k] for k in isel)
        meds[b_] = np.median(vals)
    return float(np.percentile(meds, q))


def _pava(y):
    y = list(y)
    w = [1.0] * len(y)
    i = 0
    while i < len(y) - 1:
        if y[i] > y[i + 1] + 1e-12:
            m = (y[i] * w[i] + y[i + 1] * w[i + 1]) / (w[i] + w[i + 1])
            y[i] = y[i + 1] = m
            w[i] = w[i + 1] = w[i] + w[i + 1]
            i = max(i - 1, 0)
        else:
            i += 1
    return np.array(y)


def _q90_time(nus, rhos, q, target):
    """Optical-time crossing with per-dwell loads (ridge policy: rho(nu))."""
    qf = _pava(q)
    lt = np.log(np.array(nus) * np.array(rhos))
    if qf[-1] < target:
        return None
    return float(np.exp(np.interp(target, qf, lt)))


def m1_analyze_r17(curves, n_images=24, boot_B=10000):
    """The three R17 verdicts from per-image curve data.

    curves: dict image -> {"family": str,
        "safe": (nus, rhos, psnr_by_nu), "ridge": (nus, rhos, psnr_by_nu),
        "q060_terminal": float, "ridge_terminal": float,
        "cert_cells": [bool, ...]}  (means over seeds already applied; ITT:
    a failed cell enters as nonpositive dQ / censored curve)."""
    families = {img: c["family"] for img, c in curves.items()}
    dQ = {img: (c["ridge_terminal"] - c["q060_terminal"]
                if np.isfinite(c["ridge_terminal"])
                and np.isfinite(c["q060_terminal"]) else 0.0)
          for img, c in curves.items()}
    med = float(np.median(list(dQ.values())))
    lb = _nested_boot_lb(dQ, families, B=boot_B, seed_tag=13)
    n_pos = sum(1 for v_ in dQ.values() if v_ > 0)
    operating = bool(med >= 1.0 and lb > 0 and n_pos >= math.ceil(
        0.75 * n_images))
    S = {}
    for img, c in curves.items():
        nus_s, rhos_s, qs = c["safe"]
        nus_f, rhos_f, qf = c["ridge"]
        R = qs[-1] - qs[0]
        tgt = qs[0] + 0.9 * R
        Ts = _q90_time(nus_s, rhos_s, qs, tgt)
        Tf = _q90_time(nus_f, rhos_f, qf, tgt)
        S[img] = (Ts / Tf) if (Ts and Tf) else 0.0    # ITT nonpositive
    medS = float(np.median(list(S.values())))
    lbS = _nested_boot_lb(S, families, B=boot_B, seed_tag=14)
    nS = sum(1 for v_ in S.values() if v_ > 1)
    speed = bool(medS >= 3.0 and lbS > 1 and nS >= math.ceil(0.75 * n_images))
    cert_flags = [bool(f) for c in curves.values()
                  for f in c.get("cert_cells", [])]
    cert = bool(cert_flags and all(cert_flags))
    return {"RIDGE_OPERATING_PASS": operating,
            "RIDGE_SPEED_PASS": speed,
            "DOSE_SAFE_CERT_PASS": cert,
            "median_dQ_dB": med, "dQ_LB2.5": lb, "n_dQ_pos": n_pos,
            "median_S": medS, "S_LB2.5": lbS, "n_S_gt1": nS,
            "n_cert_cells": len(cert_flags),
            "n_cert_pass": sum(cert_flags)}


# ---- resumable local sweep ------------------------------------------------ #
def _resume_key_row(r):
    return tuple(str(r[c]) for c in RESUME_COLS)


def ridge_disclosure(image, seed, nu, imageset="m1_dev"):
    blank = {c: "" for c in DISC_COLS}
    try:
        cal = ridge_scat32_calibration(image, seed, imageset)
    except (RuntimeError, FileNotFoundError):
        return blank
    i = list(cal["nu_grid"]).index(float(nu))
    rt = ridge_targets([nu])["%g" % nu]
    rows, sha = deployed_scat32()
    ach = float(cal["achieved"][i])
    rel = rows @ cal["xhat"]
    loads = ach * rel / float(rel.mean())
    out = dict(blank)
    out.update({
        "rho_R_production": rt["rho_R_production"],
        "ridge_clip_reason": rt["ridge_clip_reason"],
        "requested_mean_load": float(cal["requested"][i]),
        "achieved_mean_load": ach,
        "RIDGE_GUARD_CLIPPED": int(cal["guard_clipped"][i]),
        "incident_dose_ratio": ach / MODE_LOADS["SCAT32-060"],
        "detected_count_ratio": float(
            (nu * loads / (1 + loads)).mean() / (nu * 0.60 / 1.60)),
        "load_q5": float(np.percentile(loads, 5)),
        "load_q50": float(np.percentile(loads, 50)),
        "load_q95": float(np.percentile(loads, 95)),
        "load_max": float(loads.max()),
        "physical_peak": float((rows * (ach / float(rel.mean()))).max()),
        "ceiling_fraction_pred": float(np.mean(
            [v4.p_ceil_exact(l_, int(nu)) for l_ in
             np.percentile(loads, [5, 50, 95])])),
        "deployed_sha": sha[:16]})
    return out


def run_arm_grid(arm, out_dir=OUT_DIR, images=None, seeds=None,
                 imageset="m1"):
    from campaign import run_cell
    _confirmatory_guard(imageset)
    os.makedirs(out_dir, exist_ok=True)
    C0 = frozen_C0()
    cells = cells_for_arm(arm, C0, images=images, seeds=seeds,
                          imageset=imageset)
    path = os.path.join(out_dir, "M1_%s_rows.csv" % arm.replace("-", "_"))
    errlog = os.path.join(out_dir, "M1_%s_errors.log" % arm.replace("-", "_"))
    done = set()
    header_present = os.path.exists(path) and os.path.getsize(path) > 0
    if header_present:
        with open(path, newline="") as f:
            rd = csv.DictReader(f)
            for r in rd:
                done.add(_resume_key_row(r))
    t0 = time.time()
    with open(path, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=M1_FIELDNAMES, extrasaction="ignore")
        if not header_present:
            w.writeheader()
            f.flush()
        for ci, cell in enumerate(cells):
            try:
                rows = run_cell(cell)
            except Exception:
                with open(errlog, "a") as ef:
                    ef.write("[%s] %s\n%s\n" % (
                        time.strftime("%Y-%m-%d %H:%M:%S"), cell["cell_id"],
                        traceback.format_exc()))
                continue
            disc = (ridge_disclosure(cell["images"][0], cell["seed"],
                                     cell["nu"], imageset)
                    if arm == "RIDGE-SCAT32" else {c: "" for c in DISC_COLS})
            n_new = 0
            for r in rows:
                r.update(disc)
                r["m1_arm"] = arm
                r["stage"] = cell["stage"]
                r["cell_id"] = cell["cell_id"]
                r["descriptive_context"] = int(arm in CONTEXT_ARMS)
                k_ = _resume_key_row(r)
                if k_ in done:
                    continue
                w.writerow(r)
                done.add(k_)
                n_new += 1
            f.flush()
            os.fsync(f.fileno())
            print("[m1:%s] cell %d/%d %s -> +%d rows (%.0fs)"
                  % (arm, ci + 1, len(cells), cell["cell_id"], n_new,
                     time.time() - t0), flush=True)
    return path


# ---- CLI ------------------------------------------------------------------ #
def main(argv=None):
    ap = argparse.ArgumentParser(description="M1 runner (R17 architecture).")
    ap.add_argument("--arm", type=str, default=None, help="|".join(ARMS_ALL))
    ap.add_argument("--imageset", type=str, default="m1")
    ap.add_argument("--images", type=str, default=None)
    ap.add_argument("--seeds", type=str, default=None)
    ap.add_argument("--ridge-tables", action="store_true",
                    help="precompute the nine-dwell ridge target table")
    a = ap.parse_args(sys.argv[1:] if argv is None else argv)
    if a.ridge_tables:
        ridge_targets(NU_FULL)
        return 0
    if a.arm:
        if a.arm in RETIRED_ARMS:
            raise SystemExit("arm %r RETIRED_BY_R17" % a.arm)
        if a.arm not in ARMS_ALL:
            raise SystemExit("unknown arm %r" % a.arm)
        run_arm_grid(a.arm,
                     images=a.images.split(",") if a.images else None,
                     seeds=[int(s) for s in a.seeds.split(",")]
                     if a.seeds else None,
                     imageset=a.imageset)
        return 0
    ap.print_help()
    return 2


if __name__ == "__main__":
    sys.exit(main())
