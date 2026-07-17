"""ROUND63 campaign cell runner — the run_cell(cell) contract consumed by
shard_runner.py (spec §3, §5).

A cell is a dict of grid coordinates:
  side (16/32/64/128), pattern ('bern50'|'hadpair'|'gam4'), rho_bar, nu,
  M (signed measurements), seed, arms (list of solver arm names),
  images ('pilot8' | 'all' | explicit list),
  det: optional dict of Detector overrides for the TRUE simulator
       (dark_frac, p_ap, ap_tau, paralyzable, tau_jitter_cv, start_mode, guard),
  est: optional dict of ESTIMATOR-side mismatches
       (tau_err e.g. +0.10, dark_known bool, assume_paralyzable bool),
  use_lpips (default False; PSNR/SSIM/rad-NRMSE always),
  fista_iters (default 200), tau (default 50e-9), sigma_b (default 0.0),
  dev (default False; True -> S1 development image set, disjoint from S2),
  select_rule ('discrepancy' default | 'legacy'), select_iter (default 60),
  sel_K (folds, default 5), gof_mode (default 'refit') — see run_cell.

lam_TV selection (spec D2 §4): for iterative arms (solvers._ITER_ARMS) under the
default 'discrepancy' rule, run_cell calls select_eta.select_eta_and_fit (the
six-step cross-fitted common renewal discrepancy rule) instead of run_arm's
deprecated own-NLL fallback. 'legacy' keeps every arm on the run_arm path.
Non-iterative arms (GI/DGI/EXACT) always use run_arm.

Physical accounting per cell: T = nu*tau; optical integration time =
n_physical_rows * T (pattern switching overhead t_switch is an analysis-time
parameter, not baked into rows). Flux is set by Phi = rho_bar/(tau*E[u]) with
E[u] = 1 analytic (unit-mean patterns x sum-normalized truth).

Returns list[dict] rows: one per (image, arm) with PSNR/SSIM/LPIPS (main
flux-matched protocol), radiometric NRMSE without rescaling (physical-scale
arms only; NaN for direction-only arms), lam_tv, runtime. Extra columns:
PSNR_rad (radiometric, non-rescaled PSNR — spec D2 §4 PRIMARY metric; ''
for direction-only arms), select_runtime_s (seconds inside select_eta_and_fit,
'' when the discrepancy rule did not run), MODEL_FAIL (bool for a discrepancy-
selected iterative arm, '' otherwise), eta_star (chosen dimensionless TV level,
'' otherwise). Existing columns are unchanged so shard merging stays stable.
RNG streams: rng_for(seed, 63, 3, ...) for buckets — disjoint from pattern
(63,1,*) and image (63,2,*) streams.
"""
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
from gi_core import metrics as MET
from gi_core.utils import rng_for
from images63 import build_image_set, build_dev_image_set
from patterns import make_patterns
from physics import Detector, simulate_counts
from solvers import ArmContext, run_arm, _ITER_ARMS
import select_eta

PILOT8 = ["stl_00", "stl_01", "stl_02", "stl_03",
          "text", "usaf_bars", "fine_lines", "low_contrast"]
PHYSICAL_ARMS = {"RQL", "POISSON-LIN", "SAT-POISSON", "PRECORRECT", "QMLE",
                 "EXACT"}
# RQL (the sigma_b=0 main production arm) reconstructs physical arrival rates
# like the other count-likelihood arms, so its radiometric NRMSE is meaningful
# and must be reported (it was previously omitted -> blank rad_nrmse for the
# flagship arm; fixed round-3).

_IMG_CACHE = {}
_PAT_CACHE = {}


def _images(side, spec, dev=False):
    """Build (and cache) the image set for a cell. dev=True routes to the S1
    DEVELOPMENT set (images63.build_dev_image_set: STL-10 TRAIN split + dev_*
    structural targets), which is DISJOINT from the S2 confirmatory set — the S1
    pilot passes dev=True + an explicit list of dev_* names so it runs the same
    run_cell code path as S2 without ever touching a confirmatory image."""
    key = ("dev" if dev else "conf", side)
    if key not in _IMG_CACHE:
        _IMG_CACHE[key] = build_dev_image_set(side) if dev else build_image_set(side)
    imgs = _IMG_CACHE[key]
    if spec == "pilot8":
        return {k: imgs[k] for k in PILOT8}
    if spec == "all":
        return dict(imgs)
    return {k: imgs[k] for k in spec}


def _patterns(kind, M, n, seed):
    key = (kind, M, n, seed)
    if key not in _PAT_CACHE:
        _PAT_CACHE.clear()  # hold one pattern set (memory)
        _PAT_CACHE[key] = make_patterns(kind, M, n, seed)
    return _PAT_CACHE[key]


def run_cell(cell):
    side = int(cell["side"])
    n = side * side
    tau = float(cell.get("tau", 50e-9))
    nu = float(cell["nu"])
    T = nu * tau
    rho = float(cell["rho_bar"])
    seed = int(cell["seed"])
    M = int(cell["M"])
    sigma_b = float(cell.get("sigma_b", 0.0))
    det_kw = dict(cell.get("det", {}))
    dark_frac = float(det_kw.pop("dark_frac", 0.0))
    Phi = rho / tau  # E[u] = 1 analytic
    det_true = Detector(tau=tau, dark=dark_frac * Phi, **det_kw)
    est_kw = dict(cell.get("est", {}))
    tau_est = tau * (1.0 + float(est_kw.get("tau_err", 0.0)))
    dark_est = det_true.dark if est_kw.get("dark_known", True) else 0.0
    det_est = Detector(tau=tau_est, dark=dark_est,
                       paralyzable=bool(est_kw.get("assume_paralyzable", False)),
                       start_mode=det_true.start_mode)

    pat = _patterns(cell["pattern"], M, n, seed)
    A, meta = pat["A"], pat["meta"]
    imgs = _images(side, cell.get("images", "pilot8"),
                   dev=bool(cell.get("dev", False)))
    # production lam_TV selection: "discrepancy" (default; six-step cross-fitted
    # common renewal discrepancy rule in select_eta) or "legacy" (the deprecated
    # own-NLL rule inside run_arm, retained for ablation / back-compat). Legacy
    # keeps ALL arms on the run_arm path.
    select_rule = cell.get("select_rule", "discrepancy")
    sel_K = int(cell.get("sel_K", 5))
    gof_mode = cell.get("gof_mode", "refit")
    rows = []
    for img_name, x in imgs.items():
        u = A @ x
        kind_id = {"bern50": 1, "hadpair": 2, "gam4": 3}[cell["pattern"]]
        rng_noise = rng_for(seed, 63, 3, kind_id, int(rho * 1000), int(nu), M,
                            list(imgs).index(img_name))
        b, N = simulate_counts(u, Phi, T, det_true, rng_noise, sigma_b=sigma_b)
        ctx = ArmContext(Phi=Phi, det=det_est, T=T, side=side,
                         sigma_b=sigma_b, n_iter=int(cell.get("fista_iters", 200)),
                         select_iter=int(cell.get("select_iter", 60)),
                         pattern_kind=cell["pattern"], meta=meta)
        for arm in cell["arms"]:
            t0 = time.time()
            use_select = (select_rule == "discrepancy" and arm in _ITER_ARMS)
            if use_select:
                # production rule: six-step discrepancy selection + refit
                # (spec D2 §4). One forward-compat shim: if this build of
                # select_eta_and_fit predates the gof_mode kwarg, retry without.
                t_sel = time.time()
                try:
                    xh, info = select_eta.select_eta_and_fit(
                        arm, A, b, ctx, K=sel_K, seed=seed, gof_mode=gof_mode)
                except TypeError as e:
                    if "gof_mode" not in str(e):
                        raise
                    xh, info = select_eta.select_eta_and_fit(
                        arm, A, b, ctx, K=sel_K, seed=seed)
                sel_dt = time.time() - t_sel
                model_fail = bool(info.get("MODEL_FAIL", False))
                es = info.get("eta_star", None)
                eta_out = round(float(es), 6) if es is not None else ""
                sel_out = round(sel_dt, 3)
            else:
                xh, info = run_arm(arm, A, b, ctx)
                model_fail, eta_out, sel_out = "", "", ""
            lam_tv = float(info.get("lam_tv", np.nan)) if isinstance(info, dict) \
                else np.nan
            m = MET.main_metrics(xh, x, side,
                                 with_lpips=bool(cell.get("use_lpips", False)))
            if arm in PHYSICAL_ARMS:
                xp = np.maximum(xh, 0.0)
                rad = float(np.linalg.norm(xp - x) / np.linalg.norm(x))
                flux_dev = float(xp.sum() / x.sum() - 1.0)
                # PSNR_rad = radiometric PSNR, NO rescaling (spec D2 §4 PRIMARY
                # metric): data_range = truth max, MSE on the nonneg-clipped
                # physical-scale reconstruction vs the sum-normalized truth.
                mse_rad = float(np.mean((xp - x) ** 2))
                psnr_rad = (10.0 * np.log10(float(x.max()) ** 2 / mse_rad)
                            if mse_rad > 0 else np.inf)
                psnr_rad_out = round(psnr_rad, 4) if np.isfinite(psnr_rad) else "inf"
            else:
                rad, flux_dev = np.nan, np.nan
                psnr_rad_out = ""
            rows.append({
                "side": side, "pattern": cell["pattern"], "rho_bar": rho,
                "nu": nu, "M": M, "seed": seed, "image": img_name, "arm": arm,
                "PSNR": round(m["PSNR"], 4), "SSIM": round(m["SSIM"], 5),
                "LPIPS": (round(m["LPIPS"], 5) if np.isfinite(m["LPIPS"]) else ""),
                "rad_nrmse": (round(rad, 5) if np.isfinite(rad) else ""),
                "flux_dev": (round(flux_dev, 5) if np.isfinite(flux_dev) else ""),
                "lam_tv": (round(float(lam_tv), 6) if np.isfinite(lam_tv) else ""),
                "mean_counts": round(float(N.mean()), 2),
                "optical_time_s": meta.get("n_physical_rows", M) * T,
                "dark_frac": dark_frac, "tau_err": est_kw.get("tau_err", 0.0),
                "runtime_s": round(time.time() - t0, 2),
                "select_runtime_s": sel_out,
                "MODEL_FAIL": model_fail,
                "eta_star": eta_out,
                "PSNR_rad": psnr_rad_out,
            })
    return rows
