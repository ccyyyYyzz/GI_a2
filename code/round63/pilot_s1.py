"""ROUND63 S1 DEV PILOT (local, development-only) — GPT round-6 ruling.

S1 FREEZE CLAUSE (spec D2 §1, verbatim):
  "S1 is exploratory and development-only. S1 images and seeds are disjoint
   from all S2 confirmatory images and seeds. S1 results shall not enter
   confirmatory confidence intervals or main quantitative tables."

Accordingly this script runs ONLY on the S1 development image set
(images63.build_dev_image_set: STL-10 TRAIN split + dev_* structural targets),
disjoint from the S2 confirmatory set (STL-10 TEST split). It never touches a
confirmatory image or seed; outputs are for protocol development only.

ROUND-6 two-pass protocol (docs/ROUND63_GPT_ROUND6_RULING.md §1.6, §4, §5):

  --pass-a  C0 CALIBRATION (once).  RQL arm only; every dev (rho,nu,seed,image)
            cell is reconstructed TWICE with the two endpoint c weights
            c_force = {0.25, 0.50} (via campaign.run_cell's c_force override).
            The endpoint-oracle regret objective (§1.6) chooses the frozen
            concentration threshold C0*. Writes c_calibration.csv,
            c_threshold_regret.json, calibration_manifest.json; writes
            code/round63/C0_FROZEN.json ONLY when --freeze-c0 is also given.

  --pass-b  FINAL PROTOCOL VERIFICATION.  Normal 5-arm sweep at the frozen C0*
            (select_eta.frozen_C0(); errors out if unfrozen), then the round-6
            §4 frozen floor-censoring analysis.

Geometry (round-6 §2/§3 S2-A, shrunk to the dev set): side=64, pattern="bern50",
M=4096 (M/n=1); dwell scan nu in {5,10,20,50,100,200,500,1000,2000}; two
operating points rho_bar in {0.05, 0.6}; tau=50 ns; sigma_b=0; dark=0; active
start. Arms: RQL, POISSON-LIN, SAT-POISSON, PRECORRECT (iterative, analytic
lam_TV) + GI (display reference; non-iterative). Dev image subset
dev_stl_00..03 + dev_text + dev_fine_lines; seeds {0,1}.

Every cell is run through campaign.run_cell — the SAME code path as S2 — so the
pilot also exercises the production analytic lam_TV rule (select_eta) and its
per-cell timing is directly transferable to the shard budget.

Run:
  # tiny pass-A slice (no freeze): prints the would-be C0 choice
  python code/round63/pilot_s1.py --pass-a --nu-list 100,500 --rho-list 0.05,0.6 \
         --images 2 --seeds 0 --side 32 --M 512
  # full pass A + freeze the threshold
  python code/round63/pilot_s1.py --pass-a --freeze-c0
  # final verification at the frozen C0
  python code/round63/pilot_s1.py --pass-b
  # dev/smoke of the normal path with an explicit C0 override
  python code/round63/pilot_s1.py --smoke --C0 4.0
  # rebuild summaries from an existing CSV
  python code/round63/pilot_s1.py --pass-b --analyze-only
"""
import argparse
import csv
import hashlib
import json
import os
import sys
import time
import traceback
from datetime import datetime, timezone
from types import SimpleNamespace

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
from campaign import run_cell
from gi_core.utils import rng_for
from images63 import DEV_IMG_ROOT
from select_eta import frozen_C0, C0_FILE, C0_CANDIDATES

ROOT = os.path.dirname(os.path.dirname(HERE))

# ---- frozen S1 dev configuration (round-6 §2/§3/§5) ----------------------- #
TAU = 50e-9
PATTERN = "bern50"
DEFAULT_IMAGES = ["dev_stl_00", "dev_stl_01", "dev_stl_02", "dev_stl_03",
                  "dev_text", "dev_fine_lines"]
DEFAULT_ARMS = ["RQL", "POISSON-LIN", "SAT-POISSON", "PRECORRECT", "GI"]
DEFAULT_NU = [5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0, 1000.0, 2000.0]
DEFAULT_RHO = [0.05, 0.6]
DEFAULT_SEEDS = [0, 1]
DEFAULT_M = 4096                                    # round-6 §2: M/n = 1 at 64^2

PRIMARY_ARM = "RQL"                                 # round-6 §4 primary gate S_j
SECONDARY_TTQ = ["POISSON-LIN", "SAT-POISSON", "PRECORRECT"]

# the two endpoint c weights the Pass-A calibration reconstructs at every cell
C_ENDPOINTS = (0.25, 0.50)

# stable CSV schema = EXACTLY campaign.run_cell's current emission order
# (round-6: eta_star dropped; the analytic-rule ledger + audit columns added).
FIELDNAMES = ["side", "pattern", "rho_bar", "nu", "M", "seed", "image", "arm",
              "PSNR", "SSIM", "LPIPS", "rad_nrmse", "flux_dev", "lam_tv",
              "mean_counts", "optical_time_s", "dark_frac", "tau_err",
              "runtime_s", "select_runtime_s", "PSNR_rad",
              "C_hat", "S_N2", "V0_exact", "muprime0_exact", "omega_A",
              "c_used", "sigma_grad_arm", "lambda_tv_arm",
              "audit_status", "d_ratio", "q_d", "q_mean", "leak_suspect"]
# Pass A adds one column marking which endpoint c produced the row.
PASSA_FIELDNAMES = FIELDNAMES + ["c_force"]

_KEY_COLS = ("side", "pattern", "rho_bar", "nu", "M", "seed", "image", "arm")

# gate bookkeeping (round-6 §4): statuses excluded from the "S_gate>1" tally.
# C = both-left-censored (S_gate fixed 1.0, not a measured 1x); E = fast never
# reaches (failure, S_gate 0); F = analysis failure (S_gate 0). D naturally has
# S_gate<1 so it is not >1 anyway (it counts as a non-positive result).
_GT1_EXCLUDED = {"BOTH_LEFT_CENSORED", "FAST_RIGHT_CENSORED", "ANALYSIS_FAILURE"}
_TIE_DB = 0.02                                      # round-6 §1.6 C0 tie window


# ---- small coercion / key helpers ----------------------------------------- #
def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


def _psnr_rad(r):
    """PSNR_rad of a row as a float, or None when absent/unusable."""
    v = r.get("PSNR_rad", "")
    if v in ("", None):
        return None
    if v == "inf":
        return float("inf")
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _rkey(d):
    """Canonical resume key. str() of a native run_cell value equals the string
    csv.DictWriter wrote and DictReader read back, so out-rows and CSV-rows and
    _cell_keys all agree."""
    return tuple(str(d[k]) for k in _KEY_COLS)


def _cell_keys(cell):
    """Expected per-(image,arm) resume keys for a cell, mirroring run_cell's own
    int()/float() coercions so the strings match _rkey exactly."""
    side = str(int(cell["side"]))
    pat = str(cell["pattern"])
    rho = str(float(cell["rho_bar"]))
    nu = str(float(cell["nu"]))
    M = str(int(cell["M"]))
    seed = str(int(cell["seed"]))
    return [(side, pat, rho, nu, M, seed, img, arm)
            for img in cell["images"] for arm in cell["arms"]]


def _cf_str(c_force):
    return str(float(c_force))


def _rkey_cf(d):
    return _rkey(d) + (_cf_str(_f(d.get("c_force", ""))),)


def _cell_keys_cf(cell, c_force):
    tag = _cf_str(c_force)
    return [k + (tag,) for k in _cell_keys(cell)]


def _sha256(path):
    if not os.path.exists(path):
        return None
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def _c0str(c):
    return "inf" if (isinstance(c, float) and np.isinf(c)) else ("%g" % c)


# ---- configuration --------------------------------------------------------- #
def _floats(s):
    return [float(x) for x in s.replace(" ", "").split(",") if x != ""]


def _ints(s):
    return [int(x) for x in s.replace(" ", "").split(",") if x != ""]


def _strs(s):
    return [x.strip() for x in s.split(",") if x.strip()]


def parse_args(argv):
    ap = argparse.ArgumentParser(
        description="ROUND63 S1 dev pilot (development-only; round-6 ruling).")
    ap.add_argument("--pass-a", action="store_true", dest="pass_a",
                    help="C0 calibration pass (RQL only, c={0.25,0.50} endpoints;"
                         " round-6 §1.6/§5)")
    ap.add_argument("--pass-b", action="store_true", dest="pass_b",
                    help="final verification pass (5 arms at the frozen C0)")
    ap.add_argument("--freeze-c0", action="store_true", dest="freeze_c0",
                    help="pass A only: write code/round63/C0_FROZEN.json with the"
                         " chosen C0* (without it, only print the would-be choice)")
    ap.add_argument("--C0", type=float, default=None,
                    help="explicit concentration threshold for dev/smoke runs "
                         "(overrides the frozen file; ignored in --pass-a/--pass-b)")
    ap.add_argument("--nu-list", type=str, default=None,
                    help="comma dwell list (default 5,10,20,50,100,200,500,1000,2000)")
    ap.add_argument("--rho-list", type=str, default=None,
                    help="comma rho_bar list (default 0.05,0.6)")
    ap.add_argument("--arms", type=str, default=None,
                    help="comma arm list (default RQL,POISSON-LIN,SAT-POISSON,"
                         "PRECORRECT,GI)")
    ap.add_argument("--images", type=int, default=None,
                    help="use the first N of the default dev image list")
    ap.add_argument("--seeds", type=str, default=None,
                    help="comma seed list (default 0,1)")
    ap.add_argument("--side", type=int, default=None, help="default 64")
    ap.add_argument("--M", type=int, default=None, help="signed measurements "
                    "(default 4096 = M/n 1 at side 64)")
    ap.add_argument("--select-iter", type=int, default=60,
                    help="fista iters inside lam_TV selection (default 60)")
    ap.add_argument("--fista-iters", type=int, default=200,
                    help="final-refit fista iters (default 200)")
    ap.add_argument("--smoke", action="store_true",
                    help="fast preset: side=32,M=512,2 images,seeds{0},"
                         "nu{100,500},arms RQL+POISSON-LIN+GI (own out dir)")
    ap.add_argument("--out", type=str, default=None,
                    help="output dir (default results/round63_s1[_passA|_passB])")
    ap.add_argument("--analyze-only", action="store_true",
                    help="skip the sweep, (re)build summaries from existing CSV")
    ap.add_argument("--no-analysis", action="store_true",
                    help="sweep only (parallel worker mode)")
    return ap.parse_args(argv)


def build_cfg(a):
    smoke = bool(a.smoke)
    mode = "pass_a" if a.pass_a else ("pass_b" if a.pass_b else "normal")
    side = a.side if a.side is not None else (32 if smoke else 64)
    M = a.M if a.M is not None else (512 if smoke else DEFAULT_M)
    n_img = a.images if a.images is not None else (2 if smoke else len(DEFAULT_IMAGES))
    images = DEFAULT_IMAGES[:max(1, n_img)]
    seeds = _ints(a.seeds) if a.seeds else ([0] if smoke else list(DEFAULT_SEEDS))
    nu_list = _floats(a.nu_list) if a.nu_list else ([100.0, 500.0] if smoke
                                                    else list(DEFAULT_NU))
    rho_list = _floats(a.rho_list) if a.rho_list else list(DEFAULT_RHO)
    # Pass A is RQL-only regardless of --arms (the endpoint regret needs only RQL).
    if mode == "pass_a":
        arms = [PRIMARY_ARM]
    else:
        arms = _strs(a.arms) if a.arms else (["RQL", "POISSON-LIN", "GI"] if smoke
                                             else list(DEFAULT_ARMS))
    if a.out:
        out = a.out
    elif mode == "pass_a":
        out = os.path.join(ROOT, "results", "round63_s1_passA")
    elif mode == "pass_b":
        out = os.path.join(ROOT, "results", "round63_s1_passB")
    else:
        out = os.path.join(ROOT, "results",
                           "round63_s1_smoke" if smoke else "round63_s1")
    return SimpleNamespace(
        side=side, M=M, pattern=PATTERN, tau=TAU, images=images, seeds=seeds,
        nu_list=nu_list, rho_list=rho_list, arms=arms,
        select_iter=a.select_iter, fista_iters=a.fista_iters,
        out=out, smoke=smoke, mode=mode, freeze_c0=bool(a.freeze_c0),
        C0=a.C0)


def cfg_dict(cfg):
    return {"side": cfg.side, "M": cfg.M, "pattern": cfg.pattern, "tau": cfg.tau,
            "nu_list": cfg.nu_list, "rho_list": cfg.rho_list, "arms": cfg.arms,
            "images": cfg.images, "seeds": cfg.seeds,
            "select_iter": cfg.select_iter, "fista_iters": cfg.fista_iters,
            "smoke": cfg.smoke, "mode": cfg.mode,
            "C0": (None if cfg.C0 is None else float(cfg.C0))}


def _base_cell(cfg, rho, nu, seed, arms):
    return dict(
        side=cfg.side, pattern=cfg.pattern, rho_bar=float(rho), nu=float(nu),
        M=cfg.M, seed=int(seed), arms=list(arms), images=list(cfg.images),
        dev=True, tau=cfg.tau, sigma_b=0.0, fista_iters=cfg.fista_iters,
        select_iter=cfg.select_iter, select_rule="discrepancy", use_lpips=False)


def cells(cfg):
    """Normal / Pass-B cell list: one cell per (rho_bar, nu, seed). The frozen
    C0 is threaded to run_cell -> select_eta.select_eta_and_fit."""
    cs = []
    for rho in cfg.rho_list:
        for nu in cfg.nu_list:
            for seed in cfg.seeds:
                cell = _base_cell(cfg, rho, nu, seed, cfg.arms)
                cell["C0"] = cfg.C0
                cs.append(cell)
    return cs


def cells_passa(cfg, c_force):
    """Pass-A cell list at one endpoint c: RQL only, c forced (no C0 needed —
    select_eta's c_force branch bypasses the concentration threshold)."""
    cs = []
    for rho in cfg.rho_list:
        for nu in cfg.nu_list:
            for seed in cfg.seeds:
                cell = _base_cell(cfg, rho, nu, seed, [PRIMARY_ARM])
                cell["c_force"] = float(c_force)
                # Pass A needs only the two endpoint PSNR_rad values; the
                # descriptive audit (~40s/cell-image at M=4096) is Pass-B work
                cell["audit"] = False
                cs.append(cell)
    return cs


# ---- normal / Pass-B sweep (incremental, resumable) ----------------------- #
def sweep(cfg):
    os.makedirs(cfg.out, exist_ok=True)
    path = os.path.join(cfg.out, "pilot_rows.csv")
    errlog = os.path.join(cfg.out, "pilot_errors.log")

    done = set()
    header_present = os.path.exists(path) and os.path.getsize(path) > 0
    if header_present:
        with open(path, newline="") as f:
            rd = csv.DictReader(f)
            existing = rd.fieldnames
            for r in rd:
                done.add(_rkey(r))
        if existing is not None and list(existing) != FIELDNAMES:
            raise SystemExit(
                "[s1] existing %s has an incompatible header\n  found:  %s\n"
                "  expect: %s\nremove/rename it to re-run." % (path, existing,
                                                              FIELDNAMES))

    cell_list = cells(cfg)
    n_cells = len(cell_list)
    warned_schema = False
    t0 = time.time()
    with open(path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES, extrasaction="ignore")
        if not header_present:
            writer.writeheader()
            f.flush()
        for ci, cell in enumerate(cell_list):
            keys = _cell_keys(cell)
            if all(k in done for k in keys):
                print("[s1] SKIP %d/%d rho=%.3g nu=%g seed=%d (all %d rows present)"
                      % (ci + 1, n_cells, cell["rho_bar"], cell["nu"],
                         cell["seed"], len(keys)), flush=True)
                continue
            tc = time.time()
            try:
                rows = run_cell(cell)
            except Exception:
                tb = traceback.format_exc()
                with open(errlog, "a") as ef:
                    ef.write("[%s] cell rho=%.4g nu=%g seed=%d side=%d M=%d\n%s\n"
                             % (time.strftime("%Y-%m-%d %H:%M:%S"),
                                cell["rho_bar"], cell["nu"], cell["seed"],
                                cell["side"], cell["M"], tb))
                print("[s1] ERROR %d/%d rho=%.3g nu=%g seed=%d -> logged to %s, "
                      "continuing" % (ci + 1, n_cells, cell["rho_bar"],
                                      cell["nu"], cell["seed"], errlog), flush=True)
                continue
            if rows and not warned_schema and set(rows[0].keys()) != set(FIELDNAMES):
                warned_schema = True
                print("[s1] WARNING run_cell row schema differs from FIELDNAMES; "
                      "extra keys dropped: %s"
                      % (set(rows[0].keys()) - set(FIELDNAMES)), flush=True)
            n_new = 0
            for r in rows:
                k = _rkey(r)
                if k in done:
                    continue
                writer.writerow(r)
                done.add(k)
                n_new += 1
            f.flush()
            os.fsync(f.fileno())
            print("[s1] cell %d/%d rho=%.3g nu=%g seed=%d arms=%d imgs=%d -> "
                  "+%d rows (%.1fs, total %.0fs)"
                  % (ci + 1, n_cells, cell["rho_bar"], cell["nu"], cell["seed"],
                     len(cell["arms"]), len(cell["images"]), n_new,
                     time.time() - tc, time.time() - t0), flush=True)
    print("[s1] SWEEP DONE wall=%.0fs rows=%s" % (time.time() - t0, path),
          flush=True)
    return path


# ---- Pass-A sweep (two endpoint c per cell, resumable) -------------------- #
def sweep_passa(cfg):
    os.makedirs(cfg.out, exist_ok=True)
    path = os.path.join(cfg.out, "pass_a_rows.csv")
    errlog = os.path.join(cfg.out, "pass_a_errors.log")

    done = set()
    header_present = os.path.exists(path) and os.path.getsize(path) > 0
    if header_present:
        with open(path, newline="") as f:
            rd = csv.DictReader(f)
            existing = rd.fieldnames
            for r in rd:
                done.add(_rkey_cf(r))
        if existing is not None and list(existing) != PASSA_FIELDNAMES:
            raise SystemExit(
                "[s1:A] existing %s has an incompatible header\n  found:  %s\n"
                "  expect: %s\nremove/rename it to re-run." % (path, existing,
                                                              PASSA_FIELDNAMES))

    plan = [(cf, cell) for cf in C_ENDPOINTS for cell in cells_passa(cfg, cf)]
    n_cells = len(plan)
    t0 = time.time()
    with open(path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=PASSA_FIELDNAMES,
                                extrasaction="ignore")
        if not header_present:
            writer.writeheader()
            f.flush()
        for ci, (cf, cell) in enumerate(plan):
            keys = _cell_keys_cf(cell, cf)
            if all(k in done for k in keys):
                print("[s1:A] SKIP %d/%d c=%.2f rho=%.3g nu=%g seed=%d (present)"
                      % (ci + 1, n_cells, cf, cell["rho_bar"], cell["nu"],
                         cell["seed"]), flush=True)
                continue
            tc = time.time()
            try:
                rows = run_cell(cell)
            except Exception:
                tb = traceback.format_exc()
                with open(errlog, "a") as ef:
                    ef.write("[%s] c=%.2f rho=%.4g nu=%g seed=%d side=%d M=%d\n%s\n"
                             % (time.strftime("%Y-%m-%d %H:%M:%S"), cf,
                                cell["rho_bar"], cell["nu"], cell["seed"],
                                cell["side"], cell["M"], tb))
                print("[s1:A] ERROR %d/%d c=%.2f rho=%.3g nu=%g seed=%d -> %s"
                      % (ci + 1, n_cells, cf, cell["rho_bar"], cell["nu"],
                         cell["seed"], errlog), flush=True)
                continue
            n_new = 0
            for r in rows:
                r["c_force"] = float(cf)
                k = _rkey_cf(r)
                if k in done:
                    continue
                writer.writerow(r)
                done.add(k)
                n_new += 1
            f.flush()
            os.fsync(f.fileno())
            print("[s1:A] cell %d/%d c=%.2f rho=%.3g nu=%g seed=%d -> +%d rows "
                  "(%.1fs, total %.0fs)"
                  % (ci + 1, n_cells, cf, cell["rho_bar"], cell["nu"],
                     cell["seed"], n_new, time.time() - tc, time.time() - t0),
                  flush=True)
    print("[s1:A] SWEEP DONE wall=%.0fs rows=%s" % (time.time() - t0, path),
          flush=True)
    return path


# ---- Pass-A analysis: endpoint-oracle C0 calibration (round-6 §1.6) -------- #
def analyze_passa(cfg, freeze):
    path = os.path.join(cfg.out, "pass_a_rows.csv")
    if not os.path.exists(path):
        raise SystemExit("[s1:A] no rows CSV at %s (run --pass-a sweep first)"
                         % path)
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise SystemExit("[s1:A] rows CSV %s is empty" % path)

    # index: (image, rho, nu, seed) -> {c_force: {"P":.., "C_hat":..}}
    cells_data = {}
    for r in rows:
        if r.get("arm") != PRIMARY_ARM:
            continue
        cf = round(_f(r.get("c_force", "")), 3)
        key = (r["image"], _f(r["rho_bar"]), _f(r["nu"]), int(_f(r["seed"])))
        cells_data.setdefault(key, {})[cf] = {
            "P": _psnr_rad(r), "C_hat": _f(r.get("C_hat", ""))}

    images = sorted({k[0] for k in cells_data})
    lo, hi = round(C_ENDPOINTS[0], 3), round(C_ENDPOINTS[1], 3)

    # R_j(C0) per image + the wide per-(cell,C0) audit rows
    R = {im: {} for im in images}
    calib_rows = []
    for C0 in C0_CANDIDATES:
        c0s = _c0str(C0)
        for im in images:
            regs = []
            for key, d in cells_data.items():
                if key[0] != im or lo not in d or hi not in d:
                    continue
                P25, P50 = d[lo]["P"], d[hi]["P"]
                chat = d[lo].get("C_hat")
                if chat is None or not np.isfinite(chat):
                    chat = d[hi].get("C_hat")
                if P25 is None or P50 is None or chat is None or not np.isfinite(chat):
                    continue
                c_pick = 0.50 if chat <= C0 else 0.25
                P_pick = P50 if c_pick == 0.50 else P25
                P_max = max(P25, P50)
                reg = P_max - P_pick
                regs.append(reg)
                calib_rows.append({
                    "image": im, "rho_bar": key[1], "nu": key[2], "seed": key[3],
                    "C0": c0s, "C_hat": round(float(chat), 6),
                    "c_pick": c_pick, "P_25": round(float(P25), 5),
                    "P_50": round(float(P50), 5), "P_max": round(float(P_max), 5),
                    "P_pick": round(float(P_pick), 5), "regret": round(float(reg), 5)})
            R[im][C0] = float(np.mean(regs)) if regs else float("nan")

    # J(C0) = 90th percentile over images (numpy linear interpolation)
    J = {}
    for C0 in C0_CANDIDATES:
        vals = [R[im][C0] for im in images if np.isfinite(R[im][C0])]
        J[C0] = float(np.quantile(vals, 0.90)) if vals else float("nan")

    finite_J = {C0: j for C0, j in J.items() if np.isfinite(j)}
    if not finite_J:
        raise SystemExit("[s1:A] no finite J(C0) — pass-A data insufficient")
    J_min = min(finite_J.values())
    # tie rule: among candidates within J_min + 0.02 dB, take the LARGEST C0
    # (candidates are ascending; the last qualifying one). inf -> global c=0.50.
    within = [C0 for C0 in C0_CANDIDATES
              if C0 in finite_J and finite_J[C0] <= J_min + _TIE_DB]
    C0_star = max(within)

    # ---- write artifacts
    os.makedirs(cfg.out, exist_ok=True)
    ccsv = os.path.join(cfg.out, "c_calibration.csv")
    cc_fields = ["image", "rho_bar", "nu", "seed", "C0", "C_hat", "c_pick",
                 "P_25", "P_50", "P_max", "P_pick", "regret"]
    with open(ccsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cc_fields)
        w.writeheader()
        for row in calib_rows:
            w.writerow(row)

    regret_json = {
        "candidates": [_c0str(c) for c in C0_CANDIDATES],
        "objective": "R_j = mean_{rho,nu,seed}[max(P.25,P.50) - P_pick]; "
                     "J(C0) = Q0.90 over dev images; argmin J, tie<=0.02dB -> "
                     "larger C0 (round-6 §1.6)",
        "R_j": {im: {_c0str(C0): R[im][C0] for C0 in C0_CANDIDATES}
                for im in images},
        "J": {_c0str(C0): J[C0] for C0 in C0_CANDIDATES},
        "J_min": J_min, "tie_window_dB": _TIE_DB,
        "within_tie": [_c0str(c) for c in within],
        "C0_star": _c0str(C0_star), "C0_star_value": (
            None if np.isinf(C0_star) else float(C0_star)),
        "n_dev_images": len(images), "images": images}
    with open(os.path.join(cfg.out, "c_threshold_regret.json"), "w") as f:
        json.dump(regret_json, f, indent=2, default=_json_default)

    dev_sha_file = os.path.join(DEV_IMG_ROOT, str(cfg.side), "sha256.txt")
    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "mode": "pass_a", "side": cfg.side, "M": cfg.M, "pattern": cfg.pattern,
        "tau": cfg.tau, "nu_list": cfg.nu_list, "rho_list": cfg.rho_list,
        "images": cfg.images, "seeds": cfg.seeds, "arm": PRIMARY_ARM,
        "c_endpoints": list(C_ENDPOINTS), "C0_star": _c0str(C0_star),
        "sha256": {
            "pilot_s1.py": _sha256(os.path.join(HERE, "pilot_s1.py")),
            "select_eta.py": _sha256(os.path.join(HERE, "select_eta.py")),
            "physics.py": _sha256(os.path.join(HERE, "physics.py")),
            "campaign.py": _sha256(os.path.join(HERE, "campaign.py")),
            "dev_images_sha_file": _sha256(dev_sha_file)},
        "dev_images_sha_file_path": dev_sha_file}
    with open(os.path.join(cfg.out, "calibration_manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2, default=_json_default)

    # ---- stdout digest
    print("[s1:A] CALIBRATION  images=%d  cells=%d" % (len(images), len(cells_data)),
          flush=True)
    for C0 in C0_CANDIDATES:
        mark = "  <== C0*" if C0 == C0_star else ""
        print("       C0=%-4s  J(C0)=%s%s"
              % (_c0str(C0), ("%.4f" % J[C0]) if np.isfinite(J[C0]) else "n/a",
                 mark), flush=True)
    print("[s1:A] C0* = %s  (J_min=%.4f, tie window %.2f dB, within=%s)"
          % (_c0str(C0_star), J_min, _TIE_DB,
             ",".join(_c0str(c) for c in within)), flush=True)

    if freeze:
        payload = {"C0": (None if np.isinf(C0_star) else float(C0_star)),
                   "frozen_utc": datetime.now(timezone.utc).isoformat(),
                   "source": "passA"}
        if np.isinf(C0_star):
            # inf -> global c=0.50; select_eta.frozen_C0() reads a float, and
            # Python json round-trips Infinity, so store the literal.
            payload["C0"] = float("inf")
        with open(C0_FILE, "w") as f:
            json.dump(payload, f, indent=2)
        print("[s1:A] FROZE %s -> C0=%s" % (C0_FILE, _c0str(C0_star)), flush=True)
    else:
        print("[s1:A] (dry run) --freeze-c0 NOT given; C0_FROZEN.json NOT written."
              "  Would freeze C0=%s at %s" % (_c0str(C0_star), C0_FILE),
              flush=True)
    return {"C0_star": _c0str(C0_star), "J": {_c0str(c): J[c] for c in C0_CANDIDATES},
            "frozen": bool(freeze)}


# ---- round-6 §4 frozen floor-censoring analysis --------------------------- #
def _t_reach_interior(T, P, Q):
    """Interior crossing time by linear interpolation in log10(T) on the
    monotonized curve P (assumes P[0] < Q <= P.max(), so k>=1)."""
    k = int(np.nonzero(P >= Q)[0][0])
    logT = np.log10(T)
    if P[k] == P[k - 1]:
        lt = logT[k]
    else:
        lt = logT[k - 1] + (Q - P[k - 1]) * (logT[k] - logT[k - 1]) \
            / (P[k] - P[k - 1])
    return float(10.0 ** lt)


def _mono_curve(seed_map, seeds, T_by_nu):
    """(T asc, cummax P) seed-mean PSNR_rad curve over the given seed multiset,
    or None. Cumulative max = isotonic-lite monotonization (round-5 convention)."""
    Ts, Ps = [], []
    for nu in sorted(seed_map):
        vals = [seed_map[nu][s] for s in seeds
                if s in seed_map[nu] and seed_map[nu][s] is not None]
        if not vals:
            continue
        Ts.append(T_by_nu[nu])
        Ps.append(float(np.mean(vals)))
    if not Ts:
        return None
    T = np.asarray(Ts, dtype=np.float64)
    P = np.asarray(Ps, dtype=np.float64)
    order = np.argsort(T)
    T, P = T[order], P[order]
    return T, np.maximum.accumulate(P)


def _reach_state(curve, Q):
    """('RIGHT'|'FLOOR'|'INTERIOR', T_c). RIGHT = never reaches Q; FLOOR =
    already at/above Q at the shortest T (left-censored); INTERIOR = crosses
    strictly inside the grid."""
    if curve is None:
        return "RIGHT", float("nan")
    T, P = curve
    if float(P.max()) < Q:
        return "RIGHT", float("nan")
    if float(P[0]) >= Q:
        return "FLOOR", float(T[0])
    return "INTERIOR", _t_reach_interior(T, P, Q)


def _judge_image(entry, seeds, T_by_nu, T_min, nu_ref, q_override=None):
    """Classify one image's (safe,fast) curves per round-6 §4 and return
    (status, S_gate, T_s, T_f, Qstar). Q*_j is the RQL-derived target
    (seed-mean PSNR_rad of the safe curve at nu_ref, minus 1 dB); q_override
    supplies that same RQL target when judging the SECONDARY arms (round-6 §4:
    one shared target for every arm)."""
    safe, fast = entry["safe"], entry["fast"]
    if q_override is not None:
        Q = float(q_override)
    else:
        if nu_ref not in safe:
            return "ANALYSIS_FAILURE", 0.0, None, None, None
        qvals = [safe[nu_ref][s] for s in seeds
                 if s in safe[nu_ref] and safe[nu_ref][s] is not None]
        if not qvals:
            return "ANALYSIS_FAILURE", 0.0, None, None, None
        Q = float(np.mean(qvals)) - 1.0

    cs = _mono_curve(safe, seeds, T_by_nu)
    cf = _mono_curve(fast, seeds, T_by_nu)
    st_s, T_s = _reach_state(cs, Q)
    st_f, T_f = _reach_state(cf, Q)

    # F: safe fails / data missing (should be impossible by Q* construction).
    if st_s == "RIGHT" or cs is None:
        return "ANALYSIS_FAILURE", 0.0, T_s, T_f, Q
    # E: fast never reaches Q* by nu_max -> failure (S_gate 0).
    if st_f == "RIGHT" or cf is None:
        return "FAST_RIGHT_CENSORED", 0.0, T_s, T_f, Q
    # A/B/C/D over {FLOOR, INTERIOR} x {FLOOR, INTERIOR}
    if st_s == "INTERIOR" and st_f == "INTERIOR":
        s = (T_s / T_f) if T_f > 0 else 0.0
        return "NORMAL", float(s), T_s, T_f, Q
    if st_s == "INTERIOR" and st_f == "FLOOR":
        s = (T_s / T_min) if T_min > 0 else 0.0        # lower bound
        return "FAST_LEFT_CENSORED", float(s), T_s, T_f, Q
    if st_s == "FLOOR" and st_f == "INTERIOR":
        s = (T_min / T_f) if T_f > 0 else 0.0          # < 1: non-positive
        return "SAFE_LEFT_FAST_INTERIOR", float(s), T_s, T_f, Q
    # both FLOOR: unresolved below floor, S_gate fixed to 1.0
    return "BOTH_LEFT_CENSORED", 1.0, T_s, T_f, Q


def _gate_stats(items):
    """items: list of (status, S_gate). Returns (median over ALL, count S_gate>1
    excluding C/E/F per round-6 §4)."""
    S = [sg for (_, sg) in items]
    median = float(np.median(S)) if S else float("nan")
    n_gt1 = sum(1 for (st, sg) in items if st not in _GT1_EXCLUDED and sg > 1.0)
    return median, n_gt1


def _collect(rows, arm, rho_safe, rho_fast, images):
    """Per-image seed maps for one arm:
       data[image] = {"safe": {nu:{seed:P}}, "fast": {nu:{seed:P}}}.
    Also T_by_nu and per-image seed lists (censoring is judged, not frozen,
    so raw per-seed values are retained for bootstrap resampling)."""
    data = {im: {"safe": {}, "fast": {}} for im in images}
    T_by_nu = {}
    seeds_by_image = {im: set() for im in images}
    for r in rows:
        if r.get("arm") != arm:
            continue
        im = r.get("image")
        if im not in data:
            continue
        p = _psnr_rad(r)
        if p is None:
            continue
        rho = _f(r["rho_bar"])
        nu = _f(r["nu"])
        seed = int(_f(r["seed"]))
        T_by_nu[nu] = _f(r["optical_time_s"])
        if abs(rho - rho_safe) < 1e-12:
            which = "safe"
        elif abs(rho - rho_fast) < 1e-12:
            which = "fast"
        else:
            continue
        data[im][which].setdefault(nu, {})[seed] = p
        seeds_by_image[im].add(seed)
    seeds_by_image = {im: sorted(s) for im, s in seeds_by_image.items()}
    return data, T_by_nu, seeds_by_image


def _bootstrap_lower(data, seeds_by_image, T_by_nu, T_min, nu_ref, images, B=2000):
    """Image-level stratified bootstrap of the median S_gate (round-6 §4): every
    draw resamples seeds within image, re-forms both curves, RE-JUDGES the
    censoring class, regenerates S_gate, then resamples images. Returns the 2.5%
    lower bound of the median-S_gate distribution."""
    if not images:
        return float("nan"), np.array([])
    rng = rng_for(0, 63, 6, 1)
    n = len(images)
    meds = np.empty(B, dtype=np.float64)
    for d in range(B):
        per_img = []
        for im in images:
            sl = seeds_by_image[im]
            if sl:
                drawn = [int(s) for s in rng.choice(sl, size=len(sl), replace=True)]
            else:
                drawn = []
            st, sg, _, _, _ = _judge_image(data[im], drawn, T_by_nu, T_min, nu_ref)
            per_img.append(sg)
        idx = rng.integers(0, n, size=n)
        meds[d] = float(np.median([per_img[i] for i in idx]))
    return float(np.quantile(meds, 0.025)), meds


def censoring_analysis(rows, cfg):
    rho_safe = min(cfg.rho_list)
    rho_fast = max(cfg.rho_list)
    images = [im for im in cfg.images
              if any(r.get("image") == im and r.get("arm") == PRIMARY_ARM
                     for r in rows)]
    data, T_by_nu, seeds_by_image = _collect(rows, PRIMARY_ARM, rho_safe,
                                             rho_fast, images)
    if not T_by_nu or not images:
        return {"rho_safe": rho_safe, "rho_fast": rho_fast, "images": images,
                "per_image": [], "S_median": None, "n_Sgate_gt1": 0,
                "n_images": len(images), "note": "no PSNR_rad rows for the "
                "primary arm"}
    nu_ref = max(T_by_nu)
    T_min = float(min(T_by_nu.values()))

    per_image = []
    warnings = []
    for im in images:
        st, sg, T_s, T_f, Q = _judge_image(data[im], seeds_by_image[im],
                                           T_by_nu, T_min, nu_ref)
        per_image.append({"image": im, "status": st, "S_gate": sg,
                          "T_s": T_s, "T_f": T_f, "Qstar": Q,
                          "n_seeds": len(seeds_by_image[im])})
        if st == "ANALYSIS_FAILURE":
            w = ("[s1] LOUD WARNING cell-completeness: image %s -> "
                 "ANALYSIS_FAILURE (safe arm never reaches Q* or data missing; "
                 "round-6 §4-F). Check the cell grid." % im)
            warnings.append(w)
            print(w, flush=True)

    items = [(d["status"], d["S_gate"]) for d in per_image]
    median, n_gt1 = _gate_stats(items)
    boot_lower, _ = _bootstrap_lower(data, seeds_by_image, T_by_nu, T_min,
                                     nu_ref, images, B=2000)

    # secondary arms: point-estimate median only (descriptive). They are judged
    # against the SAME RQL-derived Q*_j per image (round-6 §4: one shared target).
    qstar_by_image = {d["image"]: d["Qstar"] for d in per_image}
    secondary = {}
    for arm in [a for a in SECONDARY_TTQ if a in cfg.arms]:
        d2, T2, s2 = _collect(rows, arm, rho_safe, rho_fast, images)
        if not T2:
            continue
        nr2, tm2 = max(T2), float(min(T2.values()))
        it2 = []
        for im in images:
            q = qstar_by_image.get(im)
            if q is None:
                continue
            st, sg, _, _, _ = _judge_image(d2[im], s2[im], T2, tm2, nr2,
                                           q_override=q)
            it2.append((st, sg))
        if it2:
            m2, g2 = _gate_stats(it2)
            secondary[arm] = {"S_median": m2, "n_Sgate_gt1": g2}

    n_status = {}
    for d in per_image:
        n_status[d["status"]] = n_status.get(d["status"], 0) + 1

    return {"rho_safe": rho_safe, "rho_fast": rho_fast, "nu_ref": nu_ref,
            "T_min": T_min, "T_by_nu": {("%g" % k): v for k, v in
                                        sorted(T_by_nu.items())},
            "images": images, "per_image": per_image,
            "S_median": median, "n_Sgate_gt1": n_gt1, "n_images": len(images),
            "status_counts": n_status,
            "bootstrap_B": 2000, "bootstrap_median_lower_2p5": boot_lower,
            "secondary_arms": secondary,
            "reference_gate": {"median>=3": (median is not None and median >= 3.0),
                               "bootstrap_lower>1": (np.isfinite(boot_lower)
                                                     and boot_lower > 1.0),
                               "n_gt1": n_gt1, "n_images": len(images)},
            "warnings": warnings}


# ---- descriptive audit summary (round-5: no binary gate) ------------------ #
def audit_table(rows):
    """Aggregate the RQL cell-level DESCRIPTIVE audit: D_ratio spread and
    LEAKAGE_SUSPECT counts by (rho, nu). Pure diagnostics — never a gate."""
    agg = {}
    for r in rows:
        if r.get("audit_status", "") != "OK":
            continue
        key = (_f(r["rho_bar"]), _f(r["nu"]))
        d = agg.setdefault(key, {"d_ratio": [], "leak": 0, "tot": 0})
        try:
            d["d_ratio"].append(float(r["d_ratio"]))
        except (TypeError, ValueError):
            pass
        d["tot"] += 1
        if str(r.get("leak_suspect", "")) == "True":
            d["leak"] += 1
    out = {}
    for k, v in sorted(agg.items()):
        dr = v["d_ratio"]
        out["%g|%g" % k] = {
            "n": v["tot"], "leak": v["leak"],
            "d_ratio_med": (sorted(dr)[len(dr) // 2] if dr else None),
            "d_ratio_max": (max(dr) if dr else None)}
    return out


# ---- runtime calibration -------------------------------------------------- #
def runtime_calibration(rows):
    agg = {}
    for r in rows:
        key = (_f(r["nu"]), r["arm"])
        d = agg.setdefault(key, {"wall": [], "sel": []})
        d["wall"].append(_f(r["runtime_s"]))
        s = r.get("select_runtime_s", "")
        d["sel"].append(_f(s) if s not in ("", None) else 0.0)
    out = {}
    for (nu, arm), d in sorted(agg.items()):
        wall = np.array(d["wall"], dtype=np.float64)
        sel = np.array(d["sel"], dtype=np.float64)
        out["%g|%s" % (nu, arm)] = {
            "n": int(wall.size),
            "wall_mean_s": float(np.nanmean(wall)),
            "wall_max_s": float(np.nanmax(wall)),
            "select_mean_s": float(np.nanmean(sel)),
            "select_max_s": float(np.nanmax(sel))}
    return out


# ---- summary writers ------------------------------------------------------- #
def _fmt_T(t):
    if t is None or (isinstance(t, float) and not np.isfinite(t)):
        return "--"
    return "%.4g" % t


def _fmt_S(s):
    if s is None or (isinstance(s, float) and not np.isfinite(s)):
        return "--"
    return "%.3f" % s


def _write_md(cfg, cen, mft, rtc, n_rows, path):
    L = []
    L.append("# ROUND63 S1 DEV PILOT SUMMARY (round-6 protocol)")
    L.append("")
    L.append("**S1 is exploratory and development-only** (spec D2 §1): dev images"
             " and seeds are disjoint from all S2 confirmatory images and seeds;"
             " these results do not enter confirmatory intervals or main tables.")
    L.append("")
    L.append("- rows: %d   images: %d   seeds: %s   arms: %s   mode: %s"
             % (n_rows, len(cfg.images), cfg.seeds, ", ".join(cfg.arms), cfg.mode))
    L.append("- side %d, %s, M %d, nu %s, rho_bar %s, tau %.0f ns, sigma_b 0, "
             "active start" % (cfg.side, cfg.pattern, cfg.M,
                               [int(x) for x in cfg.nu_list], cfg.rho_list,
                               cfg.tau * 1e9))
    L.append("- select_iter %d, fista_iters %d, C0 %s, analytic "
             "score-concentration lam_TV rule + descriptive audit"
             % (cfg.select_iter, cfg.fista_iters,
                ("None" if cfg.C0 is None else ("%g" % cfg.C0))))
    L.append("")

    # --- Table 1: frozen floor-censoring taxonomy (round-6 §4)
    L.append("## 1. Time-to-quality speedup S_gate — RQL (round-6 §4 censoring)")
    L.append("")
    L.append("T_opt = M_physical*nu*tau; T_min at the shortest nu tier. "
             "Q*_j = PSNR_rad[RQL, rho=%.3g, nu=%g] - 1 dB. S_gate per the frozen "
             "taxonomy A(normal)/B(fast-left-censored, lower bound)/"
             "C(both-left-censored, =1, unresolved below floor)/"
             "D(safe-left/fast-interior, <1)/E(fast-right-censored, fail=0)/"
             "F(analysis failure=0)."
             % (cen.get("rho_safe", float("nan")), cen.get("nu_ref", float("nan"))))
    L.append("")
    if cen.get("per_image"):
        L.append("| image | status | S_gate | T_safe | T_fast | Q*_j (dB) | seeds |")
        L.append("|---|---|---|---|---|---|---|")
        for d in cen["per_image"]:
            q = d.get("Qstar")
            L.append("| %s | %s | %s | %s | %s | %s | %d |"
                     % (d["image"], d["status"], _fmt_S(d["S_gate"]),
                        _fmt_T(d["T_s"]), _fmt_T(d["T_f"]),
                        ("%.2f" % q) if q is not None else "--",
                        d.get("n_seeds", 0)))
        L.append("| **median S_gate** | | %s | | | | |"
                 % _fmt_S(cen.get("S_median")))
        L.append("")
        sc = cen.get("status_counts", {})
        L.append("Status counts: %s"
                 % (", ".join("%s=%d" % (k, v) for k, v in sorted(sc.items()))
                    or "none"))
        L.append("")
        L.append("Gate stats (S1 dev, descriptive — NOT a confirmatory gate): "
                 "median S_gate=%s; images with S_gate>1 (excl. C/E/F)=%d/%d; "
                 "image-level stratified bootstrap (B=%d) 2.5%% lower bound of "
                 "the median = %s."
                 % (_fmt_S(cen.get("S_median")), cen.get("n_Sgate_gt1", 0),
                    cen.get("n_images", 0), cen.get("bootstrap_B", 0),
                    _fmt_S(cen.get("bootstrap_median_lower_2p5"))))
        L.append("")
        ref = cen.get("reference_gate", {})
        L.append("Reference gate check (median>=3: %s; bootstrap_lower>1: %s) — "
                 "reported for orientation only; S1 never enters the confirmatory "
                 "decision." % (ref.get("median>=3"), ref.get("bootstrap_lower>1")))
        L.append("")
        sec = cen.get("secondary_arms", {})
        if sec:
            L.append("Secondary arms (point-estimate median S_gate): %s"
                     % ", ".join("%s=%s (n>1=%d)"
                                 % (a, _fmt_S(v["S_median"]), v["n_Sgate_gt1"])
                                 for a, v in sec.items()))
            L.append("")
        if cen.get("warnings"):
            L.append("**Cell-completeness warnings:**")
            for w in cen["warnings"]:
                L.append("- %s" % w)
            L.append("")
    else:
        L.append("(no primary-arm PSNR_rad rows — nothing to classify)")
        L.append("")

    # --- Table 2: descriptive audit summary
    L.append("## 2. Descriptive measurement audit by (rho, nu) — RQL cells")
    L.append("")
    L.append("Continuous diagnostics only (round-5 ruling: no binary gate).")
    L.append("")
    if mft:
        L.append("| rho | nu | n | D_ratio med | D_ratio max | leak |")
        L.append("|---|---|---|---|---|---|")
        for k, v in mft.items():
            rho, nu = k.split("|")
            L.append("| %s | %s | %d | %s | %s | %d |"
                     % (rho, nu, v["n"],
                        ("%.3f" % v["d_ratio_med"]) if v["d_ratio_med"]
                        is not None else "n/a",
                        ("%.3f" % v["d_ratio_max"]) if v["d_ratio_max"]
                        is not None else "n/a", v["leak"]))
    else:
        L.append("(no OK-status audit rows)")
    L.append("")

    # --- Table 3: runtime calibration
    L.append("## 3. Runtime calibration per (nu, arm)  [sizes Colab shards]")
    L.append("")
    L.append("| nu | arm | n | wall mean (s) | wall max (s) | "
             "select mean (s) | select max (s) |")
    L.append("|---|---|---|---|---|---|---|")
    for k, v in rtc.items():
        nu, arm = k.split("|")
        L.append("| %s | %s | %d | %.2f | %.2f | %.2f | %.2f |"
                 % (nu, arm, v["n"], v["wall_mean_s"], v["wall_max_s"],
                    v["select_mean_s"], v["select_max_s"]))
    L.append("")
    with open(path, "w") as f:
        f.write("\n".join(L) + "\n")


def _json_default(o):
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def analyze(cfg):
    path = os.path.join(cfg.out, "pilot_rows.csv")
    if not os.path.exists(path):
        raise SystemExit("[s1] no rows CSV at %s (run the sweep first)" % path)
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise SystemExit("[s1] rows CSV %s is empty" % path)

    cen = censoring_analysis(rows, cfg)
    mft = audit_table(rows)
    rtc = runtime_calibration(rows)
    summary = {"config": cfg_dict(cfg), "n_rows": len(rows),
               "censoring_analysis": cen, "model_fail_rates": mft,
               "runtime_calibration": rtc}

    with open(os.path.join(cfg.out, "pilot_summary.json"), "w") as f:
        json.dump(summary, f, indent=2, default=_json_default)
    with open(os.path.join(cfg.out, "runtime_calibration.json"), "w") as f:
        json.dump(rtc, f, indent=2, default=_json_default)
    md = os.path.join(cfg.out, "PILOT_SUMMARY.md")
    _write_md(cfg, cen, mft, rtc, len(rows), md)

    # stdout digest
    print("[s1] ANALYSIS  rows=%d  nu_ref=%s  rho %g->%g"
          % (len(rows), _fmt_T(cen.get("nu_ref")), cen.get("rho_safe", float("nan")),
             cen.get("rho_fast", float("nan"))), flush=True)
    print("       S_median[RQL]=%s  n(S_gate>1 excl C/E/F)=%d/%d  "
          "bootstrap 2.5%% lower=%s"
          % (_fmt_S(cen.get("S_median")), cen.get("n_Sgate_gt1", 0),
             cen.get("n_images", 0), _fmt_S(cen.get("bootstrap_median_lower_2p5"))),
          flush=True)
    sc = cen.get("status_counts", {})
    print("       status: %s"
          % (", ".join("%s=%d" % (k, v) for k, v in sorted(sc.items())) or "none"),
          flush=True)
    n_leak = sum(v["leak"] for v in mft.values())
    n_tot = sum(v["n"] for v in mft.values())
    drs = [v["d_ratio_max"] for v in mft.values() if v["d_ratio_max"]]
    print("       audit: %d OK cells, leak_suspect=%d, max D_ratio=%s"
          % (n_tot, n_leak, ("%.3f" % max(drs)) if drs else "n/a"), flush=True)
    if rtc:
        wmax = max(v["wall_max_s"] for v in rtc.values())
        print("       max per-fit wall = %.2fs (see runtime_calibration.json)"
              % wmax, flush=True)
    print("[s1] wrote %s , pilot_summary.json , runtime_calibration.json"
          % md, flush=True)
    return summary


# ---- driver ---------------------------------------------------------------- #
def _resolve_c0(cfg, need_sweep):
    """Set cfg.C0 for normal / pass-b runs (pass-a needs no C0)."""
    if cfg.mode == "pass_b":
        fc = frozen_C0()
        if need_sweep and fc is None:
            raise SystemExit(
                "[s1] pass B needs a frozen C0 but %s is missing — run "
                "`--pass-a --freeze-c0` first." % C0_FILE)
        cfg.C0 = fc
    elif cfg.mode == "normal":
        if cfg.C0 is None:
            cfg.C0 = frozen_C0()
        if need_sweep and cfg.C0 is None:
            raise SystemExit(
                "[s1] no C0 available: pass --C0 X for a dev/smoke run, or "
                "freeze one via `--pass-a --freeze-c0`.")


def main(argv=None):
    a = parse_args(sys.argv[1:] if argv is None else argv)
    if a.pass_a and a.pass_b:
        raise SystemExit("[s1] choose at most one of --pass-a / --pass-b")
    cfg = build_cfg(a)
    os.makedirs(cfg.out, exist_ok=True)
    need_sweep = not a.analyze_only

    if cfg.mode == "pass_a":
        print("[s1:A] out=%s cells=%d (%d rho x %d nu x %d seed) x 2 endpoints"
              % (cfg.out, len(cfg.rho_list) * len(cfg.nu_list) * len(cfg.seeds),
                 len(cfg.rho_list), len(cfg.nu_list), len(cfg.seeds)), flush=True)
        if need_sweep:
            sweep_passa(cfg)
        if not a.no_analysis:
            analyze_passa(cfg, freeze=cfg.freeze_c0)
        print("[s1] DONE", flush=True)
        return 0

    _resolve_c0(cfg, need_sweep)
    print("[s1] out=%s mode=%s smoke=%s C0=%s cells=%d (%d rho x %d nu x %d seed)"
          % (cfg.out, cfg.mode, cfg.smoke,
             ("None" if cfg.C0 is None else ("%g" % cfg.C0)),
             len(cfg.rho_list) * len(cfg.nu_list) * len(cfg.seeds),
             len(cfg.rho_list), len(cfg.nu_list), len(cfg.seeds)), flush=True)
    if need_sweep:
        sweep(cfg)
    if not a.no_analysis:
        analyze(cfg)
    print("[s1] DONE", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
