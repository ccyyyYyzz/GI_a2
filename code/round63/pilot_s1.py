"""ROUND63 S1 DEV PILOT (local, development-only).

SWEEP / Pass-A machinery = GPT round-6 ruling; the ANALYSIS layer is the GPT
round-7 FINAL endpoint (docs/ROUND63_GPT_ROUND7_RULING.md §2-4): per-image
seed-mean + equal-weight PAVA isotonic curves, the Q90 = Qiso_safe(T_min) +
0.90*(Qiso_safe(T_max)-Qiso_safe(T_min)) target (retires the old
Q_safe(nu_max) - 1 dB target), the augmented censoring taxonomy, a 10,000-draw
nested family-stratified bootstrap, and the fixed-budget quality-gain secondary.
The analysis functions (pava/judge_image/fixed_budget_gain/bootstrap_endpoints/
endpoint_analysis) are a clean module-level API reused verbatim for S2.

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

PRIMARY_ARM = "RQL"                                 # round-7 §2 primary gate S_gate

# the two endpoint c weights the Pass-A calibration reconstructs at every cell
C_ENDPOINTS = (0.25, 0.50)

# ---- round-7 FINAL endpoint constants (docs/ROUND63_GPT_ROUND7_RULING.md) --- #
# NOTE: the ANALYSIS layer below is the round-7 Q90 endpoint (retires the old
# Q* = Q_safe(nu_max) - 1 dB target); the SWEEP / Pass-A machinery is unchanged
# from round-6. This same analysis code is reused verbatim for S2 confirmatory.
R_MIN_DB = 0.50            # §2.3 safe-range informativeness floor (dB)
Q90_FRAC = 0.90           # §2.3 fraction of the safe dynamic range to recover
NU_BUDGET = 2000.0        # §3 fixed-budget terminal dwell tier
MEDIAN_MIN = 3.0          # §2.5 primary gate: median S_gate >= 3
N_GT1_MIN = 18            # §2.5 primary gate: >=18 images with S_gate>1 (of 24)
DQ_MEDIAN_MIN = 1.0       # §3  secondary: median DeltaQ >= 1.0 dB
BOOT_B_DEFAULT = 10000    # §2.5 nested family-stratified replicate count
BOOT_STREAM = (0, 63, 6, 2)  # §2.5 rng_for(0,63,6,2) bootstrap stream

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

# gate bookkeeping (round-7 §2.5): NON-POSITIVE statuses excluded from the
# "S_gate>1" tally. BOTH_LEFT_CENSORED (S_gate fixed 1.0, unresolved below floor)
# and SAFE_RANGE_UNINFORMATIVE (S_gate 1.0, safe curve too flat to test) are
# both pinned at 1.0 and never counted positive; FAST_RIGHT_CENSORED / ANALYSIS
# _FAILURE sit at 0. SAFE_LEFT_FAST_INTERIOR naturally has S_gate<1 (not >1).
_NOT_POSITIVE = {"BOTH_LEFT_CENSORED", "FAST_RIGHT_CENSORED", "ANALYSIS_FAILURE",
                 "SAFE_RANGE_UNINFORMATIVE"}
_TIE_DB = 0.02                                      # round-6 §1.6 C0 tie window

# Family map (round-7 §2.5): prefer detail24.family_of when that module is
# importable; otherwise parse detail_<family>_<i>. detail24.py is authored
# concurrently, so the import is best-effort and must never hard-fail here.
try:                                                # noqa: SIM105
    from detail24 import family_of as _family_of_detail
except Exception:                                   # pragma: no cover - fallback
    _family_of_detail = None


def family_of(image):
    """Family label for a DETAIL-24 image name detail_<family>_<i>.
    Delegates to detail24.family_of when present, else name-prefix parsing."""
    if _family_of_detail is not None:
        try:
            return _family_of_detail(image)
        except Exception:
            pass
    parts = str(image).split("_")
    if len(parts) >= 3 and parts[0] == "detail":
        return "_".join(parts[1:-1])
    return str(image)


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
    ap.add_argument("--cohort", type=str, default="dev",
                    help="analysis cohort label (default dev). Only 'detail' "
                         "claims a positive DETAIL_SPEED_PASS gate; every other "
                         "label (dev/natural) reports the same stats with the "
                         "gate left null (round-7 §4).")
    ap.add_argument("--boot-B", type=int, default=BOOT_B_DEFAULT, dest="boot_B",
                    help="nested family-stratified bootstrap replicates "
                         "(default %d; lower it for smokes)" % BOOT_B_DEFAULT)
    ap.add_argument("--selftest-analysis", action="store_true",
                    dest="selftest_analysis",
                    help="run the PAVA + censoring-classifier structural self "
                         "test and exit (no data needed)")
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
        C0=a.C0, cohort=str(getattr(a, "cohort", "dev")),
        boot_B=int(getattr(a, "boot_B", BOOT_B_DEFAULT)))


def cfg_dict(cfg):
    return {"side": cfg.side, "M": cfg.M, "pattern": cfg.pattern, "tau": cfg.tau,
            "nu_list": cfg.nu_list, "rho_list": cfg.rho_list, "arms": cfg.arms,
            "images": cfg.images, "seeds": cfg.seeds,
            "select_iter": cfg.select_iter, "fista_iters": cfg.fista_iters,
            "smoke": cfg.smoke, "mode": cfg.mode,
            "cohort": getattr(cfg, "cohort", "dev"),
            "boot_B": int(getattr(cfg, "boot_B", BOOT_B_DEFAULT)),
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


# ---- round-7 §2 endpoint: Q90 target + PAVA isotonic + censoring ----------- #
def pava(y):
    """Ordinary equal-weight pool-adjacent-violators isotonic regression:
    the least-squares nondecreasing fit to y in its given order (§2.2). Uses
    only the ordering (equal weights), not the x spacing. Returns a float64
    array the same length/order as y."""
    y = np.asarray(y, dtype=np.float64)
    means, weights, counts = [], [], []
    for val in y:
        m, w, c = float(val), 1.0, 1
        while means and means[-1] > m:                 # pool the violated block
            pm, pw, pc = means.pop(), weights.pop(), counts.pop()
            m = (pm * pw + m * w) / (pw + w)
            w += pw
            c += pc
        means.append(m)
        weights.append(w)
        counts.append(c)
    out = np.empty(y.size, dtype=np.float64)
    i = 0
    for m, c in zip(means, counts):
        out[i:i + c] = m
        i += c
    return out


def _t_reach_interior(T, P, Q):
    """Interior crossing time by linear interpolation in log10(T) on the
    monotone curve P (assumes P[0] < Q <= P.max(), so k>=1)."""
    k = int(np.nonzero(P >= Q)[0][0])
    logT = np.log10(T)
    if P[k] == P[k - 1]:
        lt = logT[k]
    else:
        lt = logT[k - 1] + (Q - P[k - 1]) * (logT[k] - logT[k - 1]) \
            / (P[k] - P[k - 1])
    return float(10.0 ** lt)


def _seed_mean_curve(seed_map, seeds, T_by_nu):
    """Raw seed-mean PSNR_rad curve over the given seed MULTISET (duplicates
    weight naturally), sorted by ascending T_opt. Returns (T, P_raw) or None."""
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
    return T[order], P[order]


def _iso_curve(raw):
    """(T, P_raw) -> (T, P_iso) with equal-weight PAVA of P_raw in T order
    (§2.2: nondecreasing Q~(log T_opt)). None passes through."""
    if raw is None:
        return None
    T, P = raw
    return T, pava(P)


def _reach_state(curve, Q):
    """('RIGHT'|'FLOOR'|'INTERIOR', T_c). RIGHT = never reaches Q; FLOOR =
    already at/above Q at the shortest T (left-censored); INTERIOR = crosses
    strictly inside the grid. Operates on the isotonic curve."""
    if curve is None:
        return "RIGHT", float("nan")
    T, P = curve
    if float(P.max()) < Q:
        return "RIGHT", float("nan")
    if float(P[0]) >= Q:
        return "FLOOR", float(T[0])
    return "INTERIOR", _t_reach_interior(T, P, Q)


def _target(safe_iso):
    """Round-7 §2.3 per-image target from the SAFE isotonic curve.
    Returns (status, R_j, Q90, Qmin, Qmax): status is None on the normal path,
    'SAFE_RANGE_UNINFORMATIVE' when R_j<0.50 dB, 'ANALYSIS_FAILURE' when the
    safe curve is missing."""
    if safe_iso is None:
        return "ANALYSIS_FAILURE", None, None, None, None
    _, P = safe_iso
    Qmin = float(P[0])                 # Qiso_safe(T_min)
    Qmax = float(P[-1])                # Qiso_safe(T_max)
    R = Qmax - Qmin                    # observable safe-side dynamic range
    if R < R_MIN_DB:
        return "SAFE_RANGE_UNINFORMATIVE", R, None, Qmin, Qmax
    return None, R, Qmin + Q90_FRAC * R, Qmin, Qmax


def _classify(safe_iso, fast_iso, Q, T_min):
    """Round-7 §2.4 reach-state classifier of the (safe,fast) isotonic curves
    against threshold Q. Returns (status, S_gate, T_s, T_f). Kept separate from
    the target so each censoring class is directly unit-testable."""
    st_s, T_s = _reach_state(safe_iso, Q)
    st_f, T_f = _reach_state(fast_iso, Q)
    # missing/incomplete safe curve: analysis-completeness failure.
    if st_s == "RIGHT":
        return "ANALYSIS_FAILURE", 0.0, T_s, T_f
    # fast never reaches Q by nu_max -> failure (S_gate 0).
    if st_f == "RIGHT":
        return "FAST_RIGHT_CENSORED", 0.0, T_s, T_f
    if st_s == "INTERIOR" and st_f == "INTERIOR":
        s = (T_s / T_f) if T_f > 0 else 0.0
        return "NORMAL", float(s), T_s, T_f
    if st_s == "INTERIOR" and st_f == "FLOOR":
        s = (T_s / T_min) if T_min > 0 else 0.0        # >= lower bound
        return "FAST_LEFT_CENSORED", float(s), T_s, T_f
    if st_s == "FLOOR" and st_f == "INTERIOR":
        s = (T_min / T_f) if T_f > 0 else 0.0          # < 1: non-positive
        return "SAFE_LEFT_FAST_INTERIOR", float(s), T_s, T_f
    # both FLOOR: unresolved below floor, S_gate fixed to 1.0.
    return "BOTH_LEFT_CENSORED", 1.0, T_s, T_f


def judge_image(entry, seeds, T_by_nu, T_min):
    """Full round-7 §2.2-2.4 per-image judge (reused verbatim by S2 and by the
    bootstrap). Builds seed-mean + isotonic curves for both operating points
    over the seed multiset, derives the Q90 target from the safe curve, and
    classifies. Returns a dict with status, S_gate, T_s, T_f, R_j, Q90, Qmin,
    Qmax and the raw/isotonic curves (curves are kept for reporting only)."""
    safe_raw = _seed_mean_curve(entry["safe"], seeds, T_by_nu)
    fast_raw = _seed_mean_curve(entry["fast"], seeds, T_by_nu)
    safe_iso = _iso_curve(safe_raw)
    fast_iso = _iso_curve(fast_raw)
    tstatus, R, Q90, Qmin, Qmax = _target(safe_iso)
    base = {"R_j": R, "Q90": Q90, "Qmin": Qmin, "Qmax": Qmax,
            "safe_raw": safe_raw, "fast_raw": fast_raw,
            "safe_iso": safe_iso, "fast_iso": fast_iso}
    if tstatus == "ANALYSIS_FAILURE":
        return {"status": "ANALYSIS_FAILURE", "S_gate": 0.0,
                "T_s": None, "T_f": None, **base}
    if tstatus == "SAFE_RANGE_UNINFORMATIVE":
        return {"status": "SAFE_RANGE_UNINFORMATIVE", "S_gate": 1.0,
                "T_s": None, "T_f": None, **base}
    status, S_gate, T_s, T_f = _classify(safe_iso, fast_iso, Q90, T_min)
    return {"status": status, "S_gate": float(S_gate),
            "T_s": T_s, "T_f": T_f, **base}


def fixed_budget_gain(entry, seeds, nu_budget):
    """Round-7 §3 fixed-budget quality gain DeltaQ_j = RAW seed-mean PSNR_rad
    (NOT isotonic) at the terminal dwell nu_budget: fast minus safe. None when
    that grid point is missing for either operating point."""
    def raw_at(side_map):
        d = side_map.get(nu_budget)
        if not d:
            return None
        vals = [d[s] for s in seeds if s in d and d[s] is not None]
        return float(np.mean(vals)) if vals else None
    qs = raw_at(entry["safe"])
    qf = raw_at(entry["fast"])
    if qs is None or qf is None:
        return None
    return float(qf - qs)


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


def _ser_curve(c):
    """(T, P) numpy tuple -> JSON-friendly {"T":[...], "P":[...]} (or None).
    Used to preserve BOTH the raw seed-mean and the isotonic curve (§2.2)."""
    if c is None:
        return None
    T, P = c
    return {"T": [float(x) for x in T], "P": [float(x) for x in P]}


def _cohort_family(image, cohort):
    """DETAIL-24 uses the family_of(...) grouping; every other cohort treats
    each image as its own family (§2.5 degenerate stratification)."""
    return family_of(image) if cohort == "detail" else str(image)


def _build_strata(images, cohort):
    """Bootstrap strata (§2.5). DETAIL-24: {family: [its images]} so each draw
    resamples the images WITHIN each of the six families. Any non-detail cohort
    (dev/natural): one stratum of all images => a plain image bootstrap."""
    if cohort == "detail":
        strata = {}
        for im in images:
            strata.setdefault(family_of(im), []).append(im)
        return strata
    return {"__all__": list(images)}


def bootstrap_endpoints(images, strata, seeds_by_image, stat_fns, B):
    """Round-7 §2.5 nested family-stratified bootstrap of the cohort median,
    computed for EVERY statistic in stat_fns off the SAME resamples (a paired
    bootstrap: identical seed/image draws, different per-image statistics).
    Each replicate: (1) resample each image's seeds with replacement and
    recompute its statistic (which REBUILDS curves/target/censoring); (2) within
    each stratum resample its images with replacement; (3) take the cohort
    median. Deterministic stream rng_for(0,63,6,2). Returns
    {name: {"LB2.5":.., "n_finite":.., "meds": ndarray}}."""
    names = list(stat_fns)
    meds = {nm: np.empty(B, dtype=np.float64) for nm in names}
    if not images:
        for nm in names:
            meds[nm][:] = np.nan
        return {nm: {"LB2.5": float("nan"), "n_finite": 0, "meds": meds[nm]}
                for nm in names}
    rng = rng_for(*BOOT_STREAM)
    imgs_sorted = sorted(images)
    strata_sorted = sorted(strata.items())
    for d in range(B):
        per = {nm: {} for nm in names}
        for im in imgs_sorted:                          # (1) seed resample
            sl = seeds_by_image.get(im, [])
            if sl:
                drawn = [int(s) for s in rng.choice(sl, size=len(sl),
                                                    replace=True)]
            else:
                drawn = []
            for nm in names:
                per[nm][im] = stat_fns[nm](im, drawn)
        picked = []                                     # (2) image resample
        for _, fam_imgs in strata_sorted:
            k = len(fam_imgs)
            if k == 0:
                continue
            idx = rng.integers(0, k, size=k)
            picked.extend(fam_imgs[i] for i in idx)
        for nm in names:                                # (3) cohort median
            vals = [per[nm][im] for im in picked
                    if per[nm][im] is not None
                    and np.isfinite(per[nm][im])]
            meds[nm][d] = float(np.median(vals)) if vals else np.nan
    out = {}
    for nm in names:
        arr = meds[nm]
        finite = arr[np.isfinite(arr)]
        out[nm] = {"LB2.5": (float(np.quantile(finite, 0.025))
                             if finite.size else float("nan")),
                   "n_finite": int(finite.size), "meds": arr}
    return out


def endpoint_analysis(rows, cfg):
    """Round-7 FINAL endpoint analysis for one cohort. cfg.cohort selects the
    family stratification and whether a positive gate is claimed:
    cohort=='detail' -> DETAIL_SPEED_PASS true/false; every other cohort
    (dev/natural) reports the identical statistics but leaves the gate null
    (§4 no positive pass gate)."""
    cohort = getattr(cfg, "cohort", "dev")
    boot_B = int(getattr(cfg, "boot_B", BOOT_B_DEFAULT))
    rho_safe = min(cfg.rho_list)
    rho_fast = max(cfg.rho_list)
    # v-bump 2026-07-18 (outcome-blind defect fix, freeze clause "demonstrable
    # implementation defects"): the image list was pinned to cfg.images (the
    # DEV default names), so a confirmatory CSV whose cohort images differ
    # (detail_*/stl_*) yielded ZERO analyzable images — structurally visible
    # as an empty endpoint block, independent of any result value. The
    # confirmatory analysis derives its image set from the ROWS themselves;
    # an explicit cfg.images intersection is kept when it is non-empty (the
    # dev pilot path, unchanged behavior).
    row_images = sorted({r.get("image") for r in rows
                         if r.get("arm") == PRIMARY_ARM and r.get("image")})
    images = [im for im in cfg.images if im in row_images]
    if not images:
        images = row_images
    data, T_by_nu, seeds_by_image = _collect(rows, PRIMARY_ARM, rho_safe,
                                             rho_fast, images)
    if not T_by_nu or not images:
        return {"cohort": cohort, "rho_safe": rho_safe, "rho_fast": rho_fast,
                "images": images, "per_image": [], "S_median": None,
                "n_Sgate_gt1": 0, "n_images": len(images),
                "primary_gate": {"DETAIL_SPEED_PASS": None},
                "fixed_budget_quality_gain": {"pass": None},
                "note": "no PSNR_rad rows for the primary arm"}
    T_min = float(min(T_by_nu.values()))
    nus = set(T_by_nu)
    nu_budget = NU_BUDGET if NU_BUDGET in nus else max(nus)

    # ---- point estimates (actual seeds, no resampling)
    per_image = []
    warnings = []
    for im in images:
        j = judge_image(data[im], seeds_by_image[im], T_by_nu, T_min)
        dq = fixed_budget_gain(data[im], seeds_by_image[im], nu_budget)
        per_image.append({
            "image": im, "family": _cohort_family(im, cohort),
            "status": j["status"], "S_gate": j["S_gate"],
            "T_s": j["T_s"], "T_f": j["T_f"], "R_j": j["R_j"],
            "Q90_j": j["Q90"], "Qmin": j["Qmin"], "Qmax": j["Qmax"],
            "DeltaQ_budget": dq, "n_seeds": len(seeds_by_image[im]),
            # §2.2: keep the raw seed-mean curve alongside the isotonic one.
            "curves": {"safe_raw": _ser_curve(j["safe_raw"]),
                       "safe_iso": _ser_curve(j["safe_iso"]),
                       "fast_raw": _ser_curve(j["fast_raw"]),
                       "fast_iso": _ser_curve(j["fast_iso"])}})
        if j["status"] == "ANALYSIS_FAILURE":
            w = ("[s1] LOUD WARNING analysis-completeness: image %s -> "
                 "ANALYSIS_FAILURE (safe curve missing/incomplete; round-7 "
                 "§2.4). No image is deleted; check the cell grid." % im)
            warnings.append(w)
            print(w, flush=True)

    S = [d["S_gate"] for d in per_image]
    S_median = float(np.median(S)) if S else float("nan")
    n_gt1 = sum(1 for d in per_image
                if d["status"] not in _NOT_POSITIVE and d["S_gate"] > 1.0)

    DQ = [d["DeltaQ_budget"] for d in per_image if d["DeltaQ_budget"] is not None]
    DQ_median = float(np.median(DQ)) if DQ else float("nan")
    n_dq_pos = sum(1 for d in per_image
                   if d["DeltaQ_budget"] is not None and d["DeltaQ_budget"] > 0)
    n_dq_missing = sum(1 for d in per_image if d["DeltaQ_budget"] is None)

    # ---- nested family-stratified bootstrap (both endpoints, paired draws)
    strata = _build_strata(images, cohort)
    stat_fns = {
        "S_gate": (lambda im, sd:
                   judge_image(data[im], sd, T_by_nu, T_min)["S_gate"]),
        "DeltaQ": (lambda im, sd:
                   fixed_budget_gain(data[im], sd, nu_budget)),
    }
    boot = bootstrap_endpoints(images, strata, seeds_by_image, stat_fns, boot_B)
    LB_S = boot["S_gate"]["LB2.5"]
    LB_DQ = boot["DeltaQ"]["LB2.5"]

    is_detail = (cohort == "detail")
    primary_ok = bool(S_median >= MEDIAN_MIN and np.isfinite(LB_S)
                      and LB_S > 1.0 and n_gt1 >= N_GT1_MIN)
    secondary_ok = bool(DQ_median >= DQ_MEDIAN_MIN and np.isfinite(LB_DQ)
                        and LB_DQ > 0.0 and n_dq_pos >= N_GT1_MIN)

    n_status = {}
    for d in per_image:
        n_status[d["status"]] = n_status.get(d["status"], 0) + 1

    return {
        "cohort": cohort, "rho_safe": rho_safe, "rho_fast": rho_fast,
        "nu_budget": nu_budget, "T_min": T_min,
        "T_by_nu": {("%g" % k): v for k, v in sorted(T_by_nu.items())},
        "images": images, "per_image": per_image, "status_counts": n_status,
        "n_images": len(images),
        "S_median": S_median, "n_Sgate_gt1": n_gt1,
        "bootstrap": {
            "B": boot_B, "stream": list(BOOT_STREAM), "stratified_by": (
                "family" if is_detail else "image (degenerate)"),
            "S_gate_LB2.5": LB_S, "S_gate_n_finite": boot["S_gate"]["n_finite"],
            "DeltaQ_LB2.5": LB_DQ, "DeltaQ_n_finite": boot["DeltaQ"]["n_finite"]},
        "primary_gate": {
            "endpoint": "time_to_recover_90pct_safe_range (Q90)",
            "median_S_gate": S_median, "median>=3": bool(S_median >= MEDIAN_MIN),
            "bootstrap_LB2.5": LB_S, "LB>1": bool(np.isfinite(LB_S) and LB_S > 1.0),
            "n_Sgate_gt1": n_gt1, "n_gt1>=%d" % N_GT1_MIN: bool(n_gt1 >= N_GT1_MIN),
            "n_images": len(images),
            "DETAIL_SPEED_PASS": (primary_ok if is_detail else None)},
        "fixed_budget_quality_gain": {
            "endpoint": "DeltaQ_budget = RAW seed-mean PSNR_rad(fast)-(safe) @ nu=%g"
                        % nu_budget,
            "median_DeltaQ_dB": DQ_median,
            "median>=1.0": bool(DQ_median >= DQ_MEDIAN_MIN),
            "bootstrap_LB2.5": LB_DQ, "LB>0": bool(np.isfinite(LB_DQ) and LB_DQ > 0.0),
            "n_DeltaQ_pos": n_dq_pos, "n_pos>=%d" % N_GT1_MIN: bool(n_dq_pos >= N_GT1_MIN),
            "n_missing": n_dq_missing,
            "pass": (secondary_ok if is_detail else None)},
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


def _fmt_dB(v):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return "--"
    return "%.2f" % v


def _write_md(cfg, cen, mft, rtc, n_rows, path):
    cohort = cen.get("cohort", getattr(cfg, "cohort", "dev"))
    pg = cen.get("primary_gate", {})
    fb = cen.get("fixed_budget_quality_gain", {})
    bt = cen.get("bootstrap", {})
    L = []
    L.append("# ROUND63 S1 PILOT SUMMARY (round-7 Q90 endpoint, cohort=%s)"
             % cohort)
    L.append("")
    L.append("**S1 is exploratory and development-only** (spec D2 §1): dev images"
             " and seeds are disjoint from all S2 confirmatory images and seeds;"
             " these results do not enter confirmatory intervals or main tables."
             " The analysis layer here is the round-7 FINAL endpoint reused"
             " verbatim for S2.")
    L.append("")
    L.append("- rows: %d   images: %d   seeds: %s   arms: %s   mode: %s   "
             "cohort: %s"
             % (n_rows, len(cfg.images), cfg.seeds, ", ".join(cfg.arms),
                cfg.mode, cohort))
    L.append("- side %d, %s, M %d, nu %s, rho_bar %s, tau %.0f ns, sigma_b 0, "
             "active start" % (cfg.side, cfg.pattern, cfg.M,
                               [int(x) for x in cfg.nu_list], cfg.rho_list,
                               cfg.tau * 1e9))
    L.append("- select_iter %d, fista_iters %d, C0 %s, analytic "
             "score-concentration lam_TV rule + descriptive audit"
             % (cfg.select_iter, cfg.fista_iters,
                ("None" if cfg.C0 is None else ("%g" % cfg.C0))))
    L.append("")

    # --- Table 1: round-7 Q90 endpoint + censoring taxonomy (§2.2-2.5)
    L.append("## 1. Acquisition-speed endpoint S_gate — RQL (round-7 §2 Q90)")
    L.append("")
    L.append("Per (image, rho): seed-mean PSNR_rad at each nu -> equal-weight "
             "PAVA isotonic fit vs log T_opt (raw curve kept in the JSON). "
             "Safe = rho=%.3g; fast = rho=%.3g. R_j = Qiso_safe(T_max) - "
             "Qiso_safe(T_min); if R_j<%.2f dB the safe range is uninformative "
             "(S_gate=1, not positive); else Q90_j = Qiso_safe(T_min) + "
             "%.2f*R_j and crossing times come by log-T interpolation on the "
             "isotonic curves. Censoring: NORMAL(=T_s/T_f)/FAST_LEFT_CENSORED"
             "(>=T_s/T_min)/BOTH_LEFT_CENSORED(=1, unresolved)/"
             "SAFE_LEFT_FAST_INTERIOR(<1)/FAST_RIGHT_CENSORED(fail=0)/"
             "SAFE_RANGE_UNINFORMATIVE(=1)/ANALYSIS_FAILURE(0). No image deleted."
             % (cen.get("rho_safe", float("nan")),
                cen.get("rho_fast", float("nan")), R_MIN_DB, Q90_FRAC))
    L.append("")
    if cen.get("per_image"):
        L.append("| image | family | status | S_gate | T_safe | T_fast | "
                 "R_j (dB) | Q90_j (dB) | DeltaQ_budget (dB) | seeds |")
        L.append("|---|---|---|---|---|---|---|---|---|---|")
        for d in cen["per_image"]:
            L.append("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %d |"
                     % (d["image"], d.get("family", "--"), d["status"],
                        _fmt_S(d["S_gate"]), _fmt_T(d["T_s"]), _fmt_T(d["T_f"]),
                        _fmt_dB(d.get("R_j")), _fmt_dB(d.get("Q90_j")),
                        _fmt_dB(d.get("DeltaQ_budget")), d.get("n_seeds", 0)))
        L.append("| **median S_gate** | | | %s | | | | | | |"
                 % _fmt_S(cen.get("S_median")))
        L.append("")
        sc = cen.get("status_counts", {})
        L.append("Status counts: %s"
                 % (", ".join("%s=%d" % (k, v) for k, v in sorted(sc.items()))
                    or "none"))
        L.append("")

        # primary gate block
        pass_val = pg.get("DETAIL_SPEED_PASS", None)
        pass_str = ("null (non-detail cohort: no positive gate, §4)"
                    if pass_val is None else str(pass_val))
        L.append("### Primary gate (§2.5) — DETAIL_SPEED_PASS = %s" % pass_str)
        L.append("")
        L.append("- median S_gate = %s  (>=%.0f: %s)"
                 % (_fmt_S(cen.get("S_median")), MEDIAN_MIN, pg.get("median>=3")))
        L.append("- bootstrap 2.5%% LB of median S_gate = %s  (>1: %s)"
                 % (_fmt_S(bt.get("S_gate_LB2.5")), pg.get("LB>1")))
        L.append("- images with S_gate>1 (excl. non-positive classes) = %d/%d  "
                 "(>=%d: %s)"
                 % (cen.get("n_Sgate_gt1", 0), cen.get("n_images", 0),
                    N_GT1_MIN, pg.get("n_gt1>=%d" % N_GT1_MIN)))
        L.append("")

        # secondary block
        L.append("### Fixed-budget quality gain (§3, key secondary)")
        L.append("")
        L.append("DeltaQ_budget = RAW seed-mean PSNR_rad(fast) - (safe) at "
                 "nu=%g (NOT isotonic). Descriptive success: median>=%.1f dB "
                 "AND LB2.5>0 AND >=%d images with DeltaQ>0."
                 % (cen.get("nu_budget", NU_BUDGET), DQ_MEDIAN_MIN, N_GT1_MIN))
        L.append("")
        L.append("- median DeltaQ = %s dB  (>=%.1f: %s)"
                 % (_fmt_dB(fb.get("median_DeltaQ_dB")), DQ_MEDIAN_MIN,
                    fb.get("median>=1.0")))
        L.append("- bootstrap 2.5%% LB of median DeltaQ = %s  (>0: %s)"
                 % (_fmt_dB(bt.get("DeltaQ_LB2.5")), fb.get("LB>0")))
        L.append("- images with DeltaQ>0 = %d/%d  (>=%d: %s)   pass = %s"
                 % (fb.get("n_DeltaQ_pos", 0), cen.get("n_images", 0),
                    N_GT1_MIN, fb.get("n_pos>=%d" % N_GT1_MIN),
                    ("null (non-detail cohort)" if fb.get("pass") is None
                     else fb.get("pass"))))
        L.append("")

        # bootstrap diagnostics
        L.append("Bootstrap: B=%s draws, stream rng_for%s, stratified by %s; "
                 "finite draws S_gate=%s / DeltaQ=%s."
                 % (bt.get("B"), tuple(bt.get("stream", BOOT_STREAM)),
                    bt.get("stratified_by"), bt.get("S_gate_n_finite"),
                    bt.get("DeltaQ_n_finite")))
        L.append("")
        if cen.get("warnings"):
            L.append("**Analysis-completeness warnings:**")
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

    cen = endpoint_analysis(rows, cfg)
    mft = audit_table(rows)
    rtc = runtime_calibration(rows)
    summary = {"config": cfg_dict(cfg), "n_rows": len(rows),
               "endpoint_analysis": cen, "model_fail_rates": mft,
               "runtime_calibration": rtc}

    with open(os.path.join(cfg.out, "pilot_summary.json"), "w") as f:
        json.dump(summary, f, indent=2, default=_json_default)
    with open(os.path.join(cfg.out, "runtime_calibration.json"), "w") as f:
        json.dump(rtc, f, indent=2, default=_json_default)
    md = os.path.join(cfg.out, "PILOT_SUMMARY.md")
    _write_md(cfg, cen, mft, rtc, len(rows), md)

    # stdout digest
    pg = cen.get("primary_gate", {})
    fb = cen.get("fixed_budget_quality_gain", {})
    bt = cen.get("bootstrap", {})
    print("[s1] ANALYSIS  rows=%d  cohort=%s  nu_budget=%s  rho %g->%g"
          % (len(rows), cen.get("cohort"), _fmt_T(cen.get("nu_budget")),
             cen.get("rho_safe", float("nan")),
             cen.get("rho_fast", float("nan"))), flush=True)
    print("       primary Q90: median S_gate=%s  n(S_gate>1)=%d/%d  "
          "boot LB2.5=%s  DETAIL_SPEED_PASS=%s"
          % (_fmt_S(cen.get("S_median")), cen.get("n_Sgate_gt1", 0),
             cen.get("n_images", 0), _fmt_S(bt.get("S_gate_LB2.5")),
             pg.get("DETAIL_SPEED_PASS")), flush=True)
    print("       fixed-budget: median DeltaQ=%s dB  n(DeltaQ>0)=%d/%d  "
          "boot LB2.5=%s  pass=%s"
          % (_fmt_dB(fb.get("median_DeltaQ_dB")), fb.get("n_DeltaQ_pos", 0),
             cen.get("n_images", 0), _fmt_dB(bt.get("DeltaQ_LB2.5")),
             fb.get("pass")), flush=True)
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


# ---- structural self-test (round-7 analysis; no data needed) --------------- #
def _selftest():
    """Structural asserts for PAVA, the censoring classifier (one synthetic
    fixture per class incl. SAFE_RANGE_UNINFORMATIVE), the Q90 target, the
    fixed-budget gain and the family-stratified bootstrap. Verifies mechanics,
    NOT science. Run via `pilot_s1.py --selftest-analysis`."""
    n_ok = 0

    # --- PAVA hand cases (equal weight) -----------------------------------
    assert np.allclose(pava([1, 3, 2, 4]), [1, 2.5, 2.5, 4]), "pava violator pool"
    assert np.allclose(pava([3, 2, 1]), [2, 2, 2]), "pava full decreasing"
    assert np.allclose(pava([1, 2, 3]), [1, 2, 3]), "pava already isotonic"
    assert np.allclose(pava([5]), [5]), "pava singleton"
    assert pava([]).size == 0, "pava empty"
    n_ok += 1

    # --- reach-state censoring classifier (Q=1.5, T grid 1..16) -----------
    T = np.array([1., 2., 4., 8., 16.])
    T_min = 1.0
    Q = 1.5
    def iso(vals):
        return (T, pava(np.asarray(vals, dtype=float)))
    interior = iso([0, 1, 2, 3, 4])          # crosses 1.5 in the interior
    floor = iso([2, 3, 4, 5, 6])             # already >= 1.5 at T_min
    never = iso([0, 0.1, 0.2, 0.3, 0.4])     # max < 1.5

    st, sg, _, _ = _classify(interior, iso([0, 2, 3, 4, 5]), Q, T_min)
    assert st == "NORMAL" and np.isfinite(sg) and sg > 0, "NORMAL"
    st, sg, _, _ = _classify(interior, floor, Q, T_min)
    assert st == "FAST_LEFT_CENSORED" and sg > 0, "FAST_LEFT_CENSORED"
    st, sg, _, _ = _classify(floor, floor, Q, T_min)
    assert st == "BOTH_LEFT_CENSORED" and sg == 1.0, "BOTH_LEFT_CENSORED"
    st, sg, _, _ = _classify(floor, interior, Q, T_min)
    assert st == "SAFE_LEFT_FAST_INTERIOR" and 0 <= sg < 1, "SAFE_LEFT_FAST_INTERIOR"
    st, sg, _, _ = _classify(interior, never, Q, T_min)
    assert st == "FAST_RIGHT_CENSORED" and sg == 0.0, "FAST_RIGHT_CENSORED"
    st, sg, _, _ = _classify(never, interior, Q, T_min)
    assert st == "ANALYSIS_FAILURE" and sg == 0.0, "ANALYSIS_FAILURE(safe RIGHT)"
    st, sg, _, _ = _classify(None, interior, Q, T_min)
    assert st == "ANALYSIS_FAILURE" and sg == 0.0, "ANALYSIS_FAILURE(safe None)"
    n_ok += 1

    # --- judge_image: target construction + the two pinned statuses -------
    Tn = {5.: 1., 10.: 2., 20.: 4., 50.: 8., 100.: 16.}
    def cv(vals):                            # 2-seed curve, identical seeds
        return {nu: {0: v, 1: v} for nu, v in zip(sorted(Tn), vals)}

    flat = {"safe": cv([10, 10, 10, 10, 10]), "fast": cv([20, 20, 20, 20, 20])}
    j = judge_image(flat, [0, 1], Tn, 1.0)
    assert j["status"] == "SAFE_RANGE_UNINFORMATIVE" and j["S_gate"] == 1.0, \
        "SAFE_RANGE_UNINFORMATIVE (R_j<0.5)"
    assert abs(j["R_j"]) < 1e-9 and j["Q90"] is None, "uninformative R_j/Q90"

    inform = {"safe": cv([10, 10.5, 11, 11.5, 12]),
              "fast": cv([11, 11.9, 12.5, 13, 13.5])}
    j = judge_image(inform, [0, 1], Tn, 1.0)
    assert abs(j["R_j"] - 2.0) < 1e-9, "R_j = Qiso(Tmax)-Qiso(Tmin)"
    assert abs(j["Q90"] - (10.0 + 0.9 * 2.0)) < 1e-9, "Q90 = Qmin + 0.9 R_j"
    assert j["status"] in ("NORMAL", "FAST_LEFT_CENSORED") and j["S_gate"] > 0, \
        "informative image classifies positive"
    # raw + isotonic curves are both retained
    assert j["safe_raw"] is not None and j["safe_iso"] is not None, "curves kept"
    n_ok += 1

    # --- fixed-budget gain (raw seed-mean, fast - safe) -------------------
    dq = fixed_budget_gain({"safe": cv([0, 0, 0, 0, 10]),
                            "fast": cv([0, 0, 0, 0, 14])}, [0, 1], 100.0)
    assert abs(dq - 4.0) < 1e-9, "DeltaQ = raw fast-safe at nu_budget"
    assert fixed_budget_gain({"safe": {}, "fast": cv([0]*5)}, [0], 100.0) is None, \
        "DeltaQ None when a grid point is missing"
    n_ok += 1

    # --- family_of parsing + strata -------------------------------------
    assert family_of("detail_glyph_3") == "glyph", "family_of simple"
    assert family_of("detail_line_pair_1") == "line_pair", "family_of multiword"
    strata_d = _build_strata(["detail_glyph_0", "detail_glyph_1",
                              "detail_maze_0"], "detail")
    assert set(strata_d) == {"glyph", "maze"}, "detail strata by family"
    strata_dev = _build_strata(["dev_text", "dev_fine_lines"], "dev")
    assert set(strata_dev) == {"__all__"}, "non-detail = one plain stratum"
    n_ok += 1

    # --- bootstrap machinery runs end-to-end + is deterministic ----------
    imgs = {"a": inform, "b": {"safe": cv([9, 9.6, 10.2, 10.8, 11.4]),
                               "fast": cv([10, 11, 11.8, 12.4, 13])}}
    names = list(imgs)
    sbi = {im: [0, 1] for im in names}
    fns = {"S_gate": lambda im, sd: judge_image(imgs[im], sd, Tn, 1.0)["S_gate"],
           "DeltaQ": lambda im, sd: fixed_budget_gain(imgs[im], sd, 100.0)}
    b1 = bootstrap_endpoints(names, _build_strata(names, "dev"), sbi, fns, 64)
    b2 = bootstrap_endpoints(names, _build_strata(names, "dev"), sbi, fns, 64)
    assert b1["S_gate"]["meds"].shape == (64,), "bootstrap draw count"
    assert np.isfinite(b1["S_gate"]["LB2.5"]), "bootstrap S_gate LB finite"
    assert np.allclose(b1["S_gate"]["meds"], b2["S_gate"]["meds"]), \
        "bootstrap deterministic (fixed rng_for stream)"
    n_ok += 1

    print("[s1] SELFTEST-ANALYSIS PASS (%d groups; PAVA + 7 censoring classes "
          "+ Q90 target + fixed-budget + family strata + bootstrap)" % n_ok,
          flush=True)
    return 0


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
    if getattr(a, "selftest_analysis", False):
        return _selftest()
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
