"""ROUND63 S1 DEV PILOT (local, development-only) — spec D2 §7.

S1 FREEZE CLAUSE (spec D2 §1, verbatim):
  "S1 is exploratory and development-only. S1 images and seeds are disjoint
   from all S2 confirmatory images and seeds. S1 results shall not enter
   confirmatory confidence intervals or main quantitative tables."

Accordingly this script runs ONLY on the S1 development image set
(images63.build_dev_image_set: STL-10 TRAIN split + dev_* structural targets),
which is disjoint from the S2 confirmatory set (STL-10 TEST split). It NEVER
touches a confirmatory image or seed, and its outputs are for protocol
development (signal verification + shard-budget calibration) only.

Two jobs (spec D2 §7):
  (a) verify the SPEED signal at the pre-registered high-flux point rho_bar=0.6
      against the rho_bar=0.05 safe reference, on DEV images (time-to-quality);
  (b) calibrate per-cell wall time (incl. the lam_TV selection overhead) to size
      the Colab shards for S2.

Geometry (spec D2 §3 S2-A, shrunk to the dev set): side=64, pattern="bern50",
FIXED M=2048 (M/n=0.5); dwell scan nu in {20,50,100,200,500,1000,2000}; two
operating points rho_bar in {0.05, 0.6}; tau=50 ns; sigma_b=0; dark=0; active
start. Arms: RQL, POISSON-LIN, SAT-POISSON, PRECORRECT (iterative, discrepancy-
selected) + GI (display reference; non-iterative). Dev image subset
dev_stl_00..03 + dev_text + dev_fine_lines; seeds {0,1}.

Every cell is run through campaign.run_cell — the SAME code path as S2 — so the
pilot also exercises the production discrepancy lam_TV selection (spec D2 §4)
and its per-cell timing is directly transferable to the shard budget.

Run:
  D:/Anacondar/anaconda3/python.exe code/round63/pilot_s1.py --smoke   # fast
  D:/Anacondar/anaconda3/python.exe code/round63/pilot_s1.py           # full dev
  D:/Anacondar/anaconda3/python.exe code/round63/pilot_s1.py --analyze-only
"""
import argparse
import csv
import json
import os
import sys
import time
import traceback
from types import SimpleNamespace

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
from campaign import run_cell

ROOT = os.path.dirname(os.path.dirname(HERE))

# ---- frozen S1 dev configuration (spec D2 §3/§7) --------------------------- #
TAU = 50e-9
PATTERN = "bern50"
DEFAULT_IMAGES = ["dev_stl_00", "dev_stl_01", "dev_stl_02", "dev_stl_03",
                  "dev_text", "dev_fine_lines"]
DEFAULT_ARMS = ["RQL", "POISSON-LIN", "SAT-POISSON", "PRECORRECT", "GI"]
DEFAULT_NU = [20.0, 50.0, 100.0, 200.0, 500.0, 1000.0, 2000.0]
DEFAULT_RHO = [0.05, 0.6]
DEFAULT_SEEDS = [0, 1]

PRIMARY_ARM = "RQL"                                    # spec D2 §0 primary S_j
SECONDARY_TTQ = ["POISSON-LIN", "SAT-POISSON", "PRECORRECT"]   # same Q*_j target

# stable CSV schema = exactly what campaign.run_cell emits (spec D2 §4 columns)
FIELDNAMES = ["side", "pattern", "rho_bar", "nu", "M", "seed", "image", "arm",
              "PSNR", "SSIM", "LPIPS", "rad_nrmse", "flux_dev", "lam_tv",
              "mean_counts", "optical_time_s", "dark_frac", "tau_err",
              "runtime_s", "select_runtime_s", "MODEL_FAIL", "eta_star",
              "PSNR_rad", "gof_status", "gof_p", "leak_suspect"]

_KEY_COLS = ("side", "pattern", "rho_bar", "nu", "M", "seed", "image", "arm")


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
        x = float(v)
    except (TypeError, ValueError):
        return None
    return x


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


# ---- configuration --------------------------------------------------------- #
def _floats(s):
    return [float(x) for x in s.replace(" ", "").split(",") if x != ""]


def _ints(s):
    return [int(x) for x in s.replace(" ", "").split(",") if x != ""]


def _strs(s):
    return [x.strip() for x in s.split(",") if x.strip()]


def parse_args(argv):
    ap = argparse.ArgumentParser(
        description="ROUND63 S1 dev pilot (development-only; spec D2 §7).")
    ap.add_argument("--nu-list", type=str, default=None,
                    help="comma dwell list (default 20,50,100,200,500,1000,2000)")
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
                    "(default 2048 = M/n 0.5 at side 64)")
    ap.add_argument("--select-iter", type=int, default=60,
                    help="fista iters inside lam_TV selection (default 60)")
    ap.add_argument("--fista-iters", type=int, default=200,
                    help="final-refit fista iters (default 200)")
    ap.add_argument("--gof-mode", type=str, default="refit",
                    help="select_eta gof mode (default refit)")
    ap.add_argument("--smoke", action="store_true",
                    help="fast preset: side=32,M=512,2 images,seeds{0},"
                         "nu{100,500},arms RQL+POISSON-LIN+GI (own out dir)")
    ap.add_argument("--out", type=str, default=None,
                    help="output dir (default results/round63_s1)")
    ap.add_argument("--analyze-only", action="store_true",
                    help="skip the sweep, (re)build summaries from existing CSV")
    ap.add_argument("--no-analysis", action="store_true",
                    help="sweep only (parallel worker mode: a single-cell "
                         "worker CSV cannot support the S_j analysis, which "
                         "needs the full grid — run_s1_parallel.py merges the "
                         "workers and analyzes once)")
    return ap.parse_args(argv)


def build_cfg(a):
    smoke = bool(a.smoke)
    side = a.side if a.side is not None else (32 if smoke else 64)
    M = a.M if a.M is not None else (512 if smoke else 2048)
    n_img = a.images if a.images is not None else (2 if smoke else len(DEFAULT_IMAGES))
    images = DEFAULT_IMAGES[:max(1, n_img)]
    seeds = _ints(a.seeds) if a.seeds else ([0] if smoke else list(DEFAULT_SEEDS))
    nu_list = _floats(a.nu_list) if a.nu_list else ([100.0, 500.0] if smoke
                                                    else list(DEFAULT_NU))
    rho_list = _floats(a.rho_list) if a.rho_list else list(DEFAULT_RHO)
    arms = _strs(a.arms) if a.arms else (["RQL", "POISSON-LIN", "GI"] if smoke
                                         else list(DEFAULT_ARMS))
    if a.out:
        out = a.out
    else:
        out = os.path.join(ROOT, "results",
                           "round63_s1_smoke" if smoke else "round63_s1")
    return SimpleNamespace(
        side=side, M=M, pattern=PATTERN, tau=TAU, images=images, seeds=seeds,
        nu_list=nu_list, rho_list=rho_list, arms=arms,
        select_iter=a.select_iter, fista_iters=a.fista_iters,
        gof_mode=a.gof_mode, out=out, smoke=smoke)


def cfg_dict(cfg):
    return {"side": cfg.side, "M": cfg.M, "pattern": cfg.pattern, "tau": cfg.tau,
            "nu_list": cfg.nu_list, "rho_list": cfg.rho_list, "arms": cfg.arms,
            "images": cfg.images, "seeds": cfg.seeds,
            "select_iter": cfg.select_iter, "fista_iters": cfg.fista_iters,
            "gof_mode": cfg.gof_mode, "smoke": cfg.smoke}


def cells(cfg):
    """The cell list: one cell per (rho_bar, nu, seed); run_cell loops the image
    and arm axes internally. dev=True routes run_cell to the dev image set."""
    cs = []
    for rho in cfg.rho_list:
        for nu in cfg.nu_list:
            for seed in cfg.seeds:
                cs.append(dict(
                    side=cfg.side, pattern=cfg.pattern, rho_bar=float(rho),
                    nu=float(nu), M=cfg.M, seed=int(seed), arms=list(cfg.arms),
                    images=list(cfg.images), dev=True, tau=cfg.tau, sigma_b=0.0,
                    fista_iters=cfg.fista_iters, select_iter=cfg.select_iter,
                    gof_mode=cfg.gof_mode, select_rule="discrepancy", sel_K=5,
                    use_lpips=False))
    return cs


# ---- sweep (incremental, resumable) --------------------------------------- #
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


# ---- analysis 1: time-to-quality ------------------------------------------ #
def _seed_mean_at(rows, arm, rho, image, nu):
    ps = []
    for r in rows:
        if r["arm"] != arm or r["image"] != image:
            continue
        if abs(_f(r["rho_bar"]) - rho) > 1e-12 or abs(_f(r["nu"]) - nu) > 1e-9:
            continue
        p = _psnr_rad(r)
        if p is not None:
            ps.append(p)
    return float(np.mean(ps)) if ps else None


def _seed_mean_curve(rows, arm, rho, image):
    """(T_opt asc, P_raw, P_monotone, nu asc) seed-mean PSNR_rad vs T_opt, or
    None. Monotone-non-decreasing enforced via cumulative max (isotonic-lite)."""
    by_nu = {}
    for r in rows:
        if r["arm"] != arm or r["image"] != image:
            continue
        if abs(_f(r["rho_bar"]) - rho) > 1e-12:
            continue
        p = _psnr_rad(r)
        if p is None:
            continue
        d = by_nu.setdefault(_f(r["nu"]), {"p": [], "T": _f(r["optical_time_s"])})
        d["p"].append(p)
    if not by_nu:
        return None
    nus = sorted(by_nu)
    T = np.array([by_nu[nu]["T"] for nu in nus], dtype=np.float64)
    P = np.array([float(np.mean(by_nu[nu]["p"])) for nu in nus], dtype=np.float64)
    order = np.argsort(T)
    T, P = T[order], P[order]
    return T, P, np.maximum.accumulate(P), np.array(nus, dtype=np.float64)[order]


def _t_reach(T, P, Q):
    """First T reaching quality Q by linear interpolation in log10(T) on the
    monotonized curve P. Returns (T_reach, censored)."""
    ge = np.nonzero(P >= Q)[0]
    if ge.size == 0:
        return float("nan"), True                         # never reaches -> censored
    k = int(ge[0])
    if k == 0:
        return float(T[0]), False                         # already met at shortest T
    logT = np.log10(T)
    if P[k] == P[k - 1]:
        lt = logT[k]
    else:
        lt = logT[k - 1] + (Q - P[k - 1]) * (logT[k] - logT[k - 1]) \
            / (P[k] - P[k - 1])
    return float(10.0 ** lt), False


def time_to_quality(rows, cfg):
    images = [im for im in cfg.images if any(r["image"] == im for r in rows)]
    nus_all = sorted({_f(r["nu"]) for r in rows if np.isfinite(_f(r["nu"]))})
    nu_ref = max(nus_all) if nus_all else float("nan")    # spec: RQL@0.05@nu_max
    rho_safe = min(cfg.rho_list)
    rho_fast = max(cfg.rho_list)

    Qstar = {}
    for im in images:
        q = _seed_mean_at(rows, PRIMARY_ARM, rho_safe, im, nu_ref)
        Qstar[im] = (q - 1.0) if q is not None else None

    arms = [PRIMARY_ARM] + [a for a in SECONDARY_TTQ if a in cfg.arms]
    S, censored, T_reach = {}, {}, {}
    for arm in arms:
        S[arm], censored[arm], T_reach[arm] = {}, {}, {}
        for im in images:
            Q = Qstar.get(im)
            c_safe = _seed_mean_curve(rows, arm, rho_safe, im)
            c_fast = _seed_mean_curve(rows, arm, rho_fast, im)
            if Q is None or c_safe is None or c_fast is None:
                S[arm][im], censored[arm][im] = None, True
                T_reach[arm][im] = {"T_safe": None, "T_fast": None}
                continue
            Ts, cs = _t_reach(c_safe[0], c_safe[2], Q)
            Tf, cf = _t_reach(c_fast[0], c_fast[2], Q)
            cen = bool(cs or cf)
            sj = (Ts / Tf) if (not cen and Tf > 0 and np.isfinite(Ts)
                               and np.isfinite(Tf)) else None
            S[arm][im] = sj
            censored[arm][im] = cen
            T_reach[arm][im] = {"T_safe": (Ts if np.isfinite(Ts) else None),
                                "T_fast": (Tf if np.isfinite(Tf) else None)}

    medians = {}
    for arm in arms:
        vals = [v for v in S[arm].values() if v is not None]
        medians[arm] = float(np.median(vals)) if vals else None
    cen_counts = {arm: int(sum(1 for v in censored[arm].values() if v))
                  for arm in arms}
    return {"nu_ref": nu_ref, "rho_safe": rho_safe, "rho_fast": rho_fast,
            "Qstar_j": Qstar, "S_j": S, "S_median": medians,
            "censored": censored, "censored_counts": cen_counts,
            "T_reach": T_reach, "images": images, "arms": arms}


# ---- analysis 2: MODEL_FAIL incidence ------------------------------------- #
def model_fail_table(rows):
    agg = {}
    for r in rows:
        mf = r.get("MODEL_FAIL", "")
        if mf in ("", None):
            continue                                      # non-selected arm
        key = (_f(r["rho_bar"]), _f(r["nu"]), r["arm"])
        d = agg.setdefault(key, {"fail": 0, "tot": 0})
        d["tot"] += 1
        if str(mf) == "True":
            d["fail"] += 1
    return {"%g|%g|%s" % (k[0], k[1], k[2]):
            {"fail": v["fail"], "tot": v["tot"],
             "rate": (v["fail"] / v["tot"] if v["tot"] else 0.0)}
            for k, v in sorted(agg.items())}


# ---- analysis 3: runtime calibration -------------------------------------- #
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
def _fmt_S(s, cen):
    if s is not None:
        return "%.3f" % s
    return "cens" if cen else "--"


def _write_md(cfg, ttq, mft, rtc, n_rows, path):
    L = []
    L.append("# ROUND63 S1 DEV PILOT SUMMARY")
    L.append("")
    L.append("**S1 is exploratory and development-only** (spec D2 §1): dev images"
             " and seeds are disjoint from all S2 confirmatory images and seeds;"
             " these results do not enter confirmatory intervals or main tables.")
    L.append("")
    L.append("- rows: %d   images: %d   seeds: %s   arms: %s"
             % (n_rows, len(cfg.images), cfg.seeds, ", ".join(cfg.arms)))
    L.append("- side %d, %s, M %d, nu %s, rho_bar %s, tau %.0f ns, sigma_b 0, "
             "active start" % (cfg.side, cfg.pattern, cfg.M,
                               [int(x) for x in cfg.nu_list], cfg.rho_list,
                               cfg.tau * 1e9))
    L.append("- select_iter %d, fista_iters %d, gof_mode %s, discrepancy "
             "lam_TV rule (spec D2 §4)"
             % (cfg.select_iter, cfg.fista_iters, cfg.gof_mode))
    L.append("")

    # --- Table 1: time-to-quality speedup
    L.append("## 1. Time-to-quality speedup S_j (dev)")
    L.append("")
    L.append("S_j = T_opt(rho=%.3g, Q*_j) / T_opt(rho=%.3g, Q*_j), "
             "T_opt = M_physical*T. Q*_j = PSNR_rad[RQL, rho=%.3g, nu=%g] - 1 dB "
             "(same target for every arm). 'cens' = target never reached; "
             "'--' = missing."
             % (ttq["rho_safe"], ttq["rho_fast"], ttq["rho_safe"], ttq["nu_ref"]))
    L.append("")
    arms = ttq["arms"]
    hdr = "| image | Q*_j (dB) | " + " | ".join("%s S_j" % a for a in arms) + " |"
    sep = "|" + "---|" * (2 + len(arms))
    L.append(hdr)
    L.append(sep)
    for im in ttq["images"]:
        q = ttq["Qstar_j"].get(im)
        qs = "%.2f" % q if q is not None else "--"
        cells_ = " | ".join(_fmt_S(ttq["S_j"][a].get(im),
                                   ttq["censored"][a].get(im, False))
                            for a in arms)
        L.append("| %s | %s | %s |" % (im, qs, cells_))
    med = " | ".join(("%.3f" % ttq["S_median"][a]) if ttq["S_median"][a]
                     is not None else "--" for a in arms)
    L.append("| **median** | | %s |" % med)
    L.append("")
    cc = ", ".join("%s=%d" % (a, ttq["censored_counts"][a]) for a in arms)
    L.append("Censored (image count per arm): %s" % cc)
    L.append("")

    # --- Table 2: MODEL_FAIL incidence
    L.append("## 2. MODEL_FAIL incidence by (rho, nu, arm)")
    L.append("")
    if mft:
        L.append("| rho | nu | arm | fail | total | rate |")
        L.append("|---|---|---|---|---|---|")
        for k, v in mft.items():
            rho, nu, arm = k.split("|")
            L.append("| %s | %s | %s | %d | %d | %.3f |"
                     % (rho, nu, arm, v["fail"], v["tot"], v["rate"]))
    else:
        L.append("(no discrepancy-selected rows — MODEL_FAIL not applicable)")
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

    ttq = time_to_quality(rows, cfg)
    mft = model_fail_table(rows)
    rtc = runtime_calibration(rows)
    summary = {"config": cfg_dict(cfg), "n_rows": len(rows),
               "time_to_quality": ttq, "model_fail_rates": mft,
               "runtime_calibration": rtc}

    with open(os.path.join(cfg.out, "pilot_summary.json"), "w") as f:
        json.dump(summary, f, indent=2, default=_json_default)
    with open(os.path.join(cfg.out, "runtime_calibration.json"), "w") as f:
        json.dump(rtc, f, indent=2, default=_json_default)
    md = os.path.join(cfg.out, "PILOT_SUMMARY.md")
    _write_md(cfg, ttq, mft, rtc, len(rows), md)

    # stdout digest
    print("[s1] ANALYSIS  rows=%d  nu_ref=%g  rho %g->%g"
          % (len(rows), ttq["nu_ref"], ttq["rho_safe"], ttq["rho_fast"]),
          flush=True)
    for a in ttq["arms"]:
        m = ttq["S_median"][a]
        print("       S_median[%-11s] = %s   censored=%d/%d"
              % (a, ("%.3f" % m) if m is not None else "n/a",
                 ttq["censored_counts"][a], len(ttq["images"])), flush=True)
    n_fail = sum(v["fail"] for v in mft.values())
    n_tot = sum(v["tot"] for v in mft.values())
    print("       MODEL_FAIL total = %d/%d selected rows" % (n_fail, n_tot),
          flush=True)
    if rtc:
        wmax = max(v["wall_max_s"] for v in rtc.values())
        print("       max per-fit wall = %.2fs (see runtime_calibration.json)"
              % wmax, flush=True)
    print("[s1] wrote %s , pilot_summary.json , runtime_calibration.json"
          % md, flush=True)
    return summary


# ---- driver ---------------------------------------------------------------- #
def main(argv=None):
    a = parse_args(sys.argv[1:] if argv is None else argv)
    cfg = build_cfg(a)
    os.makedirs(cfg.out, exist_ok=True)
    print("[s1] out=%s smoke=%s cells=%d (%d rho x %d nu x %d seed)"
          % (cfg.out, cfg.smoke, len(cfg.rho_list) * len(cfg.nu_list)
             * len(cfg.seeds), len(cfg.rho_list), len(cfg.nu_list),
             len(cfg.seeds)), flush=True)
    if not a.analyze_only:
        sweep(cfg)
    if not a.no_analysis:
        analyze(cfg)
    print("[s1] DONE", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
