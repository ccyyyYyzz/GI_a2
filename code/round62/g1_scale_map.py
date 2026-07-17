"""ROUND62 Part 2 (T1-G1): head-space scale map — spec §2. Runs only after G0 PASS.

Maps the semiparametric head-space (rMAVE / OPG single-index oracle class) over
resolution n and frame budget M, against the honest whitening/isotron baselines,
to test whether the Phase-A MAVE advantage survives scaling.

Grid (spec §2): n in {256, 576, 1024} (sides 16/24/32; n=4096 is an OPTIONAL
extension that MUST NOT run in this harness — see --with-4096 guard) x
M in {20000, 100000} x illum {GAM4, GAM8, CORR-LOGN} x link {DT30, SAT30}
x 3 seeds x 8 images. Photons s = 1e4 throughout (the G0 / Phase-A A3 operating
point; fixed here so the map isolates the n x M axes).

Preregistered oracle class (anti-garden-of-forking-paths):
  members = {OPG, rMAVE} x bandwidth-scale {0.5, 1.0, 2.0} = 6 configs,
  ridge = 1e-4, rMAVE iterations capped at 20. rmave_single_index(return_opg=True)
  yields BOTH members from one call per bandwidth. Anchor/window params vary with
  n (logged to the output JSON); n=1024 uses a subsampled OPG (anchor 300 <= 5000).

Head-space discipline (no truth anywhere in fitting or selection):
  - each member direction is FIT on the train split (first 90% of frames) only;
  - member + bandwidth selection uses HELD-OUT FRAME MSE only: isotonic (PAV,
    increasing) link fit on the train frames, evaluated on the held-out last 10%;
  - the selected member is recorded as SEMI-ORACLE (every config is also recorded
    individually so the map is auditable).

Baselines (same protocols as Phase B, on the FULL M frames): WHITEN-LW
(LedoitWolf + Cholesky solve of the centered correlation), L-ISOTRON (frozen
protocol reimplemented locally for arbitrary (n, M): eta = 1/tr(Sigma_LW), PAV
link, <=200 iters, tol 1e-6 on 10% held-out MSE, 3 inits {WHITEN-LW, SIR-10,
random} chosen by held-out MSE), WHITEN-OR (analytic family covariance).

Metrics: MAIN protocol only (flux-matched PSNR is the gate metric; SSIM recorded).
LPIPS is skipped entirely (NaN) — AlexNet-LPIPS is undefined below 32px and is
not needed for a PSNR gate (noted in g1_gates.json).

Patterns are generated ONCE per (family, seed, n) at M=100000 and the M=20000
cell reuses the first 20000 rows (nested-prefix frozen reuse). Bucket noise is
drawn per (seed, fam, n, link, M, image). float64 throughout; NO clamping of
patterns; NO score / density arms (closed-basin rule, ROUND59 terminal).

Checkpointing mirrors phase_b.py: META_JSON is the source of truth for cell
completion (written only AFTER a cell's rows are appended), and sanitize drops
orphan rows from a crash mid-append so a restart never double-appends or biases
the pivot means. Budget guard: 10 h wall, then abort with report.
"""
import csv
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from gi_core import config as C
from gi_core import links as L
from gi_core import metrics as MET
from gi_core.families import make_family
from gi_core.images import load_truths
from gi_core.mave import rmave_single_index
from gi_core.utils import CholOp, pav_fit_predict, rng_for

ROOT = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(ROOT, "results", "round62_g1")

R62 = 62
PART = 2                       # ROUND62 Part 2 stream tag (matches g0's part-1 tag)
S = 1e4                        # photons — the G0 / Phase-A A3 operating point
NS = [256, 576, 1024]          # sides 16 / 24 / 32
SIDES = {256: 16, 576: 24, 1024: 32}
MS = [20000, 100000]
FAMS = ["GAM4", "GAM8", "CORR-LOGN"]
LINKS = ["DT30", "SAT30"]
SEEDS = C.SEEDS_EXT            # [0, 1, 2]
BWS = [0.5, 1.0, 2.0]
RMAVE_RIDGE = 1e-4
RMAVE_NITER = 20
N_JOBS = 8
BUDGET_WALL_S = 10 * 3600

# per-n anchor / window params (opg_n_anchor <= 5000 spec bound); logged to JSON
ANCHOR_PARAMS = {
    256:  dict(n_anchor=150, n_loc=400, opg_n_anchor=150, opg_n_loc=1500),
    576:  dict(n_anchor=200, n_loc=400, opg_n_anchor=200, opg_n_loc=2000),
    1024: dict(n_anchor=250, n_loc=400, opg_n_anchor=300, opg_n_loc=2000),
}

CSV = os.path.join(OUT, "g1_scale_map.csv")
META_JSON = os.path.join(OUT, "g1_run_meta.json")
GATES_JSON = os.path.join(OUT, "g1_gates.json")
HEADER = ["n", "M", "illum", "link", "seed", "image", "method",
          "PSNR", "SSIM", "chosen_member"]

# ------- gate constants (spec §2 G1 gate, preregistered) -------
GATE_N = 1024
GATE_M = 100000
GATE_LINK = "DT30"
GATE_ILLUMS = ["GAM4", "GAM8"]
GATE_DB = 1.0
GATE_MIN_POS = 6


# ======================================================================
# checkpoint / sanitize (META-as-truth, mirrors phase_b.py)
# ======================================================================
def cell_key(n, fam, seed, M, link):
    return "%d|%s|%d|%d|%s" % (int(n), fam, int(seed), int(M), link)


def sanitize_csv(path, meta):
    """Drop rows whose cell is not marked complete in META (crash mid-append),
    so restart recomputes them cleanly and pivots are never biased."""
    if not os.path.exists(path):
        return
    complete = set(meta.get("cells", {}))
    kept, dropped = [], 0
    with open(path) as f:
        rd = csv.reader(f)
        header = next(rd)
        for r in rd:
            # row = [n, M, illum, link, seed, image, method, PSNR, SSIM, chosen]
            key = cell_key(int(r[0]), r[2], int(r[4]), int(r[1]), r[3])
            if key in complete:
                kept.append(r)
            else:
                dropped += 1
    if dropped:
        print("[sanitize] %s: dropped %d orphan rows" % (path, dropped), flush=True)
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(header)
            w.writerows(kept)


def append_rows(path, rows):
    new = not os.path.exists(path)
    with open(path, "a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(HEADER)
        w.writerows(rows)
        f.flush()
        os.fsync(f.fileno())


# ======================================================================
# semiparametric oracle members (fit on the train split only)
# ======================================================================
def _rmave_one(A, b, rng, bw, ap):
    """One rMAVE fit; returns (beta_rmave, beta_opg), both unit-norm directions."""
    return rmave_single_index(
        A, b, rng, n_anchor=ap["n_anchor"], n_loc=ap["n_loc"],
        n_iter=RMAVE_NITER, ridge=RMAVE_RIDGE, bw_scale=bw,
        opg_n_anchor=ap["opg_n_anchor"], opg_n_loc=ap["opg_n_loc"],
        return_opg=True)


def mave_pair_batch(A, B, rngs, bw, ap, n_jobs=None):
    """rMAVE + OPG directions for all K image columns (parallel over images)."""
    if n_jobs is None:
        n_jobs = N_JOBS
    try:
        from joblib import Parallel, delayed
    except ImportError:
        Parallel = None
    if Parallel is not None and n_jobs != 1:
        res = Parallel(n_jobs=n_jobs)(
            delayed(_rmave_one)(A, B[:, k], rngs[k], bw, ap)
            for k in range(B.shape[1]))
    else:
        res = [_rmave_one(A, B[:, k], rngs[k], bw, ap) for k in range(B.shape[1])]
    rmave = np.stack([r[0] for r in res], axis=1)
    opg = np.stack([r[1] for r in res], axis=1)
    return rmave, opg


# ======================================================================
# baselines (WHITEN-LW / SIR-10 / L-ISOTRON reimplemented for arbitrary (n, M))
# ======================================================================
def _center_cols(B):
    return B - B.mean(axis=0, keepdims=True)


def whiten_lw_dir(A, B, chol_lw, M):
    return chol_lw.solve(A.T @ _center_cols(B) / M)


def _slice_means(A, b, H):
    """Equal-count slice means of A ordered by b (copy of estimators._slice_means)."""
    M = b.shape[0]
    order = np.argsort(b, kind="stable")
    bounds = np.linspace(0, M, H + 1).astype(int)
    Msl = np.zeros((A.shape[1], H))
    w = np.zeros(H)
    for h in range(H):
        idx = order[bounds[h]:bounds[h + 1]]
        Msl[:, h] = A[idx].mean(axis=0)
        w[h] = idx.size / M
    return Msl, w


def _sign_fix(A, x, b):
    u = A @ x
    c = np.dot(u - u.mean(), b - b.mean())
    return x if c >= 0 else -x


def sir_local(A, B, H, chol_lw, abar):
    """SIR with LW-standardized patterns (copy of estimators.sir, no RunContext)."""
    n, K = A.shape[1], B.shape[1]
    out = np.zeros((n, K))
    for k in range(K):
        Msl, w = _slice_means(A, B[:, k], H)
        Bmat = (Msl - abar[:, None]) * np.sqrt(w)[None, :]
        Cw = chol_lw.half_solve(Bmat)          # L^{-1} B
        U, sv, _ = np.linalg.svd(Cw, full_matrices=False)
        beta = chol_lw.half_solve_T(U[:, 0])   # L^{-T} eta
        out[:, k] = _sign_fix(A, beta, B[:, k])
    return out


def l_isotron_local(A, B, chol_lw, tr_lw, abar, tr_idx, val_idx, rng,
                    max_iter=200, tol=1e-6):
    """L-Isotron frozen protocol (spec §4.8) reimplemented for arbitrary (n, M).

    Faithful port of estimators.l_isotron: eta = 1/tr(Sigma_LW), PAV link on the
    train frames, 3 inits {WHITEN-LW, SIR-10, random} scaled to the LW norm,
    best-of chosen by held-out MSE. The only difference from the RunContext
    version is the held-out split: here it is the deterministic last-10% split
    (train = first 90%) shared with the semiparametric selection, rather than the
    original random-permutation split. Frames are i.i.d. across rows, so the
    prefix/suffix split is exchangeable with a random one; truth never used."""
    A_tr, A_val = A[tr_idx], A[val_idx]
    B_tr, B_val = B[tr_idx], B[val_idx]
    eta = 1.0 / tr_lw
    n, K = A.shape[1], B.shape[1]
    M = A.shape[0]

    x_lw = whiten_lw_dir(A, B, chol_lw, M)
    x_sir = sir_local(A, B, 10, chol_lw, abar)
    x_rand = rng.standard_normal((n, K))
    x_rand /= np.linalg.norm(x_rand, axis=0, keepdims=True)
    # scale-match non-LW inits to the LW norm so eta is comparable
    for X0 in (x_sir, x_rand):
        X0 *= np.linalg.norm(x_lw, axis=0, keepdims=True) / np.maximum(
            np.linalg.norm(X0, axis=0, keepdims=True), 1e-30)

    best_val = np.full(K, np.inf)
    best_X = np.zeros((n, K))
    Mtr = A_tr.shape[0]
    for X0 in (x_lw, x_sir, x_rand):
        X = X0.copy()
        prev_val = np.full(K, np.inf)
        active = np.ones(K, dtype=bool)
        run_val = np.full(K, np.inf)
        run_best_X = X.copy()
        run_best_val = np.full(K, np.inf)
        for it in range(max_iter):
            if not active.any():
                break
            U_tr = A_tr @ X
            U_val = A_val @ X
            R = np.zeros((Mtr, K))
            for k in np.nonzero(active)[0]:
                phi_tr = pav_fit_predict(U_tr[:, k], B_tr[:, k], U_tr[:, k])
                phi_val = pav_fit_predict(U_tr[:, k], B_tr[:, k], U_val[:, k])
                R[:, k] = B_tr[:, k] - phi_tr
                run_val[k] = float(np.mean((B_val[:, k] - phi_val) ** 2))
                if run_val[k] < run_best_val[k]:
                    run_best_val[k] = run_val[k]
                    run_best_X[:, k] = X[:, k]
                if abs(prev_val[k] - run_val[k]) <= tol * max(1.0, run_val[k]):
                    active[k] = False
                prev_val[k] = run_val[k]
            X = X + eta * (A_tr.T @ R) / Mtr
        for k in range(K):
            if run_best_val[k] < best_val[k]:
                best_val[k] = run_best_val[k]
                best_X[:, k] = run_best_X[:, k]
    return best_X


# ======================================================================
# one cell = (n, fam, seed, M, link)
# ======================================================================
def process_cell(A, B, X, side, n, fam, fam_id, seed, M, link, link_id,
                 chol_lw, tr_lw, abar, cov_op):
    ap = ANCHOR_PARAMS[n]
    n_val = M // 10
    tr_idx = slice(0, M - n_val)
    val_idx = slice(M - n_val, M)
    A_tr, A_val = A[tr_idx], A[val_idx]
    B_tr, B_val = B[tr_idx], B[val_idx]

    # ---- semiparametric members (fit on the train split ONLY) ----
    member_dirs = {}   # method name -> (n, 8) directions
    for bw in BWS:
        rngs = [rng_for(seed, fam_id, R62, PART, n, link_id, M,
                        int(round(bw * 10)), 100 + k) for k in range(8)]
        rmave, opg = mave_pair_batch(A_tr, B_tr, rngs, bw, ap)
        member_dirs["RMAVE-bw%g" % bw] = rmave
        member_dirs["OPG-bw%g" % bw] = opg
    member_names = list(member_dirs.keys())   # 6 configs

    # ---- member/hyperparameter selection: held-out frame MSE only (no truth) ----
    chosen = []
    semi_dir = np.zeros((A.shape[1], 8))
    for k in range(8):
        best_mse, best_name = np.inf, None
        for name in member_names:
            beta = member_dirs[name][:, k]
            u_tr = A_tr @ beta
            u_val = A_val @ beta
            pred = pav_fit_predict(u_tr, B_tr[:, k], u_val)
            mse = float(np.mean((B_val[:, k] - pred) ** 2))
            if mse < best_mse:
                best_mse, best_name = mse, name
        chosen.append(best_name)
        semi_dir[:, k] = member_dirs[best_name][:, k]

    # ---- baselines (full M frames) ----
    corr = A.T @ _center_cols(B) / M
    x_lw = chol_lw.solve(corr)                       # WHITEN-LW
    x_or = cov_op.solve(corr)                        # WHITEN-OR (analytic)
    iso_rng = rng_for(seed, fam_id, R62, PART, n, link_id, M, 200)
    x_iso = l_isotron_local(A, B, chol_lw, tr_lw, abar, tr_idx, val_idx, iso_rng)

    # ---- collect rows (MAIN protocol; LPIPS skipped -> NaN) ----
    rows = []

    def metric_row(method, direction, k, chosen_member=""):
        m = MET.main_metrics(direction[:, k], X[:, k], side, with_lpips=False)
        return [n, int(M), fam, link, seed, C.IMAGES[k], method,
                "%.6g" % m["PSNR"], "%.6g" % m["SSIM"], chosen_member]

    for k in range(8):
        rows.append(metric_row("SEMI-ORACLE", semi_dir, k, chosen[k]))
    for name in member_names:
        for k in range(8):
            rows.append(metric_row(name, member_dirs[name], k))
    for name, direction in (("WHITEN-LW", x_lw), ("L-ISOTRON", x_iso),
                            ("WHITEN-OR", x_or)):
        for k in range(8):
            rows.append(metric_row(name, direction, k))
    return rows


# ======================================================================
# gates + headroom tables
# ======================================================================
def compute_gates(meta, aborted):
    import pandas as pd

    hyper = {
        "photons_s": S, "photons_note": "fixed at the G0/A3 operating point",
        "oracle_class": {"members": ["OPG", "rMAVE"], "bw_scales": BWS,
                         "rmave_ridge": RMAVE_RIDGE, "rmave_n_iter": RMAVE_NITER,
                         "selection": "held-out last-10% frame MSE (isotonic PAV "
                                      "link), no truth; direction fit on first-90% "
                                      "train split"},
        "anchor_params_per_n": ANCHOR_PARAMS,
        "baselines": ["WHITEN-LW", "L-ISOTRON", "WHITEN-OR"],
        "l_isotron_split": "deterministic last-10% held-out (shared with member "
                           "selection); otherwise the frozen estimators.l_isotron "
                           "protocol (eta=1/tr(Sigma_LW), PAV, <=200 iters, tol "
                           "1e-6, 3 inits best-of)",
        "lpips": "skipped entirely (NaN): AlexNet-LPIPS undefined <32px and not "
                 "needed for the PSNR gate",
        "grid": {"n": NS, "M": MS, "illum": FAMS, "link": LINKS,
                 "seeds": list(SEEDS)},
        "n4096": "OPTIONAL extension, intentionally NOT run in this harness",
    }
    gates = {"hyperparams": hyper,
             "operationalization":
                 "gate cell (n=1024, M=100000); for {GAM4xDT30, GAM8xDT30}: "
                 "per image, delta = mean_seeds(SEMI-ORACLE) - max(mean_seeds("
                 "WHITEN-LW), mean_seeds(L-ISOTRON)); combo PASS iff mean(delta) "
                 ">= +1.0 dB AND >= 6/8 images with delta > 0; G1 PASS iff either "
                 "combo passes (mean-then-max, matching g0)."}

    if not os.path.exists(CSV):
        gates["verdict"] = "NOT_RUN"
        gates["reason"] = "no results CSV"
        return gates

    df = pd.read_csv(CSV)
    if df.empty:
        gates["verdict"] = "NOT_RUN"
        gates["reason"] = "empty results CSV"
        return gates
    df["PSNR"] = df["PSNR"].astype(float)
    df["SSIM"] = df["SSIM"].astype(float)
    for c in ("n", "M", "seed"):
        df[c] = df[c].astype(int)

    def img_seed_mean(sub, method):
        """3-seed mean per image (reindexed to C.IMAGES) or None if absent."""
        s = sub[sub.method == method]
        if s.empty:
            return None
        return s.groupby("image")["PSNR"].mean().reindex(C.IMAGES)

    # ---- headroom over every completed (n, M, illum, link) ----
    headroom_table = []
    for (n, M, illum, link), sub in df.groupby(["n", "M", "illum", "link"]):
        semi = img_seed_mean(sub, "SEMI-ORACLE")
        lw = img_seed_mean(sub, "WHITEN-LW")
        iso = img_seed_mean(sub, "L-ISOTRON")
        if semi is None or lw is None or iso is None:
            continue
        base = np.maximum(lw.values, iso.values)
        delta = semi.values - base
        headroom_table.append({
            "n": int(n), "M": int(M), "illum": illum, "link": link,
            "mean_headroom_db": float(np.nanmean(delta)),
            "n_pos_images": int(np.nansum(delta > 0)),
            "M_over_n": float(M) / float(n)})

    vs_n = {}
    for row in headroom_table:
        line = "%sx%s" % (row["illum"], row["link"])
        vs_n.setdefault("M%d" % row["M"], {}).setdefault(line, []).append(
            {"n": row["n"], "headroom_db": row["mean_headroom_db"]})
    for mk in vs_n:
        for line in vs_n[mk]:
            vs_n[mk][line] = sorted(vs_n[mk][line], key=lambda d: d["n"])
    vs_M_over_n = sorted(
        [{"M_over_n": r["M_over_n"], "headroom_db": r["mean_headroom_db"],
          "n": r["n"], "M": r["M"], "illum": r["illum"], "link": r["link"]}
         for r in headroom_table], key=lambda d: d["M_over_n"])

    gates["headroom_table"] = headroom_table
    gates["headroom_vs_n"] = vs_n
    gates["headroom_vs_M_over_n"] = vs_M_over_n

    # ---- G1 gate at (n=1024, M=100000) ----
    combos = {}
    passed_combos = []
    any_present = False
    for illum in GATE_ILLUMS:
        sub = df[(df.n == GATE_N) & (df.M == GATE_M) & (df.illum == illum)
                 & (df.link == GATE_LINK)]
        semi = img_seed_mean(sub, "SEMI-ORACLE")
        lw = img_seed_mean(sub, "WHITEN-LW")
        iso = img_seed_mean(sub, "L-ISOTRON")
        tag = "%sx%s" % (illum, GATE_LINK)
        if semi is None or lw is None or iso is None:
            combos[tag] = {"status": "NOT_RUN"}
            continue
        any_present = True
        base = np.maximum(lw.values, iso.values)
        delta = semi.values - base
        mean_delta = float(np.nanmean(delta))
        n_pos = int(np.nansum(delta > 0))
        ok = bool(mean_delta >= GATE_DB and n_pos >= GATE_MIN_POS)
        combos[tag] = {
            "status": "EVALUATED",
            "mean_delta_db": mean_delta,
            "n_pos_images": n_pos,
            "per_image_delta_db": {img: float(d)
                                   for img, d in zip(C.IMAGES, delta)},
            "n_seeds": int(sub[sub.method == "SEMI-ORACLE"]["seed"].nunique()),
            "pass": ok}
        if ok:
            passed_combos.append((tag, mean_delta))

    gates["gate_cell"] = {"n": GATE_N, "M": GATE_M, "link": GATE_LINK}
    gates["combos"] = combos
    gates["aborted_over_budget"] = bool(aborted)

    if not any_present:
        gates["verdict"] = "NOT_RUN"
        gates["reason"] = ("gate cell (n=1024, M=100000, DT30) not completed"
                           + (" (budget abort)" if aborted else ""))
    elif passed_combos:
        best = max(passed_combos, key=lambda t: t[1])
        gates["verdict"] = "PASS(+%.2f dB @ n=1024, %s)" % (best[1], best[0])
        gates["next"] = "T1_PROCEED_TO_G2"
    else:
        # KILL: summarize the decay curve for the ablation section
        decay = {}
        for line, pts in vs_n.get("M%d" % GATE_M, {}).items():
            decay[line] = {str(p["n"]): round(p["headroom_db"], 3) for p in pts}
        gates["verdict"] = "KILL(head-space < +1.0 dB at n=1024; decay-vs-n below)"
        gates["decay_vs_n_at_M100000"] = decay
        gates["next"] = "T1_KILL"
    return gates


# ======================================================================
# main
# ======================================================================
def main():
    if "--with-4096" in sys.argv:
        print("n=4096 is an OPTIONAL extension and is intentionally NOT "
              "implemented in this harness (ROUND62 spec §2: run only within "
              "budget, with subsampled OPG). Refusing to run. Remove --with-4096 "
              "to execute the frozen n in {256, 576, 1024} grid.", flush=True)
        return 2

    n_jobs = N_JOBS
    if "--jobs" in sys.argv:
        n_jobs = int(sys.argv[sys.argv.index("--jobs") + 1])
    globals()["N_JOBS"] = n_jobs

    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    cpu0 = time.process_time()

    meta = {"cells": {}, "aborted_over_budget": False, "wall_accum_s": 0.0,
            "anchor_params": ANCHOR_PARAMS, "n_jobs": n_jobs}
    if os.path.exists(META_JSON):
        with open(META_JSON) as f:
            meta = json.load(f)
        meta.setdefault("cells", {})
        meta.setdefault("wall_accum_s", meta.get("wall_s", 0.0))
        meta["aborted_over_budget"] = False
    sanitize_csv(CSV, meta)
    wall_prev = float(meta["wall_accum_s"])
    done_keys = set(meta["cells"])

    truths_cache = {}
    aborted = False
    for n in NS:
        if aborted:
            break
        side = SIDES[n]
        if side not in truths_cache:
            truths = load_truths(side)
            truths_cache[side] = np.stack([truths[k] for k in C.IMAGES], axis=1)
        X = truths_cache[side]                       # (n, 8) sum-normalized truths
        for fam in FAMS:
            if aborted:
                break
            fam_id = C.FAMILY_IDS[fam]
            for seed in SEEDS:
                if aborted:
                    break
                todo = [(M, link) for M in MS for link in LINKS
                        if cell_key(n, fam, seed, M, link) not in done_keys]
                if not todo:
                    continue
                family = make_family(fam, n)
                A_full = family.sample(100000, rng_for(seed, fam_id, R62, PART, n))
                A_full = np.asarray(A_full, dtype=np.float64)   # no clamping
                u_full = A_full @ X
                cov_op = family.true_cov_op()                   # WHITEN-OR analytic
                for M in MS:
                    if aborted:
                        break
                    if all(cell_key(n, fam, seed, M, lk) in done_keys
                           for lk in LINKS):
                        continue
                    A = A_full[:M]
                    u = u_full[:M]
                    from sklearn.covariance import LedoitWolf
                    lw = LedoitWolf(assume_centered=False).fit(A)
                    chol_lw = CholOp(lw.covariance_)
                    tr_lw = float(np.trace(lw.covariance_))
                    abar = A.mean(axis=0)
                    for link in LINKS:
                        key = cell_key(n, fam, seed, M, link)
                        if key in done_keys:
                            continue
                        if wall_prev + (time.time() - t0) > BUDGET_WALL_S:
                            aborted = True
                            print("BUDGET EXCEEDED - aborting per spec §0.5",
                                  flush=True)
                            break
                        t_c = time.time()
                        link_id = C.LINK_IDS[link]
                        B = np.stack([
                            L.simulate_bucket(
                                link, u[:, k], S,
                                rng_for(seed, fam_id, R62, PART, n, link_id, M, k))
                            for k in range(8)], axis=1)
                        rows = process_cell(A, B, X, side, n, fam, fam_id, seed,
                                            M, link, link_id, chol_lw, tr_lw,
                                            abar, cov_op)
                        append_rows(CSV, rows)
                        meta["cells"][key] = round(time.time() - t_c, 1)
                        done_keys.add(key)
                        meta["wall_accum_s"] = round(wall_prev + time.time() - t0, 1)
                        with open(META_JSON, "w") as f:
                            json.dump(meta, f, indent=2)
                        print("[g1] cell n=%d %s seed%d M=%d %s done (%.1fs)"
                              % (n, fam, seed, M, link, time.time() - t_c),
                              flush=True)
                del A_full, u_full

    meta["aborted_over_budget"] = aborted
    meta["wall_s"] = round(time.time() - t0, 1)
    meta["wall_accum_s"] = round(wall_prev + time.time() - t0, 1)
    meta["process_cpu_s"] = round(time.process_time() - cpu0, 1)
    with open(META_JSON, "w") as f:
        json.dump(meta, f, indent=2)

    gates = compute_gates(meta, aborted)
    gates["runtime"] = {"wall_accum_s": meta["wall_accum_s"],
                        "process_cpu_s": meta["process_cpu_s"],
                        "aborted_over_budget": bool(aborted)}
    with open(GATES_JSON, "w") as f:
        json.dump(gates, f, indent=2)
    print("G1 verdict:", gates["verdict"],
          "| cells done: %d" % len(meta["cells"]),
          "| abort=%s" % aborted, flush=True)
    return 2 if aborted else 0


if __name__ == "__main__":
    sys.exit(main())
