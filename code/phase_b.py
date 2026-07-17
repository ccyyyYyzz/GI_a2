"""Phase B — spec §8. 64x64 (n=4096), M=20000.

Core physics table: {GAM4, GAM8, CORR-LOGN, MIX-LOGN} x {LIN, DT30, FGAIN}
  x {1e3, 1e4} x 5 seeds x 8 images x estimators §4 (no MAVE-16;
  CLUSTER-WHITEN only MIX-LOGN; RANKG optional arm included).
Extended table: same illum x {SAT30, SAT50, GAMMA07, LOG} x 1e4 x 3 seeds.

Incremental checkpointing: completed (illum, seed, link, photons) combos are
skipped on restart. Budget guard: aborts (with report) past BUDGET_WALL_S.
Gate evaluation lives in gates_b.py (runs on the CSVs).
"""
import csv
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gi_core import config as C
from gi_core import estimators as E
from gi_core import links as L
from gi_core import metrics as MET
from gi_core.families import make_family
from gi_core.images import load_truths
from gi_core.utils import rng_for

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
SIDE = 64
N = C.N_B

FAMS = ["GAM4", "GAM8", "CORR-LOGN", "MIX-LOGN"]
CORE_LINKS = ["LIN", "DT30", "FGAIN"]
EXT_LINKS = ["SAT30", "SAT50", "GAMMA07", "LOG"]
BUDGET_WALL_S = 6 * 3600

CORE_CSV = os.path.join(RESULTS, "phaseB_core_metrics.csv")
EXT_CSV = os.path.join(RESULTS, "phaseB_ext_metrics.csv")
META_JSON = os.path.join(RESULTS, "phaseB_run_meta.json")
HEADER = ["image", "illum", "link", "photons", "seed", "method",
          "PSNR", "SSIM", "LPIPS", "angerr", "pearson"]


def done_combos(path):
    done = set()
    if os.path.exists(path):
        with open(path) as f:
            for r in csv.DictReader(f):
                done.add((r["illum"], int(r["seed"]), r["link"], float(r["photons"])))
    return done


def append_rows(path, rows):
    new = not os.path.exists(path)
    with open(path, "a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(HEADER)
        w.writerows(rows)


def run_estimators(ctx, B, fam_name, link, s, true_states):
    out = {}
    out["GI"] = E.gi(ctx, B)
    out["DGI"] = E.dgi(ctx, B)
    out["CORR"] = E.corr_gi(ctx, B)
    out["SIR-10"] = E.sir(ctx, B, 10)
    out["SIR-20"] = E.sir(ctx, B, 20)
    out["WHITEN-LW"] = E.whiten_lw(ctx, B)
    out["WHITEN-OR"] = E.whiten_or(ctx, B)
    out["SCORE-OR"] = E.score_or(ctx, B)
    out["RANKG"] = E.rankg(ctx, B)
    out["L-ISOTRON"] = E.l_isotron(ctx, B)
    out["MLE-OR"] = E.mle_or(ctx, B, link if link != "FGAIN" else "LIN", s)
    acc = None
    if fam_name == "MIX-LOGN":
        out["CLUSTER-WHITEN"], acc = E.cluster_whiten(ctx, B, true_states)
    return out, acc


def combo_rows(ests, X):
    rows = []
    # batch LPIPS across methods x images
    names = list(ests.keys())
    recs, refs = [], []
    dr_list = []
    for mname in names:
        for k in range(8):
            xs = MET.flux_match(ests[mname][:, k], X[:, k])
            recs.append(xs.reshape(SIDE, SIDE))
            refs.append(X[:, k].reshape(SIDE, SIDE))
            dr_list.append(float(X[:, k].max()))
    # LPIPS uses per-pair data range; process in groups of equal dr (per image)
    lp = np.zeros(len(recs))
    recs_arr = np.stack(recs)
    refs_arr = np.stack(refs)
    for k in range(8):
        idx = np.arange(k, len(recs), 8)
        lp[idx] = MET.lpips_batch(recs_arr[idx], refs_arr[idx], dr_list[k])
    return names, lp


def process_combo(ctx, X, fam_name, link, s, seed, true_states):
    from skimage.metrics import structural_similarity

    u = ctx.A @ X
    B = np.stack([
        L.simulate_bucket(link, u[:, k], s,
                          rng_for(seed, C.FAMILY_IDS[fam_name], C.STREAM_NOISE,
                                  C.LINK_IDS[link], int(s), k))
        for k in range(8)], axis=1)
    ests, acc = run_estimators(ctx, B, fam_name, link, s, true_states)
    names, lp = combo_rows(ests, X)
    rows = []
    j = 0
    for mname in names:
        for k, img in enumerate(C.IMAGES):
            x_true = X[:, k]
            xhat = ests[mname][:, k]
            xs = MET.flux_match(xhat, x_true)
            dr = float(x_true.max())
            from gi_core.utils import psnr as _psnr

            p = float(_psnr(xs, x_true, dr))
            ss = float(structural_similarity(
                x_true.reshape(SIDE, SIDE), xs.reshape(SIDE, SIDE), data_range=dr))
            d = MET.diagnostic_metrics(xhat, x_true, SIDE)
            rows.append([img, fam_name, link, s, seed, mname,
                         "%.6g" % p, "%.6g" % ss, "%.6g" % lp[j],
                         "%.6g" % d["ANGERR_DEG"], "%.6g" % d["PEARSON"]])
            j += 1
    return rows, acc


def main():
    os.makedirs(RESULTS, exist_ok=True)
    t0 = time.time()
    cpu0 = time.process_time()
    truths = load_truths(SIDE)
    X = np.stack([truths[k] for k in C.IMAGES], axis=1)
    done_core = done_combos(CORE_CSV)
    done_ext = done_combos(EXT_CSV)
    meta = {"combos": {}, "cluster_acc": {}, "aborted_over_budget": False}
    if os.path.exists(META_JSON):
        with open(META_JSON) as f:
            meta = json.load(f)

    for fam_name in FAMS:
        for seed in C.SEEDS_B:
            todo_core = [(lk, s) for lk in CORE_LINKS for s in C.PHOTONS_B
                         if (fam_name, seed, lk, s) not in done_core]
            todo_ext = []
            if seed in C.SEEDS_EXT:
                todo_ext = [(lk, 1e4) for lk in EXT_LINKS
                            if (fam_name, seed, lk, 1e4) not in done_ext]
            if not todo_core and not todo_ext:
                continue
            t_ctx = time.time()
            fam = make_family(fam_name, N)
            ctx = E.RunContext(fam, seed, with_rankg=True)
            true_states = getattr(fam, "last_states", None)
            print("[ctx] %s seed %d built in %.1fs (dropped=%d)"
                  % (fam_name, seed, time.time() - t_ctx, ctx.n_dropped), flush=True)
            for (lk, s), is_ext in ([(c, False) for c in todo_core]
                                    + [(c, True) for c in todo_ext]):
                if time.time() - t0 > BUDGET_WALL_S:
                    meta["aborted_over_budget"] = True
                    print("BUDGET EXCEEDED - aborting per spec §0.7", flush=True)
                    break
                t_c = time.time()
                rows, acc = process_combo(ctx, X, fam_name, lk, s, seed, true_states)
                append_rows(EXT_CSV if is_ext else CORE_CSV, rows)
                key = "%s|%d|%s|%g" % (fam_name, seed, lk, s)
                meta["combos"][key] = round(time.time() - t_c, 1)
                if acc is not None:
                    meta["cluster_acc"]["%s|%d" % (fam_name, seed)] = acc
                with open(META_JSON, "w") as f:
                    json.dump(meta, f, indent=2)
                print("[combo] %s seed%d %s s=%g: %.1fs" %
                      (fam_name, seed, lk, s, time.time() - t_c), flush=True)
            else:
                del ctx
                continue
            break
        else:
            continue
        break

    meta["wall_s"] = round(time.time() - t0, 1)
    meta["process_cpu_s"] = round(time.process_time() - cpu0, 1)
    with open(META_JSON, "w") as f:
        json.dump(meta, f, indent=2)
    print("PhaseB compute done: wall %.1fs cpu %.1fs abort=%s"
          % (meta["wall_s"], meta["process_cpu_s"], meta["aborted_over_budget"]),
          flush=True)
    return 2 if meta["aborted_over_budget"] else 0


if __name__ == "__main__":
    sys.exit(main())
