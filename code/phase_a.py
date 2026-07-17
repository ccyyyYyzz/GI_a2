"""Phase A — spec §7. 16x16 (n=256), 3 seeds.

A1: Stein validity phase diagram (STOP-level), witness with Mw=2e5.
A2: linear-regime honesty check (STOP only for bias; literal gate reported
    along with cosine reading and bias-vs-variance decomposition).
A3: nonlinear head-space screening (KILL-level).

Outputs: results/phaseA_witness.csv, results/phaseA_gates.json,
results/phaseA_gam1_exhibits/, results/phaseA_a2_metrics.csv,
results/phaseA_a3_metrics.csv.
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
from gi_core import witness as W
from gi_core.families import make_family
from gi_core.images import load_truths
from gi_core.utils import pearson, rng_for

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
SIDE = 16
N = C.N_A

A1_FAMILIES = ["GAUSS", "GAM1", "GAM2", "GAM3", "GAM4", "GAM8", "CORR-LOGN", "MIX-LOGN"]
A1_GATE_FAMILIES = ["GAM3", "GAM4", "GAM8", "CORR-LOGN", "GAUSS"]  # MIX informational
PREFIXES = [10000, 100000, 200000]


def cosine(a, b):
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na == 0 or nb == 0:
        return 0.0
    return float(np.dot(a, b) / (na * nb))


def run_a1():
    P, ptypes = W.build_probes(N)
    psums = P.sum(axis=0)
    pnorms = np.linalg.norm(P, axis=0)
    zero_sum = np.abs(psums) <= 1e-8 * pnorms
    rows = []
    fam_stats = {}
    for name in A1_FAMILIES:
        fam = make_family(name, N)
        rng = rng_for(0, C.FAMILY_IDS[name], C.STREAM_WITNESS)
        A = fam.sample(C.MW, rng)
        S, n_bad = fam.score(A)
        valid = np.isfinite(S).all(axis=1)
        if not valid.all():
            A, S = A[valid], S[valid]
        est = W.witness_estimates(A, S, P, centered=True, prefixes=PREFIXES)
        errs = {m: W.witness_err(est[m], P) for m in PREFIXES}
        err = errs[PREFIXES[-1]]
        if name in ("GAM3", "GAM4", "GAM8"):
            k = int(name[3:])
            q = np.full(64, W.q_gamma(k, N, A.shape[0]))
        else:
            q = W.interp_se(A, S, P, centered=True)
        for i in range(64):
            rows.append([name, ptypes[i], i, "%.6g" % err[i], "%.6g" % q[i]])
        stats = {
            "n_frames_used": int(A.shape[0]),
            "n_dropped_nonfinite_score": int(C.MW - A.shape[0]),
            "median_err": float(np.median(err)),
            "p90_err": float(np.percentile(err, 90)),
            "median_err_by_Mw": {str(m): float(np.median(errs[m])) for m in PREFIXES},
            "p90_over_median": float(np.percentile(err, 90) / np.median(err)),
        }
        if name == "GAM1":
            est_unc = W.witness_estimates(A, S, P, centered=False)[A.shape[0]]
            se_unc = W.interp_se(A, S, P, centered=False)
            nz = ~zero_sum
            cos_const = np.array([
                cosine(est_unc[:, i], psums[i] * np.ones(N)) for i in np.nonzero(nz)[0]])
            zs_ratio = (np.linalg.norm(est_unc[:, zero_sum], axis=0)
                        / (se_unc[zero_sum] * pnorms[zero_sum]))
            stats["gam1"] = {
                "n_nonzero_sum_probes": int(nz.sum()),
                "uncentered_cosine_with_const_min": float(cos_const.min()),
                "uncentered_cosine_with_const_median": float(np.median(cos_const)),
                "zero_sum_norm_over_SE_median": float(np.median(zs_ratio)),
                "zero_sum_norm_over_SE_p90": float(np.percentile(zs_ratio, 90)),
                "centered_err_median": float(np.median(err)),
            }
        fam_stats[name] = stats

    with open(os.path.join(RESULTS, "phaseA_witness.csv"), "w", newline="") as f:
        wcsv = csv.writer(f)
        wcsv.writerow(["family", "probe_type", "v_idx", "err", "q_theory"])
        wcsv.writerows(rows)

    gates = {"operationalizations": {
        "pass_gate": "per-probe ratio err/q; PASS iff median(ratio)<=1.5 and P90(ratio)<=2.0",
        "gam1_gate": ("(i) uncentered vhat on nonzero-sum probes: min cosine with (sum v)*1 > 0.99; "
                      "(ii) zero-sum probes: ||vhat_unc||/(SE*||v||) median<=1.5, P90<=2.0; "
                      "(iii) centered err median in [0.9, 1.1] (err ~ 1 = exact-zero collapse)"),
        "gam2_gate": ("(a) P90/median err ratio >= 3x GAM4's; "
                      "(b) log-log slope of median err vs Mw in {1e4,1e5,2e5} > -0.35"),
    }}
    stop = False
    for name in A1_GATE_FAMILIES:
        st = fam_stats[name]
        # per-probe ratio gate
        errv = None
        with open(os.path.join(RESULTS, "phaseA_witness.csv")) as f:
            rd = list(csv.reader(f))[1:]
        errv = np.array([float(r[3]) for r in rd if r[0] == name])
        qv = np.array([float(r[4]) for r in rd if r[0] == name])
        ratio = errv / qv
        ok = bool(np.median(ratio) <= C.A1_MEDIAN_FACTOR
                  and np.percentile(ratio, 90) <= C.A1_P90_FACTOR)
        gates[name] = {"median_ratio": float(np.median(ratio)),
                       "p90_ratio": float(np.percentile(ratio, 90)), "pass": ok}
        stop = stop or (not ok)
    g1 = fam_stats["GAM1"]["gam1"]
    gam1_ok = bool(g1["uncentered_cosine_with_const_min"] > 0.99
                   and g1["zero_sum_norm_over_SE_median"] <= 1.5
                   and g1["zero_sum_norm_over_SE_p90"] <= 2.0
                   and 0.9 <= g1["centered_err_median"] <= 1.1)
    gates["GAM1"] = {"failure_mode_confirmed": gam1_ok, **g1}
    if not gam1_ok:
        stop = True  # unexpected PASS of identity = boundary theory wrong -> STOP
    med_by_m = fam_stats["GAM2"]["median_err_by_Mw"]
    xs = np.log(np.array(PREFIXES, dtype=float))
    ys = np.log(np.array([med_by_m[str(m)] for m in PREFIXES]))
    slope = float(np.polyfit(xs, ys, 1)[0])
    ratio_gam2 = fam_stats["GAM2"]["p90_over_median"]
    ratio_gam4 = fam_stats["GAM4"]["p90_over_median"]
    sig_a = bool(ratio_gam2 >= 3.0 * ratio_gam4)
    sig_b = bool(slope > -0.35)
    gates["GAM2"] = {
        "p90_over_median": ratio_gam2, "gam4_p90_over_median": ratio_gam4,
        "tail_ratio_signature": sig_a, "loglog_slope_median_err": slope,
        "slope_signature": sig_b,
        "instability_signature_confirmed": bool(sig_a and sig_b),
        "anomaly": not (sig_a and sig_b),
    }
    gates["MIX-LOGN_informational"] = fam_stats["MIX-LOGN"]
    gates["family_stats"] = fam_stats
    gates["STOP"] = stop
    return gates


def gam1_exhibits():
    """Reconstruction-level GAM1 exhibits: uncentered constant collapse vs
    centered exact-zero collapse (LIN x s=1e4, seed 0)."""
    out_dir = os.path.join(RESULTS, "phaseA_gam1_exhibits")
    os.makedirs(out_dir, exist_ok=True)
    truths = load_truths(SIDE)
    fam = make_family("GAM1", N)
    A = fam.sample(C.M, rng_for(0, C.FAMILY_IDS["GAM1"], C.STREAM_PATTERN))
    S, _ = fam.score(A)
    info = {}
    from imageio.v2 import imwrite

    for img_name in ["cat", "airplane"]:
        x = truths[img_name]
        u = A @ x
        b = L.simulate_bucket("LIN", u, 1e4,
                              rng_for(0, C.FAMILY_IDS["GAM1"], C.STREAM_NOISE, 1))
        bc = b - b.mean()
        x_cent = -(S.T @ bc) / (C.M - 1)
        x_unc = -(S.T @ b) / C.M
        info[img_name] = {
            "centered_max_abs": float(np.abs(x_cent).max()),
            "centered_is_exact_zero_map": bool(np.abs(x_cent).max() == 0.0),
            "uncentered_std_over_mean": float(x_unc.std() / abs(x_unc.mean())),
            "uncentered_mean": float(x_unc.mean()),
            "bucket_mean": float(b.mean()),
        }
        def norm_png(v):
            v = v.reshape(SIDE, SIDE)
            rng_ = v.max() - v.min()
            if rng_ == 0:
                return np.zeros((SIDE, SIDE), dtype=np.uint8)
            return np.round(255 * (v - v.min()) / rng_).astype(np.uint8)

        imwrite(os.path.join(out_dir, "%s_truth.png" % img_name), norm_png(x))
        imwrite(os.path.join(out_dir, "%s_uncentered_const_collapse.png" % img_name),
                norm_png(x_unc))
        imwrite(os.path.join(out_dir, "%s_centered_zero_collapse.png" % img_name),
                norm_png(x_cent))
    with open(os.path.join(out_dir, "exhibit_values.json"), "w") as f:
        json.dump(info, f, indent=2)
    return info


def run_a2():
    truths = load_truths(SIDE)
    X = np.stack([truths[k] for k in C.IMAGES], axis=1)
    recs = {"WHITEN-OR": [], "SCORE-OR": []}
    rows = []
    for seed in C.SEEDS_A:
        fam = make_family("GAM4", N)
        ctx = E.RunContext(fam, seed, with_rankg=False)
        u = ctx.A @ X
        B = np.stack([
            L.simulate_bucket("LIN", u[:, k], 1e4,
                              rng_for(seed, C.FAMILY_IDS["GAM4"], C.STREAM_NOISE,
                                      C.LINK_IDS["LIN"], 4, k))
            for k in range(8)], axis=1)
        xw = E.whiten_or(ctx, B)
        xs = E.score_or(ctx, B)
        recs["WHITEN-OR"].append(xw)
        recs["SCORE-OR"].append(xs)
        for mname, xr in (("WHITEN-OR", xw), ("SCORE-OR", xs)):
            for k, img in enumerate(C.IMAGES):
                m = MET.main_metrics(xr[:, k], X[:, k], SIDE)
                mse = float(np.mean((MET.flux_match(xr[:, k], X[:, k]) - X[:, k]) ** 2))
                rows.append([img, "GAM4", "LIN", 1e4, seed, mname,
                             m["PSNR"], m["SSIM"], m["LPIPS"], mse])
    with open(os.path.join(RESULTS, "phaseA_a2_metrics.csv"), "w", newline="") as f:
        wcsv = csv.writer(f)
        wcsv.writerow(["image", "illum", "link", "photons", "seed", "method",
                       "PSNR", "SSIM", "LPIPS", "MSE_flux"])
        wcsv.writerows(rows)

    gates = {"per_method": {}}
    stop = False
    for mname in ("WHITEN-OR", "SCORE-OR"):
        stack = np.stack(recs[mname], axis=0)     # (3, n, 8)
        mean_rec = stack.mean(axis=0)
        pear = np.array([pearson(mean_rec[:, k], X[:, k]) for k in range(8)])
        cosv = np.array([cosine(mean_rec[:, k], X[:, k]) for k in range(8)])
        # bias-vs-variance decomposition: chi2 of mean residual vs seed-scatter SE
        resid = mean_rec - X                       # (n, 8)
        se = stack.std(axis=0, ddof=1) / np.sqrt(stack.shape[0])
        z = resid / np.maximum(se, 1e-30)
        chi2_per_img = (z ** 2).mean(axis=0)       # ~1 if variance-consistent
        gates["per_method"][mname] = {
            "pearson_mean_recon_per_image": pear.tolist(),
            "pearson_min": float(pear.min()), "pearson_mean": float(pear.mean()),
            "cosine_mean_recon_per_image": cosv.tolist(),
            "cosine_min": float(cosv.min()),
            "mean_resid_chi2_per_image": chi2_per_img.tolist(),
            "literal_gate_pass_all_images": bool((pear >= 0.995).all()),
        }
        if not (pear >= 0.995).all():
            stop = True
    # efficiency comparison (report-only)
    import pandas as pd

    df = pd.read_csv(os.path.join(RESULTS, "phaseA_a2_metrics.csv"))
    mpsnr = df.groupby("method")["PSNR"].mean().to_dict()
    gates["mean_PSNR"] = mpsnr
    gates["LINEAR_EFFICIENCY_ONLY"] = bool(
        mpsnr.get("SCORE-OR", -1) > mpsnr.get("WHITEN-OR", -1))
    gates["STOP_literal_pearson"] = stop
    # Bias-based STOP decision (registered in code BEFORE Phase A ran, see note):
    # STOP iff cosine(mean recon, truth) < 0.995 for any image/method (the
    # calibrated scale-invariant reading) OR chi2 of mean residual vs seed
    # scatter shows systematic bias (chi2 > 2 on any image).
    bias_stop = False
    for mname in ("WHITEN-OR", "SCORE-OR"):
        g = gates["per_method"][mname]
        if g["cosine_min"] < 0.995:
            bias_stop = True
        if max(g["mean_resid_chi2_per_image"]) > 2.0:
            bias_stop = True
    gates["STOP"] = bias_stop
    gates["note"] = (
        "Dual-reading protocol, registered pre-data: the literal Pearson>=0.995 "
        "reading is provably unpassable for SCORE-OR by finite-M variance alone "
        "(needs centered-energy fraction >=0.85, impossible for sum-normalized "
        "nonneg natural images at n=256, M=2e4, 3 seeds; SCORE-OR linear-regime "
        "variance = k/(k-2) x WHITEN-OR = 2x at k=4), hence it cannot "
        "operationalize 'bias'. The 0.995 level matches the scale-invariant "
        "cosine reading (predicted ~0.9989 SCORE / ~0.9994 WHITEN vs gate "
        "0.995). Per spec header 'STOP only for bias', STOP fires on bias "
        "evidence: cosine<0.995 or chi2 inconsistency. Both readings reported; "
        "adjudication note in REPORT.md.")
    return gates


def _a3_estimators(ctx, B, fam_name, seed):
    from joblib import Parallel, delayed

    from gi_core.mave import rmave_single_index

    out = {}
    out["SCORE-OR"] = E.score_or(ctx, B)
    out["WHITEN-OR"] = E.whiten_or(ctx, B)
    out["WHITEN-LW"] = E.whiten_lw(ctx, B)
    out["SIR-10"] = E.sir(ctx, B, 10)
    out["SIR-20"] = E.sir(ctx, B, 20)
    out["L-ISOTRON"] = E.l_isotron(ctx, B)
    out["GI"] = E.gi(ctx, B)
    out["MLE-OR"] = E.mle_or(ctx, B, "DT30", 1e4)
    seeds_int = [int(ctx.est_rng.integers(2 ** 31)) for _ in range(B.shape[1])]
    mave_cols = Parallel(n_jobs=8)(
        delayed(rmave_single_index)(ctx.A, B[:, k], np.random.default_rng(seeds_int[k]))
        for k in range(B.shape[1]))
    out["MAVE-16"] = np.stack(mave_cols, axis=1)
    return out


def run_a3():
    truths = load_truths(SIDE)
    X = np.stack([truths[k] for k in C.IMAGES], axis=1)
    rows = []
    for fam_name in ["GAM4", "CORR-LOGN"]:
        for seed in C.SEEDS_A:
            fam = make_family(fam_name, N)
            ctx = E.RunContext(fam, seed, with_rankg=False)
            u = ctx.A @ X
            B = np.stack([
                L.simulate_bucket("DT30", u[:, k], 1e4,
                                  rng_for(seed, C.FAMILY_IDS[fam_name], C.STREAM_NOISE,
                                          C.LINK_IDS["DT30"], 4, k))
                for k in range(8)], axis=1)
            ests = _a3_estimators(ctx, B, fam_name, seed)
            for mname, xr in ests.items():
                for k, img in enumerate(C.IMAGES):
                    m = MET.main_metrics(xr[:, k], X[:, k], SIDE)
                    d = MET.diagnostic_metrics(xr[:, k], X[:, k], SIDE)
                    rows.append([img, fam_name, "DT30", 1e4, seed, mname,
                                 m["PSNR"], m["SSIM"], m["LPIPS"],
                                 d["ANGERR_DEG"], d["PEARSON"]])
            print("A3 done: %s seed %d" % (fam_name, seed), flush=True)
    with open(os.path.join(RESULTS, "phaseA_a3_metrics.csv"), "w", newline="") as f:
        wcsv = csv.writer(f)
        wcsv.writerow(["image", "illum", "link", "photons", "seed", "method",
                       "PSNR", "SSIM", "LPIPS", "angerr", "pearson"])
        wcsv.writerows(rows)

    import pandas as pd

    df = pd.read_csv(os.path.join(RESULTS, "phaseA_a3_metrics.csv"))
    baselines = ["WHITEN-OR", "WHITEN-LW", "SIR-10", "SIR-20", "L-ISOTRON", "MAVE-16"]
    gates = {"combos": {}, "operationalization":
             "A3 PASS iff >=1 of the 2 combos passes (screening semantics); both reported"}
    any_pass = False
    for fam_name in ["GAM4", "CORR-LOGN"]:
        sub = df[df.illum == fam_name]
        pm = sub.groupby(["method", "image"])["PSNR"].mean().unstack()  # seed-mean
        base_best = pm.loc[baselines].max(axis=0)
        delta = pm.loc["SCORE-OR"] - base_best
        combo_pass = bool(delta.mean() >= C.A3_PSNR_DB and (delta > 0).sum() >= C.A3_MIN_POS)
        gates["combos"][fam_name] = {
            "delta_psnr_per_image": {k: float(v) for k, v in delta.items()},
            "mean_delta_psnr": float(delta.mean()),
            "n_images_positive": int((delta > 0).sum()),
            "best_baseline_per_image": {k: str(pm.loc[baselines, k].idxmax())
                                        for k in pm.columns},
            "pass": combo_pass,
        }
        any_pass = any_pass or combo_pass
    gates["A3_PASS"] = any_pass
    return gates


def main():
    os.makedirs(RESULTS, exist_ok=True)
    t0 = time.time()
    gates = {}
    print("=== A1 ===", flush=True)
    gates["A1"] = run_a1()
    gates["A1"]["gam1_exhibits"] = gam1_exhibits()
    print("A1 STOP:", gates["A1"]["STOP"], flush=True)
    print("=== A2 ===", flush=True)
    gates["A2"] = run_a2()
    print("A2 STOP (bias-based):", gates["A2"]["STOP"],
          "| literal-pearson flag:", gates["A2"]["STOP_literal_pearson"], flush=True)
    a_stop = gates["A1"]["STOP"] or gates["A2"]["STOP"]
    if a_stop:
        print("STOP-level failure -> aborting subsequent phases per spec.", flush=True)
    else:
        print("=== A3 ===", flush=True)
        gates["A3"] = run_a3()
        print("A3 PASS:", gates["A3"]["A3_PASS"], flush=True)
    gates["runtime_s"] = round(time.time() - t0, 1)
    with open(os.path.join(RESULTS, "phaseA_gates.json"), "w") as f:
        json.dump(gates, f, indent=2)
    print("PhaseA total %.1fs" % gates["runtime_s"], flush=True)
    return 1 if a_stop else 0


if __name__ == "__main__":
    sys.exit(main())
