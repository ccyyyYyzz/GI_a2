"""ROUND62 Part 1 (T1-G0): MAVE anomaly audit — spec §1. 16x16, GAM4, 3 seeds.

Audits the Phase A MAVE-16 = 30.19 dB (+6.9 dB over strongest whitening):
(a) permutation null (hard gate) + variance-matched pure-noise arm;
(b) M-scaling diagnostic (soft);
(c) MLE redo with Gaussian-renewal likelihood + exact GammaCDF surface check;
(d) bandwidth x ridge sensitivity + GAM4 x LIN gain-source classification.

Reproduction discipline: unpermuted arm reuses the EXACT Phase A streams
(patterns, DT30/LIN noise, MAVE anchor rng, Isotron rng), so the audited
number is bit-identical to the anomaly. Permutations/noise/new arms draw from
fresh round62 streams (stream tag 62). Truth never touches MAVE; all
hyperparameters logged. No score/density arms anywhere (closed-basin rule).
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
from gi_core import estimators as E
from gi_core import links as L
from gi_core import metrics as MET
from gi_core.families import make_family
from gi_core.images import load_truths
from gi_core.mave import rmave_single_index
from gi_core.utils import rng_for
from mle_renewal import likelihood_surface_check, mle_renewal

ROOT = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(ROOT, "results", "round62_g0")
SIDE, N, FAM, S = 16, C.N_A, "GAM4", 1e4
FAM_ID = C.FAMILY_IDS[FAM]
R62 = 62
BUDGET_WALL_S = 3 * 3600  # dual-reported; wall-binding per ROUND59 precedent

SENS_BW = [0.5, 1.0, 2.0]
SENS_RIDGE = [0.1, 1.0, 10.0]


def psnr_main(xhat, x_true):
    return MET.main_metrics(xhat, x_true, SIDE, with_lpips=False)["PSNR"]


def mave_batch(A, B, anchor_rngs, n_jobs=8, **kw):
    from joblib import Parallel, delayed

    cols = Parallel(n_jobs=n_jobs)(
        delayed(rmave_single_index)(A, B[:, k], anchor_rngs[k], **kw)
        for k in range(B.shape[1]))
    return np.stack(cols, axis=1)


def a3_anchor_rngs(seed):
    return [rng_for(seed, FAM_ID, C.STREAM_ESTIMATOR, 7, k) for k in range(8)]


def main():
    t0 = time.time()
    cpu0 = time.process_time()
    os.makedirs(OUT, exist_ok=True)
    truths = load_truths(SIDE)
    X = np.stack([truths[k] for k in C.IMAGES], axis=1)
    rows = []  # (arm, config, seed, M, image, PSNR)
    audit = {"hyperparams": {
        "mave_default": dict(n_anchor=150, n_loc=400, n_iter=10, ridge=1e-4,
                             opg_n_anchor=150, opg_n_loc=1500, opg_ridge=1e-3,
                             bw_rule="rmave h=1.06*std(u)*M^-0.2; opg h2=median(d2_local)"),
        "sensitivity_grid": {"bw_scale": SENS_BW, "ridge_scale": SENS_RIDGE},
        "selection_rule": "no truth anywhere; sensitivity worst config = lowest "
                          "8-image mean PSNR on seed 0, then rerun on seeds 1,2"}}

    per_seed = {}
    aborted = False
    for seed in C.SEEDS_A:
        if time.time() - t0 > BUDGET_WALL_S:
            aborted = True
            print("G0 BUDGET EXCEEDED - aborting per spec", flush=True)
            break
        fam = make_family(FAM, N)
        ctx = E.RunContext(fam, seed, with_rankg=False)
        u = ctx.A @ X
        B_dt = np.stack([
            L.simulate_bucket("DT30", u[:, k], S,
                              rng_for(seed, FAM_ID, C.STREAM_NOISE,
                                      C.LINK_IDS["DT30"], 4, k))
            for k in range(8)], axis=1)
        B_lin = np.stack([
            L.simulate_bucket("LIN", u[:, k], S,
                              rng_for(seed, FAM_ID, C.STREAM_NOISE,
                                      C.LINK_IDS["LIN"], 4, k))
            for k in range(8)], axis=1)
        # baselines (same protocols as Phase A/B)
        x_lw = E.whiten_lw(ctx, B_dt)
        x_iso = E.l_isotron(ctx, B_dt, rng=rng_for(seed, FAM_ID, C.STREAM_ESTIMATOR,
                                                   C.LINK_IDS["DT30"], 10000, 1))
        x_lw_lin = E.whiten_lw(ctx, B_lin)

        # (a) unpermuted (exact Phase A reproduction) + permutations + noise
        x_mave = mave_batch(ctx.A, B_dt, a3_anchor_rngs(seed))
        for j in range(3):
            perm = rng_for(seed, FAM_ID, R62, 1, 10 + j).permutation(ctx.M)
            x_perm = mave_batch(ctx.A, B_dt[perm], a3_anchor_rngs(seed))
            for k, img in enumerate(C.IMAGES):
                rows.append(["mave_perm%d" % j, "default", seed, ctx.M, img,
                             psnr_main(x_perm[:, k], X[:, k])])
        B_noise = np.stack([
            rng_for(seed, FAM_ID, R62, 1, 20, k).normal(
                B_dt[:, k].mean(), B_dt[:, k].std(), ctx.M)
            for k in range(8)], axis=1)
        x_noise = mave_batch(ctx.A, B_noise, a3_anchor_rngs(seed))

        # (c) MLE-renewal, multi-start
        starts = {"WHITEN-LW": x_lw, "L-ISOTRON": x_iso, "MAVE": x_mave}
        x_mle, mle_diag = mle_renewal(ctx.A, B_dt, S, starts)
        audit.setdefault("mle_renewal", {})["seed%d" % seed] = mle_diag

        # (d) LIN arm MAVE (gain-source classification)
        x_mave_lin = mave_batch(
            ctx.A, B_lin, [rng_for(seed, FAM_ID, R62, 1, 5, k) for k in range(8)])

        for k, img in enumerate(C.IMAGES):
            rows.append(["mave_unperm", "default", seed, ctx.M, img,
                         psnr_main(x_mave[:, k], X[:, k])])
            rows.append(["mave_noise", "default", seed, ctx.M, img,
                         psnr_main(x_noise[:, k], X[:, k])])
            rows.append(["whiten_lw", "-", seed, ctx.M, img,
                         psnr_main(x_lw[:, k], X[:, k])])
            rows.append(["l_isotron", "-", seed, ctx.M, img,
                         psnr_main(x_iso[:, k], X[:, k])])
            rows.append(["mle_renewal", "best_of", seed, ctx.M, img,
                         psnr_main(x_mle[:, k], X[:, k])])
            rows.append(["mave_lin", "default", seed, ctx.M, img,
                         psnr_main(x_mave_lin[:, k], X[:, k])])
            rows.append(["whiten_lw_lin", "-", seed, ctx.M, img,
                         psnr_main(x_lw_lin[:, k], X[:, k])])

        # (b) M-scaling prefixes
        for Msub in (5000, 10000):
            Ap, Bp = ctx.A[:Msub], B_dt[:Msub]
            x_m = mave_batch(Ap, Bp,
                             [rng_for(seed, FAM_ID, R62, 1, 30, Msub, k)
                              for k in range(8)])
            from sklearn.covariance import LedoitWolf

            from gi_core.utils import CholOp

            lw = LedoitWolf(assume_centered=False).fit(Ap)
            cho = CholOp(lw.covariance_)
            Bc = Bp - Bp.mean(axis=0, keepdims=True)
            x_lwp = cho.solve(Ap.T @ Bc / Msub)
            for k, img in enumerate(C.IMAGES):
                rows.append(["mave_unperm", "default", seed, Msub, img,
                             psnr_main(x_m[:, k], X[:, k])])
                rows.append(["whiten_lw", "-", seed, Msub, img,
                             psnr_main(x_lwp[:, k], X[:, k])])

        # (c) exact likelihood surface check, seed 0 only
        if seed == 0:
            audit["likelihood_surface_check"] = likelihood_surface_check(
                u[:, 0], S, rng_for(0, FAM_ID, R62, 1, 40))
        per_seed[seed] = {"ctx": ctx, "B_dt": B_dt}
        print("[g0] seed %d arms done (%.1fs)" % (seed, time.time() - t0), flush=True)

    if aborted or len(per_seed) < 3:
        # partial run: write what exists, flag loudly, no gate evaluation
        with open(os.path.join(OUT, "g0_metrics.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["arm", "config", "seed", "M", "image", "PSNR"])
            w.writerows(rows)
        audit["budget_aborted"] = True
        audit["verdict"] = "BUDGET_ABORT_PARTIAL"
        audit["runtime"] = {"wall_s": round(time.time() - t0, 1),
                            "process_cpu_s": round(time.process_time() - cpu0, 1)}
        with open(os.path.join(OUT, "g0_audit.json"), "w") as f:
            json.dump(audit, f, indent=2)
        print("G0 verdict: BUDGET_ABORT_PARTIAL", flush=True)
        return 2

    # (d) sensitivity grid on seed 0; worst config on seeds 1, 2
    ctx0, B0 = per_seed[0]["ctx"], per_seed[0]["B_dt"]
    config_means = {}
    for bw in SENS_BW:
        for rs in SENS_RIDGE:
            if time.time() - t0 > BUDGET_WALL_S:
                aborted = True
                break
            tag = "bw%g_r%g" % (bw, rs)
            if bw == 1.0 and rs == 1.0:
                vals = [r[5] for r in rows if r[0] == "mave_unperm"
                        and r[2] == 0 and r[3] == C.M]
            else:
                x_c = mave_batch(ctx0.A, B0, a3_anchor_rngs(0),
                                 bw_scale=bw, ridge_scale=rs)
                vals = [psnr_main(x_c[:, k], X[:, k]) for k in range(8)]
                for k, img in enumerate(C.IMAGES):
                    rows.append(["mave_sens", tag, 0, C.M, img, vals[k]])
            config_means[tag] = float(np.mean(vals))
            print("[g0] sens %s mean %.2f dB (%.1fs)"
                  % (tag, config_means[tag], time.time() - t0), flush=True)
        if aborted:
            break
    if aborted:
        with open(os.path.join(OUT, "g0_metrics.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["arm", "config", "seed", "M", "image", "PSNR"])
            w.writerows(rows)
        audit["budget_aborted"] = True
        audit["verdict"] = "BUDGET_ABORT_PARTIAL(sensitivity_incomplete)"
        audit["sensitivity_partial"] = config_means
        audit["runtime"] = {"wall_s": round(time.time() - t0, 1),
                            "process_cpu_s": round(time.process_time() - cpu0, 1)}
        with open(os.path.join(OUT, "g0_audit.json"), "w") as f:
            json.dump(audit, f, indent=2)
        print("G0 verdict: BUDGET_ABORT_PARTIAL(sensitivity_incomplete)", flush=True)
        return 2
    worst_tag = min(config_means, key=config_means.get)
    bw_w, rs_w = [float(v) for v in worst_tag.replace("bw", "").replace("r", "").split("_")]
    audit["sensitivity"] = {"config_mean_psnr_seed0": config_means,
                            "worst_config": worst_tag}
    for seed in (1, 2):
        if worst_tag == "bw1_r1":
            break
        x_c = mave_batch(per_seed[seed]["ctx"].A, per_seed[seed]["B_dt"],
                         a3_anchor_rngs(seed), bw_scale=bw_w, ridge_scale=rs_w)
        for k, img in enumerate(C.IMAGES):
            rows.append(["mave_sens", worst_tag, seed, C.M, img,
                         psnr_main(x_c[:, k], X[:, k])])

    # ---- write CSV
    with open(os.path.join(OUT, "g0_metrics.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["arm", "config", "seed", "M", "image", "PSNR"])
        w.writerows(rows)

    # ---- gates
    import pandas as pd

    df = pd.DataFrame(rows, columns=["arm", "config", "seed", "M", "image", "PSNR"])
    full = df[df.M == C.M]

    def img_mean(arm, config=None):
        s = full[full.arm == arm]
        if config is not None:
            s = s[s.config == config]
        return s.groupby("image")["PSNR"].mean()

    unperm = img_mean("mave_unperm")
    perm_all = full[full.arm.str.startswith("mave_perm")]
    perm_mean = perm_all.groupby("image")["PSNR"].mean()
    noise_mean = img_mean("mave_noise")
    gate_a = {
        "unperm_mean": float(unperm.mean()),
        "perm_mean": float(perm_mean.mean()),
        "noise_mean": float(noise_mean.mean()),
        "sep_perm": float(unperm.mean() - perm_mean.mean()),
        "sep_noise": float(unperm.mean() - noise_mean.mean()),
        "perm_max_over_perms_seeds": float(perm_all["PSNR"].max()),
        "pass": bool(unperm.mean() - perm_mean.mean() >= 15.0
                     and perm_mean.mean() <= 12.0
                     and unperm.mean() - noise_mean.mean() >= 15.0
                     and noise_mean.mean() <= 12.0),
    }

    # (b) M-scaling gains
    ms = {}
    for Msub in (5000, 10000, C.M):
        sub = df[df.M == Msub]
        g = (sub[sub.arm == "mave_unperm"].groupby("image")["PSNR"].mean()
             - sub[sub.arm == "whiten_lw"].groupby("image")["PSNR"].mean())
        ms[str(Msub)] = float(g.mean())
    gains_sorted_by_M = [ms["5000"], ms["10000"], ms[str(C.M)]]
    gate_b = {"gain_vs_M": ms,
              "smoothing_prior_signature": bool(
                  gains_sorted_by_M[0] > gains_sorted_by_M[1] > gains_sorted_by_M[2])}

    # (d) gain source
    gain_dt = float((img_mean("mave_unperm") - img_mean("whiten_lw")).mean())
    gain_lin = float((img_mean("mave_lin") - img_mean("whiten_lw_lin")).mean())
    if gate_a["pass"]:
        src = ("GENERAL_EFFICIENCY" if gain_lin >= 0.5 * gain_dt
               else "NONLINEARITY_SPECIFIC")
    else:
        src = "ARTIFACT"
    gate_d = {"gain_dt30_vs_lw": gain_dt, "gain_lin_vs_lw": gain_lin,
              "classification": src}

    # G0 overall: worst-config MAVE vs max(WHITEN-LW, L-ISOTRON, MLE-renewal)
    if worst_tag == "bw1_r1":
        worst = unperm
    else:
        worst = full[(full.arm == "mave_sens")
                     & (full.config == worst_tag)].groupby("image")["PSNR"].mean()
    base = pd.concat([img_mean("whiten_lw"), img_mean("l_isotron"),
                      img_mean("mle_renewal")], axis=1).max(axis=1)
    delta = (worst - base).reindex(C.IMAGES)
    g0_pass = bool(gate_a["pass"] and int((delta >= 3.0).sum()) >= 6)
    audit.update({
        "gate_a_permutation": gate_a,
        "gate_b_M_scaling": gate_b,
        "gate_d_gain_source": gate_d,
        "g0_overall": {
            "worst_config": worst_tag,
            "delta_worstMAVE_minus_bestBaseline_per_image":
                {k: float(v) for k, v in delta.items()},
            "n_images_ge_3dB": int((delta >= 3.0).sum()),
            "mean_delta": float(delta.mean()),
            "operationalization": ">=6/8 images with (worst-config MAVE - "
                                  "best of {WHITEN-LW, L-ISOTRON, MLE-renewal}) "
                                  ">= +3.0 dB, on 3-seed per-image means "
                                  "(worst config frozen on seed-0 mean)",
            "PASS": g0_pass,
        },
        "verdict": ("PASS" if g0_pass else
                    ("ARTIFACT" if not gate_a["pass"] else "KILL(margin<3dB)")),
        "runtime": {"wall_s": round(time.time() - t0, 1),
                    "process_cpu_s": round(time.process_time() - cpu0, 1)},
    })
    with open(os.path.join(OUT, "g0_audit.json"), "w") as f:
        json.dump(audit, f, indent=2)
    print("G0 verdict:", audit["verdict"],
          "| gate_a pass:", gate_a["pass"],
          "| unperm %.2f perm %.2f noise %.2f" %
          (gate_a["unperm_mean"], gate_a["perm_mean"], gate_a["noise_mean"]),
          "| worst-config n>=3dB: %d/8" % audit["g0_overall"]["n_images_ge_3dB"],
          flush=True)
    return 0 if g0_pass else 1


if __name__ == "__main__":
    sys.exit(main())
