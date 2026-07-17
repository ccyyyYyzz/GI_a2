"""ROUND62 Part 3 (P3): block-length crossover map — spec §3. 64x64, LIN link.

Drift nuisance: frame-wise multiplicative gain, log-AR(1):
  ln g_t - mu = rho (ln g_{t-1} - mu) + sigma_eta eps_t,
  sigma_ln = sqrt(ln(1+CV^2)), sigma_eta = sigma_ln sqrt(1-rho^2),
  mu = -sigma_ln^2/2 (stationary E[g] = 1), stationary init.
Bucket: b_t = Poisson(s g_t u_t)/s + N(0, sigma_r^2), s = 1e4.

Estimator family (block-centered correlation): for block length L, partition
frames sequentially, x_L = (1/M) sum_i (b_i - bbar_blk)(a_i - abar_blk)
  = (1/M) A^T btilde  exactly (block sums of btilde vanish).
Arms: L in {2,4,8,16,32,64} + RAW (no centering) + GLOBAL (= standard GI).
Readouts per arm: plain correlation and WHITEN-LW with Sigma_hat estimated on
block-centered patterns (RAW/GLOBAL use the global-centered LW). Secondary
ratio arm b_i/bbar_blk - 1 (NDHSI practical form), reported, not gated.

Patterns reuse the frozen Phase-B streams per (family, seed). Drift/noise from
round62 streams. Closed-basin: no score arms.
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
from gi_core import metrics as MET
from gi_core.families import make_family
from gi_core.images import load_truths
from gi_core.utils import CholOp, rng_for

ROOT = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(ROOT, "results", "round62_p3")
SIDE, N, S = 64, C.N_B, 1e4
R62 = 62
FAMS = ["GAM4", "MIX-LOGN"]
CVS = [0.02, 0.05, 0.1]
RHOS = [0.99, 0.9]
LS = [2, 4, 8, 16, 32, 64]
GATE_DB = 0.3
BUDGET_WALL_S = 4 * 3600


def block_means(v, L):
    """Per-frame block mean of v (M,) or (M,K) with sequential blocks of
    length L (last block possibly shorter)."""
    M = v.shape[0]
    out = np.empty_like(v, dtype=np.float64)
    for lo in range(0, M, L):
        sl = slice(lo, min(lo + L, M))
        out[sl] = v[sl].mean(axis=0, keepdims=True)
    return out


def sim_drift_gain(M, cv, rho, rng):
    sig_ln = np.sqrt(np.log1p(cv ** 2))
    sig_eta = sig_ln * np.sqrt(1.0 - rho ** 2)
    mu = -0.5 * sig_ln ** 2
    z = np.empty(M)
    z[0] = mu + sig_ln * rng.standard_normal()
    eps = rng.standard_normal(M)
    for t in range(1, M):
        z[t] = mu + rho * (z[t - 1] - mu) + sig_eta * eps[t]
    return np.exp(z)


def main():
    t0 = time.time()
    cpu0 = time.process_time()
    os.makedirs(OUT, exist_ok=True)
    truths = load_truths(SIDE)
    X = np.stack([truths[k] for k in C.IMAGES], axis=1)
    rows = []
    from skimage.metrics import structural_similarity

    def metrics_row(xhat, k):
        x_true = X[:, k]
        xs = MET.flux_match(xhat, x_true)
        dr = float(x_true.max())
        from gi_core.utils import psnr as _psnr

        p = float(_psnr(xs, x_true, dr))
        ss = float(structural_similarity(x_true.reshape(SIDE, SIDE),
                                         xs.reshape(SIDE, SIDE), data_range=dr))
        return p, ss

    aborted = False
    completed_units = []
    for fam_name in FAMS:
        fam_id = C.FAMILY_IDS[fam_name]
        for seed in C.SEEDS_EXT:
            if time.time() - t0 > BUDGET_WALL_S:
                aborted = True
                print("P3 BUDGET EXCEEDED - aborting per spec", flush=True)
                break
            fam = make_family(fam_name, N)
            A = fam.sample(C.M, rng_for(seed, fam_id, C.STREAM_PATTERN))
            M = A.shape[0]
            u = A @ X
            # per-L whitening operators from block-centered patterns
            from sklearn.covariance import LedoitWolf

            t_ctx = time.time()
            chol_global = CholOp(LedoitWolf(assume_centered=False).fit(A).covariance_)
            chol_L = {}
            for Lb in LS:
                At = A - block_means(A, Lb)
                lw = LedoitWolf(assume_centered=True).fit(At)
                chol_L[Lb] = CholOp(lw.covariance_)
                del At
            print("[p3] %s seed%d whitening ops built (%.1fs)"
                  % (fam_name, seed, time.time() - t_ctx), flush=True)

            for ci, cv in enumerate(CVS):
                for ri, rho in enumerate(RHOS):
                    B = np.empty((M, 8))
                    for k in range(8):
                        rng_n = rng_for(seed, fam_id, R62, 3, ci, ri, k)
                        g = sim_drift_gain(M, cv, rho, rng_n)
                        lam = np.maximum(S * g * u[:, k], 0.0)
                        B[:, k] = (rng_n.poisson(lam) / S
                                   + rng_n.normal(0.0, C.READOUT_COEF / S, M))
                    arms = {}
                    for Lb in LS:
                        Bt = B - block_means(B, Lb)
                        corr = A.T @ Bt / M
                        arms[("L%d" % Lb, "plain")] = corr
                        arms[("L%d" % Lb, "whiten")] = chol_L[Lb].solve(corr)
                        Ybar = block_means(B, Lb)
                        Yr = B / np.maximum(Ybar, 1e-300) - 1.0
                        arms[("L%d" % Lb, "ratio")] = A.T @ Yr / M
                    corr_raw = A.T @ B / M
                    arms[("RAW", "plain")] = corr_raw
                    arms[("RAW", "whiten")] = chol_global.solve(corr_raw)
                    Bg = B - B.mean(axis=0, keepdims=True)
                    corr_g = A.T @ Bg / M
                    arms[("GLOBAL", "plain")] = corr_g
                    arms[("GLOBAL", "whiten")] = chol_global.solve(corr_g)
                    for (arm, readout), xh in arms.items():
                        for k, img in enumerate(C.IMAGES):
                            p, ss = metrics_row(xh[:, k], k)
                            rows.append([fam_name, cv, rho, seed, img, arm,
                                         readout, "%.6g" % p, "%.6g" % ss])
                    print("[p3] %s seed%d cv=%g rho=%g done (%.1fs)"
                          % (fam_name, seed, cv, rho, time.time() - t0), flush=True)
            completed_units.append("%s|%d" % (fam_name, seed))
            del A, u, chol_L, chol_global
        else:
            continue
        break

    with open(os.path.join(OUT, "p3_crossover.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["illum", "cv", "rho", "seed", "image", "arm", "readout",
                    "PSNR", "SSIM"])
        w.writerows(rows)

    # ---- gate (ROUND61 preregistered): exists L in [4,64] beating BOTH L=2
    # and GLOBAL by >= +0.3 dB with 6/8 images, in >= 2 CV bins.
    # Operationalization (registered pre-run): primary slice rho = 0.99;
    # evaluated per readout family (plain / whiten), gate passes if either
    # readout exhibits the crossover; rho = 0.9 reported as sensitivity.
    import pandas as pd

    df = pd.DataFrame(rows, columns=["illum", "cv", "rho", "seed", "image",
                                     "arm", "readout", "PSNR", "SSIM"])
    df["PSNR"] = df["PSNR"].astype(float)
    df["cv"] = df["cv"].astype(float)
    df["rho"] = df["rho"].astype(float)
    gates = {"per_readout": {}, "operationalization":
             "primary slice rho=0.99, gate per readout family (plain/whiten), "
             "PASS if either; per-illum evaluation, pass if any illum family "
             "shows the crossover; ratio arm reported, never gated"}
    overall = []
    for readout in ("plain", "whiten"):
        detail = {}
        for fam_name in FAMS:
            for Lb in [4, 8, 16, 32, 64]:
                n_cv_ok = 0
                cv_detail = {}
                for cv in CVS:
                    sub = df[(df.illum == fam_name) & (df.cv == cv)
                             & (df.rho == 0.99) & (df.readout == readout)]
                    pm = sub.pivot_table(index="arm", columns="image",
                                         values="PSNR", aggfunc="mean")
                    if "L%d" % Lb not in pm.index:
                        continue
                    dL2 = pm.loc["L%d" % Lb] - pm.loc["L2"]
                    dGl = pm.loc["L%d" % Lb] - pm.loc["GLOBAL"]
                    ok = bool((dL2.mean() >= GATE_DB and (dL2 > 0).sum() >= 6)
                              and (dGl.mean() >= GATE_DB and (dGl > 0).sum() >= 6))
                    cv_detail[str(cv)] = {"mean_vs_L2": float(dL2.mean()),
                                          "mean_vs_GLOBAL": float(dGl.mean()),
                                          "n_pos_vs_L2": int((dL2 > 0).sum()),
                                          "n_pos_vs_GLOBAL": int((dGl > 0).sum()),
                                          "pass": ok}
                    n_cv_ok += ok
                detail["%s|L%d" % (fam_name, Lb)] = {
                    "cv_bins": cv_detail, "n_cv_pass": n_cv_ok,
                    "crossover": bool(n_cv_ok >= 2)}
                if n_cv_ok >= 2:
                    overall.append((readout, fam_name, Lb))
        gates["per_readout"][readout] = detail
    gates["MIDDLE_L_DOMINATES"] = bool(overall)
    gates["winning"] = [{"readout": r, "illum": f, "L": l} for r, f, l in overall]
    verdict = "MIDDLE_L_DOMINATES" if overall else "ENDPOINTS_ONLY"
    gates["budget_aborted"] = aborted
    gates["completed_units"] = completed_units
    if aborted:
        # partial grid: the preregistered gate is NOT the computed quantity
        verdict = "BUDGET_ABORT_PARTIAL(%s)" % verdict
    gates["verdict"] = verdict
    gates["runtime"] = {"wall_s": round(time.time() - t0, 1),
                        "process_cpu_s": round(time.process_time() - cpu0, 1)}
    with open(os.path.join(OUT, "p3_gates.json"), "w") as f:
        json.dump(gates, f, indent=2)
    print("P3 verdict:", gates["verdict"], "winning:", gates["winning"], flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
