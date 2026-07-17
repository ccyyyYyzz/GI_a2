"""Unit tests UT1-UT7 — spec §6. Run before any Phase; all must pass.
This is the ONLY layer where failures are implementation bugs to fix.
Writes results/unit_tests.json and the GAM2 divergence figure data.
"""
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
from gi_core.families import GammaFamily, make_family, marginal_var
from gi_core.images import load_truths
from gi_core.utils import rng_for

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")

UT_STREAM = 90  # dedicated UT stream id under the SEED0 system


def ut1():
    """GAUSS x LIN noiseless, M=1e5, n=256: SCORE-OR vs WHITEN-OR pixelwise
    relative diff <= 1e-2 (Gaussian-case identity; exact ratio M/(M-1))."""
    n, M = C.N_A, 100000
    fam = make_family("GAUSS", n)
    rng = rng_for(0, UT_STREAM, 1)
    A = fam.sample(M, rng)
    S, _ = fam.score(A)
    truths = load_truths(16)
    X = np.stack([truths[k] for k in C.IMAGES], axis=1)
    B = A @ X  # noiseless
    Bc = B - B.mean(axis=0, keepdims=True)
    x_score = -(S.T @ Bc) / (M - 1)
    x_whiten = fam.true_cov_op().solve(A.T @ Bc / M)
    denom = np.maximum(np.abs(x_whiten), np.abs(x_whiten).max(axis=0) * 1e-9)
    rel = np.abs(x_score - x_whiten) / denom
    worst = float(rel.max())
    return {"pass": bool(worst <= 1e-2), "max_pixel_rel_diff": worst}


def ut2():
    """Sample mean/variance vs analytic, rel diff <= 1%, per family (n=256).
    Plus MIX-LOGN oracle posterior state accuracy (report only)."""
    n, M = C.N_A, 100000
    out = {"families": {}, "pass": True}
    for name in ["GAUSS", "GAM1", "GAM2", "GAM3", "GAM4", "GAM8", "CORR-LOGN", "MIX-LOGN"]:
        fam = make_family(name, n)
        rng = rng_for(0, UT_STREAM, 2, C.FAMILY_IDS[name])
        A = fam.sample(M, rng)
        mu_hat = float(A.mean())
        var_hat = float(A.var())
        var_true = marginal_var(name)
        dmu = abs(mu_hat - 1.0)
        dvar = abs(var_hat - var_true) / var_true
        rec = {"mean": mu_hat, "var": var_hat, "var_analytic": var_true,
               "mean_rel_diff": dmu, "var_rel_diff": dvar,
               "pass": bool(dmu <= 0.01 and dvar <= 0.01)}
        if name == "MIX-LOGN":
            rec["oracle_state_accuracy"] = fam.oracle_state_accuracy(
                A[:20000], fam.last_states[:20000])
        out["families"][name] = rec
        out["pass"] = out["pass"] and rec["pass"]
    return out


def ut3():
    """Poisson mean calibration E[b] vs phi(U0) <= 1%; DT30: simulated mean
    count vs lambda*T/(1+lambda*tau) <= 1% (both photon levels)."""
    n, M = C.N_A, C.M
    fam = make_family("GAM4", n)
    rng = rng_for(0, UT_STREAM, 3)
    A = fam.sample(M, rng)
    x = load_truths(16)["cat"]
    u = A @ x
    out = {"links": {}, "pass": True}
    for link in ["LIN", "DT30", "FGAIN", "SAT30", "SAT50", "GAMMA07", "LOG"]:
        s = 1e4
        b = L.simulate_bucket(link, u, s, rng_for(0, UT_STREAM, 3, C.LINK_IDS[link]))
        eb = float(b.mean())
        target = float(L.phi_mean(link, u).mean())
        rel = abs(eb - target) / abs(target)
        rec = {"E_b": eb, "target": target, "rel_diff": rel, "pass": bool(rel <= 0.01)}
        out["links"][link] = rec
        out["pass"] = out["pass"] and rec["pass"]
    for s in C.PHOTONS_B:
        lam = s * u
        tau = C.DT30_TAU_COEF / s
        counts = L._dt30_counts(lam, s, rng_for(0, UT_STREAM, 33, int(s)))
        target = float((lam / (1.0 + lam * tau)).mean())
        got = float(counts.mean())
        rel = abs(got - target) / target
        rec = {"mean_count": got, "renewal_target": target, "rel_diff": rel,
               "pass": bool(rel <= 0.01)}
        out["links"]["DT30_counts_s%g" % s] = rec
        out["pass"] = out["pass"] and rec["pass"]
    return out


def ut4():
    """U0 = 1 (<= 0.01) for all families x all 8 images; link calibration
    phi(U0)/U0 vs design values (<= 0.01)."""
    n, M = C.N_A, C.M
    truths = load_truths(16)
    X = np.stack([truths[k] for k in C.IMAGES], axis=1)
    out = {"families": {}, "links": {}, "pass": True}
    for name in ["GAUSS", "GAM1", "GAM2", "GAM3", "GAM4", "GAM8", "CORR-LOGN", "MIX-LOGN"]:
        fam = make_family(name, n)
        A = fam.sample(M, rng_for(0, UT_STREAM, 4, C.FAMILY_IDS[name]))
        u0 = (A @ X).mean(axis=0)
        worst = float(np.abs(u0 - 1.0).max())
        rec = {"max_abs_dev_U0": worst, "pass": bool(worst <= 0.01)}
        out["families"][name] = rec
        out["pass"] = out["pass"] and rec["pass"]
    design = {"LIN": 1.0, "DT30": 0.7, "SAT30": 0.7, "SAT50": 0.5,
              "GAMMA07": 1.0, "LOG": 0.7, "FGAIN": 1.0}
    for link, tgt in design.items():
        got = float(L.phi_mean(link, np.array([1.0]))[0])
        rec = {"phi_at_U0": got, "design": tgt, "pass": bool(abs(got - tgt) <= 0.01)}
        out["links"][link] = rec
        out["pass"] = out["pass"] and rec["pass"]
    return out


def ut5():
    """Centered-witness analytic variance q_k^2 = k(n+1)/((k-2)Mw) vs Monte
    Carlo, ratio in [0.8, 1.25], k in {3,4,8} (n=256, Mw=1e5)."""
    n, Mw = C.N_A, 100000
    P, _ = W.build_probes(n)
    out = {"k": {}, "pass": True}
    for k in (3, 4, 8):
        fam = GammaFamily(n, k)
        rng = rng_for(0, UT_STREAM, 5, k)
        A = fam.sample(Mw, rng)
        S, _ = fam.score(A)
        vhat = W.witness_estimates(A, S, P)[Mw]
        err = W.witness_err(vhat, P)
        mc = float(np.mean(err ** 2))
        q2 = float(W.q_gamma(k, n, Mw) ** 2)
        ratio = mc / q2
        rec = {"mc_mean_err2": mc, "q2_analytic": q2, "ratio": ratio,
               "pass": bool(0.8 <= ratio <= 1.25)}
        out["k"][k] = rec
        out["pass"] = out["pass"] and rec["pass"]
    return out


def ut6():
    """GAM2 score second moment vs lower truncation eps: monotone growth,
    no saturation, eps in {1e-2..1e-5}. MC (N=4e9 scalar draws, chunked)
    + exact analytic overlay 4*E1(2eps) + (8eps-4)e^{-2eps}."""
    from scipy.special import exp1

    eps_list = [1e-2, 1e-3, 1e-4, 1e-5]
    N_total = 4_000_000_000
    chunk = 200_000_000
    rng = rng_for(0, UT_STREAM, 6)
    sums = np.zeros(len(eps_list))
    for _ in range(N_total // chunk):
        a = rng.gamma(2.0, 0.5, size=chunk)  # rate 2
        s2 = (1.0 / a - 2.0) ** 2
        for j, eps in enumerate(eps_list):
            sums[j] += s2[a > eps].sum()
    mc = sums / N_total
    analytic = [float(4 * exp1(2 * e) + (8 * e - 4) * np.exp(-2 * e)) for e in eps_list]
    inc = np.diff(mc)
    monotone = bool((inc > 0).all())
    no_sat = bool(inc[-1] >= 0.5 * inc[-2])
    return {"pass": monotone and no_sat, "eps": eps_list, "mc": mc.tolist(),
            "analytic": analytic, "monotone": monotone, "no_saturation": no_sat}


def ut7():
    """Main vs LS-diagnostic metric pipelines: both produced, key-disjoint,
    gates read only main keys (structural non-contamination check).
    Run at 64x64 — the scale at which LPIPS actually enters gates (AlexNet
    LPIPS is undefined at 16x16 and recorded as NaN there)."""
    truths = load_truths(64)
    x = truths["cat"]
    n = x.size
    rng = rng_for(0, UT_STREAM, 7)
    xhat = x + 0.05 * rng.standard_normal(n) * x.std()
    main = MET.main_metrics(xhat, x, 64, with_lpips=True)
    diag = MET.diagnostic_metrics(xhat, x, 64)
    disjoint = len(set(main) & set(diag)) == 0
    main_keys_ok = set(main) == {"PSNR", "SSIM", "LPIPS"}
    scale_dep = MET.main_metrics(3.0 * xhat, x, 64, with_lpips=False)
    flux_invariant = abs(scale_dep["PSNR"] - main["PSNR"]) < 1e-9
    lpips_finite = bool(np.isfinite(main["LPIPS"]))
    return {"pass": bool(disjoint and main_keys_ok and flux_invariant and lpips_finite),
            "main": main, "diagnostic": diag, "lpips_finite_at_64": lpips_finite,
            "key_disjoint": disjoint, "flux_scale_invariant": flux_invariant}


def main():
    os.makedirs(RESULTS, exist_ok=True)
    t0 = time.time()
    results = {}
    all_pass = True
    for name, fn in [("UT1", ut1), ("UT2", ut2), ("UT3", ut3), ("UT4", ut4),
                     ("UT5", ut5), ("UT6", ut6), ("UT7", ut7)]:
        t = time.time()
        try:
            r = fn()
        except Exception as ex:
            import traceback

            r = {"pass": False, "exception": repr(ex),
                 "traceback": traceback.format_exc()}
        r["runtime_s"] = round(time.time() - t, 2)
        results[name] = r
        all_pass = all_pass and r["pass"]
        print("%s: %s (%.1fs)" % (name, "PASS" if r["pass"] else "FAIL", r["runtime_s"]),
              flush=True)
    results["ALL_PASS"] = all_pass
    results["total_runtime_s"] = round(time.time() - t0, 2)
    with open(os.path.join(RESULTS, "unit_tests.json"), "w") as f:
        json.dump(results, f, indent=2)
    print("UT ALL:", "PASS" if all_pass else "FAIL", flush=True)
    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
