"""FROZEN R16 estimator diagnostic (self-contained; upload this ONE file).

This is the FROZEN-scale sibling of jitter_score_diag_colab.py, promoted to the
R16 section-5.3 frozen spec (docs/ROUND63_GPT_ROUND16_RULING_RAW.md 5.3, and
docs/R16_DIAGNOSTIC_RUNBOOK.md). It decides whether the replicated ~-6% end-point
low bias in the jitter-capped ridge-law peaks
(docs/R14_PREREGISTERED_PREDICTIONS.md) is a histogram/tail-floor estimator
artifact or a physical/theoretical residual.

=====================================================================
DIFF vs jitter_score_diag_colab.py (the 400k prep driver)
=====================================================================
The estimator mathematics is BYTE-IDENTICAL to the prep driver: _m_max,
simulate_multi, _binstats, fi_mi, fi_cs, fi_hist, _quad_peak, _cubic_peak are
copied verbatim. Only the run harness changed, exactly per R16 5.3:

  1. N_MC        : 400_000  ->  2_000_000        (5.3 "2,000,000 frames per load")
  2. SAVE_RAW_NL : False    ->  True             (5.3 "record (N,L) for every frame")
  3. ADDED       : N_BOOT = 200 block-bootstrap peak-location replicates
                   (5.3 "200 block-bootstrap replicates for the peak location;
                    peak extracted from a shape-preserving cubic interpolation
                    and, separately, a three-point quadratic; both reported").

The bootstrap RESAMPLES the recorded per-frame arrays and re-calls the UNCHANGED
fi_mi / fi_cs / fi_hist on the resampled frames, then re-extracts the peak with
the UNCHANGED _quad_peak / _cubic_peak. No estimator formula is touched. Because
the simulated frames are i.i.d. renewal paths (no serial dependence), the
"block bootstrap" reduces to a nonparametric frame resample (block length 1);
5.3 fixes no block length. Per load the base finite-difference step (dlog=0.01)
and base floor (alpha=0.5) feed the histogram peak bootstrap; the alpha/dlog
sensitivity SWEEPS remain point-estimate trajectories (the 5.4 rule tests
whether the histogram peak *moves* with alpha / as dlog^2, i.e. the trajectory,
not a per-alpha CI).

CONSTANT-FREE: nothing is fitted; no ridge-law constant is proposed (R16 5.4).

Colab usage (CPU high-RAM is enough; no GPU needed):
    !python jitter_score_diag_colab_frozen.py
Outputs -> /content/r16_diag/ :
    r16_diag_frozen_summary.json  (curves, point peaks, bootstrap peak CIs,
                                   tail masses, alpha/dlog sweeps)
    r16_diag_frozen_curves.npz
    raw_nl_cv{cv}.npz             (per-frame N,L; SAVE_RAW_NL)
    heartbeat.json                (updated every load for the watchdog)
"""
import json
import os
import time

import numpy as np

# ------------------------------------------------------------------ parameters
NU = 2000
TAU = 1.0
N_MC = 2_000_000              # FROZEN R16 5.3 (prep driver used 400_000)
CHUNK = 8_000
DLOGS = [0.0025, 0.005, 0.01, 0.02]
ALPHAS = [0.0, 0.05, 0.10, 0.50, 1.00]
DLOG_BASE = 0.01              # the currently-deployed audited step
ALPHA_BASE = 0.50            # the currently-deployed audited floor numerator
TAIL_THRESH = 50
OUTDIR = "/content/r16_diag"
SAVE_RAW_NL = True            # FROZEN R16 5.3: record (N,L) for every frame
N_BOOT = 200                  # FROZEN R16 5.3: 200 block-bootstrap peak replicates
BOOT_ALPHA = 0.05             # 95% percentile CI on the bootstrap peak location
BOOT_SEED0 = 916000           # base seed for the peak-location bootstrap

# R16 section 5.3 grids (the two replicated end points) + c=0 calibration
GRIDS = {
    0.02: np.linspace(8.5, 11.5, 31),
    0.20: np.linspace(1.75, 2.45, 29),
}
CAL_C0 = {"rho_exact": 22.2543, "grid": np.linspace(20.0, 24.5, 19)}
# measured histogram end points to be reproduced-or-not by the score estimator
MEASURED = {0.02: 9.62, 0.20: 2.08}
# R16 section-4 source-of-record corrected predictions
R16_PRED = {0.02: 10.22, 0.20: 2.17}


# ------------------------------------------------------------------ simulation
# (verbatim from jitter_score_diag_colab.py)
def _m_max(nu, rho, dmax):
    mean_ct = nu * rho / (1.0 + rho) * (1.0 + 4.0 * dmax)
    return int(mean_ct + 12.0 * np.sqrt(max(nu, 1.0)) + 50.0)


def simulate_multi(rho, cv, nu, n_mc, seed, mults, chunk=CHUNK):
    """One CRN draw set -> N,L at lam0 and counts at lam0*mult for each mult.

    mults must contain 1.0 (the base load). Returns (N0, L0, counts_dict) where
    counts_dict[mult] is the integer count array at lam0*mult.
    """
    lam0 = rho / TAU
    T = nu * TAU
    dmax = max(abs(np.log(m)) for m in mults) if mults else 0.0
    m_max = _m_max(nu, rho, dmax)
    rng = np.random.default_rng(seed)
    N0 = np.empty(n_mc, dtype=np.int64)
    L0 = np.empty(n_mc, dtype=np.float64)
    counts = {m: np.empty(n_mc, dtype=np.int64) for m in mults}
    if cv > 0:
        sig2 = np.log1p(cv * cv)
        mu = -0.5 * sig2
    pos = 0
    while pos < n_mc:
        nb = min(chunk, n_mc - pos)
        A = np.cumsum(rng.exponential(1.0, size=(nb, m_max)), axis=1)
        if cv > 0:
            B = TAU * np.exp(mu + np.sqrt(sig2) * rng.standard_normal((nb, m_max)))
        else:
            B = np.full((nb, m_max), TAU)
        D = np.concatenate([np.zeros((nb, 1)), np.cumsum(B, axis=1)[:, :-1]], axis=1)
        slack = T - D
        for m in mults:
            counts[m][pos:pos + nb] = (A <= (lam0 * m) * slack).sum(axis=1)
        n0 = counts[1.0][pos:pos + nb]
        if int(n0.max()) >= m_max:
            raise RuntimeError(f"count saturated m_max={m_max} at rho={rho}, cv={cv}")
        N0[pos:pos + nb] = n0
        # exact live time at lam0
        ar = np.arange(nb)
        prev = np.maximum(n0 - 1, 0)
        A_N = np.where(n0 >= 1, A[ar, prev], 0.0)
        D_Np1 = D[ar, np.minimum(n0, m_max - 1)]
        waits = A_N / lam0
        L0[pos:pos + nb] = waits + np.maximum(0.0, T - waits - D_Np1)
        pos += nb
    return N0, L0, counts


# ------------------------------------------------------------------ estimators
# (verbatim from jitter_score_diag_colab.py -- DO NOT alter this mathematics)
def _binstats(N, x):
    vals, inv, counts = np.unique(N, return_inverse=True, return_counts=True)
    s1 = np.bincount(inv, weights=x)
    s2 = np.bincount(inv, weights=x * x)
    mean = s1 / counts
    var = np.where(counts >= 2,
                   np.maximum(s2 - counts * mean * mean, 0.0) / np.maximum(counts - 1, 1),
                   0.0)
    return counts, mean, var


def fi_mi(N, L, lam, nu):
    n_mc = N.shape[0]
    counts, _, var_L = _binstats(N, L)
    use = counts >= 2
    e_var = np.sum((counts[use] / n_mc) * var_L[use])
    I = float(np.mean(N)) - lam * lam * float(e_var)
    tail = float(np.sum(counts[counts < TAIL_THRESH]) / n_mc)
    return I / nu, tail


def fi_cs(N, L, lam, nu):
    n_mc = N.shape[0]
    U = N - lam * L
    counts, mean_U, var_U = _binstats(N, U)
    use = counts >= 2
    I = float(np.sum((counts[use] / n_mc) * (mean_U[use] ** 2 - var_U[use] / counts[use])))
    return I / nu


def fi_hist(cm, c0, cp, nu, dlog, alpha):
    n_mc = c0.shape[0]
    lo = int(min(cm.min(), c0.min(), cp.min()))
    hi = int(max(cm.max(), c0.max(), cp.max())) + 1
    bins = np.arange(lo, hi + 1)
    p0 = np.histogram(c0, bins=bins)[0] / n_mc
    pp = np.histogram(cp, bins=bins)[0] / n_mc
    pm = np.histogram(cm, bins=bins)[0] / n_mc
    if alpha > 0:
        eps = alpha / n_mc
        pp = np.maximum(pp, eps)
        pm = np.maximum(pm, eps)
    else:
        # alpha=0: drop bins where either shifted pmf is empty (undefined log)
        keep = (pp > 0) & (pm > 0)
        p0, pp, pm = p0[keep], pp[keep], pm[keep]
    score = (np.log(pp) - np.log(pm)) / (2 * dlog)
    return float(np.sum(p0 * score ** 2)) / nu


# ------------------------------------------------------------------ peak finder
# (verbatim from jitter_score_diag_colab.py)
def _quad_peak(logr, J):
    k = int(np.argmax(J))
    k0 = min(max(k, 1), len(J) - 2)
    x = logr[k0 - 1:k0 + 2]
    y = J[k0 - 1:k0 + 2]
    den = (x[0] - x[1]) * (x[0] - x[2]) * (x[1] - x[2])
    a = (x[2] * (y[1] - y[0]) + x[1] * (y[0] - y[2]) + x[0] * (y[2] - y[1])) / den
    b = (x[2] ** 2 * (y[0] - y[1]) + x[1] ** 2 * (y[2] - y[0]) + x[0] ** 2 * (y[1] - y[2])) / den
    return float(np.exp(-b / (2 * a))) if a < 0 else float(np.exp(logr[k]))


def _cubic_peak(rhos, J):
    try:
        from scipy.interpolate import PchipInterpolator
    except Exception:
        return _quad_peak(np.log(rhos), J)
    lr = np.log(rhos)
    f = PchipInterpolator(lr, J)
    fine = np.linspace(lr[0], lr[-1], 4001)
    return float(np.exp(fine[int(np.argmax(f(fine)))]))


# ------------------------------------------------------------------ ADDED: R16 5.3 peak-location bootstrap
def _peak_quad(rhos, J):
    return _quad_peak(np.log(np.asarray(rhos, float)), np.asarray(J, float))


def _peak_cubic(rhos, J):
    return _cubic_peak(np.asarray(rhos, float), np.asarray(J, float))


def bootstrap_peaks(rhos, raw_N, raw_L, raw_cm, raw_cp, nu, n_boot, seed0, tag=""):
    """R16 5.3: 200 block-bootstrap peak-location replicates.

    Resample the recorded per-frame arrays (block length 1: i.i.d. renewal
    frames), recompute the FI curve with the UNCHANGED fi_mi / fi_cs / fi_hist,
    and re-extract the peak with the UNCHANGED _quad_peak / _cubic_peak. Both the
    three-point quadratic and the shape-preserving cubic peaks are reported for
    the two score estimators (MI, CS) and the base histogram estimator
    (dlog=0.01, alpha=0.5)."""
    rhos = np.asarray(rhos, float)
    nrho = len(rhos)
    rng = np.random.default_rng(seed0)
    keys = ["MI_quad", "MI_cubic", "CS_quad", "HIST_quad", "HIST_cubic"]
    reps = {k: np.empty(n_boot) for k in keys}
    t0 = time.time()
    for b in range(n_boot):
        Jmi = np.empty(nrho)
        Jcs = np.empty(nrho)
        Jh = np.empty(nrho)
        for i in range(nrho):
            Ni_full = raw_N[i]
            nmc = Ni_full.shape[0]
            idx = rng.integers(0, nmc, size=nmc)
            Ni = Ni_full[idx]
            Li = raw_L[i][idx]
            lam = float(rhos[i]) / TAU
            Jmi[i] = fi_mi(Ni, Li, lam, nu)[0]
            Jcs[i] = fi_cs(Ni, Li, lam, nu)
            Jh[i] = fi_hist(raw_cm[i][idx], Ni, raw_cp[i][idx], nu,
                            dlog=DLOG_BASE, alpha=ALPHA_BASE)
        reps["MI_quad"][b] = _peak_quad(rhos, Jmi)
        reps["MI_cubic"][b] = _peak_cubic(rhos, Jmi)
        reps["CS_quad"][b] = _peak_quad(rhos, Jcs)
        reps["HIST_quad"][b] = _peak_quad(rhos, Jh)
        reps["HIST_cubic"][b] = _peak_cubic(rhos, Jh)
        if (b + 1) % 20 == 0 or (b + 1) == n_boot:
            print(f"[boot {tag}] {b+1}/{n_boot} "
                  f"MI_cubic med={np.median(reps['MI_cubic'][:b+1]):.4f} "
                  f"({time.time()-t0:.0f}s)", flush=True)
    out = {}
    for k, v in reps.items():
        out[k] = {
            "median": float(np.median(v)),
            "mean": float(np.mean(v)),
            "lo": float(np.quantile(v, BOOT_ALPHA / 2)),
            "hi": float(np.quantile(v, 1 - BOOT_ALPHA / 2)),
            "std": float(v.std(ddof=1)),
            "reps": v.tolist(),
        }
    return out


def _write_heartbeat(payload):
    try:
        with open(os.path.join(OUTDIR, "heartbeat.json"), "w") as fh:
            json.dump(payload, fh)
    except Exception:
        pass


# ------------------------------------------------------------------ driver
def run_cv(cv, rhos, seed0):
    lam_mults = sorted({1.0} | {np.exp(s * d) for d in DLOGS for s in (-1, 1)})
    mexp = float(np.exp(-DLOG_BASE))
    pexp = float(np.exp(+DLOG_BASE))
    J_mi, J_cs, tails = [], [], []
    # J_hist[(dlog,alpha)] -> list over rhos
    J_hist = {(d, a): [] for d in DLOGS for a in ALPHAS}
    # retained per-frame arrays for the peak-location bootstrap (base dlog only)
    raw_N, raw_L, raw_cm, raw_cp = [], [], [], []
    t0 = time.time()
    for i, rho in enumerate(rhos):
        N, L, counts = simulate_multi(float(rho), cv, NU, N_MC,
                                      seed=seed0 + i, mults=lam_mults)
        lam = rho / TAU
        jmi, tail = fi_mi(N, L, lam, NU)
        J_mi.append(jmi)
        J_cs.append(fi_cs(N, L, lam, NU))
        tails.append(tail)
        for d in DLOGS:
            cm = counts[float(np.exp(-d))]
            cp = counts[float(np.exp(+d))]
            for a in ALPHAS:
                J_hist[(d, a)].append(fi_hist(cm, N, cp, NU, dlog=d, alpha=a))
        # retain raw frames (base dlog shifts) for the bootstrap; copy the two
        # base-shift count arrays so the rest of the counts dict can be freed
        raw_N.append(N)
        raw_L.append(L)
        raw_cm.append(counts[mexp].copy())
        raw_cp.append(counts[pexp].copy())
        del counts
        el = time.time() - t0
        print(f"[r16f cv={cv}] {i+1:2d}/{len(rhos)} rho={rho:7.4f} "
              f"J_MI={jmi:.5f} J_CS={J_cs[-1]:.5f} "
              f"J_HIST(base)={J_hist[(DLOG_BASE, ALPHA_BASE)][-1]:.5f} "
              f"tail<{TAIL_THRESH}={tail:.4f} ({el:.0f}s)", flush=True)
        _write_heartbeat({"stage": f"sweep cv={cv}", "load": i + 1,
                          "n_loads": len(rhos), "elapsed_s": round(el, 1),
                          "ts": time.time()})
    rhos = np.asarray(rhos, float)
    J_mi = np.asarray(J_mi)
    J_cs = np.asarray(J_cs)
    peaks = {
        "MI_quad": _quad_peak(np.log(rhos), J_mi),
        "MI_cubic": _cubic_peak(rhos, J_mi),
        "CS_quad": _quad_peak(np.log(rhos), J_cs),
        "HIST_base_quad": _quad_peak(np.log(rhos), np.asarray(J_hist[(DLOG_BASE, ALPHA_BASE)])),
        "HIST_base_cubic": _cubic_peak(rhos, np.asarray(J_hist[(DLOG_BASE, ALPHA_BASE)])),
    }
    hist_peak_vs_alpha = {a: _quad_peak(np.log(rhos), np.asarray(J_hist[(DLOG_BASE, a)]))
                          for a in ALPHAS}
    hist_peak_vs_dlog = {d: _quad_peak(np.log(rhos), np.asarray(J_hist[(d, ALPHA_BASE)]))
                         for d in DLOGS}

    # ---- ADDED: 200-replicate block-bootstrap peak-location CIs (R16 5.3) ----
    print(f"[r16f cv={cv}] point sweep done ({time.time()-t0:.0f}s); "
          f"starting {N_BOOT}-replicate peak bootstrap", flush=True)
    _write_heartbeat({"stage": f"bootstrap cv={cv}", "load": len(rhos),
                      "n_loads": len(rhos), "elapsed_s": round(time.time() - t0, 1),
                      "ts": time.time()})
    seed0_boot = BOOT_SEED0 + int(round(cv * 10000))
    peak_boot = bootstrap_peaks(rhos, raw_N, raw_L, raw_cm, raw_cp, NU,
                                N_BOOT, seed0_boot, tag=f"cv={cv}")

    # ---- SAVE_RAW_NL: record (N,L) for every frame (R16 5.3) ----
    if SAVE_RAW_NL:
        raw_path = os.path.join(OUTDIR, f"raw_nl_cv{cv}.npz")
        np.savez_compressed(
            raw_path,
            rhos=rhos,
            N=np.stack(raw_N).astype(np.int32),
            L=np.stack(raw_L).astype(np.float32),
        )
        print(f"[r16f cv={cv}] saved raw (N,L) -> {raw_path} "
              f"({os.path.getsize(raw_path)/1e6:.0f} MB)", flush=True)

    # free the retained frames before the next cv
    del raw_N, raw_L, raw_cm, raw_cp

    return {
        "cv": cv, "rhos": rhos.tolist(), "J_MI": J_mi.tolist(), "J_CS": J_cs.tolist(),
        "tail_mass": tails,
        "J_HIST": {f"dlog={d},alpha={a}": J_hist[(d, a)] for d in DLOGS for a in ALPHAS},
        "peaks": peaks,
        "peak_boot": peak_boot,
        "hist_peak_vs_alpha": hist_peak_vs_alpha,
        "hist_peak_vs_dlog": hist_peak_vs_dlog,
    }


def _dump(out):
    with open(os.path.join(OUTDIR, "r16_diag_frozen_summary.json"), "w") as fh:
        json.dump(out, fh, indent=2)


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    print(f"[r16f] FROZEN nu={NU} N_MC={N_MC} chunk={CHUNK} N_BOOT={N_BOOT} "
          f"SAVE_RAW_NL={SAVE_RAW_NL} dlogs={DLOGS} alphas={ALPHAS}", flush=True)
    print(f"[r16f] grids: cv=0.02 {len(GRIDS[0.02])}pts, cv=0.20 {len(GRIDS[0.20])}pts; "
          f"c=0 cal {len(CAL_C0['grid'])}pts around {CAL_C0['rho_exact']}", flush=True)
    out = {}
    t0 = time.time()
    for cv in (0.02, 0.20):
        out[str(cv)] = run_cv(cv, GRIDS[cv], seed0=int(cv * 1e6) + 101)
        _dump(out)   # per-cv checkpoint so a recycle keeps completed end points
    # deterministic c=0 calibration
    out["cal_c0"] = run_cv(0.0, CAL_C0["grid"], seed0=909091)
    out["cal_c0"]["rho_exact"] = CAL_C0["rho_exact"]
    _dump(out)

    np.savez_compressed(os.path.join(OUTDIR, "r16_diag_frozen_curves.npz"),
                        **{f"cv{k}_rhos": np.asarray(out[k]["rhos"]) for k in ("0.02", "0.2")
                           if k in out},
                        **{f"cv{k}_JMI": np.asarray(out[k]["J_MI"]) for k in ("0.02", "0.2")
                           if k in out},
                        **{f"cv{k}_JCS": np.asarray(out[k]["J_CS"]) for k in ("0.02", "0.2")
                           if k in out})

    print("\n[r16f] ===== peak summary (point | 95% bootstrap CI) =====", flush=True)
    for cv in (0.02, 0.20):
        r = out[str(cv)]
        pk = r["peaks"]
        pb = r["peak_boot"]
        print(f"[r16f] cv={cv}: MI peak quad/cubic = {pk['MI_quad']:.3f}/{pk['MI_cubic']:.3f} "
              f"| MI_cubic boot med={pb['MI_cubic']['median']:.3f} "
              f"CI[{pb['MI_cubic']['lo']:.3f},{pb['MI_cubic']['hi']:.3f}]", flush=True)
        print(f"[r16f] cv={cv}: CS peak quad = {pk['CS_quad']:.3f} "
              f"| CS_quad boot med={pb['CS_quad']['median']:.3f} "
              f"CI[{pb['CS_quad']['lo']:.3f},{pb['CS_quad']['hi']:.3f}]", flush=True)
        print(f"[r16f] cv={cv}: HIST base quad/cubic = {pk['HIST_base_quad']:.3f}/"
              f"{pk['HIST_base_cubic']:.3f} | HIST_cubic boot "
              f"CI[{pb['HIST_cubic']['lo']:.3f},{pb['HIST_cubic']['hi']:.3f}] "
              f"| R16 pred {R16_PRED[cv]} | measured {MEASURED[cv]}", flush=True)
        print(f"[r16f] cv={cv}: HIST peak vs alpha "
              + " ".join(f"a{a}={r['hist_peak_vs_alpha'][a]:.3f}" for a in ALPHAS), flush=True)
        print(f"[r16f] cv={cv}: HIST peak vs dlog "
              + " ".join(f"d{d}={r['hist_peak_vs_dlog'][d]:.3f}" for d in DLOGS), flush=True)
    cal = out["cal_c0"]
    cpb = cal["peak_boot"]
    print(f"[r16f] c=0 calibration: MI peak {cal['peaks']['MI_quad']:.3f} "
          f"(boot CI[{cpb['MI_quad']['lo']:.3f},{cpb['MI_quad']['hi']:.3f}], "
          f"exact {CAL_C0['rho_exact']}); HIST base {cal['peaks']['HIST_base_quad']:.3f}",
          flush=True)

    mi002 = out["0.02"]["peaks"]["MI_cubic"]
    mi020 = out["0.2"]["peaks"]["MI_cubic"]
    print(f"\n[r16f] wall={time.time()-t0:.0f}s  saved -> {OUTDIR}", flush=True)
    print(f"RESULT: MI_peak cv0.02={mi002:.3f} cv0.20={mi020:.3f} | "
          f"R16_pred 10.22/2.17 | measured_hist 9.62/2.08 | "
          f"HIST_base {out['0.02']['peaks']['HIST_base_cubic']:.3f}/"
          f"{out['0.2']['peaks']['HIST_base_cubic']:.3f} | N_MC={N_MC} N_BOOT={N_BOOT}",
          flush=True)


if __name__ == "__main__":
    main()
