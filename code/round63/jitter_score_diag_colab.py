"""Colab-ready R16 estimator diagnostic (self-contained; upload this ONE file).

Runs the R16 section-5.3 frozen-diagnostic machinery at Colab scale for the two
replicated end points (lognormal cv in {0.02, 0.20}, nu=2000). At every load on
the R16 grids it computes, from ONE set of common-random-number draws:

  * J_MI   -- missing-information score estimator   (no histogram, no finite diff)
  * J_CS   -- bias-corrected conditional-score cross-check
  * J_HIST -- the audited histogram log-pmf central-difference estimator, over the
              R16 sweeps  dlog in {0.0025,0.005,0.01,0.02}  and
              floor alpha/N with alpha in {0,0.05,0.10,0.50,1.00}.

It then locates the FI peak for each estimator (three-point quadratic in log-rho
plus a shape-preserving cubic), so the R16 section-5.4 decision rule can be
applied: does the score peak sit at the theory value (10.22 / 2.17) while the
histogram peak drifts with alpha/N (=> tail-floor artifact), or does the score
peak independently reproduce the measured 9.62 / 2.08 (=> physical residual)?

CONSTANT-FREE: nothing is fitted; no ridge-law constant is proposed. Per R16 5.4
no new constant may be fitted before this diagnostic rules.

Colab usage:
    # Runtime: CPU high-RAM is enough (no GPU needed).
    !python jitter_score_diag_colab.py
Outputs summaries to /content/r16_diag/ (JSON + NPZ) and prints progress lines
plus a final RESULT: line. Bump N_MC to 2_000_000 and set SAVE_RAW_NL=True for
the fully frozen section-5.3 run once this 400k prep passes.

Derivation of the identity and the toy closed forms lives in the local sibling
jitter_score_diag.py; this file mirrors its numerics exactly.
"""
import json
import os
import time

import numpy as np

# ------------------------------------------------------------------ parameters
NU = 2000
TAU = 1.0
N_MC = 400_000                # task-specified Colab scale; frozen run bumps to 2e6
CHUNK = 8_000
DLOGS = [0.0025, 0.005, 0.01, 0.02]
ALPHAS = [0.0, 0.05, 0.10, 0.50, 1.00]
DLOG_BASE = 0.01              # the currently-deployed audited step
ALPHA_BASE = 0.50            # the currently-deployed audited floor numerator
TAIL_THRESH = 50
OUTDIR = "/content/r16_diag"
SAVE_RAW_NL = False           # True only for the 2e6 frozen run (huge files)

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


# ------------------------------------------------------------------ driver
def run_cv(cv, rhos, seed0):
    lam_mults = sorted({1.0} | {np.exp(s * d) for d in DLOGS for s in (-1, 1)})
    J_mi, J_cs, tails = [], [], []
    # J_hist[(dlog,alpha)] -> list over rhos
    J_hist = {(d, a): [] for d in DLOGS for a in ALPHAS}
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
        print(f"[r16 cv={cv}] {i+1:2d}/{len(rhos)} rho={rho:7.4f} "
              f"J_MI={jmi:.5f} J_CS={J_cs[-1]:.5f} "
              f"J_HIST(base)={J_hist[(DLOG_BASE, ALPHA_BASE)][-1]:.5f} "
              f"tail<{TAIL_THRESH}={tail:.4f} ({time.time()-t0:.0f}s)", flush=True)
    rhos = np.asarray(rhos, float)
    J_mi = np.asarray(J_mi); J_cs = np.asarray(J_cs)
    peaks = {
        "MI_quad": _quad_peak(np.log(rhos), J_mi),
        "MI_cubic": _cubic_peak(rhos, J_mi),
        "CS_quad": _quad_peak(np.log(rhos), J_cs),
        "HIST_base_quad": _quad_peak(np.log(rhos), np.asarray(J_hist[(DLOG_BASE, ALPHA_BASE)])),
        "HIST_base_cubic": _cubic_peak(rhos, np.asarray(J_hist[(DLOG_BASE, ALPHA_BASE)])),
    }
    # how the histogram peak drifts with the floor alpha (at base dlog)
    hist_peak_vs_alpha = {a: _quad_peak(np.log(rhos), np.asarray(J_hist[(DLOG_BASE, a)]))
                          for a in ALPHAS}
    # dlog^2 dependence at base alpha
    hist_peak_vs_dlog = {d: _quad_peak(np.log(rhos), np.asarray(J_hist[(d, ALPHA_BASE)]))
                         for d in DLOGS}
    return {
        "cv": cv, "rhos": rhos.tolist(), "J_MI": J_mi.tolist(), "J_CS": J_cs.tolist(),
        "tail_mass": tails,
        "J_HIST": {f"dlog={d},alpha={a}": J_hist[(d, a)] for d in DLOGS for a in ALPHAS},
        "peaks": peaks,
        "hist_peak_vs_alpha": hist_peak_vs_alpha,
        "hist_peak_vs_dlog": hist_peak_vs_dlog,
    }


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    print(f"[r16] nu={NU} N_MC={N_MC} chunk={CHUNK} dlogs={DLOGS} alphas={ALPHAS}",
          flush=True)
    print(f"[r16] grids: cv=0.02 {len(GRIDS[0.02])}pts, cv=0.20 {len(GRIDS[0.20])}pts; "
          f"c=0 cal {len(CAL_C0['grid'])}pts around {CAL_C0['rho_exact']}", flush=True)
    out = {}
    t0 = time.time()
    for cv in (0.02, 0.20):
        out[str(cv)] = run_cv(cv, GRIDS[cv], seed0=int(cv * 1e6) + 101)
    # deterministic c=0 calibration
    out["cal_c0"] = run_cv(0.0, CAL_C0["grid"], seed0=909091)
    out["cal_c0"]["rho_exact"] = CAL_C0["rho_exact"]

    with open(os.path.join(OUTDIR, "r16_diag_summary.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    np.savez_compressed(os.path.join(OUTDIR, "r16_diag_curves.npz"),
                        **{f"cv{k}_rhos": np.asarray(out[k]["rhos"]) for k in ("0.02", "0.2")
                           if k in out},
                        **{f"cv{k}_JMI": np.asarray(out[k]["J_MI"]) for k in ("0.02", "0.2")
                           if k in out})

    print("\n[r16] ===== peak summary =====", flush=True)
    for cv in (0.02, 0.20):
        r = out[str(cv)]
        pk = r["peaks"]
        print(f"[r16] cv={cv}: MI peak (quad/cubic) = "
              f"{pk['MI_quad']:.3f}/{pk['MI_cubic']:.3f} | CS {pk['CS_quad']:.3f} | "
              f"HIST base {pk['HIST_base_quad']:.3f}/{pk['HIST_base_cubic']:.3f} | "
              f"R16 pred {R16_PRED[cv]} | measured {MEASURED[cv]}", flush=True)
        print(f"[r16] cv={cv}: HIST peak vs alpha "
              + " ".join(f"a{a}={r['hist_peak_vs_alpha'][a]:.3f}" for a in ALPHAS), flush=True)
        print(f"[r16] cv={cv}: HIST peak vs dlog "
              + " ".join(f"d{d}={r['hist_peak_vs_dlog'][d]:.3f}" for d in DLOGS), flush=True)
    cal = out["cal_c0"]
    print(f"[r16] c=0 calibration: MI peak {cal['peaks']['MI_quad']:.3f} "
          f"(exact {CAL_C0['rho_exact']}); HIST base {cal['peaks']['HIST_base_quad']:.3f}",
          flush=True)

    mi002, mi020 = out["0.02"]["peaks"]["MI_cubic"], out["0.2"]["peaks"]["MI_cubic"]
    print(f"\n[r16] wall={time.time()-t0:.0f}s  saved -> {OUTDIR}", flush=True)
    print(f"RESULT: MI_peak cv0.02={mi002:.3f} cv0.20={mi020:.3f} | "
          f"R16_pred 10.22/2.17 | measured_hist 9.62/2.08 | "
          f"HIST_base {out['0.02']['peaks']['HIST_base_cubic']:.3f}/"
          f"{out['0.2']['peaks']['HIST_base_cubic']:.3f} | N_MC={N_MC}", flush=True)


if __name__ == "__main__":
    main()
