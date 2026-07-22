#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
DLGI probe -- BAR 4 FINAL: profile-likelihood intervals in (log t_c, log CV).

R33 bar 4 needs CALIBRATED interval coverage.  Attempt 1 (parametric bootstrap,
percentile & log-space pivotal) UNDER-COVERS on the interior grid (t_c coverage
0.77-0.87, worst 0.80 at slow drift t_c=64), because the estimator's sampling
distribution is right-skewed at the slow-drift CRB floor (only ~T/t_c=32 correlation
lengths) and the bootstrap centred on the point estimate cannot reshape that skew.

Attempt 2 (this file, the authorized final method): PROFILE-LIKELIHOOD intervals in
the log parametrization, built from the WHITTLE likelihood of the per-exposure
log-gain residual series -- the standard, principled tool for OU/AR(1)+noise
time-series inference.  The OU+noise power spectrum is
    P_z(w; t_c, CV, sv2) = sigma_l^2 (1-phi^2)/(1 - 2 phi cos w + phi^2) + sv2,
    phi = exp(-1/t_c),  sigma_l^2 = ln(1+CV^2).
The Whittle log-likelihood ell_W = -1/2 sum_k [ log P_z(w_k) + I(w_k)/P_z(w_k) ] over
the M=N/2 independent Fourier frequencies (I = periodogram).  We profile the white
floor sv2 (absorbs residual scene-reconstruction contamination) and the other medium
parameter, and take the 95% LR interval  { theta : 2(ell_max - ell_prof(theta)) <=
chi2_{1,0.95}=3.841 }  in log(theta).  The LR interval is transformation-invariant,
so working in log fixes the right-skew that defeats the bootstrap.

Evaluated on a FRESH >=200-record protocol (6 scenes x 34 seeds = 204 records/cell),
interior cells only; fast-drift t_c=2 remains a declared edge (bar 7), excluded from
the claim.  Target: 95% coverage in [0.92, 0.98] on every interior cell.
Reports BOTH attempts.  Writes bar4_coverage.json.  CPU.
"""
import os, sys, json, time, platform
from datetime import datetime, timezone
import numpy as np
from scipy.optimize import minimize_scalar

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dl_common as C

T0 = time.time()
OUT = os.path.join(C.REPO, "results", "round63_next", "DUAL_LEDGER_PROBE")
def log(m): print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)

INTERIOR = [(16.0, 0.05), (16.0, 0.15), (16.0, 0.40),
            (64.0, 0.05), (64.0, 0.15), (64.0, 0.40)]
EDGE = [(2.0, 0.05), (2.0, 0.15), (2.0, 0.40)]
N_SEED = 34          # 6 scenes x 34 seeds = 204 records/cell (tight coverage SE ~1.5%)
N_BOOT = 200
CHI2_1_95 = 3.8415

# grids (log space)
LOG_TC = np.log(np.geomspace(0.5, 400.0, 58))
LOG_CV = np.log(np.geomspace(0.01, 1.2, 46))
LOG_SV2 = np.log(np.geomspace(1e-6, 5e-2, 14))

def _periodogram(zc, valid):
    y = np.where(valid > 0, zc, 0.0)
    N = len(y)
    Y = np.fft.rfft(y)
    I = (np.abs(Y) ** 2) / N               # periodogram, k=0..N//2
    w = 2.0 * np.pi * np.arange(len(I)) / N
    ku = slice(1, N // 2)                   # independent positive freqs (excl DC & Nyquist)
    return I[ku], w[ku]

def whittle_nll_grid(Ik, wk):
    """Return nll on the (LOG_TC x LOG_CV) grid, profiling sv2 (min over LOG_SV2).
    Vectorised over frequency."""
    phi = np.exp(-1.0 / np.exp(LOG_TC))                     # (Ntc,)
    cosw = np.cos(wk)                                       # (M,)
    shape = (1.0 - phi[:, None] ** 2) / (1.0 - 2.0 * phi[:, None] * cosw[None, :]
                                         + phi[:, None] ** 2)   # (Ntc, M)
    sig2 = np.log1p(np.exp(LOG_CV) ** 2)                    # sigma_l^2, (Ncv,)
    sv2 = np.exp(LOG_SV2)                                   # (Nsv,)
    nll = np.empty((len(LOG_TC), len(LOG_CV)))
    for j, s2 in enumerate(sig2):
        Pl = s2 * shape                                    # (Ntc, M)
        best = np.full(len(LOG_TC), np.inf)
        for v in sv2:
            Pz = Pl + v
            val = 0.5 * np.sum(np.log(Pz) + Ik[None, :] / Pz, axis=1)   # (Ntc,)
            best = np.minimum(best, val)
        nll[:, j] = best
    return nll

def profile_ci(zc, valid):
    """95% profile-likelihood CIs for (t_c, CV) in log space via the Whittle LR."""
    Ik, wk = _periodogram(zc, valid)
    nll = whittle_nll_grid(Ik, wk)
    nmin = nll.min()
    prof_tc = nll.min(axis=1)              # profile over cv (sv2 already profiled)
    prof_cv = nll.min(axis=0)
    def lr_interval(loggrid, prof):
        stat = 2.0 * (prof - nmin)
        below = stat <= CHI2_1_95
        if not below.any():
            return np.nan, np.nan
        idx = np.where(below)[0]
        lo_i, hi_i = idx[0], idx[-1]
        # linear-interpolate the crossing in log space for a smoother interval
        def cross(i_in, i_out):
            if i_out < 0 or i_out >= len(loggrid):
                return loggrid[i_in]
            s_in, s_out = stat[i_in], stat[i_out]
            if s_out == s_in:
                return loggrid[i_in]
            f = (CHI2_1_95 - s_in) / (s_out - s_in)
            return loggrid[i_in] + f * (loggrid[i_out] - loggrid[i_in])
        lo = cross(lo_i, lo_i - 1); hi = cross(hi_i, hi_i + 1)
        return np.exp(lo), np.exp(hi)
    tc_lo, tc_hi = lr_interval(LOG_TC, prof_tc)
    cv_lo, cv_hi = lr_interval(LOG_CV, prof_cv)
    return (tc_lo, tc_hi), (cv_lo, cv_hi)

# ---- attempt-1 conditional bootstrap (log-pivotal), for head-to-head ----
def boot_pivot_ci(x_hat, tc_hat, sig_hat, mu_hat, n_boot, seed):
    s = C.s_of_x(np.clip(x_hat, 0.0, None)); rng = np.random.default_rng(4242 + seed)
    phi = np.exp(-1.0 / max(tc_hat, 1e-3)); sd = sig_hat * np.sqrt(max(1 - phi ** 2, 1e-9))
    tcs = np.empty(n_boot); cvs = np.empty(n_boot)
    for b in range(n_boot):
        l = np.empty(C.N_EXP); l[0] = rng.standard_normal() * sig_hat
        for t in range(1, C.N_EXP):
            l[t] = phi * l[t - 1] + sd * rng.standard_normal()
        Yb = rng.poisson(np.maximum(np.exp(mu_hat + l) * C.PHI * s, 0)).astype(float)
        m = C.medium_estimate(Yb, x_hat, return_path=False)
        tcs[b] = m["tc_hat"]; cvs[b] = m["cv_hat"]
    cv_hat = C.cv_of_sigma_l(sig_hat)
    def piv(hat, star):
        star = star[np.isfinite(star) & (star > 0)]
        if len(star) < 20 or hat <= 0: return np.nan, np.nan
        lh = np.log(hat); q = np.percentile(np.log(star), [2.5, 97.5])
        return np.exp(2 * lh - q[1]), np.exp(2 * lh - q[0])
    return piv(tc_hat, tcs), piv(cv_hat, cvs)

def run_cell(tc, cv, n_seed, do_boot=True):
    sig = C.sigma_l_of_cv(cv)
    c = dict(prof_tc=0, prof_cv=0, boot_tc=0, boot_cv=0, n=0,
             wtc=[], wcv=[])
    for sc in C.SCENES:
        for seed in range(n_seed):
            Y, a_time, slot = C.simulate_record(sc, seed, tc, sig)
            r = C.joint_dual_ledger(Y, slot=slot, n_outer=2); m = r["med"]
            # time-ordered residuals for the Whittle
            zc, R, valid, mu = C._z_residuals(Y, r["x_hat"])
            order = np.argsort(slot); zc = zc[order]; valid = valid[order]
            (ptc_lo, ptc_hi), (pcv_lo, pcv_hi) = profile_ci(zc, valid)
            c["n"] += 1
            if np.isfinite(ptc_lo):
                c["prof_tc"] += int(ptc_lo <= tc <= ptc_hi); c["wtc"].append(ptc_hi - ptc_lo)
            if np.isfinite(pcv_lo):
                c["prof_cv"] += int(pcv_lo <= cv <= pcv_hi); c["wcv"].append(pcv_hi - pcv_lo)
            if do_boot:
                (btc_lo, btc_hi), (bcv_lo, bcv_hi) = boot_pivot_ci(
                    r["x_hat"], m["tc_hat"], m["sigma_l_hat"], m["mu_hat"], N_BOOT, seed)
                if np.isfinite(btc_lo):
                    c["boot_tc"] += int(btc_lo <= tc <= btc_hi)
                if np.isfinite(bcv_lo):
                    c["boot_cv"] += int(bcv_lo <= cv <= bcv_hi)
    n = max(c["n"], 1)
    return dict(tc=tc, cv=cv, n=c["n"],
                cover_profile_tc=c["prof_tc"] / n, cover_profile_cv=c["prof_cv"] / n,
                cover_boot_pivot_tc=c["boot_tc"] / n, cover_boot_pivot_cv=c["boot_cv"] / n,
                profile_width_tc_med=float(np.median(c["wtc"])) if c["wtc"] else np.nan,
                profile_width_cv_med=float(np.median(c["wcv"])) if c["wcv"] else np.nan)

def main():
    log(f"BAR-4 FINAL: profile-likelihood (Whittle, log-param) vs bootstrap-pivotal; "
        f"{N_SEED} seeds x 6 scenes = {6*N_SEED} records/cell, B={N_BOOT}")
    interior = {}; edge = {}
    for (tc, cv) in INTERIOR:
        r = run_cell(tc, cv, N_SEED, do_boot=True)
        interior[f"tc{int(tc)}_cv{int(cv*100)}"] = r
        log(f"  INTERIOR tc={tc:.0f} cv={cv}: PROFILE tc={r['cover_profile_tc']:.3f} "
            f"cv={r['cover_profile_cv']:.3f} | boot-piv tc={r['cover_boot_pivot_tc']:.3f} "
            f"cv={r['cover_boot_pivot_cv']:.3f}  (n={r['n']})")
    for (tc, cv) in EDGE:
        r = run_cell(tc, cv, N_SEED, do_boot=False)
        edge[f"tc{int(tc)}_cv{int(cv*100)}"] = r
        log(f"  EDGE tc={tc:.0f} cv={cv}: PROFILE tc={r['cover_profile_tc']:.3f} "
            f"cv={r['cover_profile_cv']:.3f} (declared edge, excluded)")
    prof = [interior[k]["cover_profile_tc"] for k in interior] + \
           [interior[k]["cover_profile_cv"] for k in interior]
    in_band = [0.92 <= v <= 0.98 for v in prof]
    verdict = "PASS" if all(in_band) else "FAIL"
    meta = dict(utc=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ"),
                python=platform.python_version(), numpy=np.__version__,
                method_attempt1="parametric bootstrap (percentile + log-space pivotal)",
                method_attempt2="Whittle profile-likelihood LR interval in (log t_c, log CV)",
                n_seed=N_SEED, records_per_cell=6 * N_SEED, n_boot=N_BOOT,
                target="95% coverage in [0.92,0.98] on all interior cells",
                runtime_s=time.time() - T0)
    blob = dict(meta=meta, interior=interior, edge_excluded=edge,
                interior_profile_coverages=prof,
                interior_all_in_92_98=bool(all(in_band)),
                bar4_verdict=verdict)
    fn = os.path.join(OUT, "bar4_coverage.json")
    json.dump(blob, open(fn, "w"), indent=2)
    log(f"BAR-4 FINAL verdict: {verdict}  (interior profile coverages in [0.92,0.98]: "
        f"{sum(in_band)}/{len(prof)})")
    log(f"saved {fn} ({time.time()-T0:.1f}s)")

if __name__ == "__main__":
    main()
