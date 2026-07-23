#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
EXPLORATORY / DEV-ONLY.  NOT preregistered.  No commits.

E2 shared machinery: NON-OU gain-path generators (domain extension beyond the
paper's validated scalar-OU domain) + model-adequacy (innovation-whiteness) test.

Gain-path classes (all length N_EXP, mean(a)=1 enforced, time order):
  ou        : scalar OU / AR(1) log-gain (the paper's validated domain; control)
  regime    : regime-switching OU -- Markov switch between a SLOW/weak and a
              FAST/strong (t_c, sigma) state; the autocovariance is a 2-exponential
              mixture, not a single exponential.
  heavytail : OU recursion with Student-t(df) innovations -> heavy-tailed increments
              (occasional jumps a single-Gaussian-OU cannot represent).
  quasiper  : OU + a sinusoidal drift (t_c-scale OU plus a coherent oscillation);
              the autocovariance oscillates, which a monotone exponential misfits.

Model-adequacy test = standardized one-step-ahead innovation whiteness (Ljung-Box).
Under the CORRECT temporal model the innovations are iid N(0,1); a mis-specified
temporal model leaves autocorrelation (Ljung-Box p small) and/or wrong variance.
"""
import os, sys
import numpy as np
from scipy.stats import chi2

REPO = "D:/GI_another"
sys.path.insert(0, os.path.join(REPO, "results", "round63_next", "DL_TOOL_EXPLORATORY", "code"))
import dl_tool_common as T
import dl_common as DL

N_EXP = DL.N_EXP
cv_of_sigma_l = T.cv_of_sigma_l


def _finalize(l):
    """demeaned-in-level: enforce mean(exp(l))=1 exactly (E[a]=1 physical gauge)."""
    a = np.exp(l)
    return a / a.mean()


def gen_ou(rng, tc, cv):
    sig = T.sigma_l_of_cv(cv)
    return T.ou_path(rng, sig, tc, N_EXP), dict(cls="ou", tc=tc, cv=cv)


def gen_regime(rng, tc_lo=6.0, tc_hi=60.0, cv_lo=0.10, cv_hi=0.45, p_switch=1/250.):
    """Two-state Markov regime switch: (fast,strong) <-> (slow,weak).  Per-step
    switch prob p_switch (mean dwell ~1/p).  AR(1) with regime-dependent phi & innov."""
    states = [dict(tc=tc_lo, sig=T.sigma_l_of_cv(cv_hi)),   # fast + strong
              dict(tc=tc_hi, sig=T.sigma_l_of_cv(cv_lo))]   # slow + weak
    s = int(rng.integers(0, 2))
    l = np.empty(N_EXP)
    st = states[s]; phi = np.exp(-1.0 / st["tc"])
    l[0] = st["sig"] * rng.standard_normal()
    reg = np.empty(N_EXP, dtype=int); reg[0] = s
    for n in range(1, N_EXP):
        if rng.random() < p_switch:
            s ^= 1; st = states[s]; phi = np.exp(-1.0 / st["tc"])
        sd = st["sig"] * np.sqrt(max(1 - phi * phi, 1e-6))
        l[n] = phi * l[n - 1] + sd * rng.standard_normal()
        reg[n] = s
    a = _finalize(l - l.mean())
    return a, dict(cls="regime", tc_lo=tc_lo, tc_hi=tc_hi, cv_lo=cv_lo, cv_hi=cv_hi)


def gen_heavytail(rng, tc, cv, df=3.0):
    """OU recursion with Student-t(df) innovations, variance-matched to the OU sd."""
    sig = T.sigma_l_of_cv(cv); phi = np.exp(-1.0 / tc)
    sd = sig * np.sqrt(max(1 - phi * phi, 1e-6))
    tstd = np.sqrt(df / (df - 2.0))                    # unit-variance scaling for t
    l = np.empty(N_EXP); l[0] = sig * rng.standard_normal()
    innov = rng.standard_t(df, size=N_EXP) / tstd
    for n in range(1, N_EXP):
        l[n] = phi * l[n - 1] + sd * innov[n]
    a = _finalize(l - l.mean())
    return a, dict(cls="heavytail", tc=tc, cv=cv, df=df)


def gen_quasiper(rng, tc, cv, period_lo=30.0, period_hi=120.0, amp_frac=0.7):
    """OU + coherent sinusoid.  amp_frac = sinusoid log-amplitude as a fraction of the
    OU log-sd (so total CV ~ combines both)."""
    sig = T.sigma_l_of_cv(cv)
    l_ou = np.log(np.maximum(T.ou_path(rng, sig, tc, N_EXP), 1e-12))
    P = float(rng.uniform(period_lo, period_hi)); ph = float(rng.uniform(0, 2 * np.pi))
    amp = amp_frac * sig
    l = l_ou + amp * np.sin(2 * np.pi * np.arange(N_EXP) / P + ph)
    a = _finalize(l - l.mean())
    return a, dict(cls="quasiper", tc=tc, cv=cv, period=P, amp_frac=amp_frac)


def sample_class(rng, cls):
    """Draw one gain path of the given class with randomized within-class params."""
    if cls == "ou":
        return gen_ou(rng, float(np.exp(rng.uniform(np.log(6), np.log(80)))),
                      float(rng.uniform(0.08, 0.50)))
    if cls == "regime":
        return gen_regime(rng)
    if cls == "heavytail":
        return gen_heavytail(rng, float(np.exp(rng.uniform(np.log(6), np.log(80)))),
                             float(rng.uniform(0.08, 0.50)), df=float(rng.choice([2.5, 3.0, 4.0])))
    if cls == "quasiper":
        return gen_quasiper(rng, float(np.exp(rng.uniform(np.log(8), np.log(60)))),
                            float(rng.uniform(0.06, 0.30)))
    raise ValueError(cls)


def integrated_tau(l):
    """Integrated autocorrelation time tau = 1 + 2 sum_{h>=1} rho(h) up to the first
    non-positive rho.  For an AR(1) tau=(1+phi)/(1-phi) ~ 2 t_c."""
    l = np.asarray(l, float); l = l - l.mean()
    v0 = float(l @ l) / len(l)
    tau = 1.0
    for h in range(1, 300):
        rh = float(l[:-h] @ l[h:]) / (len(l) - h) / max(v0, 1e-12)
        if rh <= 0:
            break
        tau += 2.0 * rh
    return float(tau)


def true_summary(a_time):
    """Class-agnostic ground-truth medium summary of a realized path:
      cv    = std(a)/mean(a)                            (fluctuation strength)
      tc    = single-exponential t_c fit to the CLEAN log-gain path (apples-to-apples
              with the mom_autocov estimator; the oracle-monitor's own readout)
      tau   = integrated autocorrelation time of the clean log-gain (model-free scale)."""
    a = np.asarray(a_time, float)
    cv = float(a.std() / max(a.mean(), 1e-12))
    l = np.log(np.maximum(a, 1e-12)); l = l - l.mean()
    tc, _ = T.mom_autocov(l, np.ones_like(l))
    if not np.isfinite(tc):
        tc = 1e3
    return dict(cv=cv, tc=float(tc), tau=integrated_tau(l))


def emp_autocov_tau(zc, valid, Hmax=60):
    """LEARNED-AUTOCOVARIANCE medium readout: noise-deconvolved integrated
    correlation time from the empirical autocovariance of the demeaned log-gain
    observations (gamma(h>=1) is Poisson-white-noise-free), NOT assuming a single
    exponential.  tau = 1 + 2 sum_{h>=1} gamma(h)/gamma(0) with gamma(0) extrapolated
    from the h=1,2 signal (removing the white obs-noise spike).  Robust for regime /
    quasi-periodic where the single-exponential t_c is biased."""
    v = zc * valid; vb = valid.astype(bool); n = len(zc)
    g = np.full(Hmax + 1, np.nan)
    for h in range(Hmax + 1):
        m = vb[:n - h] & vb[h:]
        if m.sum() > 10:
            g[h] = float(np.mean(v[:n - h][m] * v[h:][m]))
    if not (np.isfinite(g[1]) and np.isfinite(g[2]) and g[1] > 0):
        return float("nan"), float("nan")
    # extrapolate signal gamma(0) from the first two SIGNAL lags (log-linear), so the
    # white-noise spike in the raw gamma(0) is removed
    if g[2] > 0:
        slope = np.log(g[2]) - np.log(g[1])
        g0 = float(np.exp(np.log(g[1]) - slope))
    else:
        g0 = float(g[1])
    g0 = max(g0, g[1])
    tau = 1.0
    for h in range(1, Hmax + 1):
        if np.isfinite(g[h]) and g[h] > 0:
            tau += 2.0 * g[h] / g0
        else:
            break
    cv = cv_of_sigma_l(np.sqrt(max(g0, 1e-12)))
    return float(tau), float(cv)


# ------------------- record + z-residual helpers (true-scene, clean gain test) -----
def z_from_true_scene(Y, x_true, gauge="counts"):
    """zc (demeaned log-gain obs), R, valid, mu_hat from the record and the TRUE scene
    (isolates the temporal-GAIN model adequacy from reconstruction error)."""
    zc, R, valid, mu_hat = T._z_residuals_cfg(Y, x_true, gauge=gauge)
    return zc, R, valid, mu_hat


def simulate(x, a_time, rng):
    s = T.s_of_x(x)
    return rng.poisson(np.maximum(a_time * T.PHI * s, 0.0)).astype(np.float64)


# ------------------- innovation-whiteness (Ljung-Box) -------------------
def ljung_box(res, K=20):
    """Ljung-Box Q on the standardized innovations `res` over lags 1..K.  Returns
    (Q, p, dof).  Small p => residual autocorrelation => model mis-specified."""
    r = np.asarray(res, float); r = r - r.mean()
    n = len(r); v = float(r @ r) / n
    if v <= 0 or n <= K + 2:
        return float("nan"), float("nan"), K
    Q = 0.0
    for h in range(1, K + 1):
        rho = float(r[:-h] @ r[h:]) / (n - h) / v
        Q += rho * rho / (n - h)
    Q *= n * (n + 2)
    p = float(chi2.sf(Q, K))
    return float(Q), p, K


def ou_innovations(zc, R, valid, tc_hat, sig_hat):
    """Standardized one-step-ahead OU-Kalman innovations (the frozen model's own
    adequacy residual)."""
    phi = np.exp(-1.0 / max(tc_hat, 1e-3))
    return DL.kalman_innovations(zc, R, valid.astype(bool), phi, max(sig_hat, 1e-6))


def adequacy_pass(p, innov_std, p_thr=0.01, std_tol=0.25):
    """PASS = white (Ljung-Box p > p_thr) AND innovation variance ~ 1."""
    return bool((p > p_thr) and (abs(innov_std - 1.0) < std_tol))
