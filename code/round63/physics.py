"""ROUND63 detector physics — decoupled parameterization (spec §2).

Differences from the ROUND59 harness (deliberate, campaign-defining):
- tau is a FIXED physical detector constant; flux Phi, frame time T and pattern
  count M vary independently. Dimensionless axes rho = tau*E[lambda], nu = T/tau.
- Arrival rate per frame: lambda_i = Phi * <a_i, x> + d  (dark/background d
  enters BEFORE the dead time and occupies the detector).
- Bucket is the raw digital count N (photon counting has no analog readout);
  optional additive Gaussian sigma_b is an ablation knob, default 0.
- Frame-start conventions: 'active' (detector ready at t=0; first record waits
  Exp(lambda)) is the physical main arm; 'delayed' reproduces the ROUND59
  convention (tau + Exp); 'continuous' carries detector state across frames
  with a guard interval g between frames.
- Estimator-side models are named honestly: gaussian-renewal QUASI-likelihood
  (moment-matched CLT), with the exact renewal count likelihood as a
  small-scale reference (log-CDF-difference, numerically stable).

All simulation paths are exact event/interval-level Monte Carlo; the
vectorized fast paths cover the main arms, and a general event-loop path
covers afterpulsing / paralyzable / continuous ablations.
"""
import numpy as np
from dataclasses import dataclass, field


@dataclass
class Detector:
    tau: float                     # dead time [s]
    paralyzable: bool = False
    p_ap: float = 0.0              # afterpulse probability per recorded event
    ap_tau: float = 1e-6           # afterpulse delay scale [s] (exponential)
    dark: float = 0.0              # dark/background rate [1/s], pre-dead-time
    tau_jitter_cv: float = 0.0     # per-event dead-time jitter CV
    start_mode: str = "active"     # 'active' | 'delayed' | 'continuous'
    guard: float = 0.0             # inter-frame guard interval [s] (continuous)


# ----------------------------------------------------------------------
# exact simulation
# ----------------------------------------------------------------------
def _counts_nonpar_vec(lam, T, tau, rng, start_mode, chunk=2048):
    """Vectorized exact non-paralyzable renewal counts (no afterpulsing/jitter).

    Record times: active-start  S_m = (m-1)*tau + Gamma(m, lam)
                  delayed-start S_m = m*tau + Gamma(m, lam)
    Implemented by chunked interval accumulation (first interval Exp or
    tau+Exp, subsequent tau+Exp), exact for any lam including 0.
    """
    lam = np.asarray(lam, dtype=np.float64)
    out = np.empty(lam.shape[0], dtype=np.int64)
    hard_cap = int(np.floor(T / tau)) + 2
    first_extra = 0.0 if start_mode == "active" else tau
    for lo in range(0, lam.shape[0], chunk):
        lam_c = lam[lo:lo + chunk]
        mc = lam_c.shape[0]
        counts = np.zeros(mc, dtype=np.int64)
        t_acc = np.zeros(mc)
        active = lam_c > 0
        first = np.ones(mc, dtype=bool)
        while active.any():
            idx = np.nonzero(active)[0]
            lam_a = lam_c[idx]
            mean_n = lam_a * T / (1.0 + lam_a * tau) + 1.0
            k_alloc = int(np.clip(np.max(mean_n + 8.0 * np.sqrt(mean_n)),
                                  16, hard_cap))
            E = rng.exponential(1.0, size=(idx.size, k_alloc)) / lam_a[:, None]
            E[:, 0] += np.where(first[idx], first_extra, tau) - tau
            times = t_acc[idx, None] + np.cumsum(tau + E, axis=1)
            counts[idx] += (times <= T).sum(axis=1)
            still = times[:, -1] <= T
            t_acc[idx] = times[:, -1]
            first[idx] = False
            active[:] = False
            active[idx[still]] = True
            if k_alloc >= hard_cap:
                break
        out[lo:lo + chunk] = counts
    return out


def _counts_general_eventloop(lam, T, det, rng):
    """General exact event-level path: paralyzable, afterpulsing, dead-time
    jitter, continuous acquisition. Slower; used for ablation arms.

    Arrivals = Poisson(lam) primaries + branching afterpulses (each RECORDED
    event spawns an afterpulse arrival at +Exp(ap_tau) w.p. p_ap; afterpulses
    are arrivals: they can be recorded, trigger dead time, and branch again).
    Non-paralyzable: record iff t >= ready; ready = t + tau_eff.
    Paralyzable: every ARRIVAL extends busy window; record iff t >= busy-end
    at arrival, then busy-end = t + tau_eff (arrivals during busy extend it).
    """
    lam = np.asarray(lam, dtype=np.float64)
    M = lam.shape[0]
    N = np.zeros(M, dtype=np.int64)
    carry = 0.0  # detector not-ready time carried across frames (continuous)
    for i in range(M):
        lam_i = float(lam[i])
        if det.start_mode == "continuous":
            ready = max(0.0, carry)
        elif det.start_mode == "delayed":
            ready = det.tau
        else:
            ready = 0.0
        n_arr = rng.poisson(lam_i * T) if lam_i > 0 else 0
        arrivals = list(np.sort(rng.uniform(0.0, T, size=n_arr)))
        j = 0
        while j < len(arrivals):
            t = arrivals[j]
            j += 1
            tau_eff = det.tau
            if det.tau_jitter_cv > 0:
                tau_eff = max(0.0, det.tau * (1.0 + det.tau_jitter_cv
                                              * rng.standard_normal()))
            if det.paralyzable:
                recorded = t >= ready
                ready = max(ready, t + tau_eff)  # every arrival extends
            else:
                recorded = t >= ready
                if recorded:
                    ready = t + tau_eff
            if recorded:
                N[i] += 1
                if det.p_ap > 0 and rng.random() < det.p_ap:
                    t_ap = t + rng.exponential(det.ap_tau)
                    if t_ap < T:
                        # insert keeping sorted order
                        import bisect

                        bisect.insort(arrivals, t_ap)
        carry = ready - T + det.guard if det.start_mode == "continuous" else 0.0
    return N


def simulate_counts(u, Phi, T, det, rng, sigma_b=0.0):
    """u: (M,) frame energies <a_i, x>; returns bucket b (float64) and raw N.

    lambda_i = Phi*u_i + dark. Uses the fast vectorized path when exact
    (non-paralyzable, no afterpulsing/jitter, non-continuous), else the
    general event loop."""
    lam = Phi * np.maximum(np.asarray(u, dtype=np.float64), 0.0) + det.dark
    fast = (not det.paralyzable and det.p_ap == 0.0
            and det.tau_jitter_cv == 0.0 and det.start_mode in ("active", "delayed"))
    if fast:
        N = _counts_nonpar_vec(lam, T, det.tau, rng, det.start_mode)
    else:
        N = _counts_general_eventloop(lam, T, det, rng)
    b = N.astype(np.float64)
    if sigma_b > 0:
        b = b + rng.normal(0.0, sigma_b, size=b.shape)
    return b, N


# ----------------------------------------------------------------------
# estimator-side models (all in count units)
# ----------------------------------------------------------------------
def qmle_mean_var(lam, T, tau, sigma_b=0.0):
    """Gaussian-renewal quasi-likelihood moments (non-paralyzable renewal CLT):
    E[N] ~ lam*T/(1+lam*tau),  Var[N] ~ lam*T/(1+lam*tau)^3  (+ sigma_b^2)."""
    lam = np.maximum(lam, 1e-300)
    den = 1.0 + lam * tau
    mu = lam * T / den
    v = lam * T / den ** 3 + sigma_b ** 2 + 1e-12
    return mu, v


def _qmle_derivs(lam, T, tau):
    lam = np.maximum(lam, 1e-300)
    den = 1.0 + lam * tau
    dmu = T / den ** 2
    dv = T * (1.0 - 2.0 * lam * tau) / den ** 4
    return dmu, dv


def qmle_f_grad(x, A, b, Phi, det, T, sigma_b=0.0):
    """Normalized (per-frame) QMLE objective and gradient in x.
    lam = Phi*(A@x) + dark. Returns (f, grad)."""
    M = A.shape[0]
    lam = Phi * (A @ x) + det.dark
    mu, v = qmle_mean_var(lam, T, det.tau, sigma_b)
    dmu, dv = _qmle_derivs(lam, T, det.tau)
    r = b - mu
    f = float(np.sum(0.5 * np.log(v) + r * r / (2.0 * v))) / M
    dL_dlam = -dmu * r / v + 0.5 * dv * (1.0 / v - r * r / (v * v))
    g = Phi * (A.T @ dL_dlam) / M
    return f, g


def poisson_lin_f_grad(x, A, b, Phi, det, T, **_):
    """Ignore dead time entirely: N ~ Poisson(lam*T)."""
    M = A.shape[0]
    lam = np.maximum(Phi * (A @ x) + det.dark, 1e-300)
    mu = lam * T
    f = float(np.sum(mu - b * np.log(mu))) / M
    g = Phi * (A.T @ (T * (1.0 - b / mu))) / M
    return f, g


def poisson_satmean_f_grad(x, A, b, Phi, det, T, **_):
    """Mean-only correction: N ~ Poisson(lam*T/(1+lam*tau)) — corrects the
    compression but ignores sub-Poisson dispersion (key mechanism ablation)."""
    M = A.shape[0]
    lam = np.maximum(Phi * (A @ x) + det.dark, 1e-300)
    den = 1.0 + lam * det.tau
    mu = np.maximum(lam * T / den, 1e-300)
    f = float(np.sum(mu - b * np.log(mu))) / M
    dmu = T / den ** 2
    g = Phi * (A.T @ (dmu * (1.0 - b / mu))) / M
    return f, g


def precorrect_rates(b, T, tau):
    """Liu-2021-style non-paralyzable pre-correction: observed rate r=N/T ->
    corrected arrival-rate estimate lam_hat = r/(1 - r*tau), clipped at the
    saturation ceiling. Returns per-frame lam_hat (to be fed to a weighted
    linear solver downstream)."""
    r = np.maximum(b, 0.0) / T
    denom = np.maximum(1.0 - r * tau, 1.0 / (1.0 + 1e4))
    return r / denom


# ----------------------------------------------------------------------
# exact renewal count likelihood (reference; log-CDF-difference stable)
# ----------------------------------------------------------------------
def exact_logpmf(N, lam, T, tau, start_mode="active"):
    """log P(N = m) under exact non-paralyzable renewal counting.
    active:  S_m = (m-1)*tau + Gamma(m, lam) -> P(N>=m) = GammaCDF(T-(m-1)tau; m, lam)
    delayed: S_m = m*tau + Gamma(m, lam)     -> P(N>=m) = GammaCDF(T-m*tau; m, lam)
    """
    from scipy.stats import gamma as gd

    N = np.asarray(N)
    lam = np.asarray(lam, dtype=np.float64)
    off = 1.0 if start_mode == "active" else 0.0

    def p_ge(m):
        out = np.ones_like(lam)
        pos = m > 0
        tt = T - (m - off) * tau
        ok = pos & (tt > 0)
        out[pos & ~ok] = 0.0
        if ok.any():
            out[ok] = gd.cdf(tt[ok], a=m[ok], scale=1.0 / np.maximum(lam[ok], 1e-300))
        return out

    p = p_ge(N) - p_ge(N + 1)
    return np.log(np.maximum(p, 1e-300))


def exact_nll(x, A, N, Phi, det, T):
    lam = Phi * (A @ x) + det.dark
    return -float(np.sum(exact_logpmf(N, lam, T, det.tau, det.start_mode))) / A.shape[0]


# ----------------------------------------------------------------------
# Fisher information (theory ridge, spec §1)
# ----------------------------------------------------------------------
def fisher_gauss_renewal(lam, T, tau):
    """CLT-model Fisher information about lam per frame: (dmu)^2 / v."""
    lam = np.maximum(lam, 1e-300)
    den = 1.0 + lam * tau
    return (T / den ** 2) ** 2 * den ** 3 / (lam * T)


def fisher_exact(lam, T, tau, start_mode="active"):
    """Exact per-frame Fisher information of the count N about lam, computed
    from the exact pmf: I = sum_m p_m (d log p_m / d lam)^2 (central diff)."""
    lam = float(lam)
    m_max = int(np.floor(T / tau)) + 2
    m = np.arange(0, m_max + 1)
    eps = max(lam * 1e-4, 1e-8)

    def logp(l):
        return exact_logpmf(m, np.full_like(m, l, dtype=float), T, tau, start_mode)

    lp0 = logp(lam)
    d = (logp(lam + eps) - logp(lam - eps)) / (2 * eps)
    p = np.exp(lp0)
    ok = p > 1e-15
    return float(np.sum(p[ok] * d[ok] ** 2))
