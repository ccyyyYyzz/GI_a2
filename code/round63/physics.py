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
    carry = 0.0       # detector not-ready time carried across frames (continuous)
    carry_ap = []     # afterpulse arrival times carried into the next frame
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
        if carry_ap:
            import bisect

            pending, carry_ap = carry_ap, []
            for t_ap in pending:
                if 0.0 <= t_ap < T:
                    bisect.insort(arrivals, t_ap)
                elif t_ap >= T:
                    carry_ap.append(t_ap - T - det.guard)
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
                    elif det.start_mode == "continuous":
                        # afterpulse beyond the frame is carried into the next
                        # frame's event queue (GPT round-3 fix)
                        carry_ap.append(t_ap - T - det.guard)
        if det.start_mode == "continuous":
            # detector not-ready time remaining at the START of the next frame:
            # the guard interval g elapses between frames (round-3 sign fix)
            carry = max(0.0, ready - T - det.guard)
        else:
            carry = 0.0
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


LAM_FLOOR_REL = 1e-6  # model floor: lam >= 1e-6*Phi (documented; keeps the
                      # boundary x->0 numerically graceful without changing
                      # any physically reachable operating point)


def qmle_f_grad(x, A, b, Phi, det, T, sigma_b=0.0):
    """Normalized (per-frame) QMLE objective and gradient in x.
    lam = Phi*(A@x) + dark. Returns (f, grad)."""
    M = A.shape[0]
    lam = np.maximum(Phi * (A @ x) + det.dark, LAM_FLOOR_REL * Phi)
    mu, v = qmle_mean_var(lam, T, det.tau, sigma_b)
    dmu, dv = _qmle_derivs(lam, T, det.tau)
    r = b - mu
    f = float(np.sum(0.5 * np.log(v) + r * r / (2.0 * v))) / M
    dL_dlam = -dmu * r / v + 0.5 * dv * (1.0 / v - r * r / (v * v))
    g = Phi * (A.T @ dL_dlam) / M
    return f, g


def qmle_weights(lam, T, tau, sigma_b=0.0):
    """Quasi-likelihood weights w = 1/Var[N] at rates lam (frozen per IRLS round)."""
    _, v = qmle_mean_var(lam, T, tau, sigma_b)
    return 1.0 / v


def qmle_wls_f_grad(x, A, b, Phi, det, T, w, sigma_b=0.0):
    """Wedderburn quasi-score form of the renewal QMLE: weighted least squares
    with FROZEN weights w (updated between IRLS rounds), NO log-det term.

    Rationale (implementation note for the paper): the full Gaussian NLL's
    0.5*log v(lam) term has d(log v)/d(lam) proportional to (1-2*lam*tau) —
    its sign FLIPS at rho = lam*tau = 1/2, so at high load it rewards
    inflating lam to shrink the modeled variance, a perverse incentive that
    corrupts reconstructions in information-poor directions. The proper
    quasi-score U = dmu^T V^{-1} (b - mu) has no such term; variance enters
    only as weights. Dispersion information is retained through the weights'
    (1+rho)^3/lam frame weighting."""
    M = A.shape[0]
    lam = np.maximum(Phi * (A @ x) + det.dark, LAM_FLOOR_REL * Phi)
    den = 1.0 + lam * det.tau
    mu = lam * T / den
    dmu = T / den ** 2
    r = b - mu
    f = float(np.sum(w * r * r)) / (2.0 * M)
    g = -Phi * (A.T @ (w * dmu * r)) / M
    return f, g


def rql_f_grad(x, A, b, Phi, det, T, **_):
    """RQL — renewal quasi-likelihood data fidelity (GPT round-3 final form).

    Integrated Wedderburn quasi-score for the non-paralyzable renewal count
    model (sigma_b = 0 main model): the quasi-score U_lam = N/lam - (T - N*tau)
    integrates to the CONVEX per-frame objective
        Q(lam; N) = (T - N*tau)*lam - N*log(lam),
    convex in lam (linear minus log) and hence in x through the affine map
    lam = Phi*(A@x) + dark. No IRLS, no log-det, no frozen weights.
    Identifiability boundary: frames with N*tau >= T lose the linear term
    (true saturation boundary) — they are counted and reported upstream;
    the objective keeps them via the -N log lam term only."""
    M = A.shape[0]
    lam = np.maximum(Phi * (A @ x) + det.dark, LAM_FLOOR_REL * Phi)
    w_lin = np.maximum(T - b * det.tau, 0.0)
    f = float(np.sum(w_lin * lam - b * np.log(lam))) / M
    g = Phi * (A.T @ (w_lin - b / lam)) / M
    return f, g


def poisson_lin_f_grad(x, A, b, Phi, det, T, **_):
    """Ignore dead time entirely: N ~ Poisson(lam*T)."""
    M = A.shape[0]
    lam = np.maximum(Phi * (A @ x) + det.dark, LAM_FLOOR_REL * Phi)
    mu = lam * T
    f = float(np.sum(mu - b * np.log(mu))) / M
    g = Phi * (A.T @ (T * (1.0 - b / mu))) / M
    return f, g


def poisson_satmean_f_grad(x, A, b, Phi, det, T, **_):
    """Mean-only correction: N ~ Poisson(lam*T/(1+lam*tau)) — corrects the
    compression but ignores sub-Poisson dispersion (key mechanism ablation)."""
    M = A.shape[0]
    lam = np.maximum(Phi * (A @ x) + det.dark, LAM_FLOOR_REL * Phi)
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
def _log_G(m, lam, T, tau, off):
    """log P(N >= m) via the Poisson representation, numerically stable:
    P(Gamma(m, lam) <= t) = P(Pois(lam*t) >= m) = poisson.sf(m-1, lam*t).
    m: int array (>=1 handled; m<=0 -> log 1 = 0); t = T - (m - off)*tau."""
    from scipy.stats import poisson

    m = np.asarray(m)
    lam = np.asarray(lam, dtype=np.float64)
    out = np.full(np.broadcast(m, lam).shape, -np.inf)
    m_b, lam_b = np.broadcast_arrays(m, lam)
    zero = m_b <= 0
    out[zero] = 0.0
    tt = T - (m_b - off) * tau
    ok = (~zero) & (tt > 0) & (lam_b > 0)
    if ok.any():
        out[ok] = poisson.logsf(m_b[ok] - 1, lam_b[ok] * tt[ok])
    return out


def _logdiffexp(la, lb):
    """log(exp(la) - exp(lb)) for la >= lb, elementwise, stable."""
    la = np.asarray(la, dtype=np.float64)
    lb = np.asarray(lb, dtype=np.float64)
    with np.errstate(invalid="ignore"):
        d = lb - la
        out = la + np.log1p(-np.exp(np.minimum(d, -1e-300)))
    out[~np.isfinite(la)] = -np.inf
    out[d >= 0] = -np.inf
    return out


def exact_logpmf(N, lam, T, tau, start_mode="active"):
    """log P(N = m) under exact non-paralyzable renewal counting.
    active:  S_m = (m-1)*tau + Gamma(m, lam) -> P(N>=m) = GammaCDF(T-(m-1)tau; m, lam)
    delayed: S_m = m*tau + Gamma(m, lam)     -> P(N>=m) = GammaCDF(T-m*tau; m, lam)
    Publication-grade path (GPT round-3): Poisson logsf representation +
    logdiffexp — no linear-domain CDF subtraction, no 1e-300 clamps."""
    N = np.asarray(N)
    lam = np.asarray(lam, dtype=np.float64)
    off = 1.0 if start_mode == "active" else 0.0
    lG = _log_G(N, lam, T, tau, off)
    lG1 = _log_G(N + 1, lam, T, tau, off)
    return _logdiffexp(lG, lG1)


def exact_fisher_analytic(lam, T, tau, start_mode="active"):
    """Exact count-only Fisher information about theta = log(lam) per frame,
    using the ANALYTIC derivative (GPT round-3, no finite differences):
      G_m = P(Pois(z_m) >= m),  z_m = lam*(T - (m-1)tau)   [active start]
      dG_m/dtheta = e^{-z_m} z_m^m / (m-1)!  = m * PoisPMF(m; z_m)
      I_theta = sum_m (Gdot_m - Gdot_{m+1})^2 / p_m .
    Returns information about log-lam (scale-free)."""
    from scipy.stats import poisson

    off = 1.0 if start_mode == "active" else 0.0
    m_max = int(np.floor(T / tau)) + 2
    m = np.arange(0, m_max + 2)
    tt = T - (m - off) * tau
    z = lam * np.maximum(tt, 0.0)

    def gdot(mm, zz):
        out = np.zeros_like(zz)
        ok = (mm >= 1) & (zz > 0)
        if ok.any():
            out[ok] = np.exp(poisson.logpmf(mm[ok], zz[ok])
                             + np.log(mm[ok].astype(np.float64)))
        return out

    Gd = gdot(m, z)
    lG = _log_G(m, np.full_like(m, lam, dtype=float), T, tau, off)
    lG1 = _log_G(m + 1, np.full_like(m, lam, dtype=float), T, tau, off)
    lp = _logdiffexp(lG, lG1)
    num = (Gd - np.roll(Gd, -1))[:-1]
    lp = lp[:-1]
    ok = np.isfinite(lp) & (np.abs(num) > 0)
    with np.errstate(divide="ignore"):
        log_terms = 2.0 * np.log(np.abs(num[ok])) - lp[ok]
    # log-domain summation avoids 0*inf when p_m underflows while num^2 -> 0
    keep = log_terms > -745.0
    if not keep.any():
        return 0.0
    lmax = float(np.max(log_terms[keep]))
    return float(np.exp(lmax) * np.sum(np.exp(log_terms[keep] - lmax)))


def exact_nll(x, A, N, Phi, det, T):
    lam = np.maximum(Phi * (A @ x) + det.dark, LAM_FLOOR_REL * Phi)
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
