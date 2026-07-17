"""Bucket links and noise — spec §3.

Baseline noise: b = Poisson(s * phi(u)) / s + N(0, sigma_r^2), sigma_r = 2/s.
DT30 is the flagship physical arm: exact non-paralyzable dead-time renewal
simulation (NOT Poisson(sat(u)) — audit revision 12).
"""
import numpy as np

from . import config as C

ALPHA_SAT30 = 3.0 / 7.0          # phi(1)/1 = 1/(1+alpha) = 0.7
ALPHA_SAT50 = 1.0                # 0.5
ALPHA_LOG = np.expm1(0.7)        # log(1+alpha*1) = 0.7


def phi_mean(link, u):
    """Deterministic mean response phi(u) (the 'true phi' for MLE-OR / UT4)."""
    u = np.asarray(u, dtype=np.float64)
    if link == "LIN" or link == "FGAIN":
        return u
    if link == "DT30":
        return u / (1.0 + C.DT30_TAU_COEF * u)
    if link == "SAT30":
        return u / (1.0 + ALPHA_SAT30 * u)
    if link == "SAT50":
        return u / (1.0 + ALPHA_SAT50 * u)
    if link == "GAMMA07":
        return np.power(np.maximum(u, 0.0), 0.7)
    if link == "LOG":
        return np.log1p(ALPHA_LOG * u)
    raise ValueError(link)


def phi_mean_deriv(link, u):
    u = np.asarray(u, dtype=np.float64)
    if link == "LIN" or link == "FGAIN":
        return np.ones_like(u)
    if link == "DT30":
        return 1.0 / (1.0 + C.DT30_TAU_COEF * u) ** 2
    if link == "SAT30":
        return 1.0 / (1.0 + ALPHA_SAT30 * u) ** 2
    if link == "SAT50":
        return 1.0 / (1.0 + ALPHA_SAT50 * u) ** 2
    if link == "GAMMA07":
        return 0.7 * np.power(np.maximum(u, 1e-300), -0.3)
    if link == "LOG":
        return ALPHA_LOG / (1.0 + ALPHA_LOG * u)
    raise ValueError(link)


def _dt30_counts(lam, s, rng, chunk=2048):
    """Exact renewal counts for non-paralyzable dead time.

    Inter-record intervals iid = tau + Exp(lam), tau = (3/7)/s, window T=1.
    N = #{records with cumulative arrival time <= T}. Chunked over frames
    (primer §6: memory control). lam: 1-D array of rates (= s*u), one per frame.
    """
    lam = np.asarray(lam, dtype=np.float64)
    tau = C.DT30_TAU_COEF / s
    T = 1.0
    out = np.empty(lam.shape[0], dtype=np.int64)
    hard_cap = int(np.floor(T / tau)) + 2  # dead time bounds N by T/tau + 1
    for lo in range(0, lam.shape[0], chunk):
        lam_c = lam[lo:lo + chunk]
        mc = lam_c.shape[0]
        counts = np.zeros(mc, dtype=np.int64)
        t_acc = np.zeros(mc)
        active = np.ones(mc, dtype=bool)
        # expected count ~ lam*T/(1+lam*tau); allocate mean + 8 sigma, continue if needed
        while active.any():
            idx = np.nonzero(active)[0]
            lam_a = lam_c[idx]
            mean_n = lam_a * T / (1.0 + lam_a * tau) + 1.0
            k_alloc = int(np.clip(np.max(mean_n + 8.0 * np.sqrt(mean_n)), 16, hard_cap))
            with np.errstate(divide="ignore"):
                scale = np.where(lam_a > 0, 1.0 / np.maximum(lam_a, 1e-300), np.inf)
            E = rng.exponential(scale=scale[:, None], size=(idx.size, k_alloc))
            times = t_acc[idx, None] + np.cumsum(tau + E, axis=1)
            add = (times <= T).sum(axis=1)
            counts[idx] += add
            # frames whose allocated block did not cross T must continue
            still = times[:, -1] <= T
            t_acc[idx] = times[:, -1]
            active[:] = False
            active[idx[still]] = True
            if k_alloc >= hard_cap:
                break
        out[lo:lo + chunk] = counts
    return out


def simulate_bucket(link, u, s, rng):
    """b per frame for one image column. u: (M,) frame energies (>=0 up to
    numerical noise for physical families). Returns b (M,)."""
    u = np.asarray(u, dtype=np.float64)
    s = float(s)
    sigma_r = C.READOUT_COEF / s
    if link == "DT30":
        lam = np.maximum(s * u, 0.0)
        N = _dt30_counts(lam, s, rng)
        b = N / s
    elif link == "FGAIN":
        sg2 = np.log1p(C.FGAIN_CV ** 2)
        g = rng.lognormal(mean=-sg2 / 2.0, sigma=np.sqrt(sg2), size=u.shape)
        lam = np.maximum(s * g * u, 0.0)
        b = rng.poisson(lam) / s
    else:
        lam = np.maximum(s * phi_mean(link, u), 0.0)
        b = rng.poisson(lam) / s
    return b + rng.normal(0.0, sigma_r, size=u.shape)
