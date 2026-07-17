"""ROUND62 G0(c): MLE with the Gaussian-renewal approximate likelihood for DT30.

Observation model (harness convention: every interval includes the dead time,
b = N/s + readout):
  mean  m_b(u) = u / (1 + c u),                 c = 3/7  (= tau * s, T = 1)
  var   v_b(u) = u / (s (1 + c u)^3) + sigma_r^2,   sigma_r = 2/s
(renewal CLT: Var[N(T)] ~ sigma_int^2 T / mu_int^3 with mu_int = tau + 1/lam,
sigma_int^2 = 1/lam^2  ->  Var[N] = lam T/(1+lam tau)^3, Fano = (1+lam tau)^-2.)

NLL(x) = sum_i [ 0.5 log v_b(u_i) + (b_i - m_b(u_i))^2 / (2 v_b(u_i)) ],
u = A x, optimized by L-BFGS-B with x >= 0, multi-start (WHITEN-LW, L-ISOTRON,
MAVE), best-of by final NLL; gradient norms reported (spec §1c).

Exact-likelihood cross-check (1000-frame subsample): with all intervals
containing tau,  P(N >= m) = GammaCDF(T - m tau; shape=m, rate=lam)  (= P(S_m
<= T), S_m = m tau + Gamma(m, lam)),  P(N = m) = P(N>=m) - P(N>=m+1).  The
Gaussian approximation's likelihood surface is compared to the exact renewal
surface along the alpha * x_true ray.
"""
import numpy as np

C_DT = 3.0 / 7.0


def mean_var_b(u, s):
    u = np.maximum(u, 0.0)
    den = 1.0 + C_DT * u
    m = u / den
    v = u / (s * den ** 3) + (2.0 / s) ** 2
    return m, v


def _mean_var_derivs(u, s):
    u = np.maximum(u, 0.0)
    den = 1.0 + C_DT * u
    dm = 1.0 / den ** 2
    dv = (1.0 - 2.0 * C_DT * u) / (s * den ** 4)
    return dm, dv


def nll_and_grad(xflat, A, B, s, n, K):
    X = xflat.reshape(n, K)
    U = A @ X
    m, v = mean_var_b(U, s)
    dm, dv = _mean_var_derivs(U, s)
    r = B - m
    f = float(np.sum(0.5 * np.log(v) + r * r / (2.0 * v)))
    dL_du = -dm * r / v + 0.5 * dv * (1.0 / v - r * r / (v * v))
    g = (A.T @ dL_du).ravel()
    return f, g


def mle_renewal(A, B, s, starts, maxiter=300):
    """starts: dict name -> (n, K) init. Returns dict with per-image best-of
    reconstruction, per-start NLL and final projected-gradient inf-norms."""
    from scipy.optimize import minimize

    M, n = A.shape
    K = B.shape[1]
    runs = {}
    for name, X0 in starts.items():
        x0 = np.maximum(np.asarray(X0, dtype=np.float64), 0.0).ravel() + 1e-9
        res = minimize(nll_and_grad, x0, args=(A, B, s, n, K), jac=True,
                       method="L-BFGS-B", bounds=[(0.0, None)] * (n * K),
                       options={"maxiter": maxiter, "maxfun": 2 * maxiter})
        Xh = res.x.reshape(n, K)
        # per-image NLL for best-of selection
        U = A @ Xh
        m, v = mean_var_b(U, s)
        nll_img = np.sum(0.5 * np.log(v) + (B - m) ** 2 / (2.0 * v), axis=0)
        # convergence diagnostic = PROJECTED gradient under x >= 0 bounds
        # (raw gradient is legitimately large at zero-clamped coordinates)
        pg = np.where(res.x <= 0.0, np.minimum(res.jac, 0.0), res.jac)
        runs[name] = {"X": Xh, "nll_per_image": nll_img,
                      "grad_inf": float(np.max(np.abs(pg))),
                      "raw_grad_inf": float(np.max(np.abs(res.jac))),
                      "converged": bool(res.success), "n_iter": int(res.nit)}
    names = list(runs)
    nll_mat = np.stack([runs[nm]["nll_per_image"] for nm in names], axis=0)
    best_idx = np.argmin(nll_mat, axis=0)
    Xbest = np.stack([runs[names[best_idx[k]]]["X"][:, k] for k in range(K)], axis=1)
    diag = {nm: {"nll_per_image": runs[nm]["nll_per_image"].tolist(),
                 "grad_inf": runs[nm]["grad_inf"],
                 "converged": runs[nm]["converged"],
                 "n_iter": runs[nm]["n_iter"]} for nm in names}
    diag["best_start_per_image"] = [names[i] for i in best_idx]
    return Xbest, diag


def exact_loglik_counts(N, lam, tau, T=1.0):
    """Exact renewal log-likelihood of counts N_i at rates lam_i:
    log[ P(N>=m) - P(N>=m+1) ],  P(N>=m) = GammaCDF(T - m tau; m, lam)."""
    from scipy.stats import gamma as gamma_dist

    N = np.asarray(N)
    lam = np.asarray(lam, dtype=np.float64)

    def p_ge(m):
        out = np.ones_like(lam)
        pos = m > 0
        tt = T - m * tau
        ok = pos & (tt > 0)
        out[pos & ~ok] = 0.0
        if ok.any():
            out[ok] = gamma_dist.cdf(tt[ok], a=m[ok], scale=1.0 / lam[ok])
        return out

    p = p_ge(N) - p_ge(N + 1)
    return np.log(np.maximum(p, 1e-300))


def likelihood_surface_check(u, s, rng, alphas=None, n_frames=1000):
    """Draw exact renewal counts for n_frames energies u, then compare the
    exact and Gaussian-renewal log-lik surfaces along the alpha * u ray.
    Returns dict with normalized curves, argmaxes and max deviation."""
    import sys as _sys
    import os as _os

    _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.dirname(
        _os.path.abspath(__file__)))))
    from gi_core.links import _dt30_counts

    if alphas is None:
        alphas = np.linspace(0.7, 1.3, 61)
    u = np.asarray(u[:n_frames], dtype=np.float64)
    tau = C_DT / s
    N = _dt30_counts(np.maximum(s * u, 0.0), s, rng)
    L_ex, L_ga = [], []
    for a in alphas:
        lam = np.maximum(s * a * u, 1e-12)
        L_ex.append(float(np.sum(exact_loglik_counts(N, lam, tau))))
        # Gaussian-renewal on counts N (no readout in this check)
        mN = lam / (1.0 + lam * tau)
        vN = lam / (1.0 + lam * tau) ** 3
        L_ga.append(float(np.sum(-0.5 * np.log(vN) - (N - mN) ** 2 / (2.0 * vN))))
    L_ex = np.array(L_ex) - np.max(L_ex)
    L_ga = np.array(L_ga) - np.max(L_ga)
    return {
        "alphas": alphas.tolist(),
        "exact_norm": L_ex.tolist(),
        "gauss_norm": L_ga.tolist(),
        "argmax_exact": float(alphas[int(np.argmax(L_ex))]),
        "argmax_gauss": float(alphas[int(np.argmax(L_ga))]),
        "max_abs_dev_within_2units": float(np.max(np.abs(
            L_ex[L_ex > -2.0] - L_ga[L_ex > -2.0]))),
    }
