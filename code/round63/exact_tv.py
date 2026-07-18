"""EXACT+TV arm for the S4 image-level split (spec D2.1 §3 S4(ii)).

S4(ii) compares exact-likelihood and RQL reconstructions at 8x8 / 16x16 with
the SAME TV prior and the SAME frozen selection rule, so the gap isolates the
likelihood approximation, not the prior. The existing solvers.EXACT arm is a
no-TV L-BFGS-B reference and cannot serve that comparison; this module supplies
the missing smooth data term:

  f(x)      = -(1/M) sum_i log P(N_i = m_i | lambda_i(x)),
  grad f(x) = -(Phi/M) A^T s,   s_i = d/dlambda log P(N_i | lambda_i)

with the exact non-paralyzable active-start renewal PMF of physics.py:

  P(N >= m) = G_m = P(Pois(z_m) >= m),  z_m = lambda * t_m,
  t_m = T - (m-1) tau  (t_m <= 0 => G_m = 0),
  P(N = m) = G_m - G_{m+1},  G_0 = 1.

Analytic derivative (the identity verified against central differences in
physics.exact_fisher_analytic):

  dG_m/dlambda = t_m * e^{-z_m} z_m^{m-1}/(m-1)!  =  (m/lambda) Pois(m; z_m),
  score(m)     = (dG_m/dlambda - dG_{m+1}/dlambda) / (G_m - G_{m+1}).

Everything is computed in the log domain (poisson.logsf + logdiffexp for the
denominator; signed log-diff for the numerator) — no finite differences, no
truncation. The exact NLL is smooth but NOT convex in general; tv_fista/MFISTA
is then a monotone local method, which is the documented S4 protocol (both
arms share the same solver, budget, init, and selection rule; S4 runs at 8/16
px where the landscape is benign and multistart checks are cheap).

Deliberately standalone: it only READS physics.py/solvers.py. Registration as
arm "EXACT-TV" (solvers._PHYS_FG/_ITER_ARMS + campaign.PHYSICAL_ARMS) is a
two-line edit DEFERRED until the running S1 sweep finishes (mid-sweep edits to
imported modules would fork S1's provenance).
"""
import os
import sys

import numpy as np
from scipy.stats import poisson
from scipy.special import gammaln

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))

from physics import _log_G, _logdiffexp, exact_logpmf, LAM_FLOOR_REL


def _log_Gdot(m, lam, T, tau):
    """log |dG_m/dlambda| elementwise (the derivative is nonnegative), with
    the convention dG_0/dlambda = 0 (log -> -inf).

    dG_m/dlambda = t_m e^{-z} z^{m-1}/(m-1)!  (z = lam t_m, t_m = T-(m-1)tau)
    log          = log t_m - z + (m-1) log z - gammaln(m)
    Zero (log -inf) when m == 0 or t_m <= 0 or lam <= 0.
    """
    m = np.asarray(m, dtype=np.int64)
    lam = np.broadcast_to(np.asarray(lam, dtype=np.float64), m.shape)
    t = T - (m - 1) * tau
    out = np.full(m.shape, -np.inf, dtype=np.float64)
    ok = (m >= 1) & (t > 0) & (lam > 0)
    if ok.any():
        z = lam[ok] * t[ok]
        mm = m[ok].astype(np.float64)
        out[ok] = (np.log(t[ok]) - z + (mm - 1.0) * np.log(z)
                   - gammaln(mm))
    return out


def exact_score(N, lam, T, tau):
    """d/dlambda log P(N = m | lam) for the exact renewal PMF, elementwise,
    log-domain stable. Returns 0 where the pmf underflows to zero (the
    solver's line search then relies on neighboring measurements)."""
    N = np.asarray(N, dtype=np.int64)
    lam = np.asarray(lam, dtype=np.float64)
    la_num_hi = _log_Gdot(N, lam, T, tau)          # log dG_m
    la_num_lo = _log_Gdot(N + 1, lam, T, tau)      # log dG_{m+1}
    log_pmf = exact_logpmf(N, lam, T, tau)         # log (G_m - G_{m+1})

    score = np.zeros(N.shape, dtype=np.float64)
    ok = np.isfinite(log_pmf)
    if not ok.any():
        return score
    hi, lo, lp = la_num_hi[ok], la_num_lo[ok], log_pmf[ok]
    # numerator dG_m - dG_{m+1}: signed log-diff (either order can dominate:
    # dG is unimodal in m around z, so the sign flips across the distribution)
    num_log = _logdiffexp(np.maximum(hi, lo), np.minimum(hi, lo))
    sgn = np.where(hi >= lo, 1.0, -1.0)
    with np.errstate(over="ignore"):
        val = sgn * np.exp(num_log - lp)
    val[~np.isfinite(val)] = 0.0
    score[ok] = val
    return score


def exact_f_grad_factory(A, b, ctx):
    """(f, grad) of the exact-renewal NLL data term, per-frame normalized like
    the other physics data terms, for tv_fista. ctx supplies Phi, det (tau_est,
    dark_est), T. Signature parallels solvers._arm_factory outputs."""
    A = np.asarray(A, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    N = np.rint(b).astype(np.int64)                 # sigma_b = 0 main model
    M = A.shape[0]
    Phi, T, tau = ctx.Phi, ctx.T, ctx.det.tau
    dark = ctx.det.dark
    lam_floor = LAM_FLOOR_REL * Phi

    def f_grad(x):
        lam = np.maximum(Phi * (A @ x) + dark, lam_floor)
        lp = exact_logpmf(N, lam, T, tau)
        # a pmf underflow (impossible count under current lam) gets a large
        # finite penalty so the line search backs off instead of dying on inf
        lp = np.where(np.isfinite(lp), lp, -745.0)
        f = -float(np.sum(lp)) / M
        s = exact_score(N, lam, T, tau)
        g = -(Phi / M) * (A.T @ s)
        return f, g

    return f_grad


# ---------------------------------------------------------------- verification
def _selftest():
    """(1) score matches central-difference d log pmf; (2) factory gradient
    matches finite differences on a random 8x8 problem; (3) MLE sanity: the
    scalar exact-NLL argmin sits at the true lambda within noise."""
    from gi_core.utils import rng_for
    from physics import Detector, simulate_counts
    from solvers import ArmContext

    rng = rng_for(0, 63, 77)
    T, tau = 500.0, 1.0

    # (1) elementwise score vs central difference
    lam = np.array([0.05, 0.3, 0.6, 1.2, 3.0, 8.0])
    Ngrid = np.array([10, 120, 180, 260, 370, 440])
    h = 1e-6
    s_an = exact_score(Ngrid, lam, T, tau)
    s_fd = (exact_logpmf(Ngrid, lam + h, T, tau)
            - exact_logpmf(Ngrid, lam - h, T, tau)) / (2 * h)
    e1 = float(np.max(np.abs(s_an - s_fd) / np.maximum(np.abs(s_fd), 1e-12)))
    print("[exact_tv] (1) score vs central-diff  max rel err = %.3g" % e1,
          flush=True)
    assert e1 < 1e-5, "score mismatch"

    # (2) factory gradient on an 8x8 inverse problem
    side, n, M = 8, 64, 96
    x_true = rng.random(n)
    x_true /= x_true.sum()
    A = 2.0 * (rng.random((M, n)) < 0.5)
    det = Detector(tau=tau, dark=0.0, start_mode="active")
    u = A @ x_true
    Phi = 0.6 / (tau * float(u.mean()))
    b, _ = simulate_counts(u, Phi, T, det, rng)
    ctx = ArmContext(Phi=Phi, det=det, T=T, side=side)
    fg = exact_f_grad_factory(A, b, ctx)
    x0 = np.full(n, 1.0 / n)
    f0, g = fg(x0)
    d = rng.standard_normal(n)
    d /= np.linalg.norm(d)
    eps = 1e-7
    fp, _ = fg(x0 + eps * d)
    fm, _ = fg(x0 - eps * d)
    dd_fd = (fp - fm) / (2 * eps)
    dd_an = float(np.dot(g, d))
    e2 = abs(dd_an - dd_fd) / max(abs(dd_fd), 1e-12)
    print("[exact_tv] (2) directional grad: analytic %.8g vs fd %.8g "
          "(rel %.3g)" % (dd_an, dd_fd, e2), flush=True)
    assert e2 < 1e-5, "factory gradient mismatch"

    # (3) scalar MLE sanity at rho = 0.6
    lam_true = 0.6
    Ns, _ = simulate_counts(np.ones(4000), lam_true, T, det, rng)
    grid = np.linspace(0.4, 0.8, 401)
    nll = [-float(np.sum(exact_logpmf(np.rint(Ns).astype(np.int64),
                                      np.full(Ns.shape, g_), T, tau)))
           for g_ in grid]
    lam_hat = float(grid[int(np.argmin(nll))])
    print("[exact_tv] (3) scalar exact-MLE: lam_hat %.4f vs true %.4f"
          % (lam_hat, lam_true), flush=True)
    assert abs(lam_hat - lam_true) < 0.02, "scalar MLE off"

    # (4) end-to-end: EXACT+TV via tv_fista beats/matches its init on 8x8
    from solvers import tv_fista
    x_hat, info = tv_fista(fg, x0, 1e-3 * float(np.linalg.norm(g)),
                           n_iter=150, side=side)
    finite = bool(np.all(np.isfinite(x_hat)))
    drop = info["obj_hist"][0] - info["final_obj"]
    print("[exact_tv] (4) tv_fista run: finite=%s obj drop=%.4g iters=%d"
          % (finite, drop, info["n_iter"]), flush=True)
    assert finite and drop > 0
    print("[exact_tv] SELFTEST PASS", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(_selftest())
