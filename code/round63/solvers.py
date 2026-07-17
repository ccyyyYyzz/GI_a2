"""ROUND63 estimator / solver suite (spec §4 baseline arms) — built ON TOP of
physics.py (the data-term interface); physics.py is NOT modified.

Every iterative arm minimizes  f_data(x) + lam_tv * TV(x)  s.t. x >= 0 with the
SAME accelerated proximal solver (`tv_fista`), the SAME isotropic-TV prior, the
SAME x0 (clipped centered-GI recon rescaled to expected flux) and the SAME
iteration budget, so an arm-to-arm comparison isolates the DATA TERM (the
mechanism ablation of spec §4):

  POISSON-LIN   : physics.poisson_lin_f_grad      (dead time ignored)
  SAT-POISSON   : physics.poisson_satmean_f_grad  (mean compression only)
  QMLE          : physics.qmle_f_grad             (mean + sub-Poisson dispersion)
  PRECORRECT    : Liu-2021-style rate pre-correction + variance-weighted LS
  GI / DGI      : closed-form linear reference arms
  EXACT         : small-scale reference (smooth exact renewal NLL, nonneg L-BFGS-B,
                  NO TV — reference-only, see run_arm)

The data terms are per-frame normalized (physics divides the objective by M), so
`tv_fista`'s backtracking line search adapts the step to each arm's native scale
and one shared lam_tv grid (spec §4) is meaningful across arms. lam_tv itself is
chosen with the deployment-legal, truth-free held-out-predictive-NLL rule
(`select_lam_tv`): fit on the first 90% of frames, score the arm's OWN data-term
NLL on the last 10%.

float64 throughout; scipy imported lazily at the call sites that need it.
"""
import numpy as np
from dataclasses import dataclass

from physics import (qmle_f_grad, poisson_lin_f_grad, poisson_satmean_f_grad,
                     qmle_mean_var, precorrect_rates, exact_nll, Detector)


# ----------------------------------------------------------------------
# isotropic 2-D total variation: gradient / divergence (div = -grad^T)
# ----------------------------------------------------------------------
def _grad(x):
    """Forward differences with Neumann (zero) boundary. x:(m,n) -> (gx, gy)."""
    gx = np.zeros_like(x)
    gy = np.zeros_like(x)
    gx[:-1, :] = x[1:, :] - x[:-1, :]
    gy[:, :-1] = x[:, 1:] - x[:, :-1]
    return gx, gy


def _div(px, py):
    """Adjoint divergence: satisfies <grad x, p> = -<x, div p> exactly (checked
    in the unit tests). px, py:(m,n) -> (m,n)."""
    d = np.zeros_like(px)
    d[0, :] += px[0, :]
    d[1:-1, :] += px[1:-1, :] - px[:-2, :]
    d[-1, :] += -px[-2, :]
    d[:, 0] += py[:, 0]
    d[:, 1:-1] += py[:, 1:-1] - py[:, :-2]
    d[:, -1] += -py[:, -2]
    return d


def tv_value(z, side):
    """Isotropic TV of a flat image z (length side*side)."""
    gx, gy = _grad(np.asarray(z, dtype=np.float64).reshape(side, side))
    return float(np.sum(np.sqrt(gx * gx + gy * gy)))


def tv_prox(b, w, n_inner=20, nonneg=True):
    """Prox of  w*TV_iso + i_{x>=0}  at image b:
        argmin_{x (>=0)}  0.5||x - b||^2 + w * TV_iso(x).

    Beck-Teboulle (2009) fast gradient projection on the TV dual with the
    box constraint folded into the primal recovery x = P_C(b + w*div p); the
    dual gradient's Lipschitz constant is 8*w^2, hence the 1/(8w) step. The
    projection onto C is non-expansive so it does not break the rate.
    `n_inner` accelerated dual iterations (the spec's "Chambolle dual, 20 iters",
    accelerated form).
    """
    b = np.asarray(b, dtype=np.float64)
    if w <= 0:
        return np.maximum(b, 0.0) if nonneg else b.copy()

    def proj_C(z):
        return np.maximum(z, 0.0) if nonneg else z

    px = np.zeros_like(b)
    py = np.zeros_like(b)
    rx = px.copy()
    ry = py.copy()
    t = 1.0
    coef = 1.0 / (8.0 * w)
    x = proj_C(b)
    for _ in range(n_inner):
        x = proj_C(b + w * _div(rx, ry))
        gx, gy = _grad(x)
        px_n = rx + coef * gx
        py_n = ry + coef * gy
        norm = np.maximum(1.0, np.sqrt(px_n * px_n + py_n * py_n))
        px_n /= norm
        py_n /= norm
        t_n = 0.5 * (1.0 + np.sqrt(1.0 + 4.0 * t * t))
        m = (t - 1.0) / t_n
        rx = px_n + m * (px_n - px)
        ry = py_n + m * (py_n - py)
        px, py, t = px_n, py_n, t_n
    return proj_C(b + w * _div(px, py))


# ----------------------------------------------------------------------
# accelerated proximal-gradient solver (MFISTA + backtracking)
# ----------------------------------------------------------------------
def tv_fista(f_grad, x0, lam_tv, n_iter=200, L0=None, nonneg=True, side=None,
             callback=None, tol=1e-7, eta=2.0, n_inner=20):
    """Minimize  f(x) + lam_tv * TV_iso(x)  s.t. x >= 0 (if nonneg) over the
    side x side image, where `f_grad(x) -> (f, grad)` is any smooth data term
    (the physics.py terms, partially applied).

    MFISTA (monotone FISTA, Beck-Teboulle 2009) with backtracking line search on
    the smooth part: the Lipschitz estimate L only ever grows within an outer
    iteration and is reused, so the solver auto-adapts to each arm's per-frame
    normalized scale. Objective history is tracked; the loop stops on relative
    objective change < tol. Returns (x, info) with obj_hist and the gradient-
    mapping norm (the correct stationarity residual for a nonsmooth prox).
    """
    x = np.asarray(x0, dtype=np.float64).ravel().copy()
    if nonneg:
        x = np.maximum(x, 0.0)
    n = x.size
    if side is None:
        side = int(round(np.sqrt(n)))
    assert side * side == n, "x0 length %d is not a square image" % n

    def prox(v, w):
        return tv_prox(v.reshape(side, side), w, n_inner=n_inner,
                       nonneg=nonneg).ravel()

    def obj_tv(z):
        return lam_tv * tv_value(z, side)

    x_init = x.copy()
    L = float(L0) if L0 is not None else 1.0
    f_x, _ = f_grad(x)
    F_x = f_x + obj_tv(x)
    hist = [F_x]
    y = x.copy()
    t = 1.0

    for it in range(n_iter):
        f_y, g_y = f_grad(y)
        # --- backtracking line search on the smooth part
        while True:
            step = 1.0 / L
            z = prox(y - step * g_y, step * lam_tv)
            d = z - y
            f_z, _ = f_grad(z)
            Q = f_y + float(np.dot(g_y, d)) + 0.5 * L * float(np.dot(d, d))
            if f_z <= Q + 1e-12 * max(1.0, abs(f_y)) or L > 1e13:
                break
            L *= eta
        F_z = f_z + obj_tv(z)
        t_n = 0.5 * (1.0 + np.sqrt(1.0 + 4.0 * t * t))
        # --- MFISTA monotone selection: never accept an uphill step
        if F_z <= F_x:
            x_new, F_new = z, F_z
        else:
            x_new, F_new = x, F_x
        y = x_new + (t / t_n) * (z - x_new) + ((t - 1.0) / t_n) * (x_new - x)
        x, F_x, t = x_new, F_new, t_n
        hist.append(F_x)
        if callback is not None:
            callback(it, x, F_x)
        if abs(hist[-2] - hist[-1]) <= tol * max(1.0, abs(hist[-2])):
            break

    def _gmap_norm(z):
        _, gz = f_grad(z)
        return float(np.linalg.norm(
            L * (z - prox(z - (1.0 / L) * gz, (1.0 / L) * lam_tv))))

    info = {"obj_hist": hist, "n_iter": len(hist) - 1, "L": float(L),
            "final_obj": float(hist[-1]),
            "grad_map_norm": _gmap_norm(x),          # stationarity at solution
            "grad_map_norm_x0": _gmap_norm(x_init)}  # same L, for a ratio
    return x, info


# ----------------------------------------------------------------------
# closed-form / linear arms (pattern-kind aware)
# ----------------------------------------------------------------------
def gi_corr(A, b):
    """Centered ghost-imaging correlation: (1/M) A_c^T b_c. Scale-free direction
    (the metric protocol flux-matches downstream)."""
    A = np.asarray(A, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    Ac = A - A.mean(axis=0, keepdims=True)
    bc = b - b.mean(axis=0, keepdims=True)
    return Ac.T @ bc / A.shape[0]


def dgi(A, b):
    """Differential GI: normalize each bucket by its frame's total illumination
    r_i = sum_j a_ij before correlating (Ferri-2010 style). b:(M,) or (M,K)."""
    A = np.asarray(A, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    r = A.sum(axis=1)
    rbar = r.mean()
    if b.ndim == 1:
        Y = b - (b.mean() / rbar) * r
    else:
        Y = b - (b.mean(axis=0, keepdims=True) / rbar) * r[:, None]
    return A.T @ Y / A.shape[0]


def hadamard_pair_combine(A, b, meta):
    """For the "hadpair" pattern kind: difference each complementary row pair
    (a_+, a_-) into one SIGNED measurement BEFORE estimation, so the downstream
    system is linear-signed and any arm can run on it. `meta["pairs"]` is a
    (P, 2) int array of (i_plus, i_minus) frame indices.

    Returns (A_signed, b_signed) with A_signed = a_+ - a_-, b_signed = b_+ - b_-.
    For a Hadamard +/-1 row h realized as a_+ = (h+1)/2 and a_- = 1 - a_+, this
    recovers A_signed = h. Exposure/photon accounting for the two exposures stays
    with the CALLER (this only differences)."""
    pairs = np.asarray(meta["pairs"], dtype=np.int64)
    assert pairs.ndim == 2 and pairs.shape[1] == 2, "meta['pairs'] must be (P,2)"
    ip, im = pairs[:, 0], pairs[:, 1]
    A = np.asarray(A, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    A_signed = A[ip] - A[im]
    b_signed = b[ip] - b[im]
    return A_signed, b_signed


# ----------------------------------------------------------------------
# Liu-2021-style pre-correction + variance-weighted least squares
# ----------------------------------------------------------------------
def _precorrect_rates_paralyzable(b, T, tau):
    """Paralyzable-detector rate inversion: the recorded rate obeys
    r = lam*exp(-lam*tau)  (every arrival, recorded or not, re-arms the dead
    time). Invert for lam with the principal Lambert-W branch:
        r*tau = y*exp(-y),  y = lam*tau  =>  y = -W0(-r*tau),  lam = y/tau.
    r*tau is clipped to the invertible domain r <= 1/(e*tau) (the paralyzable
    throughput maximum at lam=1/tau); W0 then returns the physical low-load
    branch y in [0, 1]."""
    from scipy.special import lambertw

    r = np.maximum(np.asarray(b, dtype=np.float64), 0.0) / T
    z = np.maximum(-r * tau, -1.0 / np.e + 1e-12)   # domain clip r <= 1/(e*tau)
    y = -np.real(lambertw(z, 0))
    return y / tau


def precorrect_data_term(A, b, T, tau, Phi, dark, paralyzable=False,
                         w_floor=1e-8):
    """Build the Liu-2021-style variance-weighted-LS data term (f_grad) on
    (A, b). Returns (f_grad, aux) where aux carries lam_hat and the weights W.

      data term:  0.5 * (1/M) * sum_i W_i * (Phi*(A x)_i + dark - lam_hat_i)^2

    Weights = delta-method error propagation of lam_hat:
        Var[lam_hat] ~ Var[N] * (d lam_hat / d N)^2,
    Var[N] from physics.qmle_mean_var(lam_hat), count derivative
        non-paralyzable  lam_hat = r/(1-r*tau),  d lam_hat/dN = 1/(T (1-r*tau)^2)
        paralyzable      lam_hat = -W0(-r*tau)/tau,
                         d lam_hat/dN = lam_hat / (N (1 - lam_hat*tau))   (r=N/T)
    W_i = 1/Var[lam_hat_i] (floored), NORMALIZED to mean 1 so the residual scale
    matches the per-frame-normalized likelihood arms and one shared lam_tv grid
    applies. Saturated frames (r*tau -> 1) get huge Var -> ~0 weight, as they
    should. (The paralyzable branch reuses qmle_mean_var for Var[N] as a
    documented ablation approximation — exact paralyzable count variance is not
    modeled.)"""
    A = np.asarray(A, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    M = A.shape[0]
    r = np.maximum(b, 0.0) / T
    if paralyzable:
        lam_hat = _precorrect_rates_paralyzable(b, T, tau)
        n_safe = np.maximum(b, 1e-6)
        denom = np.maximum(1.0 - lam_hat * tau, 1e-6)
        dlam_dN = lam_hat / (n_safe * denom)
    else:
        lam_hat = precorrect_rates(b, T, tau)
        denom = np.maximum(1.0 - r * tau, 1.0 / (1.0 + 1e4))
        dlam_dN = 1.0 / (T * denom ** 2)
    _, varN = qmle_mean_var(lam_hat, T, tau)
    var_lam = np.maximum(varN * dlam_dN ** 2, 1e-300)
    W = 1.0 / var_lam
    W = np.maximum(W, w_floor * W.max())
    W = W / W.mean()

    def f_grad(x):
        pred = Phi * (A @ x) + dark
        res = pred - lam_hat
        wr = W * res
        f = 0.5 * float(np.dot(res, wr)) / M
        g = Phi * (A.T @ wr) / M
        return f, g

    return f_grad, {"lam_hat": lam_hat, "W": W}


def precorrect_wls(A, b, T, tau, Phi, dark, lam_tv, side, x0=None,
                   paralyzable=False, n_iter=200, nonneg=True, w_floor=1e-8,
                   **fista_kw):
    """Fit the pre-correction + variance-weighted-LS arm with tv_fista (same
    solver / prior as every other iterative arm)."""
    f_grad, aux = precorrect_data_term(A, b, T, tau, Phi, dark,
                                       paralyzable=paralyzable, w_floor=w_floor)
    if x0 is None:
        x0 = init_gi_flux(A, b, Phi, dark, T, tau)
    x, info = tv_fista(f_grad, x0, lam_tv, n_iter=n_iter, side=side,
                       nonneg=nonneg, **fista_kw)
    info["lam_hat_mean"] = float(np.mean(aux["lam_hat"]))
    info["weight_cv"] = float(np.std(aux["W"]) / max(np.mean(aux["W"]), 1e-30))
    return x, info


# ----------------------------------------------------------------------
# shared initializer + lam_tv selection (deployment-legal, truth-free)
# ----------------------------------------------------------------------
def init_gi_flux(A, b, Phi, dark, T, tau):
    """Shared x0 for every iterative arm: clipped centered-GI reconstruction
    rescaled so its implied signal rate Phi*mean(A x0) matches the flux estimated
    from the data (mean of the pre-corrected arrival rate, minus dark). Direction
    from GI, magnitude from the data — no truth used."""
    A = np.asarray(A, dtype=np.float64)
    x = np.maximum(gi_corr(A, b), 0.0)
    if x.sum() <= 0:
        x = np.ones(A.shape[1], dtype=np.float64)
    target = max(float(np.mean(precorrect_rates(b, T, tau))) - dark, 1e-12)
    cur = Phi * float(np.mean(A @ x))
    if cur > 0:
        x = x * (target / cur)
    return x


def select_lam_tv(f_grad_factory, A, b, x0, side, scale=None,
                  grid=(1e-4, 3e-4, 1e-3, 3e-3, 1e-2), n_iter=120, nonneg=True,
                  **fista_kw):
    """Truth-free lam_tv selection (spec §4). For each lam on `grid`*`scale`:
    fit the arm on the first 90% of frames, then score the arm's OWN data-term
    NLL on the held-out last 10% (deterministic split). Pick the minimizer.

    `f_grad_factory(A_sub, b_sub) -> f_grad` builds the arm's data term on a frame
    subset (so held-out scoring uses the identical likelihood). `scale` calibrates
    the grid; when None it defaults to the L2 norm of the arm's OWN data-term
    gradient at x0 (||grad f_data(x0)||_2), which puts the fixed multiplier grid
    into the regime where TV competes with the data misfit — arm-relative and
    problem-adaptive, no truth used. Returns (lam_best, info).
    """
    A = np.asarray(A, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    M = A.shape[0]
    n_tr = int(round(0.9 * M))
    n_tr = min(max(n_tr, 1), M - 1)
    A_tr, b_tr = A[:n_tr], b[:n_tr]
    A_val, b_val = A[n_tr:], b[n_tr:]
    fg_val = f_grad_factory(A_val, b_val)

    if scale is None:
        _, g0 = f_grad_factory(A, b)(np.asarray(x0, dtype=np.float64).ravel())
        scale = float(np.linalg.norm(g0)) + 1e-30
    lams = np.asarray(grid, dtype=np.float64) * float(scale)
    val_nll, fits = [], []
    for lam in lams:
        fg_tr = f_grad_factory(A_tr, b_tr)
        x, _ = tv_fista(fg_tr, x0, lam, n_iter=n_iter, side=side,
                        nonneg=nonneg, **fista_kw)
        v, _ = fg_val(x)
        val_nll.append(float(v))
        fits.append(x)
    k = int(np.argmin(val_nll))
    return float(lams[k]), {"grid": lams.tolist(), "val_nll": val_nll,
                            "lam_best": float(lams[k]), "x_best": fits[k],
                            "scale": float(scale)}


# ----------------------------------------------------------------------
# arm dispatcher
# ----------------------------------------------------------------------
@dataclass
class ArmContext:
    """Everything an arm needs beyond (A, b). One instance is shared across arms
    for a given (cell), guaranteeing identical x0 / budget / prior."""
    Phi: float
    det: Detector
    T: float
    side: int
    lam_tv: float = None                # None -> select_lam_tv (truth-free)
    sigma_b: float = 0.0
    n_iter: int = 200
    pattern_kind: str = "bernoulli"     # "bernoulli" | "hadpair" | "gam4" ...
    meta: dict = None                   # e.g. {"pairs": (P,2)} for hadpair
    x0: np.ndarray = None
    lam_grid_scale: float = None        # None -> ||grad f_data(x0)||_2 per arm
    select_grid: tuple = (1e-4, 3e-4, 1e-3, 3e-3, 1e-2)
    select_iter: int = 120
    exact_maxiter: int = 60


_LINEAR_ARMS = {"GI", "DGI"}
_ITER_ARMS = {"POISSON-LIN", "SAT-POISSON", "QMLE", "PRECORRECT"}
_PHYS_FG = {"QMLE": qmle_f_grad, "POISSON-LIN": poisson_lin_f_grad,
            "SAT-POISSON": poisson_satmean_f_grad}


def _phys_factory(arm_name, ctx):
    """f_grad_factory for the physics.py count-likelihood arms: partial over
    (Phi, det, T, sigma_b), leaving (A, b) to be bound per frame-subset."""
    fn = _PHYS_FG[arm_name]

    def factory(A_sub, b_sub):
        def f_grad(x):
            return fn(x, A_sub, b_sub, ctx.Phi, ctx.det, ctx.T,
                      sigma_b=ctx.sigma_b)
        return f_grad
    return factory


def _arm_factory(arm_name, ctx):
    """Unified f_grad_factory for any iterative arm (count-likelihood or the
    PRECORRECT weighted-LS), so lam_tv selection is identical across arms."""
    if arm_name == "PRECORRECT":
        def factory(A_sub, b_sub):
            fg, _ = precorrect_data_term(A_sub, b_sub, ctx.T, ctx.det.tau,
                                         ctx.Phi, ctx.det.dark,
                                         paralyzable=ctx.det.paralyzable)
            return fg
        return factory
    return _phys_factory(arm_name, ctx)


def run_arm(arm_name, A, b, ctx):
    """Dispatch one reconstruction arm. Returns (x_hat (n,), info).

    Pattern-kind awareness: for "hadpair", the linear GI/DGI arms run on the
    complementary-differenced SIGNED system (standard signed-Hadamard SPI); the
    count-likelihood arms (QMLE/POISSON-LIN/SAT-POISSON/EXACT) keep the RAW
    per-frame counts, because their renewal/Poisson models are defined on raw
    counts, not on differences of two counts (differencing would break the count
    statistics). PRECORRECT also runs on raw counts (rate inversion is per frame).
    """
    A = np.asarray(A, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    side = ctx.side

    # linear arms -------------------------------------------------------
    if arm_name in _LINEAR_ARMS:
        A_use, b_use = A, b
        if ctx.pattern_kind == "hadpair" and ctx.meta is not None:
            A_use, b_use = hadamard_pair_combine(A, b, ctx.meta)
        x = gi_corr(A_use, b_use) if arm_name == "GI" else dgi(A_use, b_use)
        return x, {"arm": arm_name, "linear": True,
                   "hadpair": ctx.pattern_kind == "hadpair"}

    # shared init -------------------------------------------------------
    x0 = ctx.x0
    if x0 is None:
        x0 = init_gi_flux(A, b, ctx.Phi, ctx.det.dark, ctx.T, ctx.det.tau)

    # EXACT reference (smooth exact NLL, nonneg L-BFGS-B, NO TV) ---------
    if arm_name == "EXACT":
        from scipy.optimize import minimize

        N = np.rint(np.maximum(b, 0.0)).astype(np.int64)
        n = A.shape[1]

        def nll(xf):
            return exact_nll(np.maximum(xf, 0.0), A, N, ctx.Phi, ctx.det, ctx.T)

        res = minimize(nll, np.maximum(x0, 0.0), method="L-BFGS-B",
                       bounds=[(0.0, None)] * n,
                       options={"maxiter": ctx.exact_maxiter,
                                "maxfun": 20 * ctx.exact_maxiter})
        return np.maximum(res.x, 0.0), {
            "arm": "EXACT", "reference_only": True, "fun": float(res.fun),
            "success": bool(res.success), "n_iter": int(res.nit),
            "note": "no-TV nonneg reference upper bound"}

    if arm_name not in _ITER_ARMS:
        raise ValueError("unknown arm %r" % (arm_name,))

    # lam_tv: truth-free held-out-predictive-NLL selection (same path all arms)
    factory = _arm_factory(arm_name, ctx)
    lam_tv = ctx.lam_tv
    sel_info = None
    if lam_tv is None:
        lam_tv, sel_info = select_lam_tv(
            factory, A, b, x0, side, scale=ctx.lam_grid_scale,
            grid=ctx.select_grid, n_iter=ctx.select_iter)

    # fit on all frames -------------------------------------------------
    if arm_name == "PRECORRECT":
        x, info = precorrect_wls(A, b, ctx.T, ctx.det.tau, ctx.Phi,
                                 ctx.det.dark, lam_tv, side, x0=x0,
                                 paralyzable=ctx.det.paralyzable,
                                 n_iter=ctx.n_iter)
    else:
        x, info = tv_fista(factory(A, b), x0, lam_tv, n_iter=ctx.n_iter,
                           side=side)

    info.update({"arm": arm_name, "lam_tv": float(lam_tv)})
    if sel_info is not None:
        info["lam_select"] = {"grid": sel_info["grid"],
                              "val_nll": sel_info["val_nll"],
                              "scale": sel_info.get("scale")}
    return x, info
