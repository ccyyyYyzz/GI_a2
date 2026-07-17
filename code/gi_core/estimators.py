"""Estimators — spec §4. All operate on the same frozen (A, S) per (family, seed)
and the same bucket matrix B (M x K images), same post-processing (metrics.py).

Reconstructions are returned as (n, K) float64. Direction-only methods (SIR,
MAVE) return signed directions; the metric protocol (nonneg truncation + flux
matching) is scale-invariant, so this is fair across methods.
"""
import numpy as np

from . import config as C
from .utils import CholOp, pav_fit_predict


class RunContext:
    """Frozen per-(family, seed) state: patterns, scores, covariance operators.

    Everything derived from A is computed once and reused across all links,
    photon levels and images (spec §0.8).
    """

    def __init__(self, family, seed, M=C.M, est_rng=None, with_rankg=True):
        self.family = family
        self.seed = seed
        fam_id = C.FAMILY_IDS[family.name]
        from .utils import rng_for

        pat_rng = rng_for(seed, fam_id, C.STREAM_PATTERN)
        self.est_rng = est_rng or rng_for(seed, fam_id, C.STREAM_ESTIMATOR)
        A = family.sample(M, pat_rng)
        S, n_bad = family.score(A)
        valid = np.isfinite(S).all(axis=1) & np.isfinite(A).all(axis=1)
        self.n_dropped = int(M - valid.sum())
        if self.n_dropped:
            A, S = A[valid], S[valid]
            if getattr(family, "last_states", None) is not None:
                family.last_states = family.last_states[valid]
        self.A = A
        self.S = S
        self.M = A.shape[0]
        self.n = A.shape[1]
        self.abar = A.mean(axis=0)
        self.r = A.sum(axis=1)
        self.rbar = self.r.mean()

        from sklearn.covariance import LedoitWolf

        lw = LedoitWolf(assume_centered=False).fit(A)
        self.lw_cov = lw.covariance_
        self.lw_shrinkage = float(lw.shrinkage_)
        self.chol_lw = CholOp(self.lw_cov)
        self.tr_lw = float(np.trace(self.lw_cov))
        self.true_cov_op = family.true_cov_op()

        # frozen 10% held-out frame split for L-Isotron (spec §4.8)
        perm = self.est_rng.permutation(self.M)
        n_val = self.M // 10
        self.val_idx = perm[:n_val]
        self.tr_idx = perm[n_val:]

        self.G_rank = self._rank_gaussianize(A) if with_rankg else None
        self._cluster_cache = None

    def _rank_gaussianize(self, A):
        from scipy.special import ndtri

        order = np.argsort(A, axis=0)
        ranks = np.empty_like(order)
        m = A.shape[0]
        ranks[order, np.arange(A.shape[1])[None, :]] = np.arange(m)[:, None]
        return ndtri((ranks + 0.5) / m)


def _center_cols(B):
    return B - B.mean(axis=0, keepdims=True)


def gi(ctx, B):
    Bc = _center_cols(B)
    return ctx.A.T @ Bc / ctx.M


def dgi(ctx, B):
    bbar = B.mean(axis=0, keepdims=True)
    Y = B - (bbar / ctx.rbar) * ctx.r[:, None]
    # sum_i Y_i = 0 per column, so centering A is a no-op
    return ctx.A.T @ Y / ctx.M


def corr_gi(ctx, B):
    M, K = B.shape
    q = M // 4
    W = np.zeros((M, K))
    for k in range(K):
        order = np.argsort(B[:, k], kind="stable")
        W[order[-q:], k] = 1.0 / q
        W[order[:q], k] = -1.0 / q
    return ctx.A.T @ W


def whiten_lw(ctx, B):
    return ctx.chol_lw.solve(gi(ctx, B))


def whiten_or(ctx, B):
    return ctx.true_cov_op.solve(gi(ctx, B))


def score_or(ctx, B, centered=True):
    """SCORE-OR with leave-one-out centering: b_i - bbar_{-i} = M/(M-1) (b_i - bbar),
    i.e. exactly proportional to plain centering (ROUND59 header note; the
    equivalence is documented in the report, not a bias fix)."""
    if centered:
        Bc = _center_cols(B)
        return -(ctx.S.T @ Bc) / (ctx.M - 1)
    return -(ctx.S.T @ B) / ctx.M


def rankg(ctx, B):
    Bc = _center_cols(B)
    return ctx.G_rank.T @ Bc / ctx.M


def _slice_means(A, b, H):
    """Equal-count slice means of A ordered by b. Returns (n, H), weights (H,)."""
    M = b.shape[0]
    order = np.argsort(b, kind="stable")
    bounds = np.linspace(0, M, H + 1).astype(int)
    Msl = np.zeros((A.shape[1], H))
    w = np.zeros(H)
    for h in range(H):
        idx = order[bounds[h]:bounds[h + 1]]
        Msl[:, h] = A[idx].mean(axis=0)
        w[h] = idx.size / M
    return Msl, w


def _sign_fix(ctx, x, b):
    u = ctx.A @ x
    c = np.dot(u - u.mean(), b - b.mean())
    return x if c >= 0 else -x


def sir(ctx, B, H):
    """SIR with LW-standardized patterns (spec §4.4): top generalized eigenvector
    of (sum_h w_h m_h m_h^T, Sigma_LW), computed via Cholesky whitening of the
    slice-mean matrix — algebraically identical to standardize-then-PCA."""
    n, K = ctx.n, B.shape[1]
    out = np.zeros((n, K))
    for k in range(K):
        Msl, w = _slice_means(ctx.A, B[:, k], H)
        Bmat = (Msl - ctx.abar[:, None]) * np.sqrt(w)[None, :]
        Cw = ctx.chol_lw.half_solve(Bmat)          # L^{-1} B
        U, sv, _ = np.linalg.svd(Cw, full_matrices=False)
        beta = ctx.chol_lw.half_solve_T(U[:, 0])   # L^{-T} eta
        out[:, k] = _sign_fix(ctx, beta, B[:, k])
    return out


def l_isotron(ctx, B, max_iter=200, tol=1e-6):
    """L-Isotron, frozen protocol (spec §4.8): eta = 1/tr(Sigma_LW), PAV link fit,
    10% held-out frames, 3 inits (WHITEN-LW, SIR-10, random) chosen by held-out
    likelihood (Gaussian, i.e. held-out MSE). Truth never used."""
    A_tr, A_val = ctx.A[ctx.tr_idx], ctx.A[ctx.val_idx]
    B_tr, B_val = B[ctx.tr_idx], B[ctx.val_idx]
    eta = 1.0 / ctx.tr_lw
    K = B.shape[1]

    x_lw = whiten_lw(ctx, B)
    x_sir = sir(ctx, B, 10)
    x_rand = ctx.est_rng.standard_normal((ctx.n, K))
    x_rand /= np.linalg.norm(x_rand, axis=0, keepdims=True)
    # scale-match non-LW inits to the LW norm so eta is comparable
    for X0 in (x_sir, x_rand):
        X0 *= np.linalg.norm(x_lw, axis=0, keepdims=True) / np.maximum(
            np.linalg.norm(X0, axis=0, keepdims=True), 1e-30)

    best_val = np.full(K, np.inf)
    best_X = np.zeros((ctx.n, K))
    Mtr = A_tr.shape[0]
    for X0 in (x_lw, x_sir, x_rand):
        X = X0.copy()
        prev_val = np.full(K, np.inf)
        active = np.ones(K, dtype=bool)
        run_val = np.full(K, np.inf)
        run_best_X = X.copy()
        run_best_val = np.full(K, np.inf)
        for it in range(max_iter):
            if not active.any():
                break
            U_tr = A_tr @ X
            U_val = A_val @ X
            R = np.zeros((Mtr, K))
            for k in np.nonzero(active)[0]:
                phi_tr = pav_fit_predict(U_tr[:, k], B_tr[:, k], U_tr[:, k])
                phi_val = pav_fit_predict(U_tr[:, k], B_tr[:, k], U_val[:, k])
                R[:, k] = B_tr[:, k] - phi_tr
                run_val[k] = float(np.mean((B_val[:, k] - phi_val) ** 2))
                if run_val[k] < run_best_val[k]:
                    run_best_val[k] = run_val[k]
                    run_best_X[:, k] = X[:, k]
                if abs(prev_val[k] - run_val[k]) <= tol * max(1.0, run_val[k]):
                    active[k] = False
                prev_val[k] = run_val[k]
            X = X + eta * (A_tr.T @ R) / Mtr
        for k in range(K):
            if run_best_val[k] < best_val[k]:
                best_val[k] = run_best_val[k]
                best_X[:, k] = run_best_X[:, k]
    return best_X


def mle_or(ctx, B, link, s, maxiter=150):
    """MLE with known true phi (info-upper-bound diagnostic; not in any gate's
    control set). Poisson NLL on clipped counts c = max(s*b, 0), L-BFGS-B with
    x >= 0 bounds, WHITEN-LW init. Images solved jointly (objective separable)."""
    from scipy.optimize import minimize

    from .links import phi_mean, phi_mean_deriv

    A = ctx.A
    M, n = A.shape
    K = B.shape[1]
    Ccounts = np.maximum(s * B, 0.0)
    X0 = np.maximum(whiten_lw(ctx, B), 0.0) + 1e-9

    def nll_grad(xflat):
        X = xflat.reshape(n, K)
        U = A @ X
        lam = np.maximum(s * phi_mean(link, U), 1e-300)
        f = float(np.sum(lam - Ccounts * np.log(lam)))
        dphi = phi_mean_deriv(link, U)
        G = s * dphi * (1.0 - Ccounts / lam)
        g = (A.T @ G).ravel()
        return f, g

    res = minimize(nll_grad, X0.ravel(), jac=True, method="L-BFGS-B",
                   bounds=[(0.0, None)] * (n * K),
                   options={"maxiter": maxiter, "maxfun": 2 * maxiter})
    return res.x.reshape(n, K)


def cluster_whiten(ctx, B, true_states=None):
    """MIX-LOGN classical fair control (spec §4.11): unsupervised k-means(2) on
    log-patterns -> per-cluster LW-whitened GI -> frame-weighted merge.
    Returns (Xhat, clustering_accuracy or None)."""
    if ctx._cluster_cache is None:
        from sklearn.cluster import KMeans

        feats = np.log(ctx.A)
        km = KMeans(n_clusters=2, n_init=4,
                    random_state=int(ctx.est_rng.integers(2 ** 31)))
        labels = km.fit_predict(feats)
        from sklearn.covariance import LedoitWolf

        ops = []
        for c in (0, 1):
            idx = np.nonzero(labels == c)[0]
            lw = LedoitWolf(assume_centered=False).fit(ctx.A[idx])
            ops.append((idx, CholOp(lw.covariance_)))
        acc = None
        if true_states is not None:
            a = float(np.mean(labels == true_states))
            acc = max(a, 1.0 - a)
        ctx._cluster_cache = (labels, ops, acc)
    labels, ops, acc = ctx._cluster_cache
    Xhat = np.zeros((ctx.n, B.shape[1]))
    for idx, op in ops:
        Ac = ctx.A[idx]
        Bc = B[idx] - B[idx].mean(axis=0, keepdims=True)
        xc = op.solve(Ac.T @ Bc / idx.size)
        Xhat += (idx.size / ctx.M) * xc
    return Xhat, acc
