"""RNG discipline, linear-algebra helpers, PAV wrapper. float64 throughout."""
import numpy as np
from numpy.random import default_rng

from . import config as C


def rng_for(seed, *stream):
    """All randomness hangs off np.random.default_rng(SEED0+seed) via SeedSequence
    spawn keys — spec §0.8 ('any new random source must live under this system')."""
    return default_rng([C.SEED0 + int(seed)] + [int(t) for t in stream])


def grid_kernel_cov(side, corr_len, marginal_var, jitter_rel=C.KERNEL_JITTER_REL):
    """Gaussian-kernel spatial covariance on a side x side grid.

    C_ij = marginal_var * exp(-||r_i - r_j||^2 / (2 corr_len^2)) + jitter.
    The jitter is part of the model definition (used consistently for sampling,
    oracle score and WHITEN-OR), needed because the Gaussian kernel is
    numerically rank-deficient.
    """
    ii, jj = np.meshgrid(np.arange(side), np.arange(side), indexing="ij")
    coords = np.stack([ii.ravel(), jj.ravel()], axis=1).astype(np.float64)
    d2 = ((coords[:, None, :] - coords[None, :, :]) ** 2).sum(-1)
    K = marginal_var * np.exp(-d2 / (2.0 * corr_len ** 2))
    K[np.diag_indices_from(K)] += jitter_rel * marginal_var
    return K


class CholOp:
    """Cached Cholesky factorization: solve, half-solve, sample transform, logdet."""

    def __init__(self, S):
        from scipy.linalg import cholesky

        self.L = cholesky(S, lower=True)  # S = L L^T
        self.n = S.shape[0]

    def solve(self, B):
        from scipy.linalg import solve_triangular

        Y = solve_triangular(self.L, B, lower=True)
        return solve_triangular(self.L.T, Y, lower=False)

    def half_solve(self, B):
        """L^{-1} B (whitening direction)."""
        from scipy.linalg import solve_triangular

        return solve_triangular(self.L, B, lower=True)

    def half_solve_T(self, B):
        """L^{-T} B (back-transform of whitened directions)."""
        from scipy.linalg import solve_triangular

        return solve_triangular(self.L.T, B, lower=False)

    def sample_mult(self, Z):
        """Z @ L^T : maps iid standard normal rows to covariance-S rows."""
        return Z @ self.L.T

    def logdet(self):
        return 2.0 * np.sum(np.log(np.diag(self.L)))


def pav_fit_predict(u_train, y_train, u_eval):
    """Isotonic (PAV) regression fit on (u_train, y_train), evaluated at u_eval
    with boundary clipping. sklearn's C implementation."""
    from sklearn.isotonic import IsotonicRegression

    iso = IsotonicRegression(increasing=True, out_of_bounds="clip")
    iso.fit(u_train, y_train)
    return iso.predict(u_eval)


def psnr(x, ref, data_range):
    mse = np.mean((x - ref) ** 2)
    if mse <= 0:
        return np.inf
    return 10.0 * np.log10(data_range ** 2 / mse)


def pearson(a, b):
    a = a.ravel().astype(np.float64)
    b = b.ravel().astype(np.float64)
    sa = a.std()
    sb = b.std()
    if sa == 0 or sb == 0:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])
