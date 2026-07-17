"""Illumination families: sampling, oracle scores, true covariances — spec §2.

All families have unit mean per pixel (U0 = E[<a,x>] = 1 for sum-normalized x).
float64 throughout. NO clamping of a anywhere (spec §2 GAM row: forbidden);
frames with nonfinite score are dropped and counted (expected count: 0).
"""
import numpy as np

from . import config as C
from .utils import CholOp, grid_kernel_cov


class Family:
    name = None

    def sample(self, m, rng):
        raise NotImplementedError

    def score(self, A):
        """Oracle score s(a) per frame, same shape as A. May return (S, n_dropped)."""
        raise NotImplementedError

    def true_cov_op(self):
        """Operator with .solve(B) for the true covariance (WHITEN-OR)."""
        raise NotImplementedError

    def mean(self):
        return 1.0  # all families are unit-mean per pixel


class DiagCovOp:
    def __init__(self, var):
        self.var = float(var)

    def solve(self, B):
        return B / self.var


class GaussFamily(Family):
    """a = 1 + L z, Sigma = Gaussian-kernel correlation (corr len 2 px, sigma=0.25).
    Score: -Sigma^{-1}(a - mu). Non-physical control arm (negative pixels possible)."""

    def __init__(self, n):
        self.name = "GAUSS"
        self.n = n
        side = int(round(np.sqrt(n)))
        assert side * side == n
        self.Sigma = grid_kernel_cov(side, C.CORR_LEN_PX, C.GAUSS_SIGMA ** 2)
        self.chol = CholOp(self.Sigma)

    def sample(self, m, rng):
        Z = rng.standard_normal((m, self.n))
        return 1.0 + self.chol.sample_mult(Z)

    def score(self, A):
        S = -self.chol.solve((A - 1.0).T).T
        return S, 0

    def true_cov_op(self):
        return self.chol

    def neg_pixel_frame_frac(self, A):
        return float(np.mean((A < 0).any(axis=1)))


class GammaFamily(Family):
    """iid Gamma(shape=k, rate=k): mean 1, var 1/k. Score s_j = (k-1)/a_j - k.
    k=1: score is identically -1 (boundary-violation exhibit arm).
    NO clamping — nonfinite-score frames dropped and counted."""

    def __init__(self, n, k):
        self.name = "GAM%d" % k
        self.n = n
        self.k = float(k)

    def sample(self, m, rng):
        return rng.gamma(shape=self.k, scale=1.0 / self.k, size=(m, self.n))

    def score(self, A):
        if self.k == 1.0:
            return -np.ones_like(A), 0
        with np.errstate(divide="ignore", over="ignore"):
            S = (self.k - 1.0) / A - self.k
        bad = ~np.isfinite(S).all(axis=1)
        n_bad = int(bad.sum())
        if n_bad:
            S[bad] = np.nan  # caller drops via valid_mask
        return S, n_bad

    def true_cov_op(self):
        return DiagCovOp(1.0 / self.k)

    def score_second_moment(self):
        """E[s_j^2] = k^2/(k-2) for k>2; infinite for k<=2."""
        if self.k > 2:
            return self.k ** 2 / (self.k - 2.0)
        return np.inf


class CorrLognFamily(Family):
    """z ~ N(m, Sigma_z) (Gaussian kernel, corr len ell, sigma_ln^2=0.25),
    a = e^z, m = -sigma_ln^2/2 (unit mean).
    Score: s(a) = -diag(1/a) [Sigma_z^{-1}(log a - m) + 1].
    True covariance: Cov_ij = exp(Sigma_z_ij) - 1 (unit mean)."""

    def __init__(self, n, corr_len=C.CORR_LEN_PX, name="CORR-LOGN"):
        self.name = name
        self.n = n
        side = int(round(np.sqrt(n)))
        assert side * side == n
        self.Sigma_z = grid_kernel_cov(side, corr_len, C.LOGN_SIGMA2)
        self.chol_z = CholOp(self.Sigma_z)
        self.m = -C.LOGN_SIGMA2 / 2.0
        self._cov_op = None

    def sample(self, m_frames, rng):
        Z = rng.standard_normal((m_frames, self.n))
        return np.exp(self.m + self.chol_z.sample_mult(Z))

    def score(self, A):
        logA = np.log(A)
        W = self.chol_z.solve((logA - self.m).T).T  # Sigma_z^{-1}(log a - m)
        S = -(W + 1.0) / A
        bad = ~np.isfinite(S).all(axis=1)
        n_bad = int(bad.sum())
        if n_bad:
            S[bad] = np.nan
        return S, n_bad

    def log_density(self, A):
        """log p(a) up to a constant shared across same-dimension components:
        -0.5 (z-m)^T Sigma_z^{-1} (z-m) - 0.5 logdet(Sigma_z) - sum(z)   [z = log a;
        the Jacobian term -sum log a is included for correctness although it
        cancels between MIX components]."""
        logA = np.log(A)
        Y = self.chol_z.half_solve((logA - self.m).T)  # L^{-1}(z-m)
        quad = np.sum(Y * Y, axis=0)
        return -0.5 * quad - 0.5 * self.chol_z.logdet() - logA.sum(axis=1)

    def true_cov(self):
        return np.expm1(self.Sigma_z)

    def true_cov_op(self):
        if self._cov_op is None:
            self._cov_op = CholOp(self.true_cov())
        return self._cov_op


class MixLognFamily(Family):
    """Two-state mixture (pi = 0.5/0.5) of CORR-LOGN with corr lens 2 px / 3 px,
    same marginal. Per-frame latent state; labels never given to estimators.
    Score: posterior-weighted sum, weights softmax in log domain (spec §2/§5)."""

    def __init__(self, n):
        self.name = "MIX-LOGN"
        self.n = n
        self.comps = [
            CorrLognFamily(n, corr_len=C.MIX_CORR_LENS[0], name="MIX-C0"),
            CorrLognFamily(n, corr_len=C.MIX_CORR_LENS[1], name="MIX-C1"),
        ]
        self.log_pi = np.log(np.array([0.5, 0.5]))
        self._cov_op = None
        self.last_states = None

    def sample(self, m, rng):
        states = rng.integers(0, 2, size=m)
        A = np.empty((m, self.n))
        for c in (0, 1):
            idx = np.nonzero(states == c)[0]
            if idx.size:
                A[idx] = self.comps[c].sample(idx.size, rng)
        self.last_states = states  # ground-truth latent states (evaluation only)
        return A

    def posterior_logw(self, A):
        lp = np.stack([f.log_density(A) for f in self.comps], axis=1)  # (m,2)
        lw = lp + self.log_pi[None, :]
        lw = lw - lw.max(axis=1, keepdims=True)
        w = np.exp(lw)
        w /= w.sum(axis=1, keepdims=True)
        return w

    def score(self, A):
        w = self.posterior_logw(A)
        S = np.zeros_like(A)
        tot_bad = 0
        for c in (0, 1):
            Sc, n_bad = self.comps[c].score(A)
            tot_bad = max(tot_bad, n_bad)
            S += w[:, c][:, None] * Sc
        bad = ~np.isfinite(S).all(axis=1)
        n_bad = int(bad.sum())
        if n_bad:
            S[bad] = np.nan
        return S, n_bad

    def true_cov(self):
        return 0.5 * self.comps[0].true_cov() + 0.5 * self.comps[1].true_cov()

    def true_cov_op(self):
        if self._cov_op is None:
            self._cov_op = CholOp(self.true_cov())
        return self._cov_op

    def oracle_state_accuracy(self, A, states):
        w = self.posterior_logw(A)
        pred = (w[:, 1] > 0.5).astype(int)
        acc = float(np.mean(pred == states))
        return max(acc, 1.0 - acc)


def make_family(name, n):
    if name == "GAUSS":
        return GaussFamily(n)
    if name.startswith("GAM"):
        return GammaFamily(n, int(name[3:]))
    if name == "CORR-LOGN":
        return CorrLognFamily(n)
    if name == "MIX-LOGN":
        return MixLognFamily(n)
    raise ValueError(name)


def marginal_var(name):
    """Analytic per-pixel variance (UT2)."""
    if name == "GAUSS":
        return C.GAUSS_SIGMA ** 2 * (1.0 + C.KERNEL_JITTER_REL)
    if name.startswith("GAM"):
        return 1.0 / int(name[3:])
    if name in ("CORR-LOGN", "MIX-LOGN"):
        return float(np.expm1(C.LOGN_SIGMA2 * (1.0 + C.KERNEL_JITTER_REL)))
    raise ValueError(name)
