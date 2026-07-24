# fog_sobi.py -- SOBI-style rotation fix for the fog-as-2nd-DMD cold start (E4a).
#
# The centered bucket fluctuation is, to first order, a linear BSS mixture:
#     b_t = B_t - mean_t = P diag(x) U z_t  =  H z_t   (H = P diag(x) U, M x d).
# If the medium coefficients z_k(t) are INDEPENDENT OU processes with DISTINCT
# correlation times tau_k, then every lagged source covariance E[z_t z_{t+l}^T] is
# diagonal with DISTINCT diagonals, so the lagged bucket covariances
#     R_b(l) = H diag(sigma_k^2 rho_k^l) H^T
# are simultaneously diagonalizable ONLY in the true source coordinates. SOBI
# recovers those coordinates (whiten with R_b(0), joint-diagonalize lags>=1), which
# fixes the d x d rotation ambiguity that the plain SVD Stage-B init cannot resolve.
# With a single tau (all rho_k equal) every R_b(l) is proportional to H H^T and the
# joint diagonalization is degenerate -- so SOBI helps ONLY for multi-timescale media,
# which is the attributable physical mechanism.
import numpy as np


def _lagged_cov(c, l):
    T = c.shape[0]
    R = (c[:T - l].T @ c[l:]) / (T - l)
    return 0.5 * (R + R.T)


def joint_diagonalize(mats, n_sweeps=200, tol=1e-10):
    """Cardoso-Souloumiac real Jacobi joint diagonalization of symmetric (d,d) matrices.
    Returns orthogonal V (d,d) that most-diagonalizes all mats simultaneously."""
    d = mats[0].shape[0]
    A = [m.copy() for m in mats]
    V = np.eye(d)
    for _ in range(n_sweeps):
        moved = 0.0
        for p in range(d - 1):
            for q in range(p + 1, d):
                g = np.array([[m[p, p] - m[q, q], m[p, q] + m[q, p]] for m in A])  # (K,2)
                G = g.T @ g
                ton, toff = G[0, 0] - G[1, 1], G[0, 1] + G[1, 0]
                theta = 0.5 * np.arctan2(toff, ton + np.sqrt(ton * ton + toff * toff) + 1e-300)
                c, s = np.cos(theta), np.sin(theta)
                if abs(s) < tol:
                    continue
                moved += abs(s)
                for m in A:
                    cp, cq = m[:, p].copy(), m[:, q].copy()
                    m[:, p] = c * cp + s * cq; m[:, q] = -s * cp + c * cq
                    rp, rq = m[p, :].copy(), m[q, :].copy()
                    m[p, :] = c * rp + s * rq; m[q, :] = -s * rp + c * rq
                vp, vq = V[:, p].copy(), V[:, q].copy()
                V[:, p] = c * vp + s * vq; V[:, q] = -s * vp + c * vq
        if moved < tol:
            break
    return V


def sobi_sources(b, d, lags=(1, 2, 3, 5, 8, 13, 21, 34)):
    """b:(T,M) centered fluctuation -> recovered source coeffs zc:(T,d) (up to perm/sign/scale)."""
    T, M = b.shape
    R0 = (b.T @ b) / T
    w, E = np.linalg.eigh(R0)
    idx = np.argsort(w)[::-1][:d]
    wd = np.maximum(w[idx], 1e-12); Ed = E[:, idx]
    Wh = (Ed / np.sqrt(wd)).T                    # (d,M) whitener
    c = b @ Wh.T                                  # (T,d) whitened
    Rs = [_lagged_cov(c, l) for l in lags if l < T - 5]
    V = joint_diagonalize(Rs)
    return c @ V                                  # (T,d) recovered sources
