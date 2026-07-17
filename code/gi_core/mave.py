"""MAVE-16 (spec §4.9): refined MAVE for a single-index model, OPG initialization.

Xia et al. (2002) rMAVE specialization to d=1, run only at 16x16 (n=256).
Anchors subsampled for tractability; ridge-regularized local linear steps.
Returns a signed direction; scale handled by the metric protocol.
"""
import numpy as np


def _opg_direction(A, b, rng, n_anchor=150, n_loc=1500, ridge=1e-3, bw_scale=1.0):
    M, n = A.shape
    anchors = rng.choice(M, size=min(n_anchor, M), replace=False)
    grads = np.zeros((anchors.size, n))
    # squared distances via ||x||^2 trick, chunked over anchors
    norms = (A * A).sum(axis=1)
    for j, ai in enumerate(anchors):
        d2 = norms - 2.0 * (A @ A[ai]) + norms[ai]
        idx = np.argpartition(d2, n_loc)[:n_loc]
        h2 = (np.median(d2[idx]) + 1e-12) * bw_scale ** 2
        w = np.exp(-d2[idx] / (2.0 * h2))
        Xc = A[idx] - A[ai]
        Xd = np.concatenate([np.ones((idx.size, 1)), Xc], axis=1)
        Xw = Xd * w[:, None]
        G = Xw.T @ Xd
        G[np.diag_indices_from(G)] += ridge * np.trace(G) / (n + 1)
        beta = np.linalg.solve(G, Xw.T @ b[idx])
        grads[j] = beta[1:]
    # weighted PCA of gradient vectors -> top eigvec of grads^T grads
    _, _, Vt = np.linalg.svd(grads, full_matrices=False)
    return Vt[0]


def rmave_single_index(A, b, rng, n_anchor=150, n_loc=400, n_iter=10, ridge=1e-4,
                       bw_scale=1.0, ridge_scale=1.0,
                       opg_n_anchor=150, opg_n_loc=1500, opg_ridge=1e-3,
                       return_opg=False):
    """Returns direction beta (n,), sign-fixed so corr(A beta, b) > 0.

    Default arguments reproduce the Phase A pipeline exactly (bw_scale and
    ridge_scale are neutral 1.0 multipliers). The extra knobs exist for the
    ROUND62 G0 sensitivity audit and G1 scale map; every non-default value
    used in a run is logged to that run's JSON."""
    M, n = A.shape
    beta = _opg_direction(A, b, rng, n_anchor=opg_n_anchor,
                          n_loc=min(opg_n_loc, M),
                          ridge=opg_ridge * ridge_scale, bw_scale=bw_scale)
    beta = beta / np.linalg.norm(beta)
    beta_opg = beta.copy()
    anchors = rng.choice(M, size=min(n_anchor, M), replace=False)
    for _ in range(n_iter):
        u = A @ beta
        H = np.zeros((n, n))
        rhs = np.zeros(n)
        # 1-d kernel bandwidth: Silverman on the current index
        h = (1.06 * u.std() * M ** (-0.2) + 1e-12) * bw_scale
        for ai in anchors:
            du = u - u[ai]
            idx = np.argpartition(np.abs(du), n_loc)[:n_loc]
            w = np.exp(-(du[idx] ** 2) / (2.0 * h * h))
            sw = w.sum()
            if sw <= 0:
                continue
            # local linear in the index: y ~ c + d * (u - u_a)
            wu = w * du[idx]
            Suu = np.dot(wu, du[idx])
            Su = wu.sum()
            Sy = np.dot(w, b[idx])
            Suy = np.dot(wu, b[idx])
            det = sw * Suu - Su * Su
            if abs(det) < 1e-30:
                continue
            c_j = (Suu * Sy - Su * Suy) / det
            d_j = (sw * Suy - Su * Sy) / det
            # accumulate beta-step normal equations:
            # min sum_i w_i (y_i - c_j - d_j beta^T (a_i - a_j))^2
            Xc = A[idx] - A[ai]
            Xw = Xc * (w * d_j ** 2)[:, None]
            H += Xw.T @ Xc
            rhs += Xc.T @ (w * d_j * (b[idx] - c_j))
        H[np.diag_indices_from(H)] += ridge * ridge_scale * np.trace(H) / n
        beta_new = np.linalg.solve(H, rhs)
        nrm = np.linalg.norm(beta_new)
        if nrm < 1e-30:
            break
        beta_new /= nrm
        if np.dot(beta_new, beta) < 0:
            beta_new = -beta_new
        delta = np.linalg.norm(beta_new - beta)
        beta = beta_new
        if delta < 1e-5:
            break
    u = A @ beta
    if np.corrcoef(u, b)[0, 1] < 0:
        beta = -beta
    if return_opg:
        if np.corrcoef(A @ beta_opg, b)[0, 1] < 0:
            beta_opg = -beta_opg
        return beta, beta_opg
    return beta
