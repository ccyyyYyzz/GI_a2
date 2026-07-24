# fog_common.py -- shared exact setup + metric for the "fog as a second DMD" probe.
# Matches the local CPU toys EXACTLY (fog_dmd_oxygen*.py, fog_noise_oxygen.py) so
# numbers are directly comparable. All heavy solves optionally run on GPU via torch.
import numpy as np


def make_scene(n_side):
    yy, xx = np.mgrid[0:n_side, 0:n_side] / n_side
    x = (np.exp(-((xx - .3) ** 2 + (yy - .4) ** 2) / .02)
         + .8 * np.exp(-((xx - .7) ** 2 + (yy - .6) ** 2) / .01))
    # bars scale with resolution so high-freq content is preserved at 64x64 too
    s = n_side // 32
    x[10 * s:22 * s, 15 * s:18 * s] += 1.0
    x[5 * s:8 * s, 5 * s:26 * s] += .7
    x = (x / x.max()).ravel()
    return x


def smooth_basis(c, n_side):
    """Orthonormal smooth spatial basis for the medium field (matches references)."""
    B = []
    for i in range(c):
        for j in range(c):
            g = np.zeros((c, c)); g[i, j] = 1.0
            up = np.kron(g, np.ones((n_side // c, n_side // c)))
            up = (up + .5 * np.roll(up, 1, 0) + .5 * np.roll(up, -1, 0)
                  + .5 * np.roll(up, 1, 1) + .5 * np.roll(up, -1, 1))
            B.append(up.ravel())
    Q, _ = np.linalg.qr(np.array(B).T)
    return Q  # N x c^2, orthonormal columns


def dct_basis(d, n_side, exclude_dc=True):
    """Orthonormal 2D DCT-II low-frequency smooth basis, d columns (N x d), DC excluded.
    Predeclared smooth basis family (ruling 6.3). Columns ordered by radial frequency."""
    n = n_side
    # 1D DCT-II orthonormal matrix (n x n): row k = basis mode k over pixels
    k = np.arange(n)[:, None]; i = np.arange(n)[None, :]
    D = np.cos(np.pi * (2 * i + 1) * k / (2 * n))
    D *= np.sqrt(2.0 / n); D[0] *= 1 / np.sqrt(2)          # orthonormal rows
    modes = []
    for p in range(n):
        for q in range(n):
            if exclude_dc and p == 0 and q == 0:
                continue
            modes.append((p * p + q * q, p, q))
    modes.sort()
    cols = []
    for _, p, q in modes[:d]:
        phi = np.outer(D[p], D[q])                          # (n,n) mode, orthonormal
        cols.append(phi.ravel())
    U = np.array(cols).T                                    # (N,d)
    # numerical re-orthonormalization (safety)
    Q, _ = np.linalg.qr(U)
    return Q


def lognormal_medium(U, T, sigma_f, tau, rng):
    """Positive mean-normalized (E[w]=1 per pixel) lognormal medium with fixed pixelwise
    field RMS sigma_f, independent of d (ruling 2.2/6.3).
    per-coeff sd = sigma_f*sqrt(N/d) so that RMS(g) = sigma_f; w_t = exp(g_t - var_n/2)."""
    N, d = U.shape
    coeff_sd = sigma_f * np.sqrt(N / d)
    if tau is None:
        Z = coeff_sd * rng.standard_normal((T, d)); rho = 0.0
    else:
        rho = np.exp(-1.0 / tau)
        Z = np.zeros((T, d)); Z[0] = coeff_sd * rng.standard_normal(d)
        a = np.sqrt(1.0 - rho ** 2) * coeff_sd
        for t in range(1, T):
            Z[t] = rho * Z[t - 1] + a * rng.standard_normal(d)
    G = Z @ U.T                                             # (T,N) log-field
    var_n = (coeff_sd ** 2) * (U ** 2).sum(axis=1)          # per-pixel var of g
    W = np.exp(G - 0.5 * var_n[None, :])                    # E[w]=1 per pixel exactly
    return W, Z, rho, coeff_sd


def tau_schedule(d, tau_coarse=64.0, tau_fine=4.0):
    """Per-component OU correlation times, log-spaced from tau_coarse (low-freq DCT modes,
    k=0) to tau_fine (high-freq, k=d-1). Finer spatial scales decorrelate faster (turbulence)."""
    return tau_coarse * (tau_fine / tau_coarse) ** (np.arange(d) / max(1, d - 1))


def lognormal_medium_mt(U, T, sigma_f, tau_arr, rng):
    """Multi-timescale positive lognormal medium: each DCT component z_k is an independent OU
    with its own tau_k, holding pixelwise RMS sigma_f fixed (per-coeff sd = sigma_f*sqrt(N/d)).
    Returns W:(T,N), Z:(T,d), rho:(d,), coeff_sd."""
    N, d = U.shape
    coeff_sd = sigma_f * np.sqrt(N / d)
    rho = np.exp(-1.0 / np.asarray(tau_arr, dtype=float))     # (d,)
    a = np.sqrt(1.0 - rho ** 2) * coeff_sd                    # (d,)
    Z = np.zeros((T, d)); Z[0] = coeff_sd * rng.standard_normal(d)
    for t in range(1, T):
        Z[t] = rho * Z[t - 1] + a * rng.standard_normal(d)
    G = Z @ U.T
    var_n = (coeff_sd ** 2) * (U ** 2).sum(axis=1)
    W = np.exp(G - 0.5 * var_n[None, :])
    return W, Z, rho, coeff_sd


def make_patterns(M, N, kind, rng):
    if kind == "gauss":
        return rng.standard_normal((M, N)) / np.sqrt(N)
    elif kind == "binary":
        return (rng.random((M, N)) < 0.5).astype(float)  # physical 0/1
    raise ValueError(kind)


def ou_coeffs(T, d_w, sig_w, tau, rng):
    """OU-correlated medium coefficients. rho=exp(-1/tau); stationary std = sig_w.
    tau=None -> iid (rho=0), matching the original toys."""
    if tau is None:
        return sig_w * rng.standard_normal((T, d_w)), 0.0
    rho = np.exp(-1.0 / tau)
    Z = np.zeros((T, d_w))
    Z[0] = sig_w * rng.standard_normal(d_w)
    a = np.sqrt(1.0 - rho ** 2) * sig_w
    for t in range(1, T):
        Z[t] = rho * Z[t - 1] + a * rng.standard_normal(d_w)
    return Z, rho


def projectors(P):
    Pi = np.linalg.pinv(P)
    rangeP = lambda v: Pi @ (P @ v)
    return Pi, rangeP


def null_metric(x_hat, x_true, rangeP, null_true, nrm0):
    """EXACT metric from the toys: oracle-gauge scale, then null-component error.
    Returns (amplitude null-err, amplitude total-err). null-NMSE = ne**2."""
    x_hat = np.asarray(x_hat, dtype=np.float64)
    s = np.dot(x_hat, x_true) / np.dot(x_hat, x_hat)
    x_hat = s * x_hat
    ne = np.linalg.norm((x_hat - rangeP(x_hat)) - null_true) / nrm0
    te = np.linalg.norm(x_hat - x_true) / np.linalg.norm(x_true)
    return float(ne), float(te)


def null_fraction(x_true, rangeP):
    """Fraction of scene energy living in ker(P) (the fixed-system null space)."""
    nt = x_true - rangeP(x_true)
    return float(np.linalg.norm(nt) ** 2 / np.linalg.norm(x_true) ** 2)
