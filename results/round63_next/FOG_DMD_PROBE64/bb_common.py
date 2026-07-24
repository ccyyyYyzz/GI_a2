# bb_common.py -- shared BEYOND-MODULATOR-BAND config + estimators for the R39 risk gates.
#
# Encapsulates the 16x16 beyond-band cell EXACTLY as in the session scratchpad scripts
# (mini_superres_duel.py, mini_fresh_cov.py, mini_cov_predict.py, run_e8.py) so every gate
# number is directly comparable to the R39 evidence chain:
#   - pattern band  : fx,fy <= PB (=3)  (the DMD/modulator pixel-limit wall)
#   - medium band   : 1 <= i+j <= 6     (fine speckle BEYOND the pattern band)
#   - lognormal OU medium, TAU=8, SIG_F=0.30, PHOT=1e5, M=24 band-limited nonneg patterns
#   - metric hi_err : rel. DCT-power error on the deliberate beyond-band frequencies
#
# Estimators reuse the PRODUCTION solvers in ../FOG_DMD_PROBE/fog_tracker.py (no code drift):
#   cov_estimate  (E5a moment matching), kalman_em (E7a marginal-likelihood MLE), solve_oracle.
#
# WRITE-SCOPE: this module is read-only config; gate scripts write only under FOG_DMD_PROBE64/.
import os
import sys

import numpy as np
import torch
from scipy.fftpack import dct, idct

# reuse the frozen production tracker (cov_estimate / kalman_em / solve_oracle)
_FOG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "FOG_DMD_PROBE")
sys.path.insert(0, os.path.abspath(_FOG))
import fog_tracker as ft  # noqa: E402

DEV = "cuda" if torch.cuda.is_available() else "cpu"

# ----------------------------------------------------------------- canonical cell
N_SIDE = 16
N = N_SIDE * N_SIDE
PB = 3                       # pattern band (modulator wall): fx,fy <= PB
MED_LO, MED_HI = 1, 6        # true medium band: MED_LO <= i+j <= MED_HI (fine speckle)
SIG_F = 0.30
TAU = 8.0
PHOT = 1e5
RHO = float(np.exp(-1.0 / TAU))
M_PAT = 24                   # band-limited pattern-bank size


def _D2(a, n=N_SIDE):
    return dct(dct(a.reshape(n, n), axis=0, norm="ortho"), axis=1, norm="ortho")


def _I2(A, n=N_SIDE):
    return idct(idct(A, axis=0, norm="ortho"), axis=1, norm="ortho").ravel()


def dct_spike(i, j, n=N_SIDE):
    """Real-space image of a single 2D-DCT frequency (i,j)."""
    g = np.zeros((n, n))
    g[i, j] = 1.0
    return _I2(g, n)


def medium_subspace(lo=MED_LO, hi=MED_HI, n=N_SIDE, radial_profile=None):
    """Orthonormal medium subspace Ub (N x db) spanning DCT modes with lo<=i+j<=hi.
    Returns (Ub, db, coeff_sd, Kg). coeff_sd = SIG_F*sqrt(N/db) keeps pixelwise field
    RMS = SIG_F independent of db (ruling normalization).
    radial_profile: optional callable r->weight applied per mode BEFORE the QR (used only
    to *generate* a non-flat true medium in G1(d); the estimator always declares flat)."""
    modes = [(i, j) for i in range(n) for j in range(n) if lo <= i + j <= hi]
    cols = [dct_spike(i, j, n) for (i, j) in modes]
    A = np.array(cols).T                                  # (N, db) raw spikes (already orthonormal)
    if radial_profile is not None:
        w = np.array([radial_profile(np.hypot(i, j)) for (i, j) in modes], dtype=float)
        w = w / np.sqrt((w ** 2).mean())                  # keep mean power ~1
        A = A * w[None, :]
    Ub = np.linalg.qr(A)[0]
    db = Ub.shape[1]
    coeff_sd = SIG_F * np.sqrt(N / db)
    Kg = (SIG_F ** 2 * N / db) * (Ub @ Ub.T)
    return Ub, db, coeff_sd, Kg


# deliberate beyond-band scene content (frequencies OUTSIDE the pattern band PB=3)
HI_FREQS = [(5, 2), (2, 5), (4, 4), (6, 1), (1, 6), (5, 5), (7, 0), (0, 7)]


def make_scene(seed=4, n=N_SIDE, pb=PB, hi_freqs=None):
    """In-band random DCT content + deliberate beyond-band spikes -> nonneg image in [0,1]."""
    hi_freqs = HI_FREQS if hi_freqs is None else hi_freqs
    C = np.zeros((n, n))
    rs = np.random.default_rng(seed)
    C[:pb + 1, :pb + 1] = rs.standard_normal((pb + 1, pb + 1))
    for (i, j) in hi_freqs:
        C[i, j] = 2.0 * rs.choice([-1, 1])
    x = _I2(C, n)
    x = x - x.min()
    x = x / x.max()
    hi_mask = np.zeros((n, n), bool)
    for (i, j) in hi_freqs:
        hi_mask[i, j] = True
    return x, hi_mask


def bl_patterns(M_, seed, pb=PB, n=N_SIDE):
    """Band-limited (fx,fy<=pb) incoherent nonneg random patterns (the DMD hardware limit)."""
    rp = np.random.default_rng(seed)
    ps = []
    for _ in range(M_):
        Cp = np.zeros((n, n))
        Cp[:pb + 1, :pb + 1] = rp.standard_normal((pb + 1, pb + 1))
        Cp[0, 0] = 0
        f = _I2(Cp, n)
        f /= np.abs(f).max()
        ps.append(0.5 + 0.45 * f)
    return np.array(ps)


def medium_field(T, seed, Ub, coeff_sd, rho=RHO, gen_sd=None):
    """Lognormal OU medium W (T x N), E[w]=1 per pixel. gen_sd (per-mode sd vector) allows a
    NON-flat radial profile in generation while the estimator still uses scalar coeff_sd."""
    db = Ub.shape[1]
    r = np.random.default_rng(seed)
    sd = coeff_sd if gen_sd is None else np.asarray(gen_sd, float)
    z = sd * r.standard_normal(db)
    var_n = (Ub ** 2 * (sd ** 2 if gen_sd is not None else coeff_sd ** 2)).sum(axis=1) \
        if gen_sd is not None else (coeff_sd ** 2) * (Ub ** 2).sum(axis=1)
    W = np.zeros((T, N))
    for t in range(T):
        g = Ub @ z
        W[t] = np.exp(g - 0.5 * var_n)
        z = rho * z + np.sqrt(1 - rho ** 2) * sd * r.standard_normal(db)
    return W


def buckets(P, W, x, seed, phot=PHOT):
    """Poisson bucket record for the multiplicative thin-screen model b=<P (.) w, x>.
    Returns (B, R_diag, scp)."""
    r = np.random.default_rng(seed)
    mu = np.einsum("mn,tn->tm", P, W * x)
    scp = phot / mu.mean()
    B = r.poisson(mu * scp) / scp
    R_diag = np.clip(mu, 1e-9, None) / scp
    return B, R_diag, scp


def hi_err(xh, x, hi_mask, n=N_SIDE):
    """Relative 2D-DCT power error on the beyond-band frequencies (1.0 = nothing recovered).
    Gauge: scalar amplitude match (exactly as the scratchpad scripts)."""
    xh = np.asarray(xh, float)
    s = np.dot(xh, x) / max(np.dot(xh, xh), 1e-12)
    E = _D2(s * xh - x, n) ** 2
    X0 = _D2(x, n) ** 2
    return float(E[hi_mask].sum() / X0[hi_mask].sum())


# ---------------------------------------------------------------- estimator wrappers
def est_moment(P, Ub, B, coeff_sd, rho=RHO, n_starts=3, steps=4000, dev=DEV, seed0=0):
    """E5a moment-matching (production cov_estimate). Returns x of the lowest-objective start."""
    rc = ft.cov_estimate(P, Ub, B, coeff_sd, rho, n_starts=n_starts, steps=steps,
                         dev=dev, seed0=seed0)
    return rc["per_start_x"][int(np.argmin(rc["objs"]))], rc


def est_mle(P, Ub, B, R_diag, rho, coeff_sd, x_init, n_em=25, dev=DEV, lam_x=1e-4):
    """E7a two-stage MLE: moment init -> Kalman EM (marginal likelihood)."""
    return ft.kalman_em(P, Ub, B, R_diag, rho, coeff_sd, lam_x, n_em=n_em, dev=dev,
                        x_init=x_init)


def est_oracle(P, W, B, dev=DEV, lam_x=1e-6):
    return ft.solve_oracle(P, W, B, np.ones_like(B), lam_x, dev)


def fresh_wall(P_fix_seed, Ub, coeff_sd, x, T, seed, phot=PHOT, pb=PB, m_pat=M_PAT):
    """Physics wall: FRESH band-limited patterns + mean route. Band-limited operators cannot
    leave their band -> beyond-band hi_err -> 1.000 (the mathematically blind first-moment arm)."""
    W = medium_field(T, 600 + seed, Ub, coeff_sd)
    AtA = np.zeros((N, N))
    Atb = np.zeros(N)
    rb = np.random.default_rng(680 + seed)
    scp = phot / np.einsum("mn,tn->tm", bl_patterns(m_pat, P_fix_seed, pb), W * x).mean()
    for t in range(T):
        Pf = bl_patterns(m_pat, 7000 + 1000 * seed + t, pb)
        mu2 = Pf @ (W[t] * x)
        B2 = rb.poisson(mu2 * scp) / scp
        AtA += Pf.T @ Pf
        Atb += Pf.T @ B2
    lam = 1e-3 * np.trace(AtA) / N
    return np.clip(np.linalg.solve(AtA + lam * np.eye(N), Atb), 0, None)
