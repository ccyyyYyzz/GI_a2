"""A1 Stein witness — spec §7-A1, primer §3 identity two.

v_hat = -(1/Mw) sum_i <a_i - mu, v> s(a_i);  err(v) = ||v_hat - v|| / ||v||.

Probes (64, frozen from the SEED0 probe stream): 32 forced-zero-mean random,
16 nonzero-sum random, 16 structured (constant, stripes, checkerboard,
low-frequency Fourier modes). All witness estimates use CENTERED probes
<a - mu, v>; GAM1 additionally evaluates the UNCENTERED estimator for its
predicted failure-mode exhibit.
"""
import numpy as np

from . import config as C
from .utils import rng_for


def build_probes(n):
    """(n, 64) probe matrix + probe_type list. Frozen: derived from SEED0 stream."""
    side = int(round(np.sqrt(n)))
    rng = rng_for(0, 0, C.STREAM_PROBES)
    probes = []
    types = []
    for _ in range(32):
        v = rng.standard_normal(n)
        v -= v.mean()          # forced zero-mean
        probes.append(v)
        types.append("zero_mean_random")
    for _ in range(16):
        v = rng.standard_normal(n)
        if abs(v.sum()) < 1e-6:
            v += 1.0 / n
        probes.append(v)
        types.append("nonzero_sum_random")
    ii, jj = np.meshgrid(np.arange(side), np.arange(side), indexing="ij")
    structured = [np.ones((side, side))]
    structured.append(np.where(ii % 2 == 0, 1.0, -1.0))            # horizontal stripes
    structured.append(np.where(jj % 2 == 0, 1.0, -1.0))            # vertical stripes
    structured.append(np.where((ii + jj) % 2 == 0, 1.0, -1.0))     # checkerboard
    for (kx, ky, fn) in [(1, 0, np.cos), (0, 1, np.cos), (1, 1, np.cos),
                         (2, 0, np.cos), (0, 2, np.cos), (2, 1, np.cos),
                         (1, 2, np.cos), (2, 2, np.cos),
                         (1, 0, np.sin), (0, 1, np.sin), (1, 1, np.sin),
                         (2, 1, np.sin)]:
        structured.append(fn(2 * np.pi * (kx * ii / side + ky * jj / side)))
    assert len(structured) == 16
    for v in structured:
        probes.append(v.ravel().astype(np.float64))
        types.append("structured")
    P = np.stack(probes, axis=1)
    return P, types


def witness_estimates(A, S, P, mu=1.0, centered=True, prefixes=None):
    """Compute v_hat for all probes. A: (Mw, n) patterns, S: (Mw, n) scores,
    P: (n, K) probes. Returns dict Mw_used -> (n, K) v_hat.
    prefixes: list of frame counts (nested prefixes of the same frozen draw,
    for the GAM2 convergence-rate check); default = [Mw]."""
    Mw = A.shape[0]
    if prefixes is None:
        prefixes = [Mw]
    out = {}
    for m in prefixes:
        Am, Sm = A[:m], S[:m]
        T = (Am - mu) @ P if centered else Am @ P   # (m, K)
        out[m] = -(Sm.T @ T) / m
    return out


def witness_err(vhat, P):
    """Relative L2 error per probe."""
    return np.linalg.norm(vhat - P, axis=0) / np.linalg.norm(P, axis=0)


def interp_se(A, S, P, mu=1.0, centered=True):
    """Per-probe interpolated SE (spec §7-A1): sqrt(tr Var_hat[<a-mu,v> s(a)] / Mw),
    relative to ||v||. Chunked accumulation of first/second moments."""
    Mw, n = A.shape
    K = P.shape[1]
    T = (A - mu) @ P if centered else A @ P        # (Mw, K)
    se2 = np.zeros(K)
    mean_acc = np.zeros((n, K))
    sq_acc = np.zeros((n, K))
    chunk = 20000
    for lo in range(0, Mw, chunk):
        Sc = S[lo:lo + chunk]
        Tc = T[lo:lo + chunk]
        for k in range(K):
            Z = Tc[:, k][:, None] * Sc            # (m, n)
            mean_acc[:, k] += Z.sum(axis=0)
            sq_acc[:, k] += (Z * Z).sum(axis=0)
    mean_acc /= Mw
    sq_acc /= Mw
    var = sq_acc - mean_acc ** 2
    se2 = var.sum(axis=0) / Mw
    return np.sqrt(se2) / np.linalg.norm(P, axis=0)


def q_gamma(k, n, Mw):
    """Analytic RMS relative witness error for iid Gamma(k, rate k), k > 2."""
    return np.sqrt(k * (n + 1) / ((k - 2.0) * Mw))
