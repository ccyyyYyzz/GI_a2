"""Shared 5%-sampling binary GI operator (built in-probe for cleanliness).

K = 51 binary masks on a 32x32 grid (P = 1024).  Sampling ratio K/P = 51/1024
= 4.98% ("5% linear system").  Each mask pixel ~ Bernoulli(0.5) (canonical
random binary GI illumination).  Deterministic from a fixed seed.
"""
import numpy as np

SIDE = 32
P = SIDE * SIDE
K = 51
SEED_OP = 651_000            # disjoint from scene seeds (650xxx) and DEV seeds


def build_operator(K=K, P=P, density=0.5, seed=SEED_OP):
    """Return M (K x P) binary mask matrix, float64 in {0,1}."""
    rng = np.random.default_rng(seed)
    M = (rng.random((K, P)) < density).astype(np.float64)
    # guarantee no all-zero mask (would be a trivial constraint)
    for i in range(K):
        if M[i].sum() == 0:
            M[i, rng.integers(P)] = 1.0
    return M


def load_scene(path):
    """Load a bridge scene .npz -> flat float64 vector in [0,1]."""
    x = np.load(path)["x"].astype(np.float64).ravel()
    return x
