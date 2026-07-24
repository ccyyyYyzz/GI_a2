# Confirm fog_common reproduces the established toy numbers before building trackers.
# Expected (from task brief): stacked rank 1024/1024; oracle null-err ~0.03-0.00;
# warm ALS (truth+15%) ~0.10-0.56; cold ALS ~1.0.
import numpy as np
from fog_common import (make_scene, smooth_basis, make_patterns, projectors, null_metric)

rng = np.random.default_rng(20260723)
n_side = 32; N = n_side * n_side; M = 64; T = 64; SIG_W = 0.30
x_true = make_scene(n_side)
P = make_patterns(M, N, "gauss", rng)  # matches fog_dmd_oxygen2 (gaussian)
Pi, rangeP = projectors(P)
null_true = x_true - rangeP(x_true); nrm0 = np.linalg.norm(null_true)

for c, lab in [(4, "d_w=16<M"), (8, "d_w=64=M"), (16, "d_w=256>M")]:
    U = smooth_basis(c, n_side); d_w = U.shape[1]
    Z = SIG_W * rng.standard_normal((T, d_w))
    W = np.clip(1.0 + Z @ U.T, 0.05, None)
    Bkt = np.einsum('mn,tn->tm', P, W * x_true)
    print(f"[{lab}]")
    # oracle
    Aeff = (P[None, :, :] * W[:, None, :]).reshape(T * M, N)
    xh = np.linalg.solve(Aeff.T @ Aeff + 1e-8 * np.eye(N), Aeff.T @ Bkt.reshape(-1))
    print(f"  stacked rank = {np.linalg.matrix_rank(Aeff)} / {N}")
    ne, te = null_metric(xh, x_true, rangeP, null_true, nrm0)
    print(f"  oracle (w known)            null-err {ne:.3f}  total-err {te:.3f}")
    # warm ALS (truth + 15%)
    z = Z + 0.15 * SIG_W * rng.standard_normal(Z.shape)
    for it in range(25):
        Wst = np.clip(1.0 + z @ U.T, 0.05, None)
        Aeff = (P[None, :, :] * Wst[:, None, :]).reshape(T * M, N)
        xh = np.linalg.solve(Aeff.T @ Aeff + 1e-8 * np.eye(N), Aeff.T @ Bkt.reshape(-1))
        G = P @ (U * xh[:, None]); GtG = G.T @ G + 1e-6 * np.eye(d_w)
        for t in range(T): z[t] = np.linalg.solve(GtG, G.T @ (Bkt[t] - P @ xh))
    ne, te = null_metric(xh, x_true, rangeP, null_true, nrm0)
    print(f"  ALS warm init (truth+15%)   null-err {ne:.3f}  total-err {te:.3f}")
    # cold ALS
    z = np.zeros((T, d_w))
    for it in range(25):
        Wst = np.clip(1.0 + z @ U.T, 0.05, None)
        Aeff = (P[None, :, :] * Wst[:, None, :]).reshape(T * M, N)
        xh = np.linalg.solve(Aeff.T @ Aeff + 1e-8 * np.eye(N), Aeff.T @ Bkt.reshape(-1))
        G = P @ (U * xh[:, None]); GtG = G.T @ G + 1e-6 * np.eye(d_w)
        for t in range(T): z[t] = np.linalg.solve(GtG, G.T @ (Bkt[t] - P @ xh))
    ne, te = null_metric(xh, x_true, rangeP, null_true, nrm0)
    print(f"  ALS cold init               null-err {ne:.3f}  total-err {te:.3f}")
