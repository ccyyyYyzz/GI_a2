# Test the SPECTRAL-ROTATION blind solver against cold-start hardness.
import time
import numpy as np
import torch
from fog_common import (make_scene, smooth_basis, make_patterns, projectors,
                        null_metric, ou_coeffs)
import fog_tracker as ft
def log(*a): print(*a, flush=True)
dev = 'cuda' if torch.cuda.is_available() else 'cpu'

n_side = 32; N = n_side * n_side; M = 64; lam_x = 1e-6
x_true = make_scene(n_side)

def run(c, sig_w, tau, T, seed=20260723):
    rng = np.random.default_rng(seed)
    P = make_patterns(M, N, "gauss", rng)
    Pi, rangeP = projectors(P)
    null_true = x_true - rangeP(x_true); nrm0 = np.linalg.norm(null_true)
    U = smooth_basis(c, n_side); d_w = U.shape[1]
    Z, rho = ou_coeffs(T, d_w, sig_w, tau, rng)
    W = np.clip(1.0 + Z @ U.T, 0.05, None)
    mu = np.einsum('mn,tn->tm', P, W * x_true)
    B = mu.copy(); wts = np.ones_like(B); Rd = np.ones_like(B) * 1e-8
    def m(x): return null_metric(x, x_true, rangeP, null_true, nrm0)[0]
    xo = ft.solve_oracle(P, W, B, wts, lam_x, dev); ne_o = m(xo)
    t0 = time.time()
    xs, res = ft.solve_spectral_rot(P, U, B, wts, sig_w, lam_x, n_als=60, n_restart=10, dev=dev)
    ne_s = m(xs); dt = time.time() - t0
    xr, _ = ft.solve_spectral_rot(P, U, B, wts, sig_w, lam_x, n_als=60, n_restart=10,
                                  dev=dev, refine_em=15, rho=rho, R_diag_np=Rd)
    ne_r = m(xr)
    log(f"c={c} d_w={d_w:3d} sig={sig_w} tau={str(tau):>4} T={T:4d} rho={rho:.3f} | "
        f"oracle {ne_o:.3f} | SPECROT {ne_s:.3f} | +EM {ne_r:.3f} | {dt:.1f}s")
    return ne_s

log("SPECTRAL-ROTATION blind solver -- cold-start crack test (clean)\n")
for tau in [16, 4, 64, None]:
    run(4, 0.30, tau, 128)
log("")
for T in [64, 128, 256, 512]:
    run(4, 0.30, 16, T)
log("")
run(8, 0.30, 16, 256)   # d_w=64=M
run(8, 0.15, 16, 256)
run(4, 0.15, 16, 256)
