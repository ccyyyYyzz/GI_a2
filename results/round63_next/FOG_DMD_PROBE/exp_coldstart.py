# Cold-start ESCAPE diagnostic: which strategy (if any) cracks blind-from-cold recovery?
# One config (the smoke config). Prints null-err for each strategy. Target <= 0.3.
import time
import numpy as np
import torch
from fog_common import (make_scene, smooth_basis, make_patterns, projectors,
                        null_metric, ou_coeffs)
import fog_tracker as ft

def log(*a): print(*a, flush=True)
dev = 'cuda' if torch.cuda.is_available() else 'cpu'

rng = np.random.default_rng(20260723)
n_side = 32; N = n_side * n_side; M = 64
x_true = make_scene(n_side)
P = make_patterns(M, N, "gauss", rng)
Pi, rangeP = projectors(P)
null_true = x_true - rangeP(x_true); nrm0 = np.linalg.norm(null_true)
c = 4; d_w = c * c; sig_w = 0.30; tau = 16; T = 128
U = smooth_basis(c, n_side)
Z, rho = ou_coeffs(T, d_w, sig_w, tau, rng)
W = np.clip(1.0 + Z @ U.T, 0.05, None)
mu = np.einsum('mn,tn->tm', P, W * x_true)
B = mu.copy(); wts = np.ones_like(B); Rd = np.ones_like(B) * 1e-8
lam_x = 1e-6
def m(x): return null_metric(x, x_true, rangeP, null_true, nrm0)
log(f"config d_w={d_w} sig={sig_w} tau={tau} rho={rho:.3f} T={T}; oracle & escape strategies\n")

x = ft.solve_oracle(P, W, B, wts, lam_x, dev); log(f"ORACLE                       null {m(x)[0]:.3f}")

# 1. spectral init Z -> what does the raw spectral W give x?
Zsp = ft.spectral_init_Z(P, U, B, sig_w, dev)
x = ft.track_em(P, U, B, wts, rho, sig_w, Rd, lam_x, n_em=1, Z_init=Zsp, dev=dev)
log(f"spectral-init (1 M-step)     null {m(x)[0]:.3f}")

# 2. spectral init -> full track-EM
for ne in [20, 60]:
    x = ft.track_em(P, U, B, wts, rho, sig_w, Rd, lam_x, n_em=ne, anneal=False, Z_init=Zsp, dev=dev)
    log(f"spectral -> track-EM n_em={ne:<3d}  null {m(x)[0]:.3f}")

# 3. spectral init -> joint Adam
x = ft.joint_adam(P, U, B, wts, rho, sig_w, lam_x, steps=6000, lr=1e-2, dev=dev, Z_init=Zsp)
log(f"spectral -> joint-Adam       null {m(x)[0]:.3f}")

# 4. cold + null kick -> track-EM (several kick sizes / seeds)
for kick in [0.05, 0.2]:
    best = 9
    for sd in range(3):
        torch.manual_seed(sd)
        x = ft.track_em(P, U, B, wts, rho, sig_w, Rd, lam_x, n_em=40, anneal=False,
                        null_kick=kick, dev=dev)
        best = min(best, m(x)[0])
    log(f"cold+null_kick={kick:<4}         null {best:.3f} (best of 3)")

# 5. many-restart joint Adam from random Z
best = 9
for sd in range(4):
    x = ft.joint_adam(P, U, B, wts, rho, sig_w, lam_x, steps=6000, lr=2e-2,
                      dev=dev, z_scale_init=sig_w, seed=sd)
    best = min(best, m(x)[0])
log(f"random-restart Adam          null {best:.3f} (best of 4)")

# 6. spectral init magnitude sweep (sign/scale of the rotated init matters?)
for sc in [0.5, 1.0, 2.0]:
    Zs = Zsp * sc
    x = ft.track_em(P, U, B, wts, rho, sig_w, Rd, lam_x, n_em=40, anneal=False, Z_init=Zs, dev=dev)
    log(f"spectral x{sc} -> track-EM     null {m(x)[0]:.3f}")
