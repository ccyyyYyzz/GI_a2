# Smoke test: one E1 config, validate oracle + tracking-EM(cold) + joint-Adam on GPU.
import sys, time
import numpy as np
import torch
from fog_common import (make_scene, smooth_basis, make_patterns, projectors,
                        null_metric, ou_coeffs)
import fog_tracker as ft

def log(*a): print(*a, flush=True)

dev = 'cuda' if torch.cuda.is_available() else 'cpu'
log("device:", dev, torch.cuda.get_device_name(0) if dev == 'cuda' else '')

rng = np.random.default_rng(20260723)
n_side = 32; N = n_side * n_side; M = 64
x_true = make_scene(n_side)
P = make_patterns(M, N, "gauss", rng)
Pi, rangeP = projectors(P)
null_true = x_true - rangeP(x_true); nrm0 = np.linalg.norm(null_true)

# config
c = 4; d_w = c * c; sig_w = 0.30; tau = 16; T = 128
U = smooth_basis(c, n_side)
Z, rho = ou_coeffs(T, d_w, sig_w, tau, rng)
W = np.clip(1.0 + Z @ U.T, 0.05, None)
mu = np.einsum('mn,tn->tm', P, W * x_true)  # clean buckets
B = mu.copy()
wts = np.ones_like(B)  # clean: unit weights
lam_x = 1e-6

log(f"\n[config] d_w={d_w} sig_w={sig_w} tau={tau} rho={rho:.4f} T={T}  ||null||={nrm0:.3f}")

t0 = time.time()
x_or = ft.solve_oracle(P, W, B, wts, lam_x, dev)
ne, te = null_metric(x_or, x_true, rangeP, null_true, nrm0)
log(f"ORACLE (w known)          null-err {ne:.3f}  total-err {te:.3f}   [{time.time()-t0:.1f}s]")

t0 = time.time()
x_em = ft.track_em(P, U, B, wts, rho, sig_w, wts * 0 + 1e-8, lam_x,
                   n_em=20, anneal=True, dev=dev)
ne, te = null_metric(x_em, x_true, rangeP, null_true, nrm0)
log(f"TRACK-EM cold (anneal)    null-err {ne:.3f}  total-err {te:.3f}   [{time.time()-t0:.1f}s]")

t0 = time.time()
x_em2 = ft.track_em(P, U, B, wts, rho, sig_w, wts * 0 + 1e-8, lam_x,
                    n_em=20, anneal=False, dev=dev)
ne, te = null_metric(x_em2, x_true, rangeP, null_true, nrm0)
log(f"TRACK-EM cold (no anneal) null-err {ne:.3f}  total-err {te:.3f}   [{time.time()-t0:.1f}s]")

t0 = time.time()
x_ad = ft.joint_adam(P, U, B, wts, rho, sig_w, lam_x, steps=4000, lr=2e-2, dev=dev)
ne, te = null_metric(x_ad, x_true, rangeP, null_true, nrm0)
log(f"JOINT-ADAM cold           null-err {ne:.3f}  total-err {te:.3f}   [{time.time()-t0:.1f}s]")
