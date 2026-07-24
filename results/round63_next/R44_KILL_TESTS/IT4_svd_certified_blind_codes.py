#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""IT4 -- SVD certified-blind codes (internal divergence finer #3, THE SLEDGEHAMMER GATE).

The mean-channel signed intensity is EXACTLY LINEAR in the code:
    D_z = |E^+|^2-|E^-|^2 = Re[U_c conj(U_s)] ,  U_c=prop(dmd_field(1)),
    U_s(s)=prop(dmd_field(s))  (linear in s).
So the beyond-band mean leak is a linear operator L_be(z1): in-band code coeffs ->
beyond-band leak coeffs.  Build L_be(z1) and the in-band response L_in(z1) for
z1 in {0,2,5,10} mm, STACK over z1, and solve the generalized eigenproblem
    min_v  v^T (Lbe^T Lbe) v / v^T (Lin^T Lin) v
whose bottom eigenvectors are code directions that are joint-mean-blind at EVERY z1
(relative leak = sqrt(eigenvalue)).  Take the bottom-64 -> a certified-blind code
subspace.  Then a covariance-retention run at z2=5mm compares detection noncentrality
lambda for M=32 codes drawn from the blind subspace vs M=32 random in-band codes.

GATE (coordinator): bottom-64 joint leak <= 1e-4 AND lambda within 3x of random
=> certified-blind codes exist (next-paper core + bench calibration step).
KILL: leak >= 1e-3 OR lambda drop > 3x => floor is code-independent; hardware
concentration (IT5) becomes the only wall-restoration route.
"""
import os
os.environ.pop("WAVE_DTYPE", None)          # float64 for the leak-operator SVD
import time, json
import numpy as np
import torch
import wave_twin as wt

t0 = time.time()
DEV = wt.DEV; N = wt.N; NM = wt.NMACRO; MP = wt.MACRO_PX; KP = wt.KP
ACT = wt.DMD_PIX; c = wt._cx; a0 = c - ACT // 2
Z1S = [0.0, 2.0, 5.0, 10.0]
CHUNK = 16

Uin = np.asarray(wt.band_modes(1, KP))          # (4096,120) in-band DC-free
Ube = np.asarray(wt.band_modes(KP + 1, 9))      # (4096,nbe) beyond-band annulus
n_in = Uin.shape[1]; n_be = Ube.shape[1]
Uin_t = torch.tensor(Uin, device=DEV, dtype=torch.float64)
Ube_t = torch.tensor(Ube, device=DEV, dtype=torch.float64)

def block_sum_64(D):
    reg = D[..., a0:a0 + ACT, a0:a0 + ACT].reshape(*D.shape[:-2], NM, MP, NM, MP)
    return reg.sum(dim=(-3, -1))

def render_codes(codes64):
    """dmd_field of a stack of 64x64 code patterns (linear), (M,N,N) complex128."""
    return torch.stack([wt.dmd_field(codes64[i]).to(torch.complex128) for i in range(codes64.shape[0])])

# ---- build stacked leak operators L_be, L_in over z1 ----
ones64 = np.ones((NM, NM))
Lbe_blocks = []; Lin_blocks = []
for z1mm in Z1S:
    z1 = z1mm * 1e-3
    Uc = wt.propagate(wt.dmd_field(ones64).to(torch.complex128), z1)     # (N,N)
    Lbe = torch.zeros(n_be, n_in, device=DEV, dtype=torch.float64)
    Lin = torch.zeros(n_in, n_in, device=DEV, dtype=torch.float64)
    for s in range(0, n_in, CHUNK):
        cols = Uin[:, s:s + CHUNK].T.reshape(-1, NM, NM)                 # (chunk,64,64)
        Us = torch.stack([wt.propagate(wt.dmd_field(cols[k]).to(torch.complex128), z1)
                          for k in range(cols.shape[0])])                # (chunk,N,N)
        D = (Uc[None] * Us.conj()).real                                  # Re[Uc conj(Us)]
        mu = block_sum_64(D).reshape(cols.shape[0], -1)                  # (chunk,4096)
        Lbe[:, s:s + CHUNK] = (mu @ Ube_t).T
        Lin[:, s:s + CHUNK] = (mu @ Uin_t).T
    Lbe_blocks.append(Lbe); Lin_blocks.append(Lin)
    print("[L z1=%4.1f] built  ||Lbe||=%.3e ||Lin||=%.3e [%.0fs]"
          % (z1mm, float(torch.linalg.norm(Lbe)), float(torch.linalg.norm(Lin)), time.time() - t0))

Lbe_stack = torch.cat(Lbe_blocks, 0)             # (4*nbe, n_in)
Lin_stack = torch.cat(Lin_blocks, 0)             # (4*n_in, n_in)

# ---- generalized eigenproblem: leak^2 / signal^2 ----
A = Lbe_stack.T @ Lbe_stack                       # (n_in,n_in)
B = Lin_stack.T @ Lin_stack
wB, UB = torch.linalg.eigh(0.5 * (B + B.T))
wB = torch.clamp(wB, min=wB.max() * 1e-12)
Bisq = UB @ torch.diag(wB.rsqrt()) @ UB.T
Mmat = Bisq @ A @ Bisq
evals, evecs = torch.linalg.eigh(0.5 * (Mmat + Mmat.T))   # ascending
rel_leak = torch.sqrt(torch.clamp(evals, min=0.0))        # relative joint leak per direction
# bottom-64 blind subspace (smallest relative leak)
n_blind = 64
blind_rel_leak = rel_leak[:n_blind]
V_blind = Bisq @ evecs[:, :n_blind]               # code-space directions (n_in, 64), in Uin-coeff space
bottom64_joint_leak = float(blind_rel_leak.max())
print("bottom-64 joint leak (max rel) = %.3e   [median %.3e]"
      % (bottom64_joint_leak, float(blind_rel_leak.median())))

# ---- covariance retention: blind codes vs random codes ----
def make_codes_from_coeffs(Ucoef, M, seed):
    """M signed codes as random combos of the given code-space directions, [-1,1]."""
    rp = np.random.default_rng(seed)
    C = Ucoef.cpu().numpy()                       # (n_in, k)
    P = rp.standard_normal((M, C.shape[1])) @ (Uin @ C).T     # (M,4096)
    P /= np.abs(P).max(axis=1, keepdims=True)
    return P.reshape(M, NM, NM)

def cov_lambda(codes64, z1mm=10.0, z2mm=5.0, n_bank=int(os.environ.get("COV_NBANK", "200")), ridge=1e-4, seed=300):
    z1 = z1mm * 1e-3; z2 = z2mm * 1e-3
    M = codes64.shape[0]
    Ep = torch.stack([wt.propagate(wt.dmd_field((1.0 + codes64[i]) * 0.5).to(torch.complex128), z1) for i in range(M)])
    Em = torch.stack([wt.propagate(wt.dmd_field((1.0 - codes64[i]) * 0.5).to(torch.complex128), z1) for i in range(M)])
    scr = wt.make_screen(50.0, 2 * np.pi, seed=seed); big = scr.shape[0]
    offs = wt.bank_offsets(n_bank, big, 300, seed=seed + 1)
    Gp = torch.empty(n_bank, M, NM, NM, device=DEV, dtype=wt.RDT)
    Gm = torch.empty(n_bank, M, NM, NM, device=DEV, dtype=wt.RDT)
    for tk in range(n_bank):
        ph = wt.screen_crop(scr, int(offs[tk, 0]), int(offs[tk, 1]))
        ep = wt.propagate(Ep * ph[None], z2); em = wt.propagate(Em * ph[None], z2)
        Gp[tk] = block_sum_64(ep.real ** 2 + ep.imag ** 2)
        Gm[tk] = block_sum_64(em.real ** 2 + em.imag ** 2)
    x = wt.witness_scene(4); xn = float(np.linalg.norm(x))
    d2 = wt.beyond_band_delta(0.02, xn, seed=7)
    xt = torch.as_tensor(x, device=DEV, dtype=wt.RDT); dt = torch.as_tensor(x + d2, device=DEV, dtype=wt.RDT)
    fp0 = torch.einsum('bmij,ij->bm', Gp, xt); fm0 = torch.einsum('bmij,ij->bm', Gm, xt)
    fp1 = torch.einsum('bmij,ij->bm', Gp, dt); fm1 = torch.einsum('bmij,ij->bm', Gm, dt)
    f0 = fp0 - fm0; f1 = fp1 - fm1
    scp = wt.PHOT / (0.5 * (fp0 + fm0)).mean().clamp(min=1e-30)
    C0 = torch.cov(f0.T); dC = torch.cov(f1.T) - C0
    shot0 = ((fp0 + fm0).mean(0)) / scp
    V0 = C0 + torch.diag(shot0) + ridge * float(torch.diag(C0).mean()) * torch.eye(M, device=DEV, dtype=wt.RDT)
    A2 = torch.linalg.inv(V0) @ dC
    return 0.5 * float(torch.trace(A2 @ A2))

blind_codes = make_codes_from_coeffs(V_blind, 32, seed=101)
rand_codes = wt.signed_codes(M=32, seed=10)
lam_blind = cov_lambda(blind_codes); print("lam_blind=%.4e [%.0fs]" % (lam_blind, time.time() - t0))
lam_rand = cov_lambda(rand_codes); print("lam_rand =%.4e [%.0fs]" % (lam_rand, time.time() - t0))
lam_ratio = lam_rand / lam_blind if lam_blind > 0 else float('inf')

# also: leak of a random code set vs blind, for reference (bench interpretability)
def meanleak_of_codes(codes64, z1mm=10.0):
    z1 = z1mm * 1e-3; M = codes64.shape[0]
    Uc = wt.propagate(wt.dmd_field(np.ones((NM, NM))).to(torch.complex128), z1)
    Us = torch.stack([wt.propagate(wt.dmd_field(codes64[i]).to(torch.complex128), z1) for i in range(M)])
    D = (Uc[None] * Us.conj()).real
    mu = block_sum_64(D).reshape(M, -1)
    be = torch.linalg.norm(mu @ Ube_t, dim=1); inb = torch.linalg.norm(mu @ Uin_t, dim=1) + 1e-30
    return float((be / inb).mean())

leak_blind = meanleak_of_codes(blind_codes); leak_rand = meanleak_of_codes(rand_codes)

# ---- verdict ----
leak_ok = bottom64_joint_leak <= 1e-4
lam_ok = (lam_ratio <= 3.0) and (lam_ratio >= 1 / 3.0)
if leak_ok and lam_ok:
    verdict = "PASS -- certified-blind code subspace exists (bottom-64 joint mean-leak <=1e-4) AND retains covariance detection power (lambda within 3x of random). Certified-blind codes = next-paper core + bench calibration step."
elif bottom64_joint_leak >= 1e-3:
    verdict = "KILL -- bottom-64 joint leak >= 1e-3: the mean-leak floor is code-independent (no blind subspace). Hardware concentration (IT5) becomes the only wall-restoration route."
elif not lam_ok:
    verdict = "KILL -- blind codes lose >3x covariance detection power (lambda ratio %.2f): blindness costs too much. IT5 hardware route only." % lam_ratio
else:
    verdict = "PARTIAL -- bottom-64 joint leak %.2e (between 1e-4 and 1e-3); report honestly." % bottom64_joint_leak

out = {"test": "IT4_svd_certified_blind_codes",
       "ref": "internal divergence finer #3 (sledgehammer gate); merged adjudication route 2 (code-space)",
       "params": dict(z1_mm=Z1S, n_in=n_in, n_be=n_be, n_blind=n_blind, cov_z2_mm=5.0, cov_n_bank=200),
       "gate": "bottom-64 joint leak <=1e-4 AND lambda within 3x of random = PASS; "
               "leak>=1e-3 or lambda drop>3x = KILL",
       "measured": dict(bottom64_joint_leak_max=bottom64_joint_leak,
                        bottom64_joint_leak_median=float(blind_rel_leak.median()),
                        full_rel_leak_min=float(rel_leak.min()), full_rel_leak_max=float(rel_leak.max()),
                        lam_blind=lam_blind, lam_random=lam_rand, lam_ratio_rand_over_blind=lam_ratio,
                        meanleak_blind_codes=leak_blind, meanleak_random_codes=leak_rand),
       "verdict": verdict, "elapsed_s": round(time.time() - t0, 1)}
with open(os.environ.get("OUT", "IT4_svd_certified_blind_codes.json"), "w") as f:
    json.dump(out, f, indent=2)
print("\nVERDICT:", verdict)
print("saved  elapsed %.1fs" % (time.time() - t0))
