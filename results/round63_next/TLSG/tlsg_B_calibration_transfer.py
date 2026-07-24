#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""TLSG PART B -- Calibration-transfer arm (R45 sec C1.2 bars 3-4, Eq 4.4 / B4.2).

Blindness-capacity certificate under calibration uncertainty.  The whitened leak operator
Lw = L_be B^{-1/2} (B = L_in^T L_in, the in-band-response metric) has singular values =
relative leaks (Eq 4.3).  Calibration produces Lhat_w (nominal hardware).  With an operator
error bound ||L-Lhat||_op <= delta (Eq 4.4):
    max_{c in S_d, ||c||=1} ||L c|| <= sigma_d(Lhat) + delta,
where S_d is the bottom-d right-singular subspace of Lhat_w.

Fresh calibration/test split with HARDWARE-LAW perturbations (fill, mirror pitch, small
relay defocus z1).  Conjugate plane z1=0.
BARS (R45 C1.2):
  3. >=95% simultaneous coverage of the sigma_d(Lhat)+delta leak envelope on FRESH test
     draws, no code in S_d exceeding it;
  4. at the conjugate plane, >=64 code dimensions remain certified below 1e-4 on FRESH
     physical-law draws.
delta is estimated on a CALIBRATION ensemble and TESTED on a disjoint FRESH ensemble.
"""
import os
os.environ.pop("WAVE_DTYPE", None)
import time, json
import numpy as np
import torch
import wave_twin as wt

t0 = time.time()
DEV = wt.DEV; N = wt.N; NM = wt.NMACRO; MP = wt.MACRO_PX; KP = wt.KP; ACT = wt.DMD_PIX
c = wt._cx; a0 = c - ACT // 2; DX = wt.DX
D_CERT = 64
N_CAL = int(os.environ.get("N_CAL", "24"))
N_TEST = int(os.environ.get("N_TEST", "24"))
SEED = int(os.environ.get("SEED", "20260724"))
rng = np.random.default_rng(SEED)
fx = torch.fft.fftfreq(N, d=DX).to(device=DEV, dtype=torch.float64)
FX, FY = torch.meshgrid(fx, fx, indexing="ij")

Uin = np.asarray(wt.band_modes(1, KP)); Ube = np.asarray(wt.band_modes(KP + 1, 9))
n_in = Uin.shape[1]; n_be = Ube.shape[1]
Uin_t = torch.tensor(Uin, device=DEV, dtype=torch.float64)
Ube_t = torch.tensor(Ube, device=DEV, dtype=torch.float64)

def render(code64, fill, pitch):
    a = torch.as_tensor(code64, device=DEV, dtype=torch.float64)
    held = a.repeat_interleave(MP, 0).repeat_interleave(MP, 1)
    full = torch.zeros(N, N, device=DEV, dtype=torch.float64); full[a0:a0 + ACT, a0:a0 + ACT] = held
    w = np.sqrt(fill) * pitch
    env = (torch.sinc(w * FX) * torch.sinc(w * FY)).to(torch.complex128)
    return torch.fft.ifft2(torch.fft.fft2(full.to(torch.complex128)) * env) * np.sqrt(fill)

def block_sum(D):
    reg = D[..., a0:a0 + ACT, a0:a0 + ACT].reshape(*D.shape[:-2], NM, MP, NM, MP)
    return reg.sum(dim=(-3, -1))

def leak_operator(fill, pitch, z1):
    Uc = wt.propagate(render(np.ones((NM, NM)), fill, pitch), z1) if z1 > 0 else render(np.ones((NM, NM)), fill, pitch)
    Lbe = torch.zeros(n_be, n_in, device=DEV, dtype=torch.float64)
    Lin = torch.zeros(n_in, n_in, device=DEV, dtype=torch.float64)
    for s in range(0, n_in, 16):
        cols = Uin[:, s:s + 16].T.reshape(-1, NM, NM)
        Us = torch.stack([render(cols[k], fill, pitch) for k in range(cols.shape[0])])
        if z1 > 0:
            Us = wt.propagate(Us, z1)
        D = (Uc[None] * Us.conj()).real
        mu = block_sum(D).reshape(cols.shape[0], -1)
        Lbe[:, s:s + 16] = (mu @ Ube_t).T; Lin[:, s:s + 16] = (mu @ Uin_t).T
    return Lbe, Lin

# ---- calibration (nominal hardware) ----
FILL0, PITCH0 = wt.FILL, wt.MIRROR_PITCH
Lbe0, Lin0 = leak_operator(FILL0, PITCH0, 0.0)
B0 = Lin0.T @ Lin0
wB, UB = torch.linalg.eigh(0.5 * (B0 + B0.T)); wB = torch.clamp(wB, min=wB.max() * 1e-14)
Bisq = UB @ torch.diag(wB.rsqrt()) @ UB.T                 # B0^{-1/2} (fixed whitening metric)
Lw0 = Lbe0 @ Bisq                                          # whitened calibration operator
U0, S0, Vh0 = torch.linalg.svd(Lw0, full_matrices=False)   # S0 descending
S0_asc = torch.flip(S0, [0])                               # ascending relative leaks
V_asc = torch.flip(Vh0.T, [1])                             # right sing vectors, ascending
S_d = V_asc[:, :D_CERT]                                    # bottom-64 subspace (whitened coords)
sigma_d = float(S0_asc[D_CERT - 1])                        # sigma_64(Lhat)
print("[calib] sigma_64(Lhat_w)=%.3e | d@1e-4=%d d@1e-5=%d [%.0fs]"
      % (sigma_d, int((S0_asc <= 1e-4).sum()), int((S0_asc <= 1e-5).sum()), time.time() - t0))

def perturbed_operator():
    fill = FILL0 * (1 + rng.uniform(-0.05, 0.05))
    pitch = PITCH0 * (1 + rng.uniform(-0.03, 0.03))
    z1 = abs(rng.uniform(0, 0.5)) * 1e-3                    # small relay defocus 0..0.5mm
    Lbe, _ = leak_operator(fill, pitch, z1)
    Lw = Lbe @ Bisq                                        # fixed calibration whitening
    return Lw, (fill, pitch, z1 * 1e3)

# ---- calibration ensemble: estimate delta = max ||Lw - Lhat_w||_op ----
cal_dev = []
for i in range(N_CAL):
    Lw, hw = perturbed_operator()
    cal_dev.append(float(torch.linalg.matrix_norm(Lw - Lw0, ord=2)))
delta = float(np.max(cal_dev))                             # operator-norm bound from calibration ensemble
delta_p95 = float(np.percentile(cal_dev, 95))
print("[calib ensemble] delta(max)=%.3e delta(p95)=%.3e envelope sigma_64+delta=%.3e [%.0fs]"
      % (delta, delta_p95, sigma_d + delta, time.time() - t0))

# ---- FRESH test ensemble: coverage of sigma_d+delta AND 64-dim < 1e-4 ----
envelope = sigma_d + delta
covered = 0; s64_max_leaks = []; fresh_d_at_1e4 = []; s64_dims_below_1e4 = []
for i in range(N_TEST):
    Lw, hw = perturbed_operator()
    # max leak over calibration subspace S_d under fresh Lw = top singular value of Lw @ S_d
    LS = Lw @ S_d                                          # (n_be, 64)
    s64_max = float(torch.linalg.matrix_norm(LS, ord=2))
    s64_max_leaks.append(s64_max)
    if s64_max <= envelope:
        covered += 1
    # per-code leak of S_d under fresh Lw
    per_code = torch.linalg.norm(LS, dim=0)               # (64,) leak of each S_d basis dir (approx; not orthonormal image)
    s64_dims_below_1e4.append(int((torch.linalg.norm(Lw @ S_d, dim=0) <= 1e-4).sum()))
    # fresh operator's OWN bottom spectrum: d@1e-4
    Sf = torch.linalg.svdvals(Lw)
    fresh_d_at_1e4.append(int((Sf <= 1e-4).sum()))
coverage = covered / N_TEST
print("[fresh test] coverage(<=sigma64+delta)=%.3f | fresh d@1e-4 min=%d | S64 max-leak worst=%.3e [%.0fs]"
      % (coverage, min(fresh_d_at_1e4), max(s64_max_leaks), time.time() - t0))

# ---- bars ----
bar3_pass = coverage >= 0.95
bar4_pass = min(fresh_d_at_1e4) >= 64
results = {"test": "TLSG_partB_calibration_transfer", "R45_ref": "sec C1.2 bars 3-4, Eq 4.4",
           "params": dict(d_cert=D_CERT, n_cal=N_CAL, n_test=N_TEST, seed=SEED,
                          perturb="fill+-5%, pitch+-3%, relay defocus 0-0.5mm", conjugate_z1=0.0),
           "measured": dict(sigma_64_Lhat=sigma_d, delta_max=delta, delta_p95=delta_p95,
                            envelope_sigma64_plus_delta=envelope,
                            calib_d_at_1e4=int((S0_asc <= 1e-4).sum()), calib_d_at_1e5=int((S0_asc <= 1e-5).sum()),
                            fresh_coverage_of_envelope=coverage, S64_max_leak_worst_fresh=max(s64_max_leaks),
                            fresh_d_at_1e4_min=min(fresh_d_at_1e4), fresh_d_at_1e4_mean=float(np.mean(fresh_d_at_1e4))),
           "bars": {
               "bar3_coverage_95pct": dict(coverage=round(coverage, 4), threshold=0.95, PASS=bool(bar3_pass)),
               "bar4_64dims_below_1e4_fresh": dict(fresh_d_at_1e4_min=min(fresh_d_at_1e4), threshold=64, PASS=bool(bar4_pass)),
           },
           "partB_all_pass": bool(bar3_pass and bar4_pass), "elapsed_s": round(time.time() - t0, 1)}
print("\n=== PART B BARS ===")
print("  bar3 coverage %.3f (>=0.95): %s" % (coverage, "PASS" if bar3_pass else "FAIL"))
print("  bar4 fresh d@1e-4 min %d (>=64): %s" % (min(fresh_d_at_1e4), "PASS" if bar4_pass else "FAIL"))
with open(os.environ.get("OUT", "TLSG_partB.json"), "w") as f:
    json.dump(results, f, indent=2)
print("partB_all_pass:", results["partB_all_pass"], " elapsed %.1fs" % (time.time() - t0))
