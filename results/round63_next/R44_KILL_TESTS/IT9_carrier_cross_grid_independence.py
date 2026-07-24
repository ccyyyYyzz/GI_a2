#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""IT9 -- Carrier-cross exactness, grid-independence confirm (internal finer #2 residual).

Before freezing the WAVE_TWIN_REPORT sec1 correction ("no 2k_p intensity-autocorrelation
edge; the leak is a carrier-code cross term D = Re[C conj(S)]"), confirm grid-independence.
The complementary difference is ALGEBRAICALLY exact:
    |E^+|^2-|E^-|^2 = |(C+S)/2|^2-|(C-S)/2|^2 = Re[C conj(S)]  (the |S|^2 term CANCELS).
Verify numerically at M=32 on the N=2048 grid AND a reimplemented N=4096 grid that
(a) D == Re[C conj(S)] to machine zero (|S|^2 autocorrelation absent -> no 2k_p edge),
(b) the beyond-2k_p leak (carrier pedestal) is the SAME at both grids (not an aliasing
    artifact).  Confirms => correction frozen with confidence; deviates => flag grid sensitivity.
"""
import os
os.environ.pop("WAVE_DTYPE", None)
import time, json
import numpy as np
import torch
import wave_twin as wt

t0 = time.time()
DEV = wt.DEV; NM = wt.NMACRO; MP = wt.MACRO_PX; KP = wt.KP; ACT = wt.DMD_PIX
DX = wt.DX; FILL = wt.FILL; PITCH = wt.MIRROR_PITCH
M = int(os.environ.get("M", "32"))
GRIDS = [2048, 4096]

U_in = torch.tensor(wt.band_modes(0, 2 * KP), device=DEV, dtype=torch.float64)     # in-band+core (<=2kp)
U_bey = torch.tensor(wt.band_modes(2 * KP + 1, 31), device=DEV, dtype=torch.float64)  # beyond 2kp
codes = wt.signed_codes(M=M, seed=10)

def render_at(code64, N):
    a0 = N // 2 - ACT // 2
    c = torch.as_tensor(code64, device=DEV, dtype=torch.float64)
    held = c.repeat_interleave(MP, 0).repeat_interleave(MP, 1)     # 512x512
    full = torch.zeros(N, N, device=DEV, dtype=torch.float64)
    full[a0:a0 + ACT, a0:a0 + ACT] = held
    fx = torch.fft.fftfreq(N, d=DX).to(device=DEV, dtype=torch.float64)
    FX, FY = torch.meshgrid(fx, fx, indexing="ij")
    w = np.sqrt(FILL) * PITCH
    env = (torch.sinc(w * FX) * torch.sinc(w * FY)).to(torch.complex128)
    F = torch.fft.fft2(full.to(torch.complex128))
    return torch.fft.ifft2(F * env) * np.sqrt(FILL)

def block_sum_64(D, N):
    a0 = N // 2 - ACT // 2
    reg = D[a0:a0 + ACT, a0:a0 + ACT].reshape(NM, MP, NM, MP)
    return reg.sum(dim=(1, 3)).reshape(-1)

results = {"test": "IT9_carrier_cross_grid_independence", "ref": "internal divergence finer #2 residual",
           "params": dict(M=M, grids=GRIDS),
           "check": "D == Re[C conj(S)] (|S|^2 cancels, no 2kp edge) + beyond-2kp leak grid-independent",
           "grids": {}}

Cones = {N: render_at(np.ones((NM, NM)), N) for N in GRIDS}    # common field per grid
for N in GRIDS:
    C = Cones[N]
    cc_err = []; leak = []
    for i in range(M):
        s = codes[i]
        Ep = render_at((1.0 + s) * 0.5, N); Em = render_at((1.0 - s) * 0.5, N)
        D = (Ep.real ** 2 + Ep.imag ** 2) - (Em.real ** 2 + Em.imag ** 2)
        S = render_at(s, N)
        carrier_cross = (C * S.conj()).real                       # Re[C conj(S)]
        cc_err.append(float(torch.linalg.norm(D - carrier_cross) / (torch.linalg.norm(D) + 1e-30)))
        mu = block_sum_64(D, N)
        leak.append(float(torch.linalg.norm(mu @ U_bey) / (torch.linalg.norm(mu @ U_in) + 1e-30)))
    results["grids"]["N%d" % N] = dict(carrier_cross_rel_err_max=max(cc_err),
                                       carrier_cross_rel_err_mean=float(np.mean(cc_err)),
                                       beyond_2kp_leak_mean=float(np.mean(leak)),
                                       beyond_2kp_leak_std=float(np.std(leak)))
    print("[N=%d] carrier-cross rel-err max=%.2e | beyond-2kp leak mean=%.4e std=%.2e [%.0fs]"
          % (N, max(cc_err), float(np.mean(leak)), float(np.std(leak)), time.time() - t0))

# ---- verdict ----
cc_exact = all(g["carrier_cross_rel_err_max"] < 1e-10 for g in results["grids"].values())
leak2048 = results["grids"]["N2048"]["beyond_2kp_leak_mean"]
leak4096 = results["grids"]["N4096"]["beyond_2kp_leak_mean"]
grid_indep = abs(leak2048 - leak4096) / (0.5 * (leak2048 + leak4096) + 1e-30) < 0.10
results["verdict"] = dict(
    carrier_cross_exact_both_grids=cc_exact,
    beyond_2kp_leak_2048=leak2048, beyond_2kp_leak_4096=leak4096,
    grid_rel_diff=abs(leak2048 - leak4096) / (0.5 * (leak2048 + leak4096) + 1e-30),
    P_carrier_cross_confirmed="PASS" if cc_exact else "FAIL",
    P_grid_independent="PASS" if grid_indep else "PARTIAL",
    overall=("CONFIRMED -- D == Re[C conj(S)] to machine zero on both grids (|S|^2 cancels, "
             "NO 2k_p edge) and the beyond-2kp carrier-pedestal leak agrees to %.1f%% between "
             "N=2048 and N=4096. The WAVE_TWIN_REPORT sec1 carrier-cross correction can be "
             "frozen with confidence." % (100 * abs(leak2048 - leak4096) / (0.5 * (leak2048 + leak4096) + 1e-30))
             if (cc_exact and grid_indep) else
             "DEVIATES -- flag grid sensitivity before quoting annulus profiles (see numbers)"))
results["elapsed_s"] = round(time.time() - t0, 1)
with open(os.environ.get("OUT", "IT9_carrier_cross_grid_independence.json"), "w") as f:
    json.dump(results, f, indent=2)
print("VERDICT:", results["verdict"]["overall"])
print("saved  elapsed %.1fs" % (time.time() - t0))
