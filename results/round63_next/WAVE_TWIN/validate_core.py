#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Core-physics validation for the wave twin: sampling, energy conservation,
speckle contrast (developed->1), and grain-size vs correlation-length trend."""
import time, json
import numpy as np
import torch
import wave_twin as wt

t0 = time.time()
print("=== sampling ===")
print(json.dumps(wt.check_sampling(), indent=2))

# ---- energy conservation of band-limited ASM over the used distances ----
print("\n=== propagation energy conservation ===")
codes = wt.signed_codes(M=1, seed=10)
a = (1.0 + codes[0]) * 0.5                       # nonneg amplitude in [0,1]
e0 = wt.dmd_field(a)
E_in = float((e0.real**2 + e0.imag**2).sum())
for z_mm in [5, 10, 20, 40]:
    e = wt.propagate(e0, z_mm*1e-3)
    E = float((e.real**2 + e.imag**2).sum())
    print(f"  z={z_mm:2d} mm  energy ratio = {E/E_in:.6f}")

# ---- speckle contrast: developed (>=2pi) -> 1 ; weak (0.3 rad) -> small ----
print("\n=== speckle contrast at object plane (z1=10mm, z2=5mm) ===")
z1, z2 = 10e-3, 5e-3
ediff = wt.propagate(e0, z1)                      # code field at diffuser
for tag, sig in [("weak_0.3rad", 0.3), ("developed_2pi", 2*np.pi), ("developed_4pi", 4*np.pi)]:
    for lc in [5, 15, 50]:
        scr = wt.make_screen(lc, sig, seed=1)
        offs = wt.bank_offsets(64, scr.shape[0], step_px=300, seed=2)
        Is = []
        for k in range(64):
            ph = wt.screen_crop(scr, int(offs[k,0]), int(offs[k,1]))
            S = wt.object_plane_intensity(ediff, ph, z2)
            Is.append(S)
        Is = torch.stack(Is)                     # (64,N,N)
        # measure in the central illuminated region only
        c = wt._cx; h = wt.DMD_PIX//2
        reg = Is[:, c-h:c+h, c-h:c+h]
        mean_t = reg.mean(0)
        std_t = reg.std(0)
        m = mean_t.mean().item()
        contrast = (std_t.mean().item())/ (m + 1e-30)
        print(f"  {tag:14s} l_c={lc:3d}um  contrast(std/mean over epochs) ~ {contrast:.3f}")

print(f"\nelapsed {time.time()-t0:.1f}s")
