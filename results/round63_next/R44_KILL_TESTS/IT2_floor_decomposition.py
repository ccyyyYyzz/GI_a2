#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""IT2 -- Floor decomposition + leak-law fit (internal divergence nextpaper #9).

The T1 z1=0 pixelation floor L_pix (~6.8e-3): is it a real DMD INSTRUMENT CONSTANT
(fill-factor + macro-staircase) or a twin grid artifact?  Vary fill (0.92->1.0), the
macro-pixel staircase (mp 8->16), and a relay pupil (3 NAs), all at z1=0; and fit the
leak law leak(z1)=L_pix + a*z1^p on an extended z1 grid.

DECISION: fill/staircase terms move the floor >20% => it DECOMPOSES -> L_pix quotable as a
DMD instrument constant + relay-NA spec; moves <20% => grid artifact -> the transfer
paragraph keeps only the bracket [6.8e-3, 8.9e-2].
"""
import os
os.environ.pop("WAVE_DTYPE", None)
import time, json
import numpy as np
import torch
import wave_twin as wt

t0 = time.time()
DEV = wt.DEV; N = wt.N; NM = wt.NMACRO; KP = wt.KP; DX = wt.DX; PITCH = wt.MIRROR_PITCH
M = int(os.environ.get("M", "12"))
U_in = torch.tensor(wt.band_modes(0, KP), device=DEV, dtype=torch.float64)
U_be = torch.tensor(wt.band_modes(KP + 1, 9), device=DEV, dtype=torch.float64)
codes = wt.signed_codes(M=M, seed=10)
fx = torch.fft.fftfreq(N, d=DX).to(device=DEV, dtype=torch.float64)
FX, FY = torch.meshgrid(fx, fx, indexing="ij")

def render(code64, mp, fill, pupil_kR=None):
    act = NM * mp; a0 = N // 2 - act // 2
    c = torch.as_tensor(code64, device=DEV, dtype=torch.float64)
    held = c.repeat_interleave(mp, 0).repeat_interleave(mp, 1)
    full = torch.zeros(N, N, device=DEV, dtype=torch.float64)
    full[a0:a0 + act, a0:a0 + act] = held
    w = np.sqrt(fill) * PITCH
    env = (torch.sinc(w * FX) * torch.sinc(w * FY)).to(torch.complex128)
    F = torch.fft.fft2(full.to(torch.complex128)) * env
    if pupil_kR is not None:
        F1 = 1.0 / (NM * mp * DX); thr = (pupil_kR + 0.5) * F1
        F = F * ((FX.abs() < thr) & (FY.abs() < thr)).to(torch.complex128)
    return torch.fft.ifft2(F) * np.sqrt(fill)

def block_sum(D, mp):
    act = NM * mp; a0 = N // 2 - act // 2
    reg = D[a0:a0 + act, a0:a0 + act].reshape(NM, mp, NM, mp)
    return reg.sum(dim=(1, 3)).reshape(-1)

def leak(mp, fill, z1, pupil_kR=None):
    vals = []
    for i in range(M):
        ap = render((1.0 + codes[i]) * 0.5, mp, fill, pupil_kR)
        am = render((1.0 - codes[i]) * 0.5, mp, fill, pupil_kR)
        if z1 > 0:
            ap = wt.propagate(ap, z1); am = wt.propagate(am, z1)
        D = (ap.real ** 2 + ap.imag ** 2) - (am.real ** 2 + am.imag ** 2)
        vals.append(block_sum(D, mp))
    J = torch.stack(vals)
    return float(torch.linalg.norm(J @ U_be) / (torch.linalg.norm(J @ U_in) + 1e-30))

results = {"test": "IT2_floor_decomposition", "ref": "internal divergence nextpaper #9",
           "params": dict(M=M, variants="fill{0.92,1.0} x mp{8,16} + pupil kR{5,6,8}"),
           "decision": "fill/staircase move floor >20% => instrument constant; <20% => grid artifact (bracket only)"}

# z1=0 floor decomposition
base = leak(8, 0.92, 0.0)
fill1 = leak(8, 1.0, 0.0)
mp16 = leak(16, 0.92, 0.0)
fill1_mp16 = leak(16, 1.0, 0.0)
pupils = {kR: leak(8, 0.92, 0.0, pupil_kR=kR) for kR in [5, 6, 8]}
results["floor_z1_0"] = dict(baseline_fill0p92_mp8=base, fill1p0_mp8=fill1, fill0p92_mp16=mp16,
                             fill1p0_mp16=fill1_mp16, pupil_kR=pupils)
d_fill = abs(fill1 - base) / base
d_macro = abs(mp16 - base) / base
print("[z1=0] base=%.4e fill1=%.4e(%.0f%%) mp16=%.4e(%.0f%%) pupil5=%.4e [%.0fs]"
      % (base, fill1, 100 * d_fill, mp16, 100 * d_macro, pupils[5], time.time() - t0))

# leak-law fit on extended z1 grid (baseline)
Z1 = [0.0, 0.5, 1.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0]
lk = [leak(8, 0.92, z * 1e-3) for z in Z1]
print("[z1 sweep] " + " ".join("z%g:%.2e" % (z, v) for z, v in zip(Z1, lk)))
# fit leak(z1) = L_pix + a*z1^p on z1>0 (subtract floor)
zpos = np.array([z for z in Z1 if z > 0]); lpos = np.array([v for z, v in zip(Z1, lk) if z > 0])
excess = np.clip(lpos - base, 1e-9, None)
p, logA = np.polyfit(np.log(zpos), np.log(excess), 1)
results["leak_law"] = dict(z1_mm=Z1, leak=lk, L_pix=base, fit_a=float(np.exp(logA)), fit_p=float(p),
                           leak_at_20mm=lk[-1])

# verdict
decomposes = (d_fill > 0.20) or (d_macro > 0.20)
results["verdict"] = dict(
    floor_move_fill=round(d_fill, 3), floor_move_macro=round(d_macro, 3),
    pupil_reduces_floor=pupils[8] < base,
    L_pix=base, bracket=[base, lk[-1]], fit_a=float(np.exp(logA)), fit_p=float(p),
    P_decomposes="PASS" if decomposes else "FAIL",
    overall=("DECOMPOSES -- the z1=0 floor moves %.0f%% (fill) / %.0f%% (macro-staircase): L_pix=%.2e "
             "is a real DMD instrument constant (fill + staircase), quotable with a relay-NA spec "
             "(pupil kR=8 -> %.2e). Leak law leak(z1)=L_pix + %.2e*z1^%.2f."
             % (100 * d_fill, 100 * d_macro, base, pupils[8], float(np.exp(logA)), p)
             if decomposes else
             "GRID ARTIFACT -- floor moves <20%% under fill(%.0f%%)/macro(%.0f%%): keep only the "
             "bracket [%.2e, %.2e] in the transfer paragraph." % (100 * d_fill, 100 * d_macro, base, lk[-1])))
results["elapsed_s"] = round(time.time() - t0, 1)
with open(os.environ.get("OUT", "IT2_floor_decomposition.json"), "w") as f:
    json.dump(results, f, indent=2)
print("VERDICT:", results["verdict"]["overall"])
print("saved  elapsed %.1fs" % (time.time() - t0))
