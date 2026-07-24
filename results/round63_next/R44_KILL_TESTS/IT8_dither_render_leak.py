#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""IT8 -- Dither-render leak (internal divergence finer #5).

A real DMD renders a grayscale amplitude by BINARY dithering at the mirror pitch (each
macro-pixel's amplitude is an on/off mirror pattern).  Ordered-dither (Bayer) the
complementary amplitudes a^+=(1+s)/2, a^-=(1-s)/2 into binary sub-macro patterns and
measure the beyond-band mean leak vs the analog (grayscale) render, at z1=0 and 10mm.

DECISION: dither leak >=1e-2 => PWM temporal dithering + exact whole-period bucket
integration becomes a hard line in the apparatus spec; <1e-3 => dithering is a non-issue.
"""
import os
os.environ.pop("WAVE_DTYPE", None)
import time, json
import numpy as np
import torch
import wave_twin as wt

t0 = time.time()
DEV = wt.DEV; N = wt.N; NM = wt.NMACRO; MP = wt.MACRO_PX; KP = wt.KP; c = wt._cx; ACT = wt.DMD_PIX
FILL = wt.FILL; PITCH = wt.MIRROR_PITCH; DX = wt.DX
M = int(os.environ.get("M", "16"))
a0 = c - ACT // 2
U_in = torch.tensor(wt.band_modes(0, KP), device=DEV, dtype=torch.float64)
U_be = torch.tensor(wt.band_modes(KP + 1, 9), device=DEV, dtype=torch.float64)
fx = torch.fft.fftfreq(N, d=DX).to(device=DEV, dtype=torch.float64)
FX, FY = torch.meshgrid(fx, fx, indexing="ij")
w = np.sqrt(FILL) * PITCH
ENV = (torch.sinc(w * FX) * torch.sinc(w * FY)).to(torch.complex128)
codes = wt.signed_codes(M=M, seed=10)

def bayer(n):
    if n == 1:
        return np.array([[0.0]])
    m = bayer(n // 2)
    return np.block([[4 * m, 4 * m + 2], [4 * m + 3, 4 * m + 1]])
B8 = torch.tensor((bayer(MP) + 0.5) / (MP * MP), device=DEV, dtype=torch.float64)   # thresholds in (0,1)
BTILE = B8.repeat(NM, NM)                                                            # 512x512

def analog_field(amp64):
    a = torch.as_tensor(amp64, device=DEV, dtype=torch.float64)
    held = a.repeat_interleave(MP, 0).repeat_interleave(MP, 1)
    full = torch.zeros(N, N, device=DEV, dtype=torch.float64); full[a0:a0 + ACT, a0:a0 + ACT] = held
    return torch.fft.ifft2(torch.fft.fft2(full.to(torch.complex128)) * ENV) * np.sqrt(FILL)

def dither_field(amp64):
    a = torch.as_tensor(amp64, device=DEV, dtype=torch.float64).clamp(0, 1)
    held = a.repeat_interleave(MP, 0).repeat_interleave(MP, 1)
    binary = (held > BTILE).to(torch.float64)                                        # ordered-dither binary
    full = torch.zeros(N, N, device=DEV, dtype=torch.float64); full[a0:a0 + ACT, a0:a0 + ACT] = binary
    return torch.fft.ifft2(torch.fft.fft2(full.to(torch.complex128)) * ENV) * np.sqrt(FILL)

def block_sum(D):
    reg = D[a0:a0 + ACT, a0:a0 + ACT].reshape(NM, MP, NM, MP)
    return reg.sum(dim=(1, 3)).reshape(-1)

def leak(render, z1):
    vals = []
    for i in range(M):
        ap = render((1.0 + codes[i]) * 0.5); am = render((1.0 - codes[i]) * 0.5)
        if z1 > 0:
            ap = wt.propagate(ap, z1); am = wt.propagate(am, z1)
        D = (ap.real ** 2 + ap.imag ** 2) - (am.real ** 2 + am.imag ** 2)
        vals.append(block_sum(D))
    J = torch.stack(vals)
    return float(torch.linalg.norm(J @ U_be) / (torch.linalg.norm(J @ U_in) + 1e-30))

results = {"test": "IT8_dither_render_leak", "ref": "internal divergence finer #5",
           "params": dict(M=M, dither="Bayer %dx%d ordered" % (MP, MP), z1_mm=[0.0, 10.0]),
           "decision": "dither leak >=1e-2 => PWM+whole-period integration hard line; <1e-3 => non-issue",
           "leak": {}}

for tag, ren in [("analog", analog_field), ("dither", dither_field)]:
    results["leak"][tag] = {"z1_0": leak(ren, 0.0), "z1_10": leak(ren, 10e-3)}
    print("[%s] leak z1=0: %.4e  z1=10mm: %.4e [%.0fs]"
          % (tag, results["leak"][tag]["z1_0"], results["leak"][tag]["z1_10"], time.time() - t0))

dither_excess0 = abs(results["leak"]["dither"]["z1_0"] - results["leak"]["analog"]["z1_0"])
dither_excess10 = abs(results["leak"]["dither"]["z1_10"] - results["leak"]["analog"]["z1_10"])
dmax = max(dither_excess0, dither_excess10, results["leak"]["dither"]["z1_0"])
results["verdict"] = dict(
    dither_leak_z1_0=results["leak"]["dither"]["z1_0"], analog_leak_z1_0=results["leak"]["analog"]["z1_0"],
    dither_excess_z1_0=dither_excess0, dither_excess_z1_10=dither_excess10,
    P_dither=("PWM_CONSTRAINT_NEEDED" if results["leak"]["dither"]["z1_0"] >= 1e-2
              else ("NON_ISSUE" if results["leak"]["dither"]["z1_0"] < 1e-3 else "MARGINAL")),
    overall=("PWM CONSTRAINT NEEDED -- ordered-dither render leaks %.2e into the beyond-band annulus "
             "at z1=0 (>=1e-2): PWM temporal dithering + exact whole-period bucket integration is a "
             "hard line in the apparatus spec." % results["leak"]["dither"]["z1_0"]
             if results["leak"]["dither"]["z1_0"] >= 1e-2 else
             ("DITHERING NON-ISSUE -- dither leak %.2e < 1e-3 at z1=0; the PWM constraint can be dropped."
              % results["leak"]["dither"]["z1_0"] if results["leak"]["dither"]["z1_0"] < 1e-3 else
              "MARGINAL -- dither leak %.2e in [1e-3,1e-2]; one numbered sentence in the bench doc, PWM advisable."
              % results["leak"]["dither"]["z1_0"])))
results["elapsed_s"] = round(time.time() - t0, 1)
with open(os.environ.get("OUT", "IT8_dither_render_leak.json"), "w") as f:
    json.dump(results, f, indent=2)
print("VERDICT:", results["verdict"]["overall"])
print("saved  elapsed %.1fs" % (time.time() - t0))
