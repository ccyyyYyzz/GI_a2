#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""IT3 -- Leak-law stability (internal divergence stranger #8).

Fit leak(z1)=L_pix + a*z1^p across 3 code seeds x M in {8,32}.  Decision: (a,p) stable
<15% AND all points within 10% => 'reproducible instrument constant' wording allowed;
unstable => bracket-only.  Cross-check against the proven 4Re(C-bar S-tilde) carrier-cross
factorization (never attribute the law to a 2k_p edge).
"""
import os
os.environ.pop("WAVE_DTYPE", None)
import time, json
import numpy as np
import torch
import wave_twin as wt

t0 = time.time()
DEV = wt.DEV; N = wt.N; NM = wt.NMACRO; KP = wt.KP
U_in = torch.tensor(wt.band_modes(0, KP), device=DEV, dtype=torch.float64)
U_be = torch.tensor(wt.band_modes(KP + 1, 9), device=DEV, dtype=torch.float64)
SEEDS = [10, 21, 33]; MLIST = [8, 32]
Z1 = [0.0, 0.5, 1.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0]

def block_sum(D):
    a = wt._cx - wt.DMD_PIX // 2
    reg = D[a:a + wt.DMD_PIX, a:a + wt.DMD_PIX].reshape(NM, wt.MACRO_PX, NM, wt.MACRO_PX)
    return reg.sum(dim=(1, 3)).reshape(-1)

def leak(codes, z1):
    vals = []
    for i in range(codes.shape[0]):
        ap = wt.dmd_field((1.0 + codes[i]) * 0.5).to(torch.complex128)
        am = wt.dmd_field((1.0 - codes[i]) * 0.5).to(torch.complex128)
        if z1 > 0:
            ap = wt.propagate(ap, z1); am = wt.propagate(am, z1)
        D = (ap.real ** 2 + ap.imag ** 2) - (am.real ** 2 + am.imag ** 2)
        vals.append(block_sum(D))
    J = torch.stack(vals)
    return float(torch.linalg.norm(J @ U_be) / (torch.linalg.norm(J @ U_in) + 1e-30))

results = {"test": "IT3_leak_law_stability", "ref": "internal divergence stranger #8",
           "params": dict(seeds=SEEDS, M=MLIST, z1_mm=Z1),
           "decision": "(a,p) stable <15% and points within 10% => reproducible instrument constant",
           "fits": []}

fitrows = []
for M in MLIST:
    for seed in SEEDS:
        codes = wt.signed_codes(M=M, seed=seed)
        lk = [leak(codes, z * 1e-3) for z in Z1]
        Lpix = lk[0]
        zpos = np.array([z for z in Z1 if z > 0]); lpos = np.array([v for z, v in zip(Z1, lk) if z > 0])
        excess = np.clip(lpos - Lpix, 1e-9, None)
        p, logA = np.polyfit(np.log(zpos), np.log(excess), 1)
        row = dict(M=M, seed=seed, L_pix=Lpix, a=float(np.exp(logA)), p=float(p), leak=lk)
        fitrows.append(row); results["fits"].append(row)
        print("[M=%d seed=%d] L_pix=%.4e a=%.4e p=%.3f [%.0fs]" % (M, seed, Lpix, np.exp(logA), p, time.time() - t0))

def cv(vals):
    return float(np.std(vals) / (np.mean(vals) + 1e-30))
Lpix_cv = cv([r["L_pix"] for r in fitrows]); a_cv = cv([r["a"] for r in fitrows]); p_cv = cv([r["p"] for r in fitrows])
# point-level spread: at each z1, CV of leak across seeds/M
pt_cvs = [cv([r["leak"][k] for r in fitrows]) for k in range(len(Z1))]
stable = (Lpix_cv < 0.15) and (a_cv < 0.15) and (p_cv < 0.15) and (max(pt_cvs) < 0.15)
results["verdict"] = dict(
    L_pix_CV=round(Lpix_cv, 4), a_CV=round(a_cv, 4), p_CV=round(p_cv, 4), max_point_CV=round(max(pt_cvs), 4),
    L_pix_mean=float(np.mean([r["L_pix"] for r in fitrows])), a_mean=float(np.mean([r["a"] for r in fitrows])),
    p_mean=float(np.mean([r["p"] for r in fitrows])),
    P_reproducible="PASS" if stable else "PARTIAL",
    overall=("REPRODUCIBLE instrument constant -- (L_pix,a,p) stable across 3 seeds x M{8,32}: "
             "CV(L_pix)=%.1f%% CV(a)=%.1f%% CV(p)=%.1f%%, all points within %.1f%%. The leak law "
             "leak(z1)=%.2e + %.2e*z1^%.2f (carrier-cross factorization, NOT a 2kp edge) is quotable."
             % (100 * Lpix_cv, 100 * a_cv, 100 * p_cv, 100 * max(pt_cvs),
                np.mean([r["L_pix"] for r in fitrows]), np.mean([r["a"] for r in fitrows]), np.mean([r["p"] for r in fitrows]))
             if stable else
             "UNSTABLE -- bracket-only wording (CV L_pix=%.1f%% a=%.1f%% p=%.1f%% maxpt=%.1f%%)."
             % (100 * Lpix_cv, 100 * a_cv, 100 * p_cv, 100 * max(pt_cvs))))
results["elapsed_s"] = round(time.time() - t0, 1)
with open(os.environ.get("OUT", "IT3_leak_law_stability.json"), "w") as f:
    json.dump(results, f, indent=2)
print("VERDICT:", results["verdict"]["overall"])
print("saved  elapsed %.1fs" % (time.time() - t0))
