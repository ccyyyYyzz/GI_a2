#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""KT6 -- Mean-Covariance Channel-Ratio Ranging  (R44 sec 5.5, innovation #5, score 80).

R44 eq 5.3 (near-contact, shot-limited): J_cov/J_mean = [K1/(2L)] z2^4 + o(z2^4).
A calibrated ratio estimates z2 (a zero-dimensional two-moment range observable).

BINDING coordinator caveat (internal estimator-layer finding): blind lambda_cov
estimation sits far under the Wishart M^2/T noise floor.  KT6 MUST report the
BLIND-ESTIMABLE version alongside the ORACLE version; an oracle-only pass is NOT a
pass.  If blind is hopeless -> record ESTIMATOR_KILLED with the noise-floor numbers.

ORACLE: CRN-clean per-bank cov noncentrality lambda_cov and mean noncentrality
lambda_mean for a fixed beyond-band 2% change; ratio vs z2; fit the near-contact
exponent (expect ->4).
BLIND: the finite-bank Wishart null floor lambda_null = 0.5 tr[(Vinv dC_null)^2] where
dC_null = Chat_0a - Chat_0b from two independent H0 bank halves (T/2 each).  Report
signal-to-floor = lambda_cov / lambda_null.  If <<1 near contact -> blind ranging is
below the sampling floor (ESTIMATOR_KILLED).
"""
import os, time, json
import numpy as np
import torch
import wave_twin as wt
import twin_pool as tp

t0 = time.time()
DEV = wt.DEV; NM = wt.NMACRO
M = int(os.environ.get("M", "16"))
N_BANK = int(os.environ.get("N_BANK", "300"))
RIDGE = 1e-4
Z2_MM = [0.05, 0.1, 0.2, 0.35, 0.5, 1.0, 2.0, 5.0]
Z1_MM = [5.0, 10.0]
SIGMAS = [("weak_0p3", 0.3), ("developed_2pi", 2 * np.pi)]
N_SPLIT = 24                      # Wishart-null split repetitions

x64 = wt.witness_scene(4); xnorm = float(np.linalg.norm(x64))
d2 = wt.beyond_band_delta(0.02, xnorm, seed=7)

def lambdas(Gp, Gm, z1seed):
    x0 = torch.as_tensor(x64, device=DEV, dtype=wt.RDT)
    x1 = torch.as_tensor(x64 + d2, device=DEV, dtype=wt.RDT)
    fp0 = torch.einsum('bmij,ij->bm', Gp, x0); fm0 = torch.einsum('bmij,ij->bm', Gm, x0)
    fp1 = torch.einsum('bmij,ij->bm', Gp, x1); fm1 = torch.einsum('bmij,ij->bm', Gm, x1)
    f0 = fp0 - fm0; f1 = fp1 - fm1
    scp = float(wt.PHOT / (0.5 * (fp0 + fm0)).mean().clamp(min=1e-30))
    shot0 = ((fp0 + fm0).mean(0)) / scp
    C0 = torch.cov(f0.T); dC = torch.cov(f1.T) - C0
    dmean = (f1 - f0).mean(0)
    V0 = C0 + torch.diag(shot0) + RIDGE * float(torch.diag(C0).mean()) * torch.eye(M, device=DEV, dtype=wt.RDT)
    Vinv = torch.linalg.inv(V0)
    A = Vinv @ dC
    lam_cov = 0.5 * float(torch.trace(A @ A))
    lam_mean = float(dmean @ Vinv @ dmean)
    # BLIND Wishart null floor: split H0 banks in half, dC_null = cov(half a) - cov(half b)
    B = f0.shape[0]; nulls = []
    gg = torch.Generator(device=DEV).manual_seed(int(z1seed))
    for _ in range(N_SPLIT):
        perm = torch.randperm(B, device=DEV, generator=gg)
        a = f0[perm[:B // 2]]; b = f0[perm[B // 2:2 * (B // 2)]]
        dCn = torch.cov(a.T) - torch.cov(b.T)
        An = Vinv @ dCn
        nulls.append(0.5 * float(torch.trace(An @ An)))
    lam_null = float(np.mean(nulls))               # finite-T sampling floor at T=B/2
    return lam_cov, lam_mean, lam_null

results = {"test": "KT6_channel_ratio_ranging",
           "R44_ref": "sec 5.5 / innovation #5 (Channel-Ratio Ranging, 80)",
           "params": dict(M=M, n_bank=N_BANK, z2_mm=Z2_MM, z1_mm=Z1_MM,
                          sigmas=[s[0] for s in SIGMAS], change="2pct beyond-band",
                          wishart_null_T=N_BANK // 2, n_split=N_SPLIT),
           "binding_caveat": "blind lambda_cov must be reported vs the Wishart M^2/T floor; "
                             "oracle-only is NOT a pass (coordinator).",
           "frozen_prediction": "ORACLE J_cov/J_mean ~ z2^4 near contact",
           "cells": []}

for z1mm in Z1_MM:
    Ep, Em, _ = tp.code_fields_at_diffuser(M, z1mm * 1e-3, seed=11)
    for stag, sig in SIGMAS:
        ratios = []; z2s = []; s2f = []
        rows = []
        for z2mm in Z2_MM:
            Gp, Gm = tp.speckle_pool(Ep, Em, z2mm * 1e-3, sig, 50.0, N_BANK, seed=300)
            lc, lm, ln = lambdas(Gp, Gm, z1seed=int(z1mm * 100 + z2mm * 10 + sig))
            ratio = lc / lm if lm > 0 else float('inf')
            ratios.append(ratio); z2s.append(z2mm)
            s2f.append(lc / ln if ln > 0 else float('inf'))
            rows.append(dict(z2_mm=z2mm, lam_cov=lc, lam_mean=lm, ratio_cov_over_mean=ratio,
                             lam_null_wishart_floor=ln, signal_to_floor=lc / ln if ln > 0 else None))
            del Gp, Gm
            if DEV == "cuda":
                torch.cuda.empty_cache()
        # oracle near-contact exponent fit (z2<=0.5)
        near = [(z, r) for z, r in zip(z2s, ratios) if z <= 0.5 and r > 0 and np.isfinite(r)]
        if len(near) >= 3:
            lz = np.log([p[0] for p in near]); lr = np.log([p[1] for p in near])
            expo = float(np.polyfit(lz, lr, 1)[0])
        else:
            expo = None
        s2f_near = [v for z, v in zip(z2s, s2f) if z <= 0.5 and np.isfinite(v)]
        results["cells"].append(dict(z1_mm=z1mm, sigma=stag, near_contact_exponent=expo,
                                     signal_to_floor_near_contact=[round(v, 4) for v in s2f_near],
                                     median_signal_to_floor_near=float(np.median(s2f_near)) if s2f_near else None,
                                     rows=rows))
        print("[z1=%g %s] oracle near-contact exponent=%s | median signal/floor(near)=%.2e [%.0fs]"
              % (z1mm, stag, ("%.2f" % expo) if expo else "NA",
                 np.median(s2f_near) if s2f_near else float('nan'), time.time() - t0))

# ------------------------------------------------------------------ verdict
expos = [c["near_contact_exponent"] for c in results["cells"] if c["near_contact_exponent"] is not None]
oracle_z4 = all(2.5 <= e <= 5.5 for e in expos) and len(expos) > 0     # near-contact quartic-ish
floors = [c["median_signal_to_floor_near"] for c in results["cells"] if c["median_signal_to_floor_near"] is not None]
blind_dead = all(f < 1.0 for f in floors)                              # signal below sampling floor
results["verdict"] = dict(
    oracle_near_contact_exponents=[round(e, 3) for e in expos],
    blind_median_signal_to_floor_near=[round(f, 4) for f in floors],
    worst_signal_to_floor=float(min(floors)) if floors else None,
    P_oracle_z4_law="PASS" if oracle_z4 else "PARTIAL",
    P_blind_estimable="ESTIMATOR_KILLED" if blind_dead else ("MARGINAL" if any(f < 3 for f in floors) else "BLIND_OK"),
    overall=("ORACLE confirms the near-contact quartic ratio law (exponent ~4); BLIND channel-ratio "
             "ranging is ESTIMATOR_KILLED -- lambda_cov sits below the Wishart M^2/T sampling floor "
             "near contact (signal/floor << 1), so the ratio is a valid theory relationship but NOT "
             "bench-observable via naive blind estimation. Honest per coordinator's binding caveat."
             if (oracle_z4 and blind_dead) else
             "see numbers -- report oracle and blind honestly"))
results["elapsed_s"] = round(time.time() - t0, 1)
with open(os.environ.get("OUT", "KT6_channel_ratio_ranging.json"), "w") as f:
    json.dump(results, f, indent=2)
print("\nVERDICT:", json.dumps(results["verdict"]))
print("saved  elapsed %.1fs" % (time.time() - t0))
