#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""IT6 -- Segmented-bucket P^2 law (internal divergence stranger #7).

Replace the single integrating bucket by P spatial segments (quadrant/multi-cell diode).
The covariance channel then carries the segment x segment structure, so per-bank detection
noncentrality should grow ~P^2.  Re-reduce the SAME twin speckle kernels (no re-propagation)
into P in {1,4,16} segments at z2=5mm and near-contact z2=1mm, change eps in {2,5}%.

MANDATORY (from the F-pixel kill): per-segment beyond-band MEAN-leak audit -- segmentation
must not flood the witness annulus with mean leak (mean d' must stay below covariance d').
SAMPLE FLOOR: the P*M covariance needs n_bank > P*M; M=16 keeps P=16 -> 256 dim well sampled.

PASS = scaling exponent >=1.5 AND P=4 buys >=6x cov noncentrality AND per-segment mean leak
below the covariance d'.  KILL = exponent<1.5, or P=4 gain<6x, or the mean leak floods the
witness.
"""
import os, time, json
import numpy as np
import torch
import wave_twin as wt
import twin_pool as tp

t0 = time.time()
DEV = wt.DEV; NM = wt.NMACRO
M = int(os.environ.get("M", "16"))
N_BANK = int(os.environ.get("N_BANK", "700"))
Z1 = 10e-3; RIDGE = 1e-4
PLIST = [1, 4, 16]
CELLS = [("mid_z2_5mm", 5.0), ("near_contact_z2_1mm", 1.0)]

x64 = wt.witness_scene(4); xnorm = float(np.linalg.norm(x64))
d2 = wt.beyond_band_delta(0.02, xnorm, seed=7)
d5 = wt.beyond_band_delta(0.05, xnorm, seed=7)
Ep, Em, _ = tp.code_fields_at_diffuser(M, Z1, seed=11)

def seg_buckets(G, xf):
    """segmented signed buckets for all P in PLIST. G:(B,M,64,64), xf:(64,64) -> dict P->(B,M*P)."""
    Gx = G * xf[None, None]
    out = {}
    for P in PLIST:
        q = int(round(P ** 0.5)); s = NM // q
        r = Gx.reshape(N_BANK, M, q, s, q, s).sum(dim=(3, 5))     # (B,M,q,q)
        out[P] = r.reshape(N_BANK, M * P)
    return out

def lam_mean_for(Gp, Gm, xd_field, scp):
    x0 = torch.as_tensor(x64, device=DEV, dtype=wt.RDT); x1 = torch.as_tensor(xd_field, device=DEV, dtype=wt.RDT)
    G = Gp - Gm
    f0 = seg_buckets(G, x0); f1 = seg_buckets(G, x1)
    # shot: per-segment flux = seg of (Gp+Gm) x
    Gp_seg0 = seg_buckets(Gp, x0); Gm_seg0 = seg_buckets(Gm, x0)
    res = {}
    for P in PLIST:
        d = M * P
        C0 = torch.cov(f0[P].T); dC = torch.cov(f1[P].T) - C0
        dmean = (f1[P] - f0[P]).mean(0)
        shot0 = ((Gp_seg0[P] + Gm_seg0[P]).mean(0)) / scp
        V0 = C0 + torch.diag(shot0) + RIDGE * float(torch.diag(C0).mean()) * torch.eye(d, device=DEV, dtype=wt.RDT)
        Vinv = torch.linalg.inv(V0); A = Vinv @ dC
        lam = 0.5 * float(torch.trace(A @ A))
        mean_lam = float(dmean @ Vinv @ dmean)
        res[P] = (lam, mean_lam)
    return res

results = {"test": "IT6_segmented_bucket_p2", "ref": "internal divergence stranger #7",
           "params": dict(M=M, n_bank=N_BANK, P=PLIST, cells=[c[0] for c in CELLS], z1_mm=10,
                          sample_floor_note="P*M<=256 < n_bank=%d (well sampled)" % N_BANK),
           "decision": "exponent>=1.5 AND P=4 buys>=6x AND per-segment mean d' below cov d' = PASS",
           "cells": []}

for tag, z2mm in CELLS:
    Gp, Gm = tp.speckle_pool(Ep, Em, z2mm * 1e-3, 2 * np.pi, 50.0, N_BANK, seed=300)
    g = torch.Generator(device=DEV).manual_seed(11)
    _, scp = tp.signed_buckets(Gp, Gm, x64, shot=True, gen=g)
    rows = []
    for lbl, xd in [("2pct", x64 + d2), ("5pct", x64 + d5)]:
        res = lam_mean_for(Gp, Gm, xd, scp)
        lam1 = res[1][0]; lam4 = res[4][0]; lam16 = res[16][0]
        expo = float(np.log(lam16 / lam1) / np.log(16.0)) if lam1 > 0 and lam16 > 0 else None
        buy4 = lam4 / lam1 if lam1 > 0 else None
        mean_over_cov = {P: (float(np.sqrt(res[P][1] / res[P][0])) if res[P][0] > 0 else None) for P in PLIST}
        rows.append(dict(change=lbl, cov_lam={str(P): res[P][0] for P in PLIST},
                         mean_lam={str(P): res[P][1] for P in PLIST},
                         mean_dp_over_cov_dp={str(P): mean_over_cov[P] for P in PLIST},
                         scaling_exponent_P=round(expo, 3) if expo else None,
                         P4_buys=round(buy4, 3) if buy4 else None))
        print("[%s %s] cov_lam P1=%.3e P4=%.3e P16=%.3e (exp=%.2f, P4buys=%.2f) mean/cov@P16=%s [%.0fs]"
              % (tag, lbl, lam1, lam4, lam16, expo if expo else float('nan'), buy4 if buy4 else float('nan'),
                 ("%.2f" % mean_over_cov[16]) if mean_over_cov[16] else "NA", time.time() - t0))
    results["cells"].append(dict(geometry=tag, z2_mm=z2mm, rows=rows))
    del Gp, Gm
    if DEV == "cuda":
        torch.cuda.empty_cache()

# ------------------------------------------------------------------ verdict
allrows = [r for cell in results["cells"] for r in cell["rows"]]
expos = [r["scaling_exponent_P"] for r in allrows if r["scaling_exponent_P"] is not None]
buys = [r["P4_buys"] for r in allrows if r["P4_buys"] is not None]
leak_ok = all((r["mean_dp_over_cov_dp"]["16"] is not None and r["mean_dp_over_cov_dp"]["16"] < 1.0) for r in allrows)
exp_ok = all(e >= 1.5 for e in expos) if expos else False
buy_ok = all(b >= 6.0 for b in buys) if buys else False
if exp_ok and buy_ok and leak_ok:
    verdict = ("PASS -- segmented bucket follows a super-linear P law (exponents %s, P=4 buys %s), "
               "per-segment mean leak stays below covariance d' -> quadrant/multi-cell diode enters "
               "the bench spec and the P^2 law seeds the next paper." % ([round(e,2) for e in expos], [round(b,1) for b in buys]))
elif not leak_ok:
    verdict = "KILL -- segmentation floods the witness annulus with mean leak (per-segment mean d' >= cov d' at P=16)."
elif not exp_ok:
    verdict = "KILL -- scaling exponent <1.5 (%s): segmentation does not deliver the P^2 gain." % [round(e,2) for e in expos]
elif not buy_ok:
    verdict = "KILL -- P=4 buys <6x (%s): quadrant diode not worth the complexity." % [round(b,2) for b in buys]
else:
    verdict = "PARTIAL -- report honestly."
results["verdict"] = dict(exponents=[round(e, 3) for e in expos], P4_buys=[round(b, 3) for b in buys],
                          per_segment_leak_below_cov=leak_ok, overall=verdict)
results["elapsed_s"] = round(time.time() - t0, 1)
with open(os.environ.get("OUT", "IT6_segmented_bucket_p2.json"), "w") as f:
    json.dump(results, f, indent=2)
print("VERDICT:", verdict)
print("saved  elapsed %.1fs" % (time.time() - t0))
