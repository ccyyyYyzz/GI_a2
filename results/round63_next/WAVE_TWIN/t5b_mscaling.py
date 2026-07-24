#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""T5b — M-code-count scaling + matched-cell comparison (comparability corrections).

Coordinator requirement 1: the T5 run used M=32 codes vs the sealed probe's M=128.
Test how the per-bank covariance noncentrality lambda scales with M by building ONE
M=128 speckle pool per geometry and evaluating lambda on code-subsets (first 32/64/128)
of the SAME pool (no re-propagation).  If lambda saturates near the in-band DOF
((k_p+1)^2-1 = 35, SCRAMBLE M_eff~13-15) then M=32 is already representative; if it
grows ~M^2 the sealed 453 must be re-projected.  Also reports the MEAN-channel lambda
scaling (mean ~ M vs cov ~ M_eff) that governs the T1 leak ratio.

Requirement 2: map each geometry via the T2 dictionary (contrast->sigma_f, grain->k_w
class kwf) to the MATCHED sealed cell and quote that cell's sealed T_det (not best cell).
"""
import time, json, os
import numpy as np
import torch
import wave_twin as wt
import twin_pool as tp

t0 = time.time()
DEV = wt.DEV; NM = wt.NMACRO
M_MAX = 128
N_BANK = int(os.environ.get("N_BANK", "400"))
Z1 = 10e-3

Ep, Em, _ = tp.code_fields_at_diffuser(M_MAX, Z1, seed=10)
x64 = wt.witness_scene(4); xnorm = float(np.linalg.norm(x64))
d2 = wt.beyond_band_delta(0.02, xnorm, seed=7)
d5 = wt.beyond_band_delta(0.05, xnorm, seed=7)

def lam_for_subset(Gp, Gm, xd, scp, Msub, ridge=1e-4):
    gp = Gp[:, :Msub]; gm = Gm[:, :Msub]
    x0t = torch.tensor(x64, device=DEV, dtype=wt.RDT); x1t = torch.tensor(xd, device=DEV, dtype=wt.RDT)
    fp0 = torch.einsum('bmij,ij->bm', gp, x0t); fm0 = torch.einsum('bmij,ij->bm', gm, x0t)
    fp1 = torch.einsum('bmij,ij->bm', gp, x1t); fm1 = torch.einsum('bmij,ij->bm', gm, x1t)
    f0 = fp0-fm0; f1 = fp1-fm1
    C0 = torch.cov(f0.T); dC = torch.cov(f1.T) - C0
    dmean = (f1-f0).mean(0)
    shot0 = ((fp0+fm0).mean(0))/scp
    V0 = C0 + torch.diag(shot0) + ridge*float(torch.diag(C0).mean())*torch.eye(Msub, device=DEV, dtype=wt.RDT)
    Vinv = torch.linalg.inv(V0); A = Vinv@dC
    lam = 0.5*float(torch.trace(A@A)); mean_lam = float(dmean@Vinv@dmean)
    return lam, mean_lam

# sealed analytic 2% map (subset) for matched-cell lookup
SEALED_2PCT = {  # (sf, shape, kwf, claim) -> T_det
    (1.0,'flat',4,1.25):1879,(1.0,'flat',4,1.5):2157,(1.0,'flat',4,1.8):2996,
    (1.0,'k^-1',4,1.25):1014,(1.0,'k^-1',4,1.5):1242,(1.0,'k^-1',4,1.8):1879,
    (1.0,'k^-2',4,1.25):513,(1.0,'k^-2',4,1.5):713,(1.0,'k^-2',4,1.8):1214,
    (0.6,'flat',4,1.25):2614,(0.6,'k^-2',4,1.25):606,
}

results = {"note":"T5b M-scaling + matched-cell. sealed best cell 453 is the BEST not the matched cell.",
           "n_bank":N_BANK,"M_values":[32,64,128],
           "in_band_dof_kp":(wt.KP+1)**2-1,"geometries":[]}

for geo_tag, z2mm, l_c, contrast, grain_um in [
        ("mid_z2_5mm_lc50", 5.0, 50.0, 1.01, 18.6),
        ("near_contact_z2_1mm_lc50", 1.0, 50.0, 0.87, 20.6)]:
    z2 = z2mm*1e-3
    Gp, Gm = tp.speckle_pool(Ep, Em, z2, 2*np.pi, l_c, N_BANK, seed=300)
    g = torch.Generator(device=DEV).manual_seed(11)
    _, scp = tp.signed_buckets(Gp, Gm, x64, shot=True, gen=g)
    rows = []
    for lbl, xd in [("2pct", x64+d2), ("5pct", x64+d5)]:
        scal = {}
        for Msub in [32, 64, 128]:
            lam, mlam = lam_for_subset(Gp, Gm, xd, scp, Msub)
            scal[Msub] = dict(cov_lam=lam, mean_lam=mlam,
                              T_det=round(25.0/lam,1) if lam>0 else None,
                              mean_over_cov_dp=round(float(np.sqrt(mlam/lam)),3) if lam>0 else None)
        rows.append(dict(change=lbl, scaling=scal,
                         lam_ratio_128_over_32=round(scal[128]["cov_lam"]/scal[32]["cov_lam"],3)))
        print(f"[{geo_tag} {lbl}] cov_lam M=32:{scal[32]['cov_lam']:.3e} "
              f"64:{scal[64]['cov_lam']:.3e} 128:{scal[128]['cov_lam']:.3e} "
              f"(x{scal[128]['cov_lam']/scal[32]['cov_lam']:.2f})  Tdet@128={scal[128]['T_det']}")
    # matched sealed cell: contrast->sf, grain->kwf (all fine grains => kwf=4)
    sf_class = min([0.3,0.6,1.0], key=lambda s:abs({0.3:0.298,0.6:0.503,1.0:0.696}[s]-min(contrast,0.696)))
    matched = {k:v for k,v in SEALED_2PCT.items() if k[0]==sf_class and k[2]==4}
    results["geometries"].append(dict(geometry=geo_tag, z2_mm=z2mm, l_c_um=l_c,
        measured_contrast=contrast, measured_grain_um=grain_um, mapped_sf_class=sf_class,
        mapped_kwf_class=4, matched_sealed_cells_2pct_Tdet=matched, rows=rows))
    del Gp, Gm; torch.cuda.empty_cache()

with open("T5b_MSCALING.json","w") as f:
    json.dump(results, f, indent=2)
print(f"\nsaved T5b_MSCALING.json  elapsed {time.time()-t0:.1f}s")
