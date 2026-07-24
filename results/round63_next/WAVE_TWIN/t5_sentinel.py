#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""T5 END-TO-END SENTINEL.

Run the frozen covariance detector (sealed_common conventions: signed buckets,
per-bank sample covariance, matched score W = V0^-1 dC V0^-1, noncentrality
d'^2 = T_eff * lambda) on TWIN-generated banks.  H0 (x) vs H1 (x + beyond-band
change at 2% and 5% energy) at near-contact (z2=1mm) and mid (z2=5mm).

Per-bank noncentrality  lambda = 1/2 * tr[(V0^-1 dC)^2]   (Gaussian covariance
discrimination; equals the sealed c^T J_B c convention).  T_det(d'=5) = 25/lambda.
Compare to the sealed statistical twin: best cell 453 banks @ 2%, MC power LCB 0.99,
D3 target TPR 0.988.  AUC by bootstrap over T_eff-bank sample covariances.
Verdict: TWIN_CONFIRMS / TWIN_DEVIATES(where).
"""
import time, json, os
import numpy as np
import torch
import wave_twin as wt
import twin_pool as tp

t0 = time.time()
DEV = wt.DEV; NM = wt.NMACRO
M = int(os.environ.get("M", "32"))
N_BANK = int(os.environ.get("N_BANK", "1000"))
Z1 = 10e-3
SEALED_BEST_TDET = 453.0
SEALED_MC_LCB = 0.99

Ep, Em, _ = tp.code_fields_at_diffuser(M, Z1, seed=10)
x64 = wt.witness_scene(4)
xnorm = float(np.linalg.norm(x64))
d2 = wt.beyond_band_delta(0.02, xnorm, seed=7)
d5 = wt.beyond_band_delta(0.05, xnorm, seed=7)

def clean_deltaC_and_V0(Gp, Gm, x0, x1, scp, ridge=1e-4):
    """CRN-clean covariance change + null covariance (medium + shot).

    Noiseless signed fluxes for x0 and x1 on the SAME diffuser banks -> the common
    medium fluctuation cancels in dC = cov(f1)-cov(f0) (variance reduction).  Null
    covariance V0 = cov(f0) + diag(shot var), shot var = (flux^+ + flux^-)/scp."""
    x0t = torch.as_tensor(x0, device=DEV, dtype=wt.RDT)
    x1t = torch.as_tensor(x1, device=DEV, dtype=wt.RDT)
    fp0 = torch.einsum('bmij,ij->bm', Gp, x0t); fm0 = torch.einsum('bmij,ij->bm', Gm, x0t)
    fp1 = torch.einsum('bmij,ij->bm', Gp, x1t); fm1 = torch.einsum('bmij,ij->bm', Gm, x1t)
    f0 = fp0 - fm0; f1 = fp1 - fm1                       # noiseless signed fluxes
    C0m = torch.cov(f0.T); C1m = torch.cov(f1.T)
    dC = C1m - C0m                                       # CRN-clean covariance change
    dmean = (f1 - f0).mean(0)                            # CRN-clean MEAN-channel shift (leak)
    shot0 = ((fp0 + fm0).mean(0))/scp                    # per-code mean signed-bucket shot var
    V0 = C0m + torch.diag(shot0) + ridge*float(torch.diag(C0m).mean())*torch.eye(M, device=DEV, dtype=wt.RDT)
    Vinv = torch.linalg.inv(V0)
    A = Vinv @ dC
    lam = 0.5*float(torch.trace(A @ A))                  # per-bank COV noncentrality (sealed c^T J_B c)
    mean_lam = float(dmean @ Vinv @ dmean)              # per-bank MEAN-channel noncentrality (the leak)
    W = Vinv @ dC @ Vinv
    return lam, W, dC, mean_lam

def bootstrap_auc(b0, b1, W, T_eff, n_mc=1500, seed=0):
    """Bootstrap AUC of the matched covariance statistic at T_eff banks."""
    g = torch.Generator(device=DEV).manual_seed(seed)
    n0 = b0.shape[0]; n1 = b1.shape[0]
    wv = W
    def stats(b):
        n = b.shape[0]
        idx = torch.randint(0, n, (n_mc, T_eff), device=DEV, generator=g)
        B = b[idx]                                     # (n_mc,T_eff,M)
        Bc = B - B.mean(1, keepdim=True)
        C = torch.einsum('kti,ktj->kij', Bc, Bc)/(T_eff-1)
        return torch.einsum('ij,kij->k', wv, C)        # (n_mc,)
    t0s = stats(b0); t1s = stats(b1)
    # empirical AUC = P(t1 > t0)
    auc = float((t1s[:,None] > t0s[None,:]).float().mean())
    dprime = float((t1s.mean()-t0s.mean())/(torch.sqrt(0.5*(t0s.var()+t1s.var()))+1e-30))
    return auc, dprime

results = {"note": "T5 end-to-end sentinel on wave-twin banks vs sealed detector.",
           "sealed_best_Tdet_2pct": SEALED_BEST_TDET, "sealed_mc_lcb": SEALED_MC_LCB,
           "sealed_D3_target_tpr": 0.988, "M": M, "n_bank": N_BANK, "geometries": []}

# developed speckle; coarse grain (l_c=50) ~ sealed best-cell kwf~1; also fine grain check
for geo_tag, z2mm, l_c in [("near_contact_z2_1mm_lc50", 1.0, 50.0),
                            ("mid_z2_5mm_lc50", 5.0, 50.0),
                            ("near_contact_z2_1mm_lc15", 1.0, 15.0)]:
    z2 = z2mm*1e-3
    Gp, Gm = tp.speckle_pool(Ep, Em, z2, 2*np.pi, l_c, N_BANK, seed=300)
    g = torch.Generator(device=DEV).manual_seed(11)
    b0, scp = tp.signed_buckets(Gp, Gm, x64, shot=True, gen=g)
    b2, _ = tp.signed_buckets(Gp, Gm, x64 + d2, shot=True, scp=scp, gen=g)
    b5, _ = tp.signed_buckets(Gp, Gm, x64 + d5, shot=True, scp=scp, gen=g)
    rows = []
    for lbl, bx, xd in [("2pct", b2, x64 + d2), ("5pct", b5, x64 + d5)]:
        lam, W, dC, mean_lam = clean_deltaC_and_V0(Gp, Gm, x64, xd, scp)
        T_det = float(25.0/lam) if lam > 0 else float('inf')
        # AUC at sealed T_det (453) and at own T_det (capped to pool)
        Teval = min(int(SEALED_BEST_TDET), N_BANK-1)
        auc453, dp453 = bootstrap_auc(b0, bx, W, Teval, seed=1)
        Town = min(max(int(T_det), 2), N_BANK-1)
        aucown, dpown = bootstrap_auc(b0, bx, W, Town, seed=2)
        mean_over_cov = float(np.sqrt(mean_lam/lam)) if lam > 0 else float('inf')  # mean d' / cov d'
        rows.append(dict(change=lbl, per_bank_cov_noncentrality=round(lam,6),
                         per_bank_mean_noncentrality=round(mean_lam,8),
                         mean_dp_over_cov_dp=round(mean_over_cov,4),
                         T_det_dprime5=round(T_det,1),
                         auc_at_T453=round(auc453,4), dprime_at_T453=round(dp453,3),
                         auc_at_Town=round(aucown,4), T_own=Town))
        print(f"[{geo_tag:26s} {lbl}] cov_lam={lam:.4e} T_det={T_det:7.1f} "
              f"AUC@453={auc453:.4f} d'@453={dp453:.2f} mean/cov_dp={mean_over_cov:.2e}")
    results["geometries"].append(dict(geometry=geo_tag, z2_mm=z2mm, l_c_um=l_c, rows=rows))
    del Gp, Gm; torch.cuda.empty_cache()

# verdict: compare 2% near-contact T_det to sealed 453 and AUC to ~1.0
nc = results["geometries"][0]["rows"][0]     # near-contact 2pct
tdet = nc["T_det_dprime5"]; auc = nc["auc_at_T453"]
if tdet < 2000 and auc > 0.95:
    verdict = f"TWIN_CONFIRMS (near-contact 2% T_det={tdet} vs sealed 453; AUC@453={auc})"
elif auc > 0.9:
    verdict = f"TWIN_CONFIRMS_QUALITATIVELY (T_det={tdet} slower than sealed 453; AUC@453={auc})"
else:
    verdict = f"TWIN_DEVIATES (near-contact 2% AUC@453={auc}, T_det={tdet})"
results["verdict"] = verdict
print("\nVERDICT:", verdict)
with open("T5_SENTINEL.json","w") as f:
    json.dump(results, f, indent=2)
print(f"saved T5_SENTINEL.json  elapsed {time.time()-t0:.1f}s")
