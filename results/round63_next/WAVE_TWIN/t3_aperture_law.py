#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""T3 APERTURE LAW, REAL CURVE.

Prediction 3 (R39 sec 7.5): covariance recovery follows the aperture B_p (+) B_w
and the usable edge moves with the measured grain size.  We measure the physical
curve: covariance-channel RESPONSE vs scene spatial frequency k, for three grain
sizes (l_c in {5,15,50} um, developed 2pi), and locate the usable edge.

Response(k) = || Cov_banks(b | x0 + eps*grating_k) - Cov_banks(b | x0) ||_F ,
signed buckets b over independent diffuser banks, common-random across k so the
finite-difference is precise.  Prediction: response flat up to k ~ k_p + k_w then
falls; the edge grows as the grain gets finer (k_w larger).
"""
import time, json, os
import numpy as np
import torch
import wave_twin as wt
import twin_pool as tp

t0 = time.time()
DEV = wt.DEV; NM = wt.NMACRO; KP = wt.KP
N_BANK = int(os.environ.get("N_BANK", "300"))
M = int(os.environ.get("M", "24"))
Z1 = 10e-3; Z2 = 3e-3                                   # near-contact working point (speckle exists)
EPS = 0.05                                              # grating perturbation amplitude (rel)

# Accessible grain-spanning diffusers (grain is illumination-diffraction-limited to
# the fine regime in this geometry, so k_w stays near the 64-grid Nyquist; the coarsest
# reachable grain just begins to close the aperture).  (l_c_um, z2_mm) -> grain.
GRAIN_SET = [(15.0, 1.0), (60.0, 1.0), (120.0, 3.0)]   # ~5, ~48, ~96 um grains
Ep, Em, _ = tp.code_fields_at_diffuser(M, Z1, seed=10)
x0 = np.full((NM, NM), 0.5)                             # uniform base object
x0v = torch.tensor(x0, device=DEV, dtype=wt.RDT)

tp_xx, tp_yy = np.meshgrid(np.arange(NM), np.arange(NM), indexing="ij")

def grating64(k, seed=0):
    """Unit-norm cosine grating at Chebyshev freq k on the 64 grid (random orientation)."""
    rng = np.random.default_rng(seed+k)
    # pick an integer (kx,ky) with max(|kx|,|ky|)=k
    if k == 0:
        g = np.ones((NM, NM))
    else:
        opts = [(kx,ky) for kx in range(-k,k+1) for ky in range(-k,k+1)
                if max(abs(kx),abs(ky))==k]
        kx,ky = opts[rng.integers(len(opts))]
        ph = 2*np.pi*(kx*tp_xx + ky*tp_yy)/NM
        g = np.cos(ph)
    g = g/ (np.linalg.norm(g)+1e-30)
    return g

def cov_response_curve(l_c, z2, seed=200):
    Gp, Gm = tp.speckle_pool(Ep, Em, z2, 2*np.pi, l_c, N_BANK, seed=seed)
    # common scp from x0
    b0, scp = tp.signed_buckets(Gp, Gm, x0, shot=True,
                                gen=torch.Generator(device=DEV).manual_seed(1))
    C0 = torch.cov(b0.T)
    xnorm = float(np.linalg.norm(x0))
    ks = list(range(0, 29))
    resp = []
    for k in ks:
        dk = grating64(k, seed=7)*(EPS*xnorm)           # perturbation, ||dk||=eps||x0||
        xk = x0 + dk
        bk, _ = tp.signed_buckets(Gp, Gm, xk, shot=True, scp=scp,
                                  gen=torch.Generator(device=DEV).manual_seed(1))
        Ck = torch.cov(bk.T)
        resp.append(float(torch.linalg.norm(Ck - C0)))
    resp = np.array(resp)
    return ks, resp

# measure grain (reuse T2 convention quickly) for edge annotation
def grain_um(l_c, z2, seed=5, n=60):
    scr = wt.make_screen(l_c, 2*np.pi, seed=seed); big=scr.shape[0]
    offs = wt.bank_offsets(n, big, 300, seed=seed+1)
    c=wt._cx; H=wt.DMD_PIX//2
    acc=torch.zeros(2*H,2*H,device=DEV,dtype=wt.RDT)
    Eflat=wt.propagate(wt.dmd_field(np.ones((NM,NM))), Z1)
    for t in range(n):
        ph=wt.screen_crop(scr,int(offs[t,0]),int(offs[t,1]))
        e=wt.propagate(Eflat*ph, z2); S=(e.real**2+e.imag**2)[c-H:c+H,c-H:c+H]
        dS=S-S.mean(); F=torch.fft.fft2(dS.to(wt.CDT))
        acc+=torch.fft.fftshift(torch.fft.ifft2(F*F.conj()).real)
    ac=(acc/n).cpu().numpy(); n2=ac.shape[0]; cc=n2//2
    yy,xx=np.mgrid[0:n2,0:n2]; r=np.sqrt((xx-cc)**2+(yy-cc)**2).astype(int)
    prof=np.bincount(r.ravel(),ac.ravel())/np.maximum(np.bincount(r.ravel()),1); prof/=prof[0]
    below=np.where(prof<0.5)[0]; hw=below[0] if len(below) else n2
    return float(2*hw*wt.DX*1e6)

results={"note":"T3 aperture law: cov response vs scene freq for 3 grain sizes; edge ~ min(k_p+k_w, Nyquist=32). "
         "grain is illumination-diffraction-limited to the fine regime here, so the aperture stays wide.",
         "kp":KP,"z1_mm":10,"eps":EPS,"n_bank":N_BANK,"M":M,"nyquist":NM//2,"curves":[]}
L_OBJ_UM = NM*wt.MACRO_PX*wt.DX*1e6
for l_c, z2mm in GRAIN_SET:
    z2 = z2mm*1e-3
    ks, resp = cov_response_curve(l_c, z2, seed=200+int(l_c))
    wg = grain_um(l_c, z2)
    k_w = min(L_OBJ_UM/max(wg,1e-6), NM/2)
    edge_pred = min(KP + k_w, NM/2)
    # measured usable edge: last k where response > 20% of the low-k plateau
    plateau = resp[1:6].mean()
    thr = 0.2*plateau
    usable = [k for k,r in zip(ks,resp) if r>thr]
    edge_meas = max(usable) if usable else 0
    results["curves"].append(dict(l_c_um=l_c, z2_mm=z2mm, grain_um=round(wg,1),
        k_w_mapped=round(float(k_w),2), edge_pred_kp_plus_kw=round(float(edge_pred),2),
        edge_measured=edge_meas, ks=ks, response=[round(float(r),6) for r in resp],
        plateau=round(float(plateau),6)))
    print(f"[l_c={l_c:5.0f}um z2={z2mm}mm grain={wg:5.1f}um k_w~{k_w:5.1f}] "
          f"edge_pred={edge_pred:5.1f} edge_meas={edge_meas}")

with open("T3_APERTURE_LAW.json","w") as f:
    json.dump(results, f, indent=2)
print(f"\nsaved T3_APERTURE_LAW.json  elapsed {time.time()-t0:.1f}s")
