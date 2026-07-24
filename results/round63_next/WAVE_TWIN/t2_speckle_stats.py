#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""T2 SPECKLE STATISTICS -> the declared-law dictionary.

At the object plane, measure speckle GRAIN SIZE (intensity-autocorrelation width)
and CONTRAST (std/mean) versus correlation length l_c in {5,15,50} um and
z2 in {0,1,5,20} mm, for the developed (2*pi) and weak (0.3 rad) diffuser.

Map each geometry to the statistical twin's (k_w, sigma_f):
  * k_w  (medium band cut-off, 64-grid Chebyshev) from grain size w_g:
        k_w ~ L_obj / w_g   (L_obj = 64*32um object extent), capped at 32 (Nyquist),
        and k_w/k_p vs the statistical levels kwf in {1,2,4}.
  * sigma_f (medium RMS contrast) from the intensity contrast; statistical
        SIG_EFF = {0.298,0.503,0.696} for sf {0.3,0.6,1.0}; developed speckle
        saturates at contrast ~ 1.0.
"""
import time, json, os
import numpy as np
import torch
import wave_twin as wt

t0 = time.time()
DEV = wt.DEV; c = wt._cx
N_EPOCH = int(os.environ.get("N_EPOCH", "150"))
Z1 = 10e-3
HALF = wt.DMD_PIX // 2                                  # central region half-width (px)
L_OBJ_UM = wt.NMACRO * wt.MACRO_PX * wt.DX * 1e6        # 2048 um
KP = wt.KP

# one representative code (flat-ish illumination fills the region); use uniform amplitude
a_flat = np.ones((wt.NMACRO, wt.NMACRO))
E_flat = wt.propagate(wt.dmd_field(a_flat), Z1)         # field at diffuser

def radial_fwhm(ac2d):
    """FWHM (in pixels) of a centered 2D autocorrelation via radial average."""
    n = ac2d.shape[0]; cc = n//2
    yy, xx = np.mgrid[0:n, 0:n]
    r = np.sqrt((xx-cc)**2 + (yy-cc)**2).astype(int)
    prof = np.bincount(r.ravel(), ac2d.ravel())/np.maximum(np.bincount(r.ravel()),1)
    prof = prof/prof[0]
    below = np.where(prof < 0.5)[0]
    if len(below)==0: return float(n)
    i = below[0]
    # linear interp between i-1 and i for half-max crossing
    if i==0: return 0.0
    x0,x1 = prof[i-1], prof[i]
    frac = (x0-0.5)/(x0-x1+1e-12)
    hwhm = (i-1)+frac
    return 2.0*hwhm

def measure(l_c, sigma_phi, z2, n_epoch=N_EPOCH, seed=3):
    scr = wt.make_screen(l_c, sigma_phi, seed=seed); big=scr.shape[0]
    offs = wt.bank_offsets(n_epoch, big, 300, seed=seed+1)
    # per-pixel temporal mean/var for contrast; averaged spatial autocorr for grain
    acc_mean = torch.zeros(2*HALF, 2*HALF, device=DEV, dtype=wt.RDT)
    acc_sq = torch.zeros_like(acc_mean)
    acc_ac = torch.zeros_like(acc_mean)
    for t in range(n_epoch):
        ph = wt.screen_crop(scr, int(offs[t,0]), int(offs[t,1]))
        e = E_flat*ph
        if z2>0: e = wt.propagate(e, z2)
        S = (e.real**2+e.imag**2)[c-HALF:c+HALF, c-HALF:c+HALF]
        acc_mean += S; acc_sq += S*S
        dS = S - S.mean()
        F = torch.fft.fft2(dS.to(wt.CDT))
        ac = torch.fft.ifft2(F*F.conj()).real
        acc_ac += torch.fft.fftshift(ac)
    mean_t = acc_mean/n_epoch
    var_t = acc_sq/n_epoch - mean_t**2
    contrast = float((var_t.clamp(min=0).sqrt().mean())/(mean_t.mean()+1e-30))
    ac = (acc_ac/n_epoch).cpu().numpy()
    fwhm_px = radial_fwhm(ac)
    w_g_um = fwhm_px*wt.DX*1e6
    return contrast, w_g_um

def map_kw(w_g_um):
    if w_g_um <= 1e-6: return None, None
    k_w = L_OBJ_UM/w_g_um
    k_w_capped = min(k_w, wt.NMACRO/2)
    return k_w, k_w_capped/KP

rows = []
SIG_EFF = {0.3:0.298, 0.6:0.503, 1.0:0.696}
for tag, sig in [("developed_2pi", 2*np.pi), ("weak_0.3rad", 0.3)]:
    for l_c in [5, 15, 50]:
        for z2mm in [0.0, 1.0, 5.0, 20.0]:
            z2 = z2mm*1e-3
            contrast, w_g = measure(l_c, sig, z2)
            k_w, kw_over_kp = map_kw(w_g)
            # nearest statistical kwf level
            kwf_level = None
            if kw_over_kp is not None:
                kwf_level = min([1,2,4], key=lambda kv: abs(kv-kw_over_kp))
            # nearest statistical sigma_f
            sf_level = min(SIG_EFF, key=lambda sv: abs(SIG_EFF[sv]-min(contrast,0.696)))
            row = dict(regime=tag, l_c_um=l_c, z2_mm=z2mm,
                       contrast=round(contrast,4), grain_um=round(w_g,2),
                       mapped_k_w=None if k_w is None else round(k_w,2),
                       mapped_kw_over_kp=None if kw_over_kp is None else round(kw_over_kp,3),
                       nearest_kwf_level=kwf_level,
                       mapped_sigma_f_contrast=round(contrast,4),
                       nearest_stat_sf=sf_level)
            rows.append(row)
            print(f"[{tag:13s} l_c={l_c:3d} z2={z2mm:4.1f}] contrast={contrast:.3f} "
                  f"grain={w_g:6.1f}um  k_w~{None if k_w is None else round(k_w,1)} "
                  f"kw/kp~{None if kw_over_kp is None else round(kw_over_kp,2)} (kwf~{kwf_level})")

out = dict(note="T2 speckle dictionary. maps wave geometry -> statistical (k_w,sigma_f).",
           L_obj_um=L_OBJ_UM, kp=KP, n_epoch=N_EPOCH, stat_SIG_EFF=SIG_EFF,
           stat_kwf_levels=[1,2,4], rows=rows)
with open("T2_SPECKLE_STATS.json","w") as f:
    json.dump(out, f, indent=2)
print(f"\nsaved T2_SPECKLE_STATS.json  elapsed {time.time()-t0:.1f}s")
