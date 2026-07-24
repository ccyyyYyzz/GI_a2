#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""T4 MULTIPLICATIVE-TO-CONVOLUTIVE boundary.

Prediction 4 (R39 sec 7.5): the thin-screen method fails as the screen becomes
nonlocal / convolutive.  The multiplicative model applies the diffuser phase
POINTWISE at the object plane; the true field applies it at the diffuser plane
and PROPAGATES z2 to the object:
    E_true(r) = prop_z2( screen . E_code_diff )(r)      [convolutive truth]
    E_mult(r) = screen(r) . prop_z2(E_code_diff)(r)     [pointwise model]
They coincide at z2=0 and decorrelate as z2 grows.  We measure the complex-field
correlation and the intensity*object correlation vs z2, and locate the boundary,
mapping the decorrelation to the alpha of the sealed D5 'convolutive_blend' axis
(levels 0.1, 0.25, 0.5), i.e. alpha = 1 - correlation.
"""
import time, json, os
import numpy as np
import torch
import wave_twin as wt

t0 = time.time()
DEV = wt.DEV; c = wt._cx; NM = wt.NMACRO; H = wt.DMD_PIX//2
Z1 = 10e-3
N_EPOCH = int(os.environ.get("N_EPOCH", "60"))

# flat illumination code field at diffuser; a witness object for the intensity*object corr
E_code_diff = wt.propagate(wt.dmd_field(np.ones((NM, NM))), Z1)
x_full = wt.object_full(wt.witness_scene(4))
reg = slice(c-H, c+H)

def corr_complex(a, b):
    """Magnitude of normalized complex correlation over the object region."""
    a = a[reg, reg].reshape(-1); b = b[reg, reg].reshape(-1)
    a = a - a.mean(); b = b - b.mean()
    num = torch.abs((a.conj()*b).sum())
    den = torch.sqrt((a.conj()*a).real.sum()*(b.conj()*b).real.sum()) + 1e-30
    return float(num/den)

def corr_real(a, b):
    a = a[reg, reg].reshape(-1); b = b[reg, reg].reshape(-1)
    a = a - a.mean(); b = b - b.mean()
    return float((a*b).sum()/(torch.sqrt((a*a).sum()*(b*b).sum())+1e-30))

results = {"note": "T4 multiplicative->convolutive. alpha=1-corr maps to D5 convolutive_blend.",
           "z1_mm": 10, "n_epoch": N_EPOCH, "regimes": []}

for tag, sig, l_c in [("developed_2pi_lc15", 2*np.pi, 15.0),
                      ("weak_0.3rad_lc15", 0.3, 15.0),
                      ("developed_2pi_lc5", 2*np.pi, 5.0)]:
    scr = wt.make_screen(l_c, sig, seed=17); big = scr.shape[0]
    offs = wt.bank_offsets(N_EPOCH, big, 300, seed=18)
    rows = []
    for z2mm in [0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]:
        z2 = z2mm*1e-3
        E_code_obj = E_code_diff if z2 == 0 else wt.propagate(E_code_diff, z2)
        fc = []; ic = []
        for t in range(N_EPOCH):
            ph = wt.screen_crop(scr, int(offs[t,0]), int(offs[t,1]))
            E_true = ph*E_code_diff if z2 == 0 else wt.propagate(ph*E_code_diff, z2)
            E_mult = ph*E_code_obj
            fc.append(corr_complex(E_true, E_mult))
            St = (E_true.real**2+E_true.imag**2)*x_full
            Sm = (E_mult.real**2+E_mult.imag**2)*x_full
            ic.append(corr_real(St, Sm))
        fcm = float(np.mean(fc)); icm = float(np.mean(ic))
        rows.append(dict(z2_mm=z2mm, field_corr=round(fcm,4), intensity_obj_corr=round(icm,4),
                         alpha_field=round(1-fcm,4), alpha_intensity=round(1-icm,4)))
    # boundaries: z2 at which alpha crosses D5 levels (interp on intensity*object corr)
    z2s = np.array([r["z2_mm"] for r in rows]); al = np.array([r["alpha_intensity"] for r in rows])
    def z_at(level):
        idx = np.where(al >= level)[0]
        if len(idx)==0: return None
        i = idx[0]
        if i==0: return float(z2s[0])
        x0,x1 = al[i-1], al[i]; z0,z1_ = z2s[i-1], z2s[i]
        return round(float(z0 + (level-x0)/(x1-x0+1e-12)*(z1_-z0)), 3)
    results["regimes"].append(dict(regime=tag, l_c_um=l_c, rows=rows,
        z2mm_alpha_0p10=z_at(0.10), z2mm_alpha_0p25=z_at(0.25), z2mm_alpha_0p50=z_at(0.50)))
    print(f"[{tag:20s}] alpha crossings: 0.10@{z_at(0.10)} 0.25@{z_at(0.25)} 0.50@{z_at(0.50)} mm")
    for r in rows:
        print(f"    z2={r['z2_mm']:5.1f}mm field_corr={r['field_corr']:.3f} "
              f"int*obj_corr={r['intensity_obj_corr']:.3f} alpha={r['alpha_intensity']:.3f}")

with open("T4_MULT_TO_CONV.json","w") as f:
    json.dump(results, f, indent=2)
print(f"\nsaved T4_MULT_TO_CONV.json  elapsed {time.time()-t0:.1f}s")
