#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""KT5 -- Etendue / N_eff collapse (R44 sec 5.4, innovation #4, score 100).

R44 eqs 5.1-5.2: speckle-contrast C_b^2 = Var(b)/E[b]^2 = 1/N_eff (effective number of
independently integrated grains).  With n photoelectrons and a Gaussian covariance
channel, visibility r=n/N_eff and I_logQ = 1/2 sum_j [r_j/(1+r_j)]^2 -> photons cease
to help for n>>N_eff.  The physical resource is the pair (etendue, independent banks).

TEST: sweep the illumination aperture diameter (>=4 diameters) at fixed screen
statistics; compute N_eff directly from the intensity-bucket contrast; sweep photon
count n; test whether the detection curves T_det(n) COLLAPSE against n/N_eff.

N_eff = E[flux]^2 / Var(flux) over banks (noiseless + exposure).  T_det(n) = 25/lambda(n),
lambda(n)=0.5 tr[(V0(n)^{-1} dC)^2] with shot variance ~ flux/n folded into V0(n).
PREDICTION: N_eff grows with aperture area; T_det(n) curves collapse onto one function
of n/N_eff.
"""
import os, time, json
import numpy as np
import torch
import wave_twin as wt

t0 = time.time()
DEV = wt.DEV; N = wt.N; NM = wt.NMACRO; c = wt._cx
ACT = wt.DMD_PIX; a0 = c - ACT // 2
M = int(os.environ.get("M", "16"))
N_BANK = int(os.environ.get("N_BANK", "300"))
Z1 = 10e-3; Z2 = 5e-3; RIDGE = 1e-4
APER = [0.40, 0.55, 0.70, 0.85, 1.00]      # illumination aperture diameter (fraction of ACT)
NPHOT = [1e3, 3e3, 1e4, 3e4, 1e5]

# disk window of diameter frac*ACT centered on the object region
yy, xx = torch.meshgrid(torch.arange(N, device=DEV, dtype=torch.float64),
                        torch.arange(N, device=DEV, dtype=torch.float64), indexing="ij")
rr = torch.sqrt((xx - c) ** 2 + (yy - c) ** 2)

def disk(frac):
    return (rr <= 0.5 * frac * ACT).to(torch.complex128)

x64 = wt.witness_scene(4); xnorm = float(np.linalg.norm(x64))
d5 = wt.beyond_band_delta(0.05, xnorm, seed=7)
S = wt.signed_codes(M=M, seed=11)

def block_sum_64(D):
    reg = D[..., a0:a0 + ACT, a0:a0 + ACT].reshape(*D.shape[:-2], NM, wt.MACRO_PX, NM, wt.MACRO_PX)
    return reg.sum(dim=(-3, -1))

def aperture_pool(frac):
    W = disk(frac)
    Ep = torch.stack([wt.propagate(wt.dmd_field((1.0 + S[i]) * 0.5).to(torch.complex128) * W, Z1) for i in range(M)])
    Em = torch.stack([wt.propagate(wt.dmd_field((1.0 - S[i]) * 0.5).to(torch.complex128) * W, Z1) for i in range(M)])
    scr = wt.make_screen(50.0, 2 * np.pi, seed=300); big = scr.shape[0]
    offs = wt.bank_offsets(N_BANK, big, 300, seed=301)
    Gp = torch.empty(N_BANK, M, NM, NM, device=DEV, dtype=wt.RDT)
    Gm = torch.empty(N_BANK, M, NM, NM, device=DEV, dtype=wt.RDT)
    for tk in range(N_BANK):
        ph = wt.screen_crop(scr, int(offs[tk, 0]), int(offs[tk, 1]))
        ep = wt.propagate(Ep * ph[None], Z2); em = wt.propagate(Em * ph[None], Z2)
        Gp[tk] = block_sum_64(ep.real ** 2 + ep.imag ** 2)
        Gm[tk] = block_sum_64(em.real ** 2 + em.imag ** 2)
    return Gp, Gm

def lam_at_n(Gp, Gm, n):
    x0 = torch.as_tensor(x64, device=DEV, dtype=wt.RDT); x1 = torch.as_tensor(x64 + d5, device=DEV, dtype=wt.RDT)
    fp0 = torch.einsum('bmij,ij->bm', Gp, x0); fm0 = torch.einsum('bmij,ij->bm', Gm, x0)
    fp1 = torch.einsum('bmij,ij->bm', Gp, x1); fm1 = torch.einsum('bmij,ij->bm', Gm, x1)
    f0 = fp0 - fm0; f1 = fp1 - fm1
    scp = n / (0.5 * (fp0 + fm0)).mean().clamp(min=1e-30)
    C0 = torch.cov(f0.T); dC = torch.cov(f1.T) - C0
    shot0 = ((fp0 + fm0).mean(0)) / scp
    V0 = C0 + torch.diag(shot0) + RIDGE * float(torch.diag(C0).mean()) * torch.eye(M, device=DEV, dtype=wt.RDT)
    A = torch.linalg.inv(V0) @ dC
    lam = 0.5 * float(torch.trace(A @ A))
    # N_eff from + exposure flux contrast
    neff = float((fp0.mean() ** 2) / (fp0.var() + 1e-30))
    return lam, neff

results = {"test": "KT5_etendue_neff_collapse",
           "R44_ref": "sec 5.4 / innovation #4 (Etendue-Information Conductance, 100)",
           "params": dict(M=M, n_bank=N_BANK, z1_mm=10, z2_mm=5, aperture_frac=APER, n_photons=NPHOT),
           "prediction": "N_eff grows with aperture area; T_det(n) collapses vs n/N_eff",
           "apertures": []}

collapse_points = []
for frac in APER:
    Gp, Gm = aperture_pool(frac)
    _, neff = lam_at_n(Gp, Gm, 1e4)
    rows = []
    for n in NPHOT:
        lam, _ = lam_at_n(Gp, Gm, n)
        tdet = 25.0 / lam if lam > 0 else float('inf')
        rows.append(dict(n=n, lam=lam, T_det=tdet, n_over_Neff=n / neff))
        collapse_points.append((n / neff, tdet, frac))
    results["apertures"].append(dict(aperture_frac=frac, N_eff=neff, rows=rows))
    print("[aper=%.2f] N_eff=%.1f | T_det(n=1e4)=%.1f  r=n/Neff at 1e4 = %.2f [%.0fs]"
          % (frac, neff, 25.0 / rows[2]["lam"], 1e4 / neff, time.time() - t0))
    del Gp, Gm
    if DEV == "cuda":
        torch.cuda.empty_cache()

# ------------------------------------------------------------------ collapse test
# group collapse points by nearest n/Neff decade and measure spread of T_det across apertures
neffs = [a["N_eff"] for a in results["apertures"]]
neff_grows = all(neffs[i] <= neffs[i + 1] * 1.05 for i in range(len(neffs) - 1))
# collapse metric: at matched r=n/Neff, T_det should agree across apertures.
# bin log10(r), compute coefficient of variation of T_det within bins that have >=2 apertures.
import collections
bins = collections.defaultdict(list)
for r, td, frac in collapse_points:
    bins[round(np.log10(r), 0)].append(td)
cvs = []
for k, v in bins.items():
    if len(v) >= 2:
        cvs.append(float(np.std(v) / (np.mean(v) + 1e-30)))
mean_cv = float(np.mean(cvs)) if cvs else None
# reference: WITHOUT the n/Neff rescaling, T_det at fixed n varies across apertures -- compute that spread
fixed_n_td = collections.defaultdict(list)
for a in results["apertures"]:
    for row in a["rows"]:
        fixed_n_td[row["n"]].append(row["T_det"])
raw_cvs = [float(np.std(v) / (np.mean(v) + 1e-30)) for v in fixed_n_td.values() if len(v) >= 2]
mean_raw_cv = float(np.mean(raw_cvs)) if raw_cvs else None

# NORMALIZED collapse: the etendue-conductance law says T_det/N_eff is a single function
# of r=n/N_eff (photon-saturated regime -> T_det ∝ N_eff, i.e. T_det/N_eff ~ const).
norm_bins = collections.defaultdict(list)
for r, td, frac in collapse_points:
    neff_f = next(a["N_eff"] for a in results["apertures"] if a["aperture_frac"] == frac)
    norm_bins[round(np.log10(r), 0)].append(td / neff_f)
norm_cvs = [float(np.std(v) / (np.mean(v) + 1e-30)) for v in norm_bins.values() if len(v) >= 2]
mean_norm_cv = float(np.mean(norm_cvs)) if norm_cvs else None
# also T_det/N_eff at each fixed n across apertures (saturated-regime constancy)
sat_cv = []
for a_n in fixed_n_td:
    vals = [row["T_det"] / a["N_eff"] for a in results["apertures"] for row in a["rows"] if row["n"] == a_n]
    if len(vals) >= 2:
        sat_cv.append(float(np.std(vals) / (np.mean(vals) + 1e-30)))
mean_sat_cv = float(np.mean(sat_cv)) if sat_cv else None

collapse_ok = (mean_norm_cv is not None) and (mean_norm_cv < 0.35)
results["verdict"] = dict(
    N_eff_vs_aperture=neffs, N_eff_grows_with_aperture=bool(neff_grows),
    raw_CV_Tdet_vs_n_over_Neff=mean_cv, raw_CV_at_fixed_n=mean_raw_cv,
    normalized_CV_Tdet_over_Neff_vs_r=mean_norm_cv, saturated_CV_Tdet_over_Neff=mean_sat_cv,
    P_Neff_grows="PASS" if neff_grows else "PARTIAL",
    P_Tdet_over_Neff_collapses="PASS" if collapse_ok else "PARTIAL",
    overall=("N_eff grows with illumination aperture (etendue: %s); the etendue-conductance "
             "law holds -- T_det/N_eff collapses to a single function of r=n/N_eff (normalized "
             "CV=%.2f; saturated-regime T_det/N_eff CV=%.2f), so the resource is (etendue, "
             "independent banks) and photons stop helping for n>>N_eff." %
             ([round(x, 1) for x in neffs], mean_norm_cv if mean_norm_cv else float('nan'),
              mean_sat_cv if mean_sat_cv else float('nan'))
             if collapse_ok else
             "N_eff grows with aperture (PASS); T_det/N_eff collapse is partial (normCV=%s, satCV=%s) "
             "-- report honestly" % (mean_norm_cv, mean_sat_cv)))
results["elapsed_s"] = round(time.time() - t0, 1)
with open(os.environ.get("OUT", "KT5_etendue_neff_collapse.json"), "w") as f:
    json.dump(results, f, indent=2)
print("\nVERDICT:", json.dumps(results["verdict"]))
print("saved  elapsed %.1fs" % (time.time() - t0))
