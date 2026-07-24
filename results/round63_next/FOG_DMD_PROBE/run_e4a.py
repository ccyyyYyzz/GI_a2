# E4a -- CLASSICAL rotation fix via multi-timescale medium + SOBI joint diagonalization.
# Best cell: M=96, d=48, sigma_f=0.15, T=128, photons {clean,1e6,1e5}.
# Compare single-tau vs multi-tau media, each WITH/WITHOUT the SOBI rotation-fix init.
# Attribution: SOBI can only fix the rotation when lagged covariances are distinctly
# diagonalizable -> requires multi-timescale (distinct tau_k). Single-tau is the control.
import time, json, numpy as np, torch
from fog_common import *
import fog_tracker as ft
import fog_sobi as sobi
def log(*a): print(*a, flush=True)
dev = 'cuda' if torch.cuda.is_available() else 'cpu'
n_side = 32; N = n_side * n_side; x_true = make_scene(n_side)
M = 96; d = 48; T = 128; sigma_f = 0.15
tau_arr = tau_schedule(d, 64.0, 4.0)          # multi-timescale schedule (coarse->fine)
log(f"device {dev}; cell M={M} d={d} sf={sigma_f} T={T}")
log(f"multi-tau schedule: tau_0={tau_arr[0]:.1f} .. tau_{d-1}={tau_arr[-1]:.2f} "
    f"(rho {np.exp(-1/tau_arr[0]):.3f}..{np.exp(-1/tau_arr[-1]):.3f})")

def build(medium, photons, seed):
    rng = np.random.default_rng(4040 + seed + (0 if photons is None else int(np.log10(photons))))
    P = make_patterns(M, N, "binary", rng); Pi, rangeP = projectors(P)
    null_true = x_true - rangeP(x_true); nrm0 = np.linalg.norm(null_true)
    U = dct_basis(d, n_side)
    if medium == 'multi':
        W, Z, rho, csd = lognormal_medium_mt(U, T, sigma_f, tau_arr, rng)
    else:  # single tau=16 (baseline temporal model)
        W, Z, rho, csd = lognormal_medium(U, T, sigma_f, 16.0, rng)
    mu = np.einsum('mn,tn->tm', P, W * x_true)
    if photons is None:
        B = mu.copy(); wts = np.ones_like(B); Rd = np.ones_like(B) * 1e-8; lam = 1e-3
    else:
        s = photons / mu.mean(); B = rng.poisson(np.clip(mu, 0, None) * s) / s
        wts = (s / np.clip(mu, 1e-9, None)); Rd = np.clip(mu, 1e-9, None) / s; lam = 3e-3
    return dict(P=P, U=U, W=W, Z=Z, B=B, wts=wts, Rd=Rd, rho=rho, csd=csd, lam=lam,
                rangeP=rangeP, null_true=null_true, nrm0=nrm0)

def nmse(x, c): return null_metric(x, x_true, c['rangeP'], c['null_true'], c['nrm0'])[0] ** 2

def sobi_quality(zc, Ztrue):
    # per-source max |corr| after SOBI (1.0 = perfectly separated true coordinate)
    zc = (zc - zc.mean(0)) / (zc.std(0) + 1e-12)
    Zt = (Ztrue - Ztrue.mean(0)) / (Ztrue.std(0) + 1e-12)
    C = np.abs(zc.T @ Zt) / zc.shape[0]         # (d,d)
    return float(C.max(axis=1).mean())

results = []
for medium in ['single', 'multi']:
    for photons in [None, 1e6, 1e5]:
        base, fix, sq, froz, truez = [], [], [], [], []
        oclist = []
        for sd in range(5):
            c = build(medium, photons, sd)
            oc = nmse(ft.solve_oracle(c['P'], c['W'], c['B'], c['wts'],
                                      1e-6 if photons is None else 1e-3, dev), c)
            oclist.append(oc)
            lamx = 1e-6 if photons is None else 1e-3
            # ATTRIBUTION: frozen KNOWN medium (linear-approx of true W) -> the ceiling if the
            # medium were not re-estimated; and true-Z-seeded refinement -> exposes the drift.
            Wl = np.clip(1.0 + c['Z'] @ c['U'].T, 0.05, None)
            froz.append(nmse(ft.solve_oracle(c['P'], Wl, c['B'], c['wts'], lamx, dev), c))
            rt = ft.staged_blind(c['P'], c['U'], c['B'], c['wts'], c['rho'], c['csd'], c['Rd'], c['lam'],
                                 dev=dev, seeds=(sd,), n_als=90, refine_em=15, sobi_z=c['Z'])
            truez.append(nmse(rt['x'], c))
            # baseline (random-seed staged blind)
            rb = ft.staged_blind(c['P'], c['U'], c['B'], c['wts'], c['rho'], c['csd'], c['Rd'], c['lam'],
                                 dev=dev, seeds=(sd,), n_als=90, refine_em=15)
            base.append(nmse(rb['x'], c))
            # SOBI rotation-fix init
            zc = sobi.sobi_sources(c['B'] - c['B'].mean(0), d)
            sq.append(sobi_quality(zc, c['Z']))
            rf = ft.staged_blind(c['P'], c['U'], c['B'], c['wts'], c['rho'], c['csd'], c['Rd'], c['lam'],
                                 dev=dev, seeds=(sd,), n_als=90, refine_em=15, sobi_z=zc)
            fix.append(nmse(rf['x'], c))
        oref = float(np.median(oclist))
        bm, fm = float(np.median(base)), float(np.median(fix))
        cap_b = (1 - bm) / (1 - oref) if (1 - oref) > 1e-6 else float('nan')
        cap_f = (1 - fm) / (1 - oref) if (1 - oref) > 1e-6 else float('nan')
        row = dict(medium=medium, photons=None if photons is None else float(photons),
                   oracle_nmse=oref, frozen_knownmedium_nmse=float(np.median(froz)),
                   trueZ_seed_nmse=float(np.median(truez)),
                   baseline_nmse_median=bm, sobifix_nmse_median=fm,
                   baseline_all=base, sobifix_all=fix, sobi_quality=float(np.median(sq)),
                   cap_baseline=cap_b, cap_sobifix=cap_f)
        results.append(row)
        log(f"[{medium:6s} ph={('clean' if photons is None else f'{photons:.0e}')}] oracle {oref:.3f} | "
            f"frozen-known {row['frozen_knownmedium_nmse']:.3f} | true-Z-seed {row['trueZ_seed_nmse']:.3f} | "
            f"baseline {bm:.3f} | SOBI {fm:.3f} | sobi-qual {row['sobi_quality']:.2f}")

with open('E4a_results.json', 'w') as f:
    json.dump(dict(cell=dict(M=M, d=d, T=T, sigma_f=sigma_f, tau_multi=list(map(float, tau_arr))),
                   results=results), f, indent=2)
best_multi = min((r for r in results if r['medium'] == 'multi'), key=lambda r: r['sobifix_nmse_median'])
passF3 = any((r['medium'] == 'multi' and r['sobifix_nmse_median'] <= 0.25 and r['cap_sobifix'] >= 0.5)
             for r in results)
log(f"\n=== E4a: best multi-tau SOBI-fix NMSE = {best_multi['sobifix_nmse_median']:.3f} "
    f"(cap {best_multi['cap_sobifix']*100:.0f}%) at ph={best_multi['photons']} ===")
log(f"ROTATION_FIX verdict: {'WORKS' if passF3 else 'FAILS'} (F3-shadow: <=0.25 & >=50%)")
