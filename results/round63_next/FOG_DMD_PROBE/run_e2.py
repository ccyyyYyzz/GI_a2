# E2 -- NOISE SCALING CURVE (regularized). null-NMSE vs mean detected photons/bucket for
# ORACLE and BLIND cold-start, at the primary M=96,d=48,sf=0.30 lognormal cell. Tikhonov on x
# + OU prior on z. Fixed-averaging baseline = 1.0 (cannot see ker P). Saves JSON + PNG.
import time, json, numpy as np, torch
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from fog_common import *
import fog_tracker as ft
def log(*a): print(*a, flush=True)
dev = 'cuda' if torch.cuda.is_available() else 'cpu'
n_side = 32; N = n_side * n_side; x_true = make_scene(n_side)
M = 96; d = 48; T = 128; sigma_f = 0.30

def nmse(x, rangeP, null_true, nrm0): return null_metric(x, x_true, rangeP, null_true, nrm0)[0] ** 2

photons_list = [1e3, 1e4, 1e5, 1e6, 1e7]
seeds = list(range(5))
curve = []
t0 = time.time()
for photons in photons_list:
    oc, bc = [], []
    for sd in seeds:
        rng = np.random.default_rng(555 + sd + int(np.log10(photons)))
        P = make_patterns(M, N, "binary", rng); Pi, rangeP = projectors(P)
        null_true = x_true - rangeP(x_true); nrm0 = np.linalg.norm(null_true)
        U = dct_basis(d, n_side)
        W, Z, rho, csd = lognormal_medium(U, T, sigma_f, None, rng)
        mu = np.einsum('mn,tn->tm', P, W * x_true)
        s = photons / mu.mean(); B = rng.poisson(np.clip(mu, 0, None) * s) / s
        wts = (s / np.clip(mu, 1e-9, None)); Rd = np.clip(mu, 1e-9, None) / s
        lam = 3e-3   # Tikhonov (regularized so low-photon degrades toward 1.0, not explode)
        oc.append(nmse(ft.solve_oracle(P, W, B, wts, lam, dev), rangeP, null_true, nrm0))
        r = ft.staged_blind(P, U, B, wts, rho, csd, Rd, lam, dev=dev, seeds=(sd,), n_als=90, refine_em=15)
        bc.append(nmse(r['x'], rangeP, null_true, nrm0))
    row = dict(photons=photons, oracle_nmse_median=float(np.median(oc)),
               oracle_nmse_all=[float(v) for v in oc],
               blind_nmse_median=float(np.median(bc)), blind_nmse_all=[float(v) for v in bc],
               baseline_nmse=1.0)
    curve.append(row)
    log(f"photons {photons:.0e} | oracle NMSE {row['oracle_nmse_median']:.3f} | "
        f"blind NMSE {row['blind_nmse_median']:.3f} | baseline 1.000")

with open('E2_noise_curve.json', 'w') as f:
    json.dump(dict(cell=dict(M=M, d=d, T=T, sigma_f=sigma_f, medium='lognormal', seeds=len(seeds)),
                   curve=curve, wall_s=time.time() - t0), f, indent=2)

# PNG
ph = [r['photons'] for r in curve]
orc = [r['oracle_nmse_median'] for r in curve]
bl = [r['blind_nmse_median'] for r in curve]
plt.figure(figsize=(6, 4.2))
plt.axhline(1.0, color='0.6', ls='--', lw=1.2, label='fixed-averaging baseline (=1.0)')
plt.plot(ph, orc, 'o-', color='#1f77b4', label='oracle (medium known)')
plt.plot(ph, bl, 's-', color='#d62728', label='blind cold-start (staged)')
plt.axhline(0.25, color='green', ls=':', lw=1.0, label='F3 bar (NMSE=0.25)')
plt.xscale('log'); plt.xlabel('mean detected photons / bucket'); plt.ylabel('null-space NMSE')
plt.title(f'E2 noise scaling  (M={M}, d={d}, T={T}, $\\sigma_f$={sigma_f}, lognormal)')
plt.ylim(-0.05, 1.15); plt.legend(fontsize=8, loc='center left'); plt.grid(alpha=0.3)
plt.tight_layout(); plt.savefig('E2_noise_curve.png', dpi=140)
log(f"saved E2_noise_curve.json + .png  [{time.time()-t0:.1f}s]")
