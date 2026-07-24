# E5 -- covariance-domain resurrection. E5a is decisive: estimate x WITHOUT ever estimating the
# medium realizations z_t, by matching the (shot-free, lag>=1) bucket covariances. No z -> nothing
# to collude with. Declared law: single-tau OU (tau=16, r_l=rho^l known), K_w=coeff_sd^2 U U^T.
import time, json, numpy as np, torch
from fog_common import *
import fog_tracker as ft
def log(*a): print(*a, flush=True)
dev = 'cuda' if torch.cuda.is_available() else 'cpu'
n = 32; N = n * n; M = 96; d = 48; sf = 0.15; TAU = 16.0
x_true = make_scene(n)

def build(T, photons, seed):
    rng = np.random.default_rng(50500 + seed + (0 if photons is None else int(np.log10(photons))))
    P = make_patterns(M, N, "binary", rng); Pi, rangeP = projectors(P)
    nt = x_true - rangeP(x_true); nrm0 = np.linalg.norm(nt)
    U = dct_basis(d, n)
    W, Z, rho, csd = lognormal_medium(U, T, sf, TAU, rng)
    mu = np.einsum('mn,tn->tm', P, W * x_true)
    if photons is None:
        B = mu.copy()
    else:
        s = photons / mu.mean(); B = rng.poisson(np.clip(mu, 0, None) * s) / s
    return dict(P=P, U=U, W=W, B=B, rho=rho, csd=csd, rangeP=rangeP, nt=nt, nrm0=nrm0)

def nm(x, c): return null_metric(x, x_true, c['rangeP'], c['nt'], c['nrm0'])[0] ** 2

def cov_run(c, n_starts=6):
    r = ft.cov_estimate(c['P'], c['U'], c['B'], c['csd'], c['rho'], n_starts=n_starts,
                        steps=5000, dev=dev)
    vals = [nm(xx, c) for xx in r['per_start_x']]
    best = r['per_start_x'][int(np.argmin(r['objs']))]      # pick by DATA objective (blind)
    return float(np.median(vals)), float(nm(best, c)), float(np.std(vals))

R = {}
t0 = time.time()

# --- exact-target injectivity ceiling (diagnostic) ---
c = build(128, None, 0)
Pt = torch.tensor(c['P'], device=dev, dtype=torch.float64); Ut = torch.tensor(c['U'], device=dev, dtype=torch.float64)
xt = torch.tensor(x_true, device=dev, dtype=torch.float64); H = Pt @ (Ut * xt[:, None])
r_exact = ft.cov_estimate(c['P'], c['U'], c['B'], c['csd'], c['rho'], n_starts=1, steps=8000, dev=dev,
                          That_override=(H @ H.t()).cpu().numpy())
R['exact_target_ceiling'] = float(nm(r_exact['per_start_x'][0], c))
log(f"exact-target ceiling (injectivity): null-NMSE {R['exact_target_ceiling']:.3f}")

# --- E5a headline cell: T=128 (spec), photons, multi-start agreement (drift test) ---
log("\n=== E5a cell T=128 (multi-start agreement = drift/collusion test) ===")
R['cell_T128'] = []
for photons in [None, 1e6, 1e5]:
    med, best, agree = [], [], []
    for sd in range(5):
        c = build(128, photons, sd); m, bnm, ag = cov_run(c)
        med.append(m); best.append(bnm); agree.append(ag)
    row = dict(T=128, photons=None if photons is None else float(photons),
               nmse_median=float(np.median(med)), nmse_best=float(np.median(best)),
               agreement_std=float(np.median(agree)))
    R['cell_T128'].append(row)
    log(f"  ph={'clean' if photons is None else f'{photons:.0e}'}: null-NMSE median {row['nmse_median']:.3f} "
        f"| multi-start agreement std {row['agreement_std']:.3f} (drift-free if small)")

# --- E5a T-scaling (clean): the statistical sample-complexity curve ---
log("\n=== E5a T-scaling (clean): sample-complexity resurrection curve ===")
R['T_scaling'] = []
for T in [128, 256, 512, 1024, 2048, 4096]:
    vals = []
    for sd in range(3):
        c = build(T, None, sd); m, bnm, ag = cov_run(c, n_starts=3); vals.append(bnm)
    row = dict(T=T, nmse_median=float(np.median(vals)))
    R['T_scaling'].append(row)
    log(f"  T={T:5d}: null-NMSE {row['nmse_median']:.3f}")

# --- E5a photon-robustness at a working T (lags>=1 are shot-free) ---
log("\n=== E5a photon-robustness at T=1024 (lags>=1 shot-free) ===")
R['photon_at_T1024'] = []
for photons in [None, 1e6, 1e5, 1e4]:
    vals = []
    for sd in range(3):
        c = build(1024, photons, sd); m, bnm, ag = cov_run(c, n_starts=3); vals.append(bnm)
    row = dict(T=1024, photons=None if photons is None else float(photons), nmse_median=float(np.median(vals)))
    R['photon_at_T1024'].append(row)
    log(f"  ph={'clean' if photons is None else f'{photons:.0e}'}: null-NMSE {row['nmse_median']:.3f}")

R['meta'] = dict(M=M, d=d, sigma_f=sf, tau=TAU, medium='lognormal single-tau OU',
                 metric='null-NMSE', wall_s=time.time() - t0)
json.dump(R, open('E5_results.json', 'w'), indent=2)
Tcurve = {r['T']: r['nmse_median'] for r in R['T_scaling']}
Tstar = next((T for T in sorted(Tcurve) if Tcurve[T] <= 0.25), None)
log(f"\n=== E5a VERDICT ===")
log(f"multi-start agreement std @T128 clean = {R['cell_T128'][0]['agreement_std']:.3f} (drift-free)")
log(f"injectivity ceiling (exact target) = {R['exact_target_ceiling']:.3f}")
log(f"crosses 0.25 bar at T ~ {Tstar} (vs pathwise T=128 which fails)")
log(f"COVARIANCE_STICKS: {'YES' if R['cell_T128'][0]['agreement_std'] < 0.05 else 'NO'}; "
    f"recovers null given T: {'YES' if Tstar else 'NO'}")
log(f"total wall {time.time()-t0:.0f}s")
