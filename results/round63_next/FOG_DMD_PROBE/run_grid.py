# E1 grid (R38-aligned): negative controls (M=64) + primary blind-capacity cells (M=96),
# fixed-sigma_f lognormal medium, physical 0/1 patterns, staged data-only cold start,
# >=5-seed medians, oracle + linear-approx-oracle + blind-cold + warm-Q diagnostic ceiling,
# + OU-prior-weakened sensitivity. Metric: null-NMSE (= amplitude null-err^2). Baseline=1.
import time, json, argparse, numpy as np, torch
from fog_common import *
import fog_tracker as ft
def log(*a): print(*a, flush=True)
dev = 'cuda' if torch.cuda.is_available() else 'cpu'
n_side = 32; N = n_side * n_side; x_true = make_scene(n_side)

def build(M, d, T, sigma_f, tau, photons, seed, medium='lognormal'):
    rng = np.random.default_rng(seed)
    P = make_patterns(M, N, "binary", rng)               # physical 0/1
    Pi, rangeP = projectors(P)
    null_true = x_true - rangeP(x_true); nrm0 = np.linalg.norm(null_true)
    if d == 0:                                            # scalar-gain control (d=1 constant)
        U = np.ones((N, 1)) / np.sqrt(N); d_eff = 1
    else:
        U = dct_basis(d, n_side); d_eff = d
    if medium == 'linear':                               # model-matched affine medium (reference)
        coeff_sd = sigma_f * np.sqrt(N / U.shape[1])
        if tau is None:
            Z = coeff_sd * rng.standard_normal((T, U.shape[1])); rho = 0.0
        else:
            rho = np.exp(-1.0 / tau); Z = np.zeros((T, U.shape[1]))
            Z[0] = coeff_sd * rng.standard_normal(U.shape[1]); a = np.sqrt(1 - rho ** 2) * coeff_sd
            for t in range(1, T): Z[t] = rho * Z[t - 1] + a * rng.standard_normal(U.shape[1])
        W = np.clip(1.0 + Z @ U.T, 0.05, None)
    else:
        W, Z, rho, coeff_sd = lognormal_medium(U, T, sigma_f, tau, rng)
    mu = np.einsum('mn,tn->tm', P, W * x_true)
    if photons is None:
        B = mu.copy(); wts = np.ones_like(B); Rd = np.ones_like(B) * 1e-8; lam = 1e-3
    else:
        s = photons / mu.mean(); B = rng.poisson(np.clip(mu, 0, None) * s) / s
        wts = (s / np.clip(mu, 1e-9, None)); Rd = np.clip(mu, 1e-9, None) / s; lam = 3e-3
    return dict(P=P, U=U, W=W, Z=Z, B=B, wts=wts, Rd=Rd, rho=rho, coeff_sd=coeff_sd, lam=lam,
                rangeP=rangeP, null_true=null_true, nrm0=nrm0, d_eff=d_eff)

def nmse(x, c): return null_metric(x, x_true, c['rangeP'], c['null_true'], c['nrm0'])[0] ** 2

def warm_Q(c, lam=1e-4, iters=120):
    P = torch.tensor(c['P'], device=dev, dtype=ft.DTYPE); U = torch.tensor(c['U'], device=dev, dtype=ft.DTYPE)
    B = torch.tensor(c['B'], device=dev, dtype=ft.DTYPE); d = c['U'].shape[1]
    Bc = B - B.mean(0, keepdim=True); _, _, Vh = torch.linalg.svd(Bc.t(), full_matrices=False); Zh = Vh[:d].t()
    Q = torch.linalg.lstsq(Zh, torch.tensor(c['Z'], device=dev, dtype=ft.DTYPE)).solution
    PtP = P.t() @ P; PtBt = P.t() @ B.t()
    for it in range(iters):
        W = torch.clamp(1.0 + (Zh @ Q) @ U.t(), min=0.05); x = torch.clamp(ft.solve_x_clean(PtP, PtBt, W, lam), min=0.0)
        H = P @ (U * x[:, None]); HtH = H.t() @ H + 1e-8 * torch.eye(d, device=dev, dtype=ft.DTYPE)
        y = B - (P @ x)[None, :]; Q = torch.linalg.solve(Zh.t() @ Zh, Zh.t() @ y @ H) @ torch.linalg.inv(HtH)
    return nmse(x.detach().cpu().numpy(), c)

def blind_cold(c, seeds, ou_scale=1.0):
    vals = []
    for sd in seeds:
        r = ft.staged_blind(c['P'], c['U'], c['B'], c['wts'], c['rho'], c['coeff_sd'], c['Rd'], c['lam'],
                            dev=dev, seeds=(sd,), n_als=90, refine_em=15, ou_prior_scale=ou_scale)
        vals.append(nmse(r['x'], c))
    return vals

def cell(tag, M, d, T, sigma_f, tau, photons, seeds, do_warm=False, ou_scale=1.0, medium='lognormal'):
    # oracle & linear-approx-oracle use median over the same seeds
    ne_or, ne_lin, bcs = [], [], []
    for sd in seeds:
        c = build(M, d, T, sigma_f, tau, photons, 20260723 + sd + M * 17 + d * 3 + int(sigma_f * 100), medium)
        ne_or.append(nmse(ft.solve_oracle(c['P'], c['W'], c['B'], c['wts'], 1e-6 if photons is None else 1e-3, dev), c))
        Wl = np.clip(1.0 + c['Z'] @ c['U'].T, 0.05, None)
        ne_lin.append(nmse(ft.solve_oracle(c['P'], Wl, c['B'], c['wts'], 1e-6 if photons is None else 1e-3, dev), c))
        r = ft.staged_blind(c['P'], c['U'], c['B'], c['wts'], c['rho'], c['coeff_sd'], c['Rd'], c['lam'],
                            dev=dev, seeds=(sd,), n_als=90, refine_em=15, ou_prior_scale=ou_scale)
        bcs.append(nmse(r['x'], c))
    wc = warm_Q(build(M, d, T, sigma_f, tau, photons, 20260723 + seeds[0] + M * 17 + d * 3 + int(sigma_f * 100), medium)) if do_warm else None
    orm = float(np.median(ne_or)); bcm = float(np.median(bcs))
    cap = (min(T, c['d_eff'] + 1)) * (M - c['d_eff']) / N
    frac = (1 - bcm) / (1 - orm) if (1 - orm) > 1e-6 else float('nan')   # oracle improvement captured
    row = dict(tag=tag, M=M, d=c['d_eff'], T=T, sigma_f=sigma_f, tau=tau, medium=medium,
               photons=None if photons is None else float(photons), chi=cap,
               oracle_nmse=orm, linapprox_oracle_nmse=float(np.median(ne_lin)),
               blind_cold_nmse_median=bcm, blind_cold_nmse_all=[float(v) for v in bcs],
               warm_Q_ceiling_nmse=wc, frac_oracle_improvement=frac)
    log(f"{tag:26s} M={M} d={c['d_eff']:3d} sf={sigma_f} tau={str(tau):>4} "
        f"ph={'clean' if photons is None else f'{photons:.0e}'} chi={cap:.2f} | "
        f"oracle {orm:.3f} linapx {row['linapprox_oracle_nmse']:.3f} | blind-cold {bcm:.3f} "
        f"(cap {frac*100:.0f}%)" + (f" | warmQ {wc:.3f}" if wc is not None else ""))
    return row

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--seeds', type=int, default=5)
    ap.add_argument('--out', default='E1_results.json'); ap.add_argument('--quick', type=int, default=0)
    a = ap.parse_args(); seeds = list(range(a.seeds)); t0 = time.time(); R = {}
    T = 128
    log("=== NEGATIVE CONTROLS (M=64: theorem predicts blind failure) ===")
    R['controls'] = []
    ctrl = [('ctrl scalar d=1', 64, 0), ('ctrl d=16', 64, 16), ('ctrl d=31', 64, 31), ('ctrl d=64', 64, 64)]
    if a.quick: ctrl = ctrl[:2]
    for tag, M, d in ctrl:
        R['controls'].append(cell(tag, M, d, T, 0.30, None, None, seeds, do_warm=(d in (16, 31))))
    log("\n=== PRIMARY BLIND-CAPACITY CELLS (M=96) ===")
    R['primary'] = []
    dlist = [48] if a.quick else [32, 48, 64]
    sflist = [0.30] if a.quick else [0.15, 0.30]
    phlist = [None, 1e5] if a.quick else [None, 1e6, 1e5, 1e4]
    for d in dlist:
        for sf in sflist:
            for ph in phlist:
                R['primary'].append(cell(f'M96 d{d} sf{sf}', 96, d, T, sf, None, ph, seeds,
                                         do_warm=(ph is None)))
    log("\n=== MODEL-MATCHED LINEAR-MEDIUM REFERENCE (isolate cold-start vs lognormal mismatch) ===")
    R['linear_reference'] = []
    for d in ([48] if a.quick else [32, 48, 64]):
        R['linear_reference'].append(cell(f'LINREF M96 d{d} sf0.30', 96, d, T, 0.30, None, None, seeds,
                                          do_warm=True, medium='linear'))
    log("\n=== OU-PRIOR SENSITIVITY (medium OU tau=16; prior weakened 4x) ===")
    R['ou_sensitivity'] = []
    for scale, lbl in [(1.0, 'ou_full'), (0.25, 'ou_weak4x')]:
        R['ou_sensitivity'].append(cell(f'M96 d48 sf0.30 OU {lbl}', 96, 48, T, 0.30, 16, 1e5, seeds, ou_scale=scale))
    R['meta'] = dict(n_side=n_side, N=N, T=T, seeds=a.seeds, wall_s=time.time() - t0, dev=str(dev),
                     metric='null-NMSE = amplitude null-err^2; baseline(fixed-avg)=1.0',
                     success_F3='blind_cold_nmse_median<=0.25 AND frac_oracle_improvement>=0.5')
    with open(a.out, 'w') as f: json.dump(R, f, indent=2)
    # verdict
    prim = R['primary']
    best = min(prim, key=lambda r: r['blind_cold_nmse_median'])
    passF3 = any((r['blind_cold_nmse_median'] <= 0.25 and r['frac_oracle_improvement'] >= 0.5) for r in prim)
    log(f"\n=== E1 VERDICT: best primary blind-cold NMSE = {best['blind_cold_nmse_median']:.3f} "
        f"at {best['tag']} ph={best['photons']} (captures {best['frac_oracle_improvement']*100:.0f}%) ===")
    log(f"F3 SUCCESS SIGNAL (<=0.25 & >=50% captured): {'MET' if passF3 else 'NOT MET'}")
    log(f"total wall {time.time()-t0:.1f}s -> {a.out}")

if __name__ == '__main__': main()
