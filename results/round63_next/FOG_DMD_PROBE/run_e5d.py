# E5d -- DECLARED-LAW MISMATCH robustness (R38 F5 kill test) for the E5a covariance estimator.
# Data generated with TRUE law; estimator uses a DECLARED (possibly wrong) law. Same cell
# M=96, d=48, sigma_f=0.15; T in {2048,1024}; 1e5 photons; 5 seeds. Report null-NMSE degradation
# vs the matched baseline. Verdict: MISMATCH_ROBUST vs EXACT_LAW_LOAD_BEARING.
import time, json, numpy as np, torch
from fog_common import *
import fog_tracker as ft
def log(*a): print(*a, flush=True)
dev = 'cuda' if torch.cuda.is_available() else 'cpu'
n = 32; N = n * n; M = 96; D_EST = 48; SF_EST = 0.15; TAU_EST = 16.0
x_true = make_scene(n)
U48 = dct_basis(D_EST, n)                              # estimator basis (also true, when matched)
RHO_EST = np.exp(-1.0 / TAU_EST)
CSD_EST = SF_EST * np.sqrt(N / D_EST)

def rotate_basis(U, eps, rng):
    Nn, dd = U.shape
    R = rng.standard_normal((Nn, dd)); R = R - U @ (U.T @ R)   # perpendicular component
    R = R / (np.linalg.norm(R, axis=0, keepdims=True) + 1e-12)
    Um = U * np.sqrt(1 - eps ** 2) + R * eps
    Ue, _ = np.linalg.qr(Um)
    ang = np.arccos(np.clip(np.linalg.svd(U.T @ Ue, compute_uv=False), -1, 1))
    return Ue, float(np.sin(ang).mean())              # mean sin(principal angle) = actual mismatch

def gen(T, photons, seed, U_true, sf_true, tau_true):
    rng = np.random.default_rng(90900 + seed + (0 if photons is None else int(np.log10(photons))))
    P = make_patterns(M, N, "binary", rng); Pi, rangeP = projectors(P)
    nt = x_true - rangeP(x_true); nrm0 = np.linalg.norm(nt)
    W, Z, rho_t, csd_t = lognormal_medium(U_true, T, sf_true, tau_true, rng)
    mu = np.einsum('mn,tn->tm', P, W * x_true)
    s = photons / mu.mean(); B = rng.poisson(np.clip(mu, 0, None) * s) / s
    return P, B, rangeP, nt, nrm0

def estimate(P, B, U_est, csd_est, rho_est, rangeP, nt, nrm0, n_starts=3):
    r = ft.cov_estimate(P, U_est, B, csd_est, rho_est, n_starts=n_starts, steps=4500, dev=dev)
    best = r['per_start_x'][int(np.argmin(r['objs']))]
    return null_metric(best, x_true, rangeP, nt, nrm0)[0] ** 2

def main():
  R = {}; t0 = time.time()
  for T in [2048, 1024]:
    log(f"\n########## T={T}, 1e5 photons ##########")
    key = f"T{T}"; R[key] = {}
    # matched baseline
    base = []
    for sd in range(5):
        P, B, rp, nt, nrm0 = gen(T, 1e5, sd, U48, SF_EST, TAU_EST)
        base.append(estimate(P, B, U48, CSD_EST, RHO_EST, rp, nt, nrm0))
    b0 = float(np.median(base)); R[key]['matched'] = b0
    log(f"MATCHED baseline null-NMSE = {b0:.3f}")

    def axis(name, cfgs):
        R[key][name] = []
        for lbl, gkw, ekw in cfgs:
            vals = []; measured = None
            for sd in range(5):
                rng = np.random.default_rng(313 + sd)
                U_true = gkw.get('U_true', U48); U_est = ekw.get('U_est', U48)
                if 'rot_eps' in ekw:
                    U_est, measured = rotate_basis(U48, ekw['rot_eps'], rng)  # estimator's rotated basis
                P, B, rp, nt, nrm0 = gen(T, 1e5, sd, gkw.get('U_true', U48),
                                         gkw.get('sf', SF_EST), gkw.get('tau', TAU_EST))
                vals.append(estimate(P, B, U_est, ekw.get('csd', CSD_EST), ekw.get('rho', RHO_EST),
                                     rp, nt, nrm0))
            m = float(np.median(vals)); deg = (m - b0) / b0 * 100
            R[key][name].append(dict(label=lbl, nmse=m, degradation_pct=deg,
                                     measured_mismatch=measured))
            log(f"  [{name}] {lbl:16s} null-NMSE {m:.3f}  (deg {deg:+.0f}%"
                + (f", angle {measured*100:.0f}%" if measured is not None else "") + ")")

    # 1. basis rotation (estimator uses rotated U; simulator keeps U48)
    axis("basis_rotation", [(f"rot~{int(e*100)}%", {}, {'rot_eps': e}) for e in (0.05, 0.10, 0.20)])
    # 2. wrong tau (simulator tau in {8,32}; estimator assumes 16)
    axis("wrong_tau", [(f"true_tau={tt}", {'tau': float(tt)}, {}) for tt in (8, 32)])
    # 3. wrong sigma_f (simulator {0.10,0.20}; estimator assumes 0.15 -> pure scale/gauge)
    axis("wrong_sigma_f", [(f"true_sf={sfx}", {'sf': sfx}, {}) for sfx in (0.10, 0.20)])
    # 4. subspace dim mismatch (simulator d in {40,56}; estimator d=48)
    dcfgs = []
    for dt in (40, 56):
        Ut = dct_basis(dt, n)
        dcfgs.append((f"true_d={dt}", {'U_true': Ut, 'sf': SF_EST}, {}))
    axis("dim_mismatch", dcfgs)

  R['meta'] = dict(cell="M=96 d=48 sf=0.15", photons=1e5, seeds=5, wall_s=time.time() - t0)
  json.dump(R, open('E5d_results.json', 'w'), indent=2)
  # verdict per F5 spirit: degradation <25% in the 10%-rotation and 2x-tau cells
  def find(k, ax, sub):
      return next((r for r in R[k][ax] if sub in r['label']), None)
  rot10 = find('T2048', 'basis_rotation', '10%')
  tau8 = find('T2048', 'wrong_tau', 'tau=8'); tau32 = find('T2048', 'wrong_tau', 'tau=32')
  worst = max(abs(rot10['degradation_pct']), abs(tau8['degradation_pct']), abs(tau32['degradation_pct']))
  robust = worst < 25
  log(f"\n=== E5d VERDICT ===")
  log(f"key cells @T=2048: rot10% deg {rot10['degradation_pct']:+.0f}% | tau8 {tau8['degradation_pct']:+.0f}% | tau32 {tau32['degradation_pct']:+.0f}%")
  log(f"{'MISMATCH_ROBUST' if robust else 'EXACT_LAW_LOAD_BEARING'} "
      f"(worst key-cell degradation {worst:.0f}% {'<' if robust else '>='} 25%)")
  log(f"total wall {time.time()-t0:.0f}s")

if __name__ == '__main__':
    main()
