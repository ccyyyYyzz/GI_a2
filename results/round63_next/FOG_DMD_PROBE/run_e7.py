# E7 -- statistical efficiency (E7a) + exposure economics (E7b). Two-stage MLE:
#   moment-matching (E5a, drift-free global basin) -> marginal-likelihood EM (kalman_em, efficient).
# Incremental JSON checkpointing (another session intermittently clears GPU python procs).
import time, json, sys, numpy as np, torch
from fog_common import *
import fog_tracker as ft
def log(*a): print(*a, flush=True)
dev = 'cuda' if torch.cuda.is_available() else 'cpu'
OUT = 'E7_results.json'
R = {}
def save():
    json.dump(R, open(OUT, 'w'), indent=2)

def make_cell(ns, M, d):
    x = make_scene(ns); N = ns * ns
    return x, N

def gen(ns, N, M, d, T, m, photons, sf, tau, seed):
    """m = patterns displayed per epoch (m<=M). Returns data for the m-subset code bank."""
    rng = np.random.default_rng(4000 + seed + T + m)
    x = make_scene(ns)
    P = make_patterns(m, N, "binary", rng)              # m x N code bank actually used
    Pi, rangeP = projectors(P)
    nt = x - rangeP(x); nrm0 = np.linalg.norm(nt)
    U = dct_basis(d, ns)
    W, Z, rho, csd = lognormal_medium(U, T, sf, tau, rng)
    mu = np.einsum('mn,tn->tm', P, W * x)
    if photons is None:
        B = mu.copy(); Rd = np.ones_like(B) * 1e-8
    else:
        s = photons / mu.mean(); B = rng.poisson(np.clip(mu, 0, None) * s) / s
        Rd = np.clip(mu, 1e-9, None) / s
    return dict(x=x, P=P, U=U, B=B, Rd=Rd, rho=rho, csd=csd, rangeP=rangeP, nt=nt, nrm0=nrm0)

def nm(x, c): return null_metric(x, c['x'], c['rangeP'], c['nt'], c['nrm0'])[0] ** 2

def moment(c, n_starts=3):
    r = ft.cov_estimate(c['P'], c['U'], c['B'], c['csd'], c['rho'], n_starts=n_starts, steps=4000, dev=dev)
    return r['per_start_x'][int(np.argmin(r['objs']))]

def mle(c, x_init, n_em=25):
    return ft.kalman_em(c['P'], c['U'], c['B'], c['Rd'], c['rho'], c['csd'], 1e-4, n_em=n_em,
                        dev=dev, x_init=x_init)

def e7a(tag, ns, M, d, Tlist, seeds, sf=0.15, tau=16.0, photons=1e5):
    R.setdefault('E7a', {})[tag] = []
    N = ns * ns
    for T in Tlist:
        mo, ml, tmo, tml = [], [], [], []
        for sd in seeds:
            c = gen(ns, N, M, d, T, M, photons, sf, tau, sd)
            t0 = time.time(); xm = moment(c); tmo.append(time.time() - t0); mo.append(nm(xm, c))
            t0 = time.time(); xL = mle(c, xm); tml.append(time.time() - t0); ml.append(nm(xL, c))
        row = dict(T=T, moment_nmse=float(np.median(mo)), mle_nmse=float(np.median(ml)),
                   sec_moment=float(np.median(tmo)), sec_mle=float(np.median(tml)))
        R['E7a'][tag].append(row); save()
        log(f"[{tag}] T={T:5d}: moment {row['moment_nmse']:.3f} | MLE {row['mle_nmse']:.3f} "
            f"| t_mom {row['sec_moment']:.1f}s t_mle {row['sec_mle']:.1f}s")

def e7b(tag, ns, M, d, total_exposures, m_list, seeds, sf=0.15, tau=16.0, photons=1e5):
    R.setdefault('E7b', {})[tag] = []
    N = ns * ns
    for m in m_list:
        T = total_exposures // m
        vals = []
        for sd in seeds:
            c = gen(ns, N, m, d, T, m, photons, sf, tau, sd)   # m patterns per epoch, T epochs
            xm = moment(c); xL = mle(c, xm); vals.append(nm(xL, c))
        cap = (min(T, d + 1)) * (m - d) / N if m > d else 0.0
        row = dict(m=m, T=T, mT=m * T, chi=cap, mle_nmse=float(np.median(vals)))
        R['E7b'][tag].append(row); save()
        log(f"[{tag} exposure] m={m:3d} T={T:6d} (mT={m*T}) chi={cap:.2f}: MLE null-NMSE {row['mle_nmse']:.3f}")

if __name__ == '__main__':
    t0 = time.time()
    log("=== E7a MINI (16x16, M=48, d=24): MLE vs moment, crossing T ===")
    e7a('mini', 16, 48, 24, [128, 256, 512, 1024], list(range(5)))
    log("\n=== E7a CONFIRM (32x32, M=96, d=48) ===")
    e7a('full32', 32, 96, 48, [512, 1024], list(range(3)))
    log("\n=== E7b exposure economics (mini d=24, fixed total exposures = 48*4096) ===")
    # m must exceed d=24 for blind capacity; sweep the m=d boundary + two viable points
    e7b('mini', 16, 48, 24, 48 * 4096, [24, 36, 48], list(range(3)))
    R['meta'] = dict(wall_s=time.time() - t0, note='two-stage MLE = moment-matching -> kalman_em EM')
    save()
    # crossing summary
    def cross(rows, key):
        below = [r['T'] for r in rows if r[key] <= 0.25]
        return min(below) if below else None
    mi = R['E7a']['mini']
    log(f"\n=== E7a VERDICT (mini) ===")
    log(f"moment crosses 0.25 at T~{cross(mi,'moment_nmse')} | MLE crosses at T~{cross(mi,'mle_nmse')}")
    log(f"total wall {time.time()-t0:.0f}s -> {OUT}")
