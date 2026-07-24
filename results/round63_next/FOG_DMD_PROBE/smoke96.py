# Smoke: staged blind solver at primary M=96 cell + M=64 negative control. Data-only.
import time, numpy as np, torch
from fog_common import *
import fog_tracker as ft
def log(*a): print(*a, flush=True)
dev = 'cuda' if torch.cuda.is_available() else 'cpu'
log("device", dev, "dtype", ft.DTYPE)
n_side=32; N=n_side*n_side; x_true = make_scene(n_side)

def cell(M, d, T, sigma_f, tau, photons, seed=1, n_seeds=3):
    rng = np.random.default_rng(777+seed+M*10+d)
    P = make_patterns(M, N, "binary", rng)          # physical 0/1
    Pi, rangeP = projectors(P)
    null_true = x_true - rangeP(x_true); nrm0 = np.linalg.norm(null_true)
    nf = null_fraction(x_true, rangeP)
    U = dct_basis(d, n_side)
    W, Z, rho, coeff_sd = lognormal_medium(U, T, sigma_f, tau, rng)
    mu = np.einsum('mn,tn->tm', P, W*x_true)
    if photons is None:
        B=mu.copy(); wts=np.ones_like(B); Rd=np.ones_like(B)*1e-8
    else:
        s=photons/mu.mean(); Y=rng.poisson(np.clip(mu,0,None)*s)/s; B=Y
        wts=(s/np.clip(mu,1e-9,None)); Rd=np.clip(mu,1e-9,None)/s
    def nm(x): return null_metric(x,x_true,rangeP,null_true,nrm0)
    xo=ft.solve_oracle(P,W,B,wts,1e-6 if photons is None else 1e-3,dev); neo=nm(xo)[0]
    t0=time.time()
    res=ft.staged_blind(P,U,B,wts,rho,coeff_sd,Rd,1e-3,dev=dev,
                        seeds=tuple(range(n_seeds)),n_als=90,refine_em=20)
    neb=nm(res['x'])[0]; dt=time.time()-t0
    cap=(min(T,d+1))*(M-d)/N
    log(f"M={M} d={d:3d} T={T} sf={sigma_f} tau={str(tau):>4} ph={'clean' if photons is None else f'{photons:.0e}'}"
        f" | chi={cap:.2f} nullfrac={nf:.2f} | oracle NMSE={neo**2:.3f} | blind NMSE={neb**2:.3f}"
        f" (amp {neb:.3f}) | {dt:.1f}s")
    return neo, neb

log("\n== primary M=96 cell (clean) ==")
cell(96, 48, 128, 0.30, None, photons=None)
log("== negative control M=64 d=16 (theorem: blind impossible) ==")
cell(64, 16, 128, 0.30, None, photons=None)
