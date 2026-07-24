# What null-NMSE is ACHIEVABLE from a warm (truth) start? Bounds the model/optimizer ceiling.
# Isolates: Q-drift vs subspace-truncation vs lognormal-mismatch. Diagnostic only (warm).
import numpy as np, torch
from fog_common import *
import fog_tracker as ft
def log(*a): print(*a, flush=True)
dev='cuda' if torch.cuda.is_available() else 'cpu'
n_side=32; N=n_side*n_side; x_true=make_scene(n_side); M=96; d=48; T=128

def setup(sigma_f, medium, seed=3):
    rng=np.random.default_rng(2024+seed)
    P=make_patterns(M,N,"binary",rng); Pi,rangeP=projectors(P)
    null_true=x_true-rangeP(x_true); nrm0=np.linalg.norm(null_true)
    U=dct_basis(d,n_side); coeff_sd=sigma_f*np.sqrt(N/d)
    if medium=='lognormal': W,Z,rho,_=lognormal_medium(U,T,sigma_f,None,rng)
    else: Z=coeff_sd*rng.standard_normal((T,d)); rho=0.0; W=np.clip(1.0+Z@U.T,0.05,None)
    return dict(P=P,U=U,W=W,Z=Z,B=np.einsum('mn,tn->tm',P,W*x_true),rangeP=rangeP,
                null_true=null_true,nrm0=nrm0,coeff_sd=coeff_sd,rho=rho)
def nm(x,s): return null_metric(x,x_true,s['rangeP'],s['null_true'],s['nrm0'])[0]**2

for medium in ['linear','lognormal']:
    s=setup(0.30, medium); P,U,W,Z,B=s['P'],s['U'],s['W'],s['Z'],s['B']
    wts=np.ones_like(B); csd=s['coeff_sd']
    log(f"\n=== medium={medium} sf=0.30 M={M} d={d} T={T} ===")
    # 1) oracle: true W
    xo=ft.solve_oracle(P,W,B,wts,1e-6,dev); log(f" oracle (true W)            NMSE {nm(xo,s):.4f}")
    # 2) linear-approx oracle: W=1+UZ_true (quantifies lognormal mismatch alone)
    Wl=np.clip(1.0+Z@U.T,0.05,None)
    xl=ft.solve_oracle(P,Wl,B,wts,1e-6,dev); log(f" linear-approx oracle       NMSE {nm(xl,s):.4f}")
    # 3) warm full-Z LINEAR Adam from Z_true (linear-model ceiling)
    xa=ft.joint_adam_reg(P,U,B,wts,0.0,csd,1e-4,steps=4000,lr=1e-2,dev=dev,Z_init=Z,nonneg_x=True)
    log(f" warm full-Z linear Adam    NMSE {nm(xa,s):.4f}")
    # 4) warm full-Z LOGNORMAL Adam from Z_true (model-matched ceiling)
    xg=ft.joint_adam_reg(P,U,B,wts,0.0,csd,1e-4,steps=4000,lr=1e-2,dev=dev,Z_init=Z,nonneg_x=True,lognormal=True)
    log(f" warm full-Z lognormal Adam NMSE {nm(xg,s):.4f}")
