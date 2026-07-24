# Decisive test: does blind SPECROT crack the null when effective diversity is high?
# Vary T and temporal correlation (tau). float64. Report oracle vs SPECROT null-err + timing.
import time, numpy as np, torch
from fog_common import *
import fog_tracker as ft
def log(*a): print(*a, flush=True)
dev = 'cuda' if torch.cuda.is_available() else 'cpu'
log("device", dev, "dtype", ft.DTYPE)
n_side=32; N=n_side*n_side; M=64; lam_x=1e-6
x_true = make_scene(n_side)

def run(c, sig_w, tau, T, n_restart=12, n_als=80):
    rng = np.random.default_rng(20260723 + T + (0 if tau is None else tau))
    P = make_patterns(M,N,"gauss",rng); Pi,rangeP = projectors(P)
    null_true = x_true-rangeP(x_true); nrm0=np.linalg.norm(null_true)
    U=smooth_basis(c,n_side); d_w=U.shape[1]
    Z,rho = ou_coeffs(T,d_w,sig_w,tau,rng); W=np.clip(1.0+Z@U.T,0.05,None)
    mu=np.einsum('mn,tn->tm',P,W*x_true); B=mu.copy(); wts=np.ones_like(B); Rd=np.ones_like(B)*1e-8
    def m(x): return null_metric(x,x_true,rangeP,null_true,nrm0)[0]
    xo=ft.solve_oracle(P,W,B,wts,lam_x,dev); neo=m(xo)
    t0=time.time()
    xs,res=ft.solve_spectral_rot(P,U,B,wts,sig_w,lam_x,n_als=n_als,n_restart=n_restart,dev=dev)
    nes=m(xs); dt=time.time()-t0
    log(f"c={c} d_w={d_w:3d} sig={sig_w} tau={str(tau):>4} T={T:4d} rho={rho:.3f} | "
        f"oracle {neo:.3f} | SPECROT {nes:.3f} resid {res:.2e} | {dt:.1f}s")
    return neo, nes

log("--- iid medium (max diversity), vary T ---")
for T in [64,128,256,512,1024]:
    run(4, 0.30, None, T)
log("--- d_w=16, tau sweep at T=512 ---")
for tau in [4,16,64]:
    run(4, 0.30, tau, 512)
log("--- d_w=64=M, iid, vary T ---")
for T in [256,512,1024]:
    run(8, 0.30, None, T)
