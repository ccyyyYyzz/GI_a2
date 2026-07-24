import time, numpy as np, torch
from fog_common import *
import fog_tracker as ft
dev = 'cuda' if torch.cuda.is_available() else 'cpu'
def log(*a): print(*a, flush=True)
n_side=32; N=n_side*n_side; M=64; lam_x=1e-6
x_true = make_scene(n_side)
rng = np.random.default_rng(20260723)
P = make_patterns(M,N,"gauss",rng); Pi,rangeP = projectors(P)
null_true = x_true-rangeP(x_true); nrm0=np.linalg.norm(null_true)
c=4; U=smooth_basis(c,n_side); d_w=U.shape[1]; sig_w=0.30; tau=16; T=128
Z,rho = ou_coeffs(T,d_w,sig_w,tau,rng); W=np.clip(1.0+Z@U.T,0.05,None)
mu=np.einsum('mn,tn->tm',P,W*x_true); B=mu.copy(); wts=np.ones_like(B); Rd=np.ones_like(B)*1e-8
def m(x): return null_metric(x,x_true,rangeP,null_true,nrm0)[0]
torch.cuda.synchronize() if dev=='cuda' else None
t0=time.time(); xo=ft.solve_oracle(P,W,B,wts,lam_x,dev); torch.cuda.synchronize() if dev=='cuda' else None
log(f"oracle    null {m(xo):.3f}  [{time.time()-t0:.2f}s]")
t0=time.time(); xs,res=ft.solve_spectral_rot(P,U,B,wts,sig_w,lam_x,n_als=60,n_restart=10,dev=dev)
torch.cuda.synchronize() if dev=='cuda' else None
log(f"SPECROT   null {m(xs):.3f}  resid {res:.3e}  [{time.time()-t0:.2f}s]")
t0=time.time(); xr,_=ft.solve_spectral_rot(P,U,B,wts,sig_w,lam_x,n_als=60,n_restart=10,dev=dev,refine_em=15,rho=rho,R_diag_np=Rd)
log(f"SPECROT+EM null {m(xr):.3f}  [{time.time()-t0:.2f}s]")
