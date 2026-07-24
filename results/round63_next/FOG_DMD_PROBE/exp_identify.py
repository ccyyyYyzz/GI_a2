# Decisive IDENTIFIABILITY probe: is blind failure non-identifiability or optimization?
# One clean config. Instrument the actual solutions & residuals at 4 points:
#   (a) ORACLE (W known)      -- data ceiling
#   (b) TRUTH-plugged x-solve at true W then read residual (== oracle) and null
#   (c) SPECROT cold          -- where blind lands + its residual & null magnitude
#   (d) SPECROT warm-started from TRUE Q -- does truth STAY a stable ALS min?
#   (e) residual of the fixed-A (range-only) x under true W and under best-blind W
import numpy as np, torch
from fog_common import *
import fog_tracker as ft
def log(*a): print(*a, flush=True)
dev = 'cuda' if torch.cuda.is_available() else 'cpu'
n_side=32; N=n_side*n_side; M=64; lam_x=1e-6
x_true = make_scene(n_side)
rng = np.random.default_rng(20260723)
P = make_patterns(M,N,"gauss",rng); Pi,rangeP = projectors(P)
null_true = x_true-rangeP(x_true); nrm0=np.linalg.norm(null_true)
c=4; U=smooth_basis(c,n_side); d_w=U.shape[1]; sig_w=0.30; tau=None; T=256
Z,rho = ou_coeffs(T,d_w,sig_w,tau,rng); W=np.clip(1.0+Z@U.T,0.05,None)
mu=np.einsum('mn,tn->tm',P,W*x_true); B=mu.copy(); wts=np.ones_like(B); Rd=np.ones_like(B)*1e-8
def full_resid(x, Wm):
    pred = np.einsum('mn,tn->tm', P, Wm*x); return float(((pred-B)**2).sum())
def diag(tag, x, Wm=None):
    xn = x - rangeP(x); ne,te = null_metric(x, x_true, rangeP, null_true, nrm0)
    r = full_resid(x, Wm) if Wm is not None else float('nan')
    log(f"  {tag:32s} null-err {ne:.3f} tot-err {te:.3f} ||xnull||={np.linalg.norm(xn):.3f} "
        f"||x||={np.linalg.norm(x):.2f} resid(givenW)={r:.3e}")

log(f"config iid T={T} d_w={d_w} sig={sig_w}  ||null_true||={nrm0:.3f} ||B||^2={np.sum(B**2):.2f}\n")

# (a) oracle
xo = ft.solve_oracle(P,W,B,wts,lam_x,dev); diag("ORACLE (true W)", xo, W)

# (c) SPECROT cold, report its recovered W too
xs, res = ft.solve_spectral_rot(P,U,B,wts,sig_w,lam_x,n_als=100,n_restart=16,dev=dev)
diag("SPECROT cold", xs, None)
log(f"     SPECROT reported data-residual(weighted-sum) = {res:.3e}  vs oracle-resid = {full_resid(xo,W):.3e}")

# (d) SPECROT warm from TRUE Q: seed ALS at the true medium. Does truth stay?
#     True Z lies in Zhat subspace: Q_true = pinv(Zhat) Z.  Run ALS from there.
b = (B - B.mean(0)); import numpy.linalg as la
_,_,Vh = la.svd(b.T, full_matrices=False); Zhat = Vh[:d_w].T   # (T,d_w)
Q_true = la.lstsq(Zhat, Z, rcond=None)[0]                       # (d_w,d_w) s.t. Zhat Q_true ~ Z
log(f"     ||Zhat Q_true - Z|| / ||Z|| = {la.norm(Zhat@Q_true - Z)/la.norm(Z):.3e} (subspace check)")
# manual ALS from Q_true (no restarts) to see stability
Pt=torch.tensor(P,device=dev,dtype=ft.DTYPE); Ut=torch.tensor(U,device=dev,dtype=ft.DTYPE)
Bt=torch.tensor(B,device=dev,dtype=ft.DTYPE); Zh=torch.tensor(Zhat,device=dev,dtype=ft.DTYPE)
PtP=Pt.t()@Pt; PtBt=Pt.t()@Bt.t()
Q=torch.tensor(Q_true,device=dev,dtype=ft.DTYPE)
for it in range(100):
    Wm=torch.clamp(1.0+(Zh@Q)@Ut.t(),min=0.05); x=ft.solve_x_clean(PtP,PtBt,Wm,lam_x)
    H=Pt@(Ut*x[:,None]); HtH=H.t()@H+1e-8*torch.eye(d_w,device=dev,dtype=ft.DTYPE)
    y=Bt-(Pt@x)[None,:]; ZtZ=Zh.t()@Zh
    Q=torch.linalg.solve(ZtZ, Zh.t()@y@H)@torch.linalg.inv(HtH)
xw=x.detach().cpu().numpy(); Ww=(1.0+(Zhat@Q.detach().cpu().numpy())@U.T)
diag("ALS warm from TRUE Q", xw, np.clip(Ww,0.05,None))

# (e) range-only x, residual under true W
x_rng = rangeP(x_true); diag("range-only x_true (null=0)", x_rng, W)
