# Isolate the M=96 blind gap: optimization vs lognormal-model-mismatch.
import time, numpy as np, torch
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
    if medium=='lognormal':
        W,Z,rho,_=lognormal_medium(U,T,sigma_f,None,rng)
    else:  # linear clipped (old model), same coeff_sd
        Z=coeff_sd*rng.standard_normal((T,d)); rho=0.0; W=np.clip(1.0+Z@U.T,0.05,None)
    mu=np.einsum('mn,tn->tm',P,W*x_true); B=mu.copy(); wts=np.ones_like(B); Rd=np.ones_like(B)*1e-8
    return dict(P=P,U=U,W=W,Z=Z,B=B,wts=wts,Rd=Rd,rho=rho,coeff_sd=coeff_sd,
                rangeP=rangeP,null_true=null_true,nrm0=nrm0,Zc=None)

def nm(x,s): return null_metric(x,x_true,s['rangeP'],s['null_true'],s['nrm0'])[0]

def warm_Q_ceiling(s, lam=1e-4, iters=150):
    # DIAGNOSTIC CEILING ONLY: init ALS at the true medium mapping in Zhat basis.
    P=torch.tensor(s['P'],device=dev,dtype=ft.DTYPE); U=torch.tensor(s['U'],device=dev,dtype=ft.DTYPE)
    B=torch.tensor(s['B'],device=dev,dtype=ft.DTYPE)
    Bc=B-B.mean(0,keepdim=True); _,_,Vh=torch.linalg.svd(Bc.t(),full_matrices=False); Zh=Vh[:d].t()
    Zt=torch.tensor(s['Z'],device=dev,dtype=ft.DTYPE)
    Q=torch.linalg.lstsq(Zh,Zt).solution
    PtP=P.t()@P; PtBt=P.t()@B.t()
    for it in range(iters):
        W=torch.clamp(1.0+(Zh@Q)@U.t(),min=0.05); x=torch.clamp(ft.solve_x_clean(PtP,PtBt,W,lam),min=0.0)
        H=P@(U*x[:,None]); HtH=H.t()@H+1e-8*torch.eye(d,device=dev,dtype=ft.DTYPE)
        y=B-(P@x)[None,:]; Q=torch.linalg.solve(Zh.t()@Zh, Zh.t()@y@H)@torch.linalg.inv(HtH)
    return nm(x.detach().cpu().numpy(), s)

for medium in ['linear','lognormal']:
    for sf in [0.15,0.30]:
        s=setup(sf,medium)
        xo=ft.solve_oracle(s['P'],s['W'],s['B'],s['wts'],1e-6,dev); neo=nm(xo,s)
        r=ft.staged_blind(s['P'],s['U'],s['B'],s['wts'],s['rho'],s['coeff_sd'],s['Rd'],1e-3,
                          dev=dev,seeds=(0,1,2,3),n_als=120,refine_em=25)
        neb=nm(r['x'],s)
        wc=warm_Q_ceiling(s)
        log(f"medium={medium:9s} sf={sf} | oracle {neo**2:.3f} | blind-cold {neb**2:.3f} (amp {neb:.3f})"
            f" | warm-Q-ALS-ceiling {wc**2:.3f}")
