# P1: moment-tower oxygen — does 3rd-order cumulant matching unlock small-M?
import numpy as np, torch
exec(open('mini_cov_predict.py').read().split('print("== 1.')[0])
SIG_F=0.30; PHOT=1e12  # strong non-Gaussianity; effectively clean (oxygen test)
TAU=8.0
sc2=SIG_F**2*N/d
Kg=sc2*(U0@U0.T)  # pixel log-field covariance
def estimate_tower(B,P,use3,iters=4000,seed=0):
    rho=np.exp(-1/TAU); M_=P.shape[0]
    Bc=B-B.mean(0,keepdims=True); T=B.shape[0]
    S1=(Bc[:-1].T@Bc[1:])/(T-1); S2=(Bc[:-2].T@Bc[2:])/(T-2)
    G=( (0.5*(S1+S1.T))/rho + (0.5*(S2+S2.T))/rho**2 )/2
    C3=np.einsum('ti,tj,tk->ijk',Bc,Bc,Bc)/T if use3 else None
    x0,_,_,_=np.linalg.lstsq(P,B.mean(0),rcond=None); x0=np.clip(x0,1e-3,None)
    tP=torch.tensor(P,device=dev,dtype=torch.float64); tK=torch.tensor(Kg,device=dev,dtype=torch.float64)
    tG=torch.tensor(G,device=dev,dtype=torch.float64)
    nG=float((G**2).sum())
    if use3:
        tC3=torch.tensor(C3,device=dev,dtype=torch.float64); nC3=float((C3**2).sum())
    xp=torch.tensor(x0,device=dev,dtype=torch.float64,requires_grad=True)
    opt=torch.optim.Adam([xp],lr=2e-2)
    for it in range(iters):
        opt.zero_grad(); xr=torch.relu(xp)
        XP=tP*xr[None,:]                     # M x N  (P diag(x))
        G_mod=XP@tK@XP.T                     # 2nd order (linear-K model)
        loss=((G_mod-tG)**2).sum()/nG
        if use3:
            Bm=tK@XP.T                        # N x M
            # T3_pqr ~ K_pq K_pr + perms  ->  C3_ijk = sum_p (Pdiagx)_ip Bm_pj Bm_pk + 2 perms
            E1=torch.einsum('ip,pj,pk->ijk',XP,Bm,Bm)
            C3_mod=E1+E1.permute(1,0,2)+E1.permute(2,1,0)
            loss=loss+((C3_mod-tC3)**2).sum()/nC3
        loss.backward(); opt.step()
        if it==3000:
            for g_ in opt.param_groups: g_['lr']=5e-3
    return np.clip(xp.detach().cpu().numpy(),0,None)
T=4096
print(f"== P1 tower test: N={N}, d={d}, sig_f={SIG_F}, T={T}, clean ==")
for M_ in [12,20,48]:
    Pm=(np.random.default_rng(3).random((M_,N))<0.5).astype(float)
    eq2=M_*(M_+1)//2; eq3=M_*(M_+1)*(M_+2)//6
    v2=[];v3=[]
    for s in range(3):
        B,_=gen(T,U0,Pm,tau=TAU,seed=300+s)
        v2.append(nullnmse(estimate_tower(B,Pm,False,seed=s),Pm))
        v3.append(nullnmse(estimate_tower(B,Pm,True ,seed=s),Pm))
    print(f"M={M_:2d} (2nd-eqs {eq2:4d}, +3rd {eq3:5d} vs N={N}): 2nd-only {np.median(v2):.3f} {[round(q,3) for q in v2]} | +3rd {np.median(v3):.3f} {[round(q,3) for q in v3]}")
