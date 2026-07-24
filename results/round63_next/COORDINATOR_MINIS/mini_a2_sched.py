# Mini A2: slot-mode chronology — ordered vs randomized schedule, covariance route
import numpy as np, torch
exec(open('mini_cov_predict.py').read().split('print("== 1.')[0])  # setup: x, U0, P_rand, N,M,d, SIG_F, dev, nullnmse
TAU_S=48.0; LAGS=list(range(1,9)); PHOT=1e5
def gen_slots(T_slots,sched,seed=0):
    r=np.random.default_rng(seed); rho=np.exp(-1/TAU_S); sc=SIG_F*np.sqrt(N/d)
    z=sc*r.standard_normal(d); B=np.zeros(T_slots)
    Z=np.zeros((T_slots,d))
    for s in range(T_slots): Z[s]=z; z=rho*z+np.sqrt(1-rho**2)*sc*r.standard_normal(d)
    W=np.exp(Z@U0.T-0.5*SIG_F**2)
    mu=np.einsum('sn,n->s',W*x[None,:]*0+W*x[None,:],np.ones(N))  # placeholder
    mu=np.array([P_rand[sched[s]]@(W[s]*x) for s in range(T_slots)])
    sc_p=PHOT/mu.mean(); B=r.poisson(mu*sc_p)/sc_p
    return B
def sched_ordered(T_slots): return np.array([s%M for s in range(T_slots)])
def sched_random(T_slots,seed=1):
    r=np.random.default_rng(seed); out=[]
    for _ in range(T_slots//M): out+=list(r.permutation(M))
    return np.array(out)
def cov_slot_estimate(B,sched,iters=3000):
    rho=np.exp(-1/TAU_S); sc2=SIG_F**2*N/d
    mean_b=np.array([B[sched==i].mean() for i in range(M)])
    db=B-mean_b[sched]
    num=np.zeros((M,M)); den=np.zeros((M,M)); cnt=np.zeros((M,M))
    T_slots=len(B)
    for l in LAGS:
        rl=rho**l
        i_=sched[:-l]; j_=sched[l:]; p=db[:-l]*db[l:]
        np.add.at(num,(i_,j_),rl*p); np.add.at(den,(i_,j_),rl*rl); np.add.at(cnt,(i_,j_),1)
    num=num+num.T; den=den+den.T; cnt=cnt+cnt.T
    mask=cnt>=5
    G=np.zeros((M,M)); G[mask]=num[mask]/np.maximum(den[mask],1e-12)
    cov_frac=mask.sum()/(M*M)
    x0,_,_,_=np.linalg.lstsq(P_rand,mean_b,rcond=None); x0=np.clip(x0,1e-3,None)
    tP=torch.tensor(P_rand,device=dev,dtype=torch.float64); tU=torch.tensor(U0,device=dev,dtype=torch.float64)
    tG=torch.tensor(G/sc2,device=dev,dtype=torch.float64); tm=torch.tensor(mask,device=dev)
    xp=torch.tensor(x0,device=dev,dtype=torch.float64,requires_grad=True)
    opt=torch.optim.Adam([xp],lr=2e-2)
    for it in range(iters):
        opt.zero_grad(); xr=torch.relu(xp)
        H=(tP*xr[None,:])@tU
        R=(H@H.T-tG); loss=(R[tm]**2).sum()
        loss.backward(); opt.step()
        if it==2000:
            for g_ in opt.param_groups: g_['lr']=5e-3
    return np.clip(xp.detach().cpu().numpy(),0,None), cov_frac
for T_slots in [12288, 49152]:
    for name,sch_fn in [("ordered",sched_ordered),("randomized",sched_random)]:
        v=[];cf=0
        for s in range(3):
            sched=sch_fn(T_slots) if name=="ordered" else sched_random(T_slots,seed=10+s)
            B=gen_slots(T_slots,sched,seed=200+s)
            xh,cf=cov_slot_estimate(B,sched)
            v.append(nullnmse(xh,P_rand))
        print(f"T_slots={T_slots:6d} {name:11s} coverage {cf:.2f}  null-NMSE median {np.median(v):.3f} (all {[round(q,3) for q in v]})")
