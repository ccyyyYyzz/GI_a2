# F4 duel: fixed-bank covariance route vs FRESH known patterns (medium as noise), equal budget
import numpy as np, torch
exec(open('mini_cov_predict.py').read().split('print("== 1.')[0])
SIG_F=0.30; TAU=8.0; PHOT=1e5
Texp=48*2048   # total exposures, both arms
rho=np.exp(-1/TAU); sc=SIG_F*np.sqrt(N/d)
def medium_path(T,seed):
    r=np.random.default_rng(seed); z=sc*r.standard_normal(d); Z=np.zeros((T,d))
    for t in range(T): Z[t]=z; z=rho*z+np.sqrt(1-rho**2)*sc*r.standard_normal(d)
    return np.exp(Z@U0.T-0.5*SIG_F**2)
res={"fixed":[], "fresh":[], "fresh_nullpart":[]}
for s in range(3):
    r=np.random.default_rng(700+s)
    # ---- ARM 1: fixed bank M=48, T=2048 epochs (moment-matching covariance route)
    W=medium_path(2048,800+s)
    mu=np.einsum('mn,tn->tm',P_rand,W*x); scp=PHOT/mu.mean()
    B=r.poisson(mu*scp)/scp
    # (reuse estimator from mini_cov_predict via run()-internals: inline compact version)
    Bc=B-B.mean(0,keepdims=True); T=2048
    S1=(Bc[:-1].T@Bc[1:])/(T-1); S2=(Bc[:-2].T@Bc[2:])/(T-2)
    G=((0.5*(S1+S1.T))/rho+(0.5*(S2+S2.T))/rho**2)/2/(SIG_F**2*N/d)
    x0,_,_,_=np.linalg.lstsq(P_rand,B.mean(0),rcond=None); x0=np.clip(x0,1e-3,None)
    tP=torch.tensor(P_rand,device=dev,dtype=torch.float64); tU=torch.tensor(U0,device=dev,dtype=torch.float64)
    tG=torch.tensor(G,device=dev,dtype=torch.float64)
    xp=torch.tensor(x0,device=dev,dtype=torch.float64,requires_grad=True)
    opt=torch.optim.Adam([xp],lr=2e-2)
    for it in range(4000):
        opt.zero_grad(); xr=torch.relu(xp)
        H=(tP*xr[None,:])@tU
        ((H@H.T-tG)**2).sum().backward(); opt.step()
        if it==3000:
            for g_ in opt.param_groups: g_['lr']=5e-3
    xf=np.clip(xp.detach().cpu().numpy(),0,None)
    sg=np.dot(xf,x)/np.dot(xf,xf); e_fixed=np.linalg.norm(sg*xf-x)/np.linalg.norm(x)
    # ---- ARM 2: fresh random binary pattern EVERY exposure; medium changes per epoch of 48 slots
    # same medium path budget: 2048 medium states, 48 fresh patterns within each state
    Pf_big_err_num=0; AtA=np.zeros((N,N)); Atb=np.zeros(N)
    W2=medium_path(2048,800+s)   # same statistics
    rb=np.random.default_rng(900+s)
    for t in range(2048):
        Pf=(rb.random((48,N))<0.5).astype(float)
        mu2=Pf@(W2[t]*x); B2=rb.poisson(mu2*scp)/scp
        AtA+=Pf.T@Pf; Atb+=Pf.T@B2
    lam=1e-3*np.trace(AtA)/N
    xg=np.linalg.solve(AtA+lam*np.eye(N),Atb); xg=np.clip(xg,0,None)
    sg2=np.dot(xg,x)/np.dot(xg,xg); e_fresh=np.linalg.norm(sg2*xg-x)/np.linalg.norm(x)
    # fresh arm error in the FIXED bank's null space (the content fixed-mean-route can never see)
    Pi=np.linalg.pinv(P_rand); proj=lambda v: v-Pi@(P_rand@v)
    e_fresh_null=np.linalg.norm(proj(sg2*xg-x))/np.linalg.norm(proj(x))
    e_fixed_null=np.linalg.norm(proj(sg*xf-x))/np.linalg.norm(proj(x))
    res["fixed"].append((e_fixed,e_fixed_null)); res["fresh"].append((e_fresh,e_fresh_null))
print("=== F4 duel (equal exposures 48x2048, equal photons 1e5, same medium stats) ===")
for k in ("fixed","fresh"):
    a=np.array(res[k])
    print(f"{k:6s}: full-scene NMSE median {np.median(a[:,0]):.3f} {[round(q,3) for q in a[:,0]]} | fixed-bank-null err median {np.median(a[:,1]):.3f}")
print("note: fresh arm has NO repeated pattern pairs -> medium covariance/statistics NOT estimable (no second product)")
