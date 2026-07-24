# F4c: crossover map — where (if anywhere) does fixed-bank covariance beat fresh patterns?
import numpy as np, torch
exec(open('mini_cov_predict.py').read().split('print("== 1.')[0])
TAU=8.0; rho=np.exp(-1/TAU)
def duel(sig_f,phot,seed):
    sc=sig_f*np.sqrt(N/d); sc2=sig_f**2*N/d
    r=np.random.default_rng(seed); rb=np.random.default_rng(seed+1000)
    T=2048
    z=sc*r.standard_normal(d); W=np.zeros((T,N))
    for t in range(T):
        W[t]=np.exp(U0@z-0.5*sig_f**2); z=rho*z+np.sqrt(1-rho**2)*sc*r.standard_normal(d)
    # fixed arm
    mu=np.einsum('mn,tn->tm',P_rand,W*x); scp=phot/mu.mean()
    B=r.poisson(np.clip(mu*scp,0,None)).astype(float)/scp
    Bc=B-B.mean(0,keepdims=True)
    S1=(Bc[:-1].T@Bc[1:])/(T-1); S2=(Bc[:-2].T@Bc[2:])/(T-2)
    G=((0.5*(S1+S1.T))/rho+(0.5*(S2+S2.T))/rho**2)/2/sc2
    x0,_,_,_=np.linalg.lstsq(P_rand,B.mean(0),rcond=None); x0=np.clip(x0,1e-3,None)
    tP=torch.tensor(P_rand,device=dev,dtype=torch.float64); tU=torch.tensor(U0,device=dev,dtype=torch.float64)
    tG=torch.tensor(G,device=dev,dtype=torch.float64)
    xp=torch.tensor(x0,device=dev,dtype=torch.float64,requires_grad=True)
    opt=torch.optim.Adam([xp],lr=2e-2)
    for it in range(3500):
        opt.zero_grad(); xr=torch.relu(xp)
        H=(tP*xr[None,:])@tU
        ((H@H.T-tG)**2).sum().backward(); opt.step()
        if it==2800:
            for g_ in opt.param_groups: g_['lr']=5e-3
    xf=np.clip(xp.detach().cpu().numpy(),0,None)
    sg=np.dot(xf,x)/max(np.dot(xf,xf),1e-12)
    e_fix=np.linalg.norm(sg*xf-x)/np.linalg.norm(x)
    # fresh arm
    AtA=np.zeros((N,N)); Atb=np.zeros(N)
    for t in range(T):
        Pf=(rb.random((48,N))<0.5).astype(float)
        mu2=Pf@(W[t]*x); B2=rb.poisson(np.clip(mu2*scp,0,None)).astype(float)/scp
        AtA+=Pf.T@Pf; Atb+=Pf.T@B2
    lam=1e-3*np.trace(AtA)/N
    xg=np.clip(np.linalg.solve(AtA+lam*np.eye(N),Atb),0,None)
    sg2=np.dot(xg,x)/max(np.dot(xg,xg),1e-12)
    e_fr=np.linalg.norm(sg2*xg-x)/np.linalg.norm(x)
    return e_fix,e_fr
print("sig_f  photons | fixed-cov  fresh | winner")
for sig_f in [0.3,0.6,1.0]:
    for phot in [1e3,1e4,1e5]:
        ef=[];eg=[]
        for s in range(2):
            a,b=duel(sig_f,phot,1234+17*s); ef.append(a); eg.append(b)
        mf,mg=np.median(ef),np.median(eg)
        print(f"{sig_f:4.1f}  {phot:8.0e} |   {mf:.3f}    {mg:.3f} | {'FIXED' if mf<mg else 'fresh'}")
