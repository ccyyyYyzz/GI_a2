# Hole check: can FRESH patterns + COVARIANCE readout also recover beyond-band content?
import numpy as np, torch
exec(open('mini_cov_predict.py').read().split('print("== 1.')[0])
from scipy.fftpack import idct,dct
D2=lambda a: dct(dct(a.reshape(n,n),axis=0,norm='ortho'),axis=1,norm='ortho')
I2=lambda A: idct(idct(A,axis=0,norm='ortho'),axis=1,norm='ortho').ravel()
SIG_F=0.30; TAU=8.0; PHOT=1e5; rho=np.exp(-1/TAU); PB=3
band=[(i,j) for i in range(7) for j in range(7) if 1<=i+j<=6]
mm=[]
for (i,j) in band:
    g=np.zeros((n,n)); g[i,j]=1.; mm.append(I2(g))
Ub=np.linalg.qr(np.array(mm).T)[0]; db=Ub.shape[1]
sc=SIG_F*np.sqrt(N/db); sc2=SIG_F**2*N/db; Kg=sc2*(Ub@Ub.T)
C=np.zeros((n,n)); rs=np.random.default_rng(4)
C[:PB+1,:PB+1]=rs.standard_normal((PB+1,PB+1))
hi=[(5,2),(2,5),(4,4),(6,1),(1,6),(5,5),(7,0),(0,7)]
for (i,j) in hi: C[i,j]=2.0*rs.choice([-1,1])
xs=I2(C); xs=xs-xs.min(); x=xs/xs.max()
hi_mask=np.zeros((n,n),bool)
for (i,j) in hi: hi_mask[i,j]=True
def bl_pats(M_,seed):
    rp=np.random.default_rng(seed); ps=[]
    for _ in range(M_):
        Cp=np.zeros((n,n)); Cp[:PB+1,:PB+1]=rp.standard_normal((PB+1,PB+1)); Cp[0,0]=0
        f=I2(Cp); f/=np.abs(f).max(); ps.append(0.5+0.45*f)
    return np.array(ps)
def hi_err(xh):
    s=np.dot(xh,x)/max(np.dot(xh,xh),1e-12); E=D2(s*xh-x)**2; X0=D2(x)**2
    return E[hi_mask].sum()/X0[hi_mask].sum()
T=2048; S=24
res=[]
for sd in range(2):
    r=np.random.default_rng(760+sd); rb=np.random.default_rng(860+sd)
    z=sc*r.standard_normal(db); W=np.zeros((T,N))
    for t in range(T): W[t]=np.exp(Ub@z-0.5*SIG_F**2); z=rho*z+np.sqrt(1-rho**2)*sc*r.standard_normal(db)
    Plist=[]; Blist=[]
    for t in range(T):
        Pf=bl_pats(S,5000+100*sd+t)
        mu=Pf@(W[t]*x); scp=PHOT/max(mu.mean(),1e-9)
        Blist.append(rb.poisson(mu*scp)/scp); Plist.append(Pf)
    # mean-route init (in-band only)
    AtA=np.zeros((N,N)); Atb=np.zeros(N)
    for t in range(T): AtA+=Plist[t].T@Plist[t]; Atb+=Plist[t].T@Blist[t]
    lam=1e-3*np.trace(AtA)/N
    x0=np.clip(np.linalg.solve(AtA+lam*np.eye(N),Atb),1e-3,None)
    tK=torch.tensor(Kg,device=dev,dtype=torch.float64)
    tP=torch.stack([torch.tensor(p,device=dev,dtype=torch.float64) for p in Plist])   # T x S x N
    dB=torch.stack([torch.tensor(Blist[t]-Plist[t]@x0,device=dev,dtype=torch.float64) for t in range(T)])
    # residual target: within-epoch outer products; model a_i^T K a_j with a=p*(x); refit db against current x each iter is complex;
    # simpler: joint objective on raw products vs model, x free:
    xp=torch.tensor(x0,device=dev,dtype=torch.float64,requires_grad=True)
    tB=torch.stack([torch.tensor(Blist[t],device=dev,dtype=torch.float64) for t in range(T)])
    opt=torch.optim.Adam([xp],lr=2e-2)
    idx=np.arange(T)
    for it in range(2500):
        opt.zero_grad(); xr=torch.relu(xp)
        sub=torch.tensor(np.random.default_rng(it).choice(T,64,replace=False),device=dev)
        Ps=tP[sub]; Bs=tB[sub]
        mean_pred=torch.einsum('tsn,n->ts',Ps,xr)
        resid=Bs-mean_pred
        A=Ps*xr[None,None,:]                     # t s n
        G_pred=torch.einsum('tsn,nm,trm->tsr',A,tK,A)
        obs=resid[:,:,None]*resid[:,None,:]
        m=~torch.eye(S,dtype=torch.bool,device=dev)
        loss=((G_pred-obs)[:,m]**2).sum()
        loss.backward(); opt.step()
        if it==1800:
            for g_ in opt.param_groups: g_['lr']=5e-3
    res.append(hi_err(np.clip(xp.detach().cpu().numpy(),0,None)))
print(f"fresh+covariance beyond-band err: {[round(q,3) for q in res]} (mean-route wall = 1.000; fixed+cov at richer budget was 0.47-0.79)")
