# THE aperture-synthesis prediction test: recoverable freqs = pattern freqs (+/-) medium band?
import numpy as np, torch
exec(open('mini_cov_predict.py').read().split('print("== 1.')[0])
SIG_F=0.30; PHOT=1e12; TAU=8.0
# medium band: DCT modes with 1 <= fx+fy <= 2 (low annulus), d small
from scipy.fftpack import idct,dct
band=[(i,j) for i in range(4) for j in range(4) if 1<=i+j<=2]
modes=[]
for (i,j) in band:
    g=np.zeros((n,n)); g[i,j]=1.
    modes.append(idct(idct(g,axis=0,norm='ortho'),axis=1,norm='ortho').ravel())
Ub=np.linalg.qr(np.array(modes).T)[0]; db=Ub.shape[1]
# patterns: cos+sin at 8 chosen freqs, nonneg
pf=[(0,3),(3,0),(2,2),(4,1),(1,4),(5,3),(3,5),(6,0)]
pats=[]
for (fi,fj) in pf:
    c=np.cos(2*np.pi*(fi*xx+fj*yy)); s_=np.sin(2*np.pi*(fi*xx+fj*yy))
    pats+=[0.5*(1+c).ravel(),0.5*(1+s_).ravel()]
Pp=np.array(pats); Mp=Pp.shape[0]
# predicted coverage: DCT-index sum set {f +/- k} for f in pf, k in band (plus f itself via DC of w? w centered -> no)
cov_mask=np.zeros((n,n),bool)
for (fi,fj) in pf:
    for (ki,kj) in band:
        for si in (fi+ki,abs(fi-ki)):
            for sj in (fj+kj,abs(fj-kj)):
                if si<n and sj<n: cov_mask[si,sj]=True
sc2=SIG_F**2*N/db
Kg=sc2*(Ub@Ub.T)
def gen2(T,seed):
    r=np.random.default_rng(seed); rho=np.exp(-1/TAU); sc=SIG_F*np.sqrt(N/db)
    z=sc*r.standard_normal(db); Z=np.zeros((T,db))
    for t in range(T): Z[t]=z; z=rho*z+np.sqrt(1-rho**2)*sc*r.standard_normal(db)
    W=np.exp(Z@Ub.T-0.5*SIG_F**2)
    mu=np.array([Pp[k%Mp]@(W[t]*x) for t,k in zip(range(T),range(T))])  # cycle patterns per slot? use epoch mode:
    return None
# epoch mode buckets (all patterns per epoch)
def genE(T,seed):
    r=np.random.default_rng(seed); rho=np.exp(-1/TAU); sc=SIG_F*np.sqrt(N/db)
    z=sc*r.standard_normal(db); Z=np.zeros((T,db))
    for t in range(T): Z[t]=z; z=rho*z+np.sqrt(1-rho**2)*sc*r.standard_normal(db)
    W=np.exp(Z@Ub.T-0.5*SIG_F**2)
    return np.einsum('mn,tn->tm',Pp,W*x)
def est(B,iters=4000):
    rho=np.exp(-1/TAU)
    Bc=B-B.mean(0,keepdims=True); T=B.shape[0]
    S1=(Bc[:-1].T@Bc[1:])/(T-1); S2=(Bc[:-2].T@Bc[2:])/(T-2)
    G=((0.5*(S1+S1.T))/rho+(0.5*(S2+S2.T))/rho**2)/2
    x0,_,_,_=np.linalg.lstsq(Pp,B.mean(0),rcond=None); x0=np.clip(x0,1e-3,None)
    tP=torch.tensor(Pp,device=dev,dtype=torch.float64); tK=torch.tensor(Kg,device=dev,dtype=torch.float64)
    tG=torch.tensor(G,device=dev,dtype=torch.float64)
    xp=torch.tensor(x0,device=dev,dtype=torch.float64,requires_grad=True)
    opt=torch.optim.Adam([xp],lr=2e-2)
    for it in range(iters):
        opt.zero_grad(); xr=torch.relu(xp)
        XP=tP*xr[None,:]
        loss=((XP@tK@XP.T-tG)**2).sum()
        loss.backward(); opt.step()
        if it==3000:
            for g_ in opt.param_groups: g_['lr']=5e-3
    return np.clip(xp.detach().cpu().numpy(),0,None)
errs_in=[];errs_out=[]
for s in range(3):
    B=genE(4096,seed=400+s)
    xh=est(B)
    sgn=np.dot(xh,x)/max(np.dot(xh,xh),1e-12); xh=sgn*xh
    E=dct(dct((xh-x).reshape(n,n),axis=0,norm='ortho'),axis=1,norm='ortho')**2
    X0=dct(dct(x.reshape(n,n),axis=0,norm='ortho'),axis=1,norm='ortho')**2
    m=cov_mask& (X0>1e-6); o=(~cov_mask)&(X0>1e-6)
    errs_in.append(E[m].sum()/X0[m].sum()); errs_out.append(E[o].sum()/X0[o].sum())
print(f"predicted-coverage fraction of plane: {cov_mask.mean():.2f}")
print(f"IN-coverage  relative error: median {np.median(errs_in):.3f}  {[round(q,3) for q in errs_in]}")
print(f"OUT-coverage relative error: median {np.median(errs_out):.3f}  {[round(q,3) for q in errs_out]}")
print("PREDICTION: in << out  (aperture law holds)  |  in ~= out (law fails)")
