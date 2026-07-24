# Aperture law test v2: BAND-LIMITED RANDOM patterns (incoherent but spectrally confined)
import numpy as np, torch
exec(open('mini_cov_predict.py').read().split('print("== 1.')[0])
from scipy.fftpack import idct,dct
SIG_F=0.30; PHOT=1e12; TAU=8.0
D2 = lambda a: dct(dct(a.reshape(n,n),axis=0,norm='ortho'),axis=1,norm='ortho')
I2 = lambda A: idct(idct(A,axis=0,norm='ortho'),axis=1,norm='ortho').ravel()
# medium band: annulus 1<=fx+fy<=2
band=[(i,j) for i in range(4) for j in range(4) if 1<=i+j<=2]
Ub=np.linalg.qr(np.array([I2((lambda g:( [g.__setitem__((i,j),1.),g][1]))(np.zeros((n,n)))) for (i,j) in band]).T)[0]
db=Ub.shape[1]; sc2=SIG_F**2*N/db; Kg=sc2*(Ub@Ub.T)
# band-limited random patterns: support fx,fy <= 5 (plus DC), random coeffs, nonneg
rngp=np.random.default_rng(9); M_=24; pats=[]
pmask=np.zeros((n,n),bool)
for i in range(6):
    for j in range(6): pmask[i,j]=True
for _ in range(M_):
    C=np.zeros((n,n)); C[pmask]=rngp.standard_normal(pmask.sum()); C[0,0]=0
    f=I2(C); f=f/np.abs(f).max()
    pats.append(0.5+0.45*f)
Pp=np.array(pats)
# predicted coverage: pattern support (<=5) Minkowski-sum medium band (<=2) => fx,fy<=7 roughly; build exact sum-set
cov=np.zeros((n,n),bool)
for pi in range(6):
    for pj in range(6):
        for (ki,kj) in band:
            for si in (pi+ki,abs(pi-ki)):
                for sj in (pj+kj,abs(pj-kj)):
                    if si<n and sj<n: cov[si,sj]=True
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
ei=[];eo=[]
for s in range(3):
    B=genE(4096,seed=500+s); xh=est(B)
    sg=np.dot(xh,x)/max(np.dot(xh,xh),1e-12); xh=sg*xh
    E=D2(xh-x)**2; X0=D2(x)**2
    m=cov&(X0>1e-6); o=(~cov)&(X0>1e-6)
    ei.append(E[m].sum()/X0[m].sum()); eo.append(E[o].sum()/X0[o].sum())
print(f"coverage fraction {cov.mean():.2f}  scene-energy in-coverage {D2(x)[cov].__pow__(2).sum()/ (D2(x)**2).sum():.2f}")
print(f"IN-coverage  rel err: median {np.median(ei):.3f} {[round(q,3) for q in ei]}")
print(f"OUT-coverage rel err: median {np.median(eo):.3f} {[round(q,3) for q in eo]}")
