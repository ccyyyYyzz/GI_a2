# Beyond-modulator-band duel: pattern band-limited (DMD limit), medium band EXCEEDS it.
import numpy as np, torch
exec(open('mini_cov_predict.py').read().split('print("== 1.')[0])
from scipy.fftpack import idct,dct
D2=lambda a: dct(dct(a.reshape(n,n),axis=0,norm='ortho'),axis=1,norm='ortho')
I2=lambda A: idct(idct(A,axis=0,norm='ortho'),axis=1,norm='ortho').ravel()
SIG_F=0.30; TAU=8.0; PHOT=1e5; rho=np.exp(-1/TAU)
PB=3   # pattern band: modulator can only make fx,fy <= 3
# medium band: 1<=i+j<=6 (fine speckle beyond modulator band)
band=[(i,j) for i in range(7) for j in range(7) if 1<=i+j<=6]
mm=[]
for (i,j) in band:
    g=np.zeros((n,n)); g[i,j]=1.; mm.append(I2(g))
Ub=np.linalg.qr(np.array(mm).T)[0]; db=Ub.shape[1]
sc=SIG_F*np.sqrt(N/db); sc2=SIG_F**2*N/db; Kg=sc2*(Ub@Ub.T)
# scene: energy INSIDE pattern band + deliberately OUTSIDE it (test content at 4<=max freq<=7)
C=np.zeros((n,n)); rs=np.random.default_rng(4)
C[:PB+1,:PB+1]=rs.standard_normal((PB+1,PB+1))
hi=[(5,2),(2,5),(4,4),(6,1),(1,6),(5,5),(7,0),(0,7)]
for (i,j) in hi: C[i,j]=2.0*rs.choice([-1,1])
xs=I2(C); xs=xs-xs.min(); xs=xs/xs.max()
x=xs  # override scene
hi_mask=np.zeros((n,n),bool)
for (i,j) in hi: hi_mask[i,j]=True
# band-limited nonneg random patterns (both arms share this hardware limit)
def bl_patterns(M_,seed):
    rp=np.random.default_rng(seed); ps=[]
    for _ in range(M_):
        Cp=np.zeros((n,n)); Cp[:PB+1,:PB+1]=rp.standard_normal((PB+1,PB+1)); Cp[0,0]=0
        f=I2(Cp); f/=np.abs(f).max(); ps.append(0.5+0.45*f)
    return np.array(ps)
M_=24; T=2048
Pfix=bl_patterns(M_,10)
def med(T,seed):
    r=np.random.default_rng(seed); z=sc*r.standard_normal(db); W=np.zeros((T,N))
    for t in range(T): W[t]=np.exp(Ub@z-0.5*SIG_F**2); z=rho*z+np.sqrt(1-rho**2)*sc*r.standard_normal(db)
    return W
def hi_err(xh):
    s=np.dot(xh,x)/max(np.dot(xh,xh),1e-12); E=D2(s*xh-x)**2; X0=D2(x)**2
    return E[hi_mask].sum()/X0[hi_mask].sum()
resF=[];resG=[]
for s in range(3):
    W=med(T,600+s); r=np.random.default_rng(650+s)
    # fixed bank + covariance
    mu=np.einsum('mn,tn->tm',Pfix,W*x); scp=PHOT/mu.mean()
    B=r.poisson(mu*scp)/scp
    Bc=B-B.mean(0,keepdims=True)
    S1=(Bc[:-1].T@Bc[1:])/(T-1); S2=(Bc[:-2].T@Bc[2:])/(T-2)
    G=((0.5*(S1+S1.T))/rho+(0.5*(S2+S2.T))/rho**2)/2
    x0,_,_,_=np.linalg.lstsq(Pfix,B.mean(0),rcond=None); x0=np.clip(x0,1e-3,None)
    tP=torch.tensor(Pfix,device=dev,dtype=torch.float64); tK=torch.tensor(Kg,device=dev,dtype=torch.float64)
    tG=torch.tensor(G,device=dev,dtype=torch.float64)
    xp=torch.tensor(x0,device=dev,dtype=torch.float64,requires_grad=True)
    opt=torch.optim.Adam([xp],lr=2e-2)
    for it in range(4000):
        opt.zero_grad(); xr=torch.relu(xp)
        XP=tP*xr[None,:]
        ((XP@tK@XP.T-tG)**2).sum().backward(); opt.step()
        if it==3000:
            for g_ in opt.param_groups: g_['lr']=5e-3
    resF.append(hi_err(np.clip(xp.detach().cpu().numpy(),0,None)))
    # fresh band-limited patterns + mean route (same hardware limit!)
    AtA=np.zeros((N,N)); Atb=np.zeros(N); rb=np.random.default_rng(680+s)
    for t in range(T):
        Pf=bl_patterns(M_,7000+1000*s+t)
        mu2=Pf@(W[t]*x); B2=rb.poisson(mu2*scp)/scp
        AtA+=Pf.T@Pf; Atb+=Pf.T@B2
    lam=1e-3*np.trace(AtA)/N
    xg=np.clip(np.linalg.solve(AtA+lam*np.eye(N),Atb),0,None)
    resG.append(hi_err(xg))
print(f"=== beyond-modulator-band content (freqs outside pattern band <= {PB}) ===")
print(f"fixed-bank covariance : rel err median {np.median(resF):.3f} {[round(q,3) for q in resF]}")
print(f"fresh patterns (mean) : rel err median {np.median(resG):.3f} {[round(q,3) for q in resG]}")
print("(1.0 = recovers nothing; the fresh arm is physics-limited to the pattern band)")
