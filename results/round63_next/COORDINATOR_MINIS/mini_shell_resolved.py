# Shell-resolved beyond-band recovery: error vs distance from the modulator band edge
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
# scene: in-band base + witnesses on shells max(i,j)=4..7
C=np.zeros((n,n)); rs=np.random.default_rng(4)
C[:PB+1,:PB+1]=rs.standard_normal((PB+1,PB+1))
shells={4:[],5:[],6:[],7:[]}
for i in range(8):
    for j in range(8):
        m=max(i,j)
        if m in shells and min(i,j)<=m:
            if (i,j) not in [(a,b) for a in range(PB+1) for b in range(PB+1)]:
                shells[m].append((i,j))
for m in shells:
    for (i,j) in shells[m]: C[i,j]=1.2*rs.choice([-1,1])
xs=I2(C); xs=xs-xs.min(); x=xs/xs.max()
def bl_pats(M_,seed):
    rp=np.random.default_rng(seed); ps=[]
    for _ in range(M_):
        Cp=np.zeros((n,n)); Cp[:PB+1,:PB+1]=rp.standard_normal((PB+1,PB+1)); Cp[0,0]=0
        f=I2(Cp); f/=np.abs(f).max(); ps.append(0.5+0.45*f)
    return np.array(ps)
Pfix=bl_pats(24,10); T=2048
def shell_errs(xh):
    s=np.dot(xh,x)/max(np.dot(xh,xh),1e-12); E=D2(s*xh-x)**2; X0=D2(x)**2
    out={}
    for m,locs in shells.items():
        num=sum(E[i,j] for (i,j) in locs); den=sum(X0[i,j] for (i,j) in locs)
        out[m]=num/max(den,1e-12)
    return out
acc={m:[] for m in shells}
for s in range(3):
    r=np.random.default_rng(600+s)
    z=sc*r.standard_normal(db); W=np.zeros((T,N))
    for t in range(T): W[t]=np.exp(Ub@z-0.5*SIG_F**2); z=rho*z+np.sqrt(1-rho**2)*sc*r.standard_normal(db)
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
    se=shell_errs(np.clip(xp.detach().cpu().numpy(),0,None))
    for m in shells: acc[m].append(se[m])
print("shell (max freq) | Delta-k from band edge | rel err (3 seeds)")
for m in sorted(shells):
    print(f"  {m}              |  {m-PB}                     | median {np.median(acc[m]):.3f}  {[round(q,3) for q in acc[m]]}")
