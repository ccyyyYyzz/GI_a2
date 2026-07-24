# Independent mini-scale (16x16) covariance-route test: verify + predict E5d/E6/R39
import numpy as np, torch
dev = 'cuda' if torch.cuda.is_available() else 'cpu'
rng = np.random.default_rng(11)
n=16; N=n*n; M=48; d=24; SIG_F=0.15; TAU=8.0; PHOT=1e5
# scene with strong null content
yy,xx=np.mgrid[0:n,0:n]/n
x=(np.exp(-((xx-.3)**2+(yy-.4)**2)/.03)+.8*np.exp(-((xx-.7)**2+(yy-.6)**2)/.02))
x[5:11,8:10]+=1.0; x[2:4,3:13]+=.7; x=(x/x.max()).ravel()
# smooth basis: 24 lowest non-DC 2D DCT modes
from scipy.fftpack import idct
modes=[]
fr=[(i,j) for i in range(6) for j in range(6) if (i,j)!=(0,0)]
fr.sort(key=lambda t:t[0]**2+t[1]**2)
for (i,j) in fr[:d]:
    g=np.zeros((n,n)); g[i,j]=1.
    m=idct(idct(g,axis=0,norm='ortho'),axis=1,norm='ortho')
    modes.append(m.ravel())
U0=np.linalg.qr(np.array(modes).T)[0]
def med_field(Z,U): 
    g=Z@U.T; return np.exp(g-0.5*SIG_F**2)
def gen(T,U,P,tau=TAU,seed=0):
    r=np.random.default_rng(seed); rho=np.exp(-1/tau)
    sc=SIG_F*np.sqrt(N/d)
    Z=np.zeros((T,d)); z=sc*r.standard_normal(d)
    for t in range(T): Z[t]=z; z=rho*z+np.sqrt(1-rho**2)*sc*r.standard_normal(d)
    W=med_field(Z,U)
    mu=np.einsum('mn,tn->tm',P,W*x)
    s=PHOT/mu.mean(); B=r.poisson(mu*s)/s
    return B,Z
def estimate(B,P,U,tau=TAU,iters=3000):
    rho=np.exp(-1/tau); sc2=(SIG_F**2*N/d)
    Bc=B-B.mean(0,keepdims=True); T=B.shape[0]
    S1=(Bc[:-1].T@Bc[1:])/(T-1); S2=(Bc[:-2].T@Bc[2:])/(T-2)
    S1=0.5*(S1+S1.T); S2=0.5*(S2+S2.T)
    G=(S1/rho+S2/rho**2)/2/sc2
    x0,_,_,_=np.linalg.lstsq(P,B.mean(0),rcond=None); x0=np.clip(x0,1e-3,None)
    tP=torch.tensor(P,device=dev,dtype=torch.float64); tU=torch.tensor(U,device=dev,dtype=torch.float64)
    tG=torch.tensor(G,device=dev,dtype=torch.float64)
    xp=torch.tensor(x0,device=dev,dtype=torch.float64,requires_grad=True)
    opt=torch.optim.Adam([xp],lr=2e-2)
    for it in range(iters):
        opt.zero_grad()
        xr=torch.relu(xp)
        H=(tP*xr[None,:])@tU
        loss=((H@H.T-tG)**2).sum()
        loss.backward(); opt.step()
        if it==2000: 
            for g_ in opt.param_groups: g_['lr']=5e-3
    return np.clip(xp.detach().cpu().numpy(),0,None)
Pi_cache={}
def nullnmse(xh,P):
    key=id(P)
    if key not in Pi_cache: Pi_cache[key]=np.linalg.pinv(P)
    Pi=Pi_cache[key]
    proj=lambda v: v-Pi@(P@v)
    s=np.dot(xh,x)/max(np.dot(xh,xh),1e-12); xs=s*xh
    return float(np.linalg.norm(proj(xs-x))**2/np.linalg.norm(proj(x))**2)
P_rand=(rng.random((M,N))<0.5).astype(float)
def run(tag,T,U_true,U_est,P,tau_true=TAU,tau_est=TAU,seeds=3):
    v=[]
    for s in range(seeds):
        B,_=gen(T,U_true,P,tau=tau_true,seed=100+s)
        xh=estimate(B,P,U_est,tau=tau_est)
        v.append(nullnmse(xh,P))
    print(f"{tag:38s} null-NMSE median {np.median(v):.3f}  (all {[round(q,3) for q in v]})")
    return np.median(v)
print("== 1. T-scaling (matched law, random patterns) ==")
for T in [256,512,1024,2048,4096]:
    run(f"T={T}",T,U0,U0,P_rand)
print("== 2. Mismatch @ T=2048 (predict E5d) ==")
def rot(U,eps,seed=5):
    r=np.random.default_rng(seed); R=r.standard_normal(U.shape)
    Up=np.linalg.qr(R-U@(U.T@R))[0][:,:U.shape[1]]
    return np.linalg.qr(np.sqrt(1-eps**2)*U+eps*Up)[0]
run("rot 10%",2048,U0,rot(U0,0.10),P_rand)
run("rot 20%",2048,U0,rot(U0,0.20),P_rand)
run("tau 2x wrong (est 16, true 8)",2048,U0,U0,P_rand,tau_est=16.0)
print("== 3. Design: Fourier-lattice overlap patterns (predict E6-A1) ==")
# cosine patterns on a frequency lattice spanning beyond U's band, nonneg
pats=[]
fl=[(i,j) for i in range(0,12,2) for j in range(0,12,2)][:M//2]
for (fi,fj) in fl:
    ph=np.cos(2*np.pi*(fi*xx+fj*yy)); pats+= [0.5*(1+ph).ravel(),0.5*(1-ph).ravel()]
P_lat=np.array(pats[:M])
for T in [512,2048]:
    run(f"lattice T={T}",T,U0,U0,P_lat)
    run(f"random  T={T} (ref)",T,U0,U0,P_rand)
