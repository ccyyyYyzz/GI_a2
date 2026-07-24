# Gate-2 oxygen: does the revealed null content survive SHOT NOISE? (oracle w, physical patterns)
import numpy as np
rng = np.random.default_rng(7)
n_side=32; N=n_side*n_side; M=64; T=64; SIG_W=0.30
yy,xx=np.mgrid[0:n_side,0:n_side]/n_side
x=(np.exp(-((xx-.3)**2+(yy-.4)**2)/.02)+.8*np.exp(-((xx-.7)**2+(yy-.6)**2)/.01))
x[10:22,15:18]+=1.0; x[5:8,5:26]+=.7; x=(x/x.max()).ravel()
P=(rng.random((M,N))<0.5).astype(float)          # physical 0/1 patterns
def sb(c):
    B=[]
    for i in range(c):
        for j in range(c):
            g=np.zeros((c,c)); g[i,j]=1.
            up=np.kron(g,np.ones((n_side//c,n_side//c)))
            up=up+.5*np.roll(up,1,0)+.5*np.roll(up,-1,0)+.5*np.roll(up,1,1)+.5*np.roll(up,-1,1)
            B.append(up.ravel())
    Q,_=np.linalg.qr(np.array(B).T); return Q
U=sb(8); d_w=U.shape[1]
Z=SIG_W*rng.standard_normal((T,d_w)); W=np.clip(1.+Z@U.T,0.05,None)
Pi_=np.linalg.pinv(P); rangeP=lambda v:Pi_@(P@v)
null_true=x-rangeP(x); nrm0=np.linalg.norm(null_true)
mu=np.einsum('mn,tn->tm',P,W*x)                   # clean buckets (T x M), mean ~ N/2*0.3 ~ 150 "units"
for counts in [1e2,1e3,1e4,1e5,1e6]:
    s=counts/mu.mean()                             # scale so mean bucket = `counts` photons
    Y=rng.poisson(mu*s)/s
    Aeff=(P[None,:,:]*W[:,None,:]).reshape(T*M,N)
    # Poisson ~ weighted LS, weights 1/var = s/mu
    wts=(s/np.clip(mu,1e-9,None)).reshape(-1)
    AtWA=Aeff.T@(Aeff*wts[:,None]); AtWb=Aeff.T@(wts*Y.reshape(-1))
    xh=np.linalg.solve(AtWA+1e-10*np.trace(AtWA)/N*np.eye(N),AtWb)
    ne=np.linalg.norm((xh-rangeP(xh))-null_true)/nrm0
    te=np.linalg.norm(xh-x)/np.linalg.norm(x)
    # fixed-A baseline at the SAME total photon budget (all T*M exposures, medium averaged)
    print(f"mean bucket counts {counts:>9.0e}: oracle null-err {ne:.3f}  total-err {te:.3f}")
print(f"(fixed-A impossibility baseline: null-err = 1.000 at ANY photon count)")
