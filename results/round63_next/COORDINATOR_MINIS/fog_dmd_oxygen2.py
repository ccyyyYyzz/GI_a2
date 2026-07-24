import numpy as np
rng = np.random.default_rng(20260723)
n_side = 32; N = n_side*n_side; M = 64; T = 64; SIG_W = 0.30
yy, xx = np.mgrid[0:n_side,0:n_side]/n_side
x_true = (np.exp(-((xx-.3)**2+(yy-.4)**2)/.02) + .8*np.exp(-((xx-.7)**2+(yy-.6)**2)/.01))
x_true[10:22,15:18]+=1.0; x_true[5:8,5:26]+=.7
x_true = (x_true/x_true.max()).ravel()
P = rng.standard_normal((M,N))/np.sqrt(N)
def smooth_basis(c):
    B=[]
    for i in range(c):
        for j in range(c):
            g=np.zeros((c,c)); g[i,j]=1.0
            up=np.kron(g,np.ones((n_side//c,n_side//c)))
            up=up+.5*np.roll(up,1,0)+.5*np.roll(up,-1,0)+.5*np.roll(up,1,1)+.5*np.roll(up,-1,1)
            B.append(up.ravel())
    Q,_=np.linalg.qr(np.array(B).T); return Q
Pi_ = np.linalg.pinv(P)
rangeP = lambda v: Pi_@(P@v)
null_true = x_true - rangeP(x_true); nrm0=np.linalg.norm(null_true)

def report(tag, x_hat):
    s = np.dot(x_hat,x_true)/np.dot(x_hat,x_hat); x_hat=s*x_hat
    ne = np.linalg.norm((x_hat-rangeP(x_hat))-null_true)/nrm0
    te = np.linalg.norm(x_hat-x_true)/np.linalg.norm(x_true)
    print(f"  {tag:34s} null-err {ne:.3f}  total-err {te:.3f}")

for c,lab in [(4,"d_w=16<M"),(8,"d_w=64=M"),(16,"d_w=256>M")]:
    U=smooth_basis(c); d_w=U.shape[1]
    Z=SIG_W*rng.standard_normal((T,d_w)); W=np.clip(1.0+Z@U.T,0.05,None)
    Bkt=np.einsum('mn,tn->tm',P,W*x_true)
    print(f"[{lab}]")
    # (a) ORACLE: medium known -> does the record CONTAIN null info?
    Aeff=(P[None,:,:]*W[:,None,:]).reshape(T*M,N)
    xh=np.linalg.solve(Aeff.T@Aeff+1e-8*np.eye(N),Aeff.T@Bkt.reshape(-1))
    print(f"  stacked rank = {np.linalg.matrix_rank(Aeff)} / {N}")
    report("oracle (w known)",xh)
    # (b) ALS from truth-adjacent init (basin test)
    z=Z+0.15*SIG_W*rng.standard_normal(Z.shape)
    for it in range(25):
        Wst=np.clip(1.0+z@U.T,0.05,None)
        Aeff=(P[None,:,:]*Wst[:,None,:]).reshape(T*M,N)
        xh=np.linalg.solve(Aeff.T@Aeff+1e-8*np.eye(N),Aeff.T@Bkt.reshape(-1))
        G=P@(U*xh[:,None]); GtG=G.T@G+1e-6*np.eye(d_w)
        for t in range(T): z[t]=np.linalg.solve(GtG,G.T@(Bkt[t]-P@xh))
    report("ALS warm init (truth+15%)",xh)
    # (c) ALS cold (as before, more iters)
    z=np.zeros((T,d_w))
    for it in range(25):
        Wst=np.clip(1.0+z@U.T,0.05,None)
        Aeff=(P[None,:,:]*Wst[:,None,:]).reshape(T*M,N)
        xh=np.linalg.solve(Aeff.T@Aeff+1e-8*np.eye(N),Aeff.T@Bkt.reshape(-1))
        G=P@(U*xh[:,None]); GtG=G.T@G+1e-6*np.eye(d_w)
        for t in range(T): z[t]=np.linalg.solve(GtG,G.T@(Bkt[t]-P@xh))
    report("ALS cold init",xh)
