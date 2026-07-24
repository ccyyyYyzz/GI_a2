# Can the FRESH-pattern arm also read the medium (t_c, sigma) using its own x-hat?
import numpy as np
exec(open('mini_cov_predict.py').read().split('print("== 1.')[0])
SIG_F=0.30; TAU=8.0; PHOT=1e5
rho=np.exp(-1/TAU); sc=SIG_F*np.sqrt(N/d)
r=np.random.default_rng(701); rb=np.random.default_rng(901)
# medium path over slots (each slot = one exposure, fresh pattern; medium evolves per SLOT block)
T=2048; S=48  # 2048 medium states x 48 fresh exposures each
z=sc*r.standard_normal(d); W=np.zeros((T,N))
for t in range(T):
    W[t]=np.exp(U0@z-0.5*SIG_F**2); z=rho*z+np.sqrt(1-rho**2)*sc*r.standard_normal(d)
AtA=np.zeros((N,N)); Atb=np.zeros(N); slots=[]
for t in range(T):
    Pf=(rb.random((S,N))<0.5).astype(float)
    mu=Pf@(W[t]*x); scp=PHOT/max(mu.mean(),1e-9); B=rb.poisson(mu*scp)/scp
    AtA+=Pf.T@Pf; Atb+=Pf.T@B; slots.append((Pf,B))
lam=1e-3*np.trace(AtA)/N
xh=np.clip(np.linalg.solve(AtA+lam*np.eye(N),Atb),0,None)
# medium readout: for lag l, E[db_i(t) db_j(t+l)] = r_l * a_i^T Kg a_j, a=p .* xh
Kg_unit=U0@U0.T   # known basis, unknown scale sc2
num={0:0.,1:0.,2:0.,4:0.}; den={0:0.,1:0.,2:0.,4:0.}
for l in num:
    for t in range(0,T-l,7):  # subsample for speed
        P1,B1=slots[t]; P2,B2=slots[t+l]
        a1=P1*xh[None,:]; a2=P2*xh[None,:]
        pred=a1@Kg_unit@a2.T          # S x S matrix of predicted (up to sc2*r_l)
        db1=B1-P1@xh; db2=B2-P2@xh
        obs=np.outer(db1,db2)
        num[l]+=np.sum(pred*obs); den[l]+=np.sum(pred*pred)
g={l:num[l]/den[l] for l in num}     # = sc2 * r_l  (shot noise inflates l=0)
r1=g[1]/g[0] if g[0]>0 else np.nan; r2=g[2]/g[1]; r4=g[4]/g[2]
import math
tau_est=[-1/math.log(max(min(rr,0.999),1e-6)) for rr in (g[1]/g[0], (g[2]/g[1]))]
tau_from_12=-1/math.log(max(min(g[2]/g[1],0.999),1e-6))
tau_from_24=-2/math.log(max(min(g[4]/g[2],0.999),1e-6))
sc2_true=SIG_F**2*N/d
print(f"true: sc2={sc2_true:.3f} tau={TAU}")
print(f"fresh-arm medium readout: sc2_hat(from l=1, corrected) = {g[1]/rho:.3f}")
print(f"tau_hat from lag ratio 2/1: {tau_from_12:.1f}   from 4/2: {tau_from_24:.1f}   (true 8.0)")
