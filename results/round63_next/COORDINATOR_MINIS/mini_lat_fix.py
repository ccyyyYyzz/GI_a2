import numpy as np
exec(open('mini_cov_predict.py').read().split('print("== 1.')[0])  # reuse setup/defs
# proper quadrature lattice: cos AND sin at each freq (skip DC sin), nonneg, full rank
pats=[]; 
fl=[(i,j) for i in range(0,10,2) for j in range(0,10,2)]
fl.sort(key=lambda t:t[0]**2+t[1]**2)
for (fi,fj) in fl:
    if (fi,fj)==(0,0): continue
    c=np.cos(2*np.pi*(fi*xx+fj*yy)); s_=np.sin(2*np.pi*(fi*xx+fj*yy))
    pats+=[0.5*(1+c).ravel(),0.5*(1+s_).ravel()]
P_lat2=np.array(pats[:M])
print("rank check:", np.linalg.matrix_rank(P_lat2), "/", M)
for T in [512,2048]:
    run(f"quad-lattice T={T}",T,U0,U0,P_lat2)
