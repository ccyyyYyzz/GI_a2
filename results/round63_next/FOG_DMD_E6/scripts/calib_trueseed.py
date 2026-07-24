# Calibrate the true-seed full-Z refinement so A0 reproduces the Sec.11c drift (~0.67).
import numpy as np, fog_e6 as fe, time
n=32; N=1024; M=96; d=48; sf=0.15; tau=16.0; T=128; S=96
x=fe.make_scene(n); U=fe.dct_basis(d,n)
rngp=np.random.default_rng(777); P=fe.make_binary_bank(M,N,rngp)
Pi,rangeP=fe.projectors(P); nt=x-rangeP(x); nrm0=np.linalg.norm(nt)
rngm=np.random.default_rng(1); W,Z,rho,csd=fe.lognormal_medium(U,T,sf,tau,rngm)
idx=fe.cartesian_idx(T,M)
rngn=np.random.default_rng(7)
for ph in [None,1e5]:
    b,R=fe.add_photons(fe.forward_clean(P,idx,W,x),ph,rngn)
    xo=fe.oracle_solve(P,idx,W,b,1e-4); orn=fe.null_metric(xo,x,rangeP,nt,nrm0)[0]
    print(f"--- photons={ph}  oracle={orn:.3f} ---")
    for oul in [0.0,0.03,0.1,0.3]:
        for outer,zst in [(30,60),(60,60)]:
            t0=time.time()
            xr=fe.pathwise_fullz_refine(P,U,idx,b,csd,rho,Z,1e-4,dev='cuda',
                                        outer=outer,zsteps=zst,zlr=6e-3,ou_lambda=oul)
            ne=fe.null_metric(xr,x,rangeP,nt,nrm0)[0]
            print(f"  oul={oul:<4} outer={outer} zsteps={zst}: TS={ne:.3f} ({time.time()-t0:.0f}s)")
