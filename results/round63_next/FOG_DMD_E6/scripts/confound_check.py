# Confound check: when a Fourier arm's true-seed STAYS, is it because the collusion is killed
# (truth is the residual minimum) or because the medium is unobservable (z can't move)?
# Instrument the full-Z refine: track null-NMSE, medium movement ||dz||/||z||, and data residual
# across outer iterations, for A0 (random, drifts) vs A1 (Fourier, ?).
import numpy as np, torch, fog_e6 as fe
n=32;N=1024;M=96;d=48;sf=0.15;tau=16.0;T=128;S=96
x=fe.make_scene(n);U=fe.dct_basis(d,n)
rngp=np.random.default_rng(777);Pb=fe.make_binary_bank(M,N,rngp)
Pf,meta=fe.make_fourier_overlap_bank(M,n,spacing=3,band=4.0)
dev=torch.device('cuda' if torch.cuda.is_available() else 'cpu')

def probe(P,idx,label):
    Pi,rangeP=fe.projectors(P);nt=x-rangeP(x);nrm0=np.linalg.norm(nt)
    rngm=np.random.default_rng(1);W,Z,rho,csd=fe.lognormal_medium(U,T,sf,tau,rngm)
    rngn=np.random.default_rng(7);b,_=fe.add_photons(fe.forward_clean(P,idx,W,x),None,rngn)
    Pt=fe._to(P,dev);Ut=fe._to(U,dev).t();idt=torch.as_tensor(idx,device=dev,dtype=torch.long)
    bt=fe._to(b,dev);var_n=(csd**2)*(fe._to(U,dev)**2).sum(1)
    Ztru=fe._to(Z,dev);Pg=Pt[idt]
    z=Ztru.clone()
    # null at truth-seed BEFORE refinement
    Wt=torch.exp(z@Ut-0.5*var_n[None,:]);x0=fe._solve_x(Pt,idt,Wt,bt,1e-4).cpu().numpy()
    print(f"[{label}] seed=truth null0={fe.null_metric(x0,x,rangeP,nt,nrm0)[0]:.3f}")
    rho_v=torch.full((d,),float(rho),device=dev,dtype=Pt.dtype)
    q=torch.clamp(1-rho_v**2,min=1e-4)*(csd**2)
    for it in range(30):
        with torch.no_grad():
            Wt=torch.exp(z@Ut-0.5*var_n[None,:]);xx=fe._solve_x(Pt,idt,Wt,bt,1e-4)
        z=z.detach().requires_grad_(True);opt=torch.optim.Adam([z],lr=6e-3)
        for _ in range(60):
            opt.zero_grad();Wt=torch.exp(z@Ut-0.5*var_n[None,:])
            pred=torch.einsum('tsn,tn->ts',Pg,Wt*xx[None,:])
            innov=z[1:]-rho_v[None,:]*z[:-1]
            loss=((pred-bt)**2).sum()+0.05*((innov**2/q[None,:]).sum()+(z[0]**2/csd**2).sum())
            loss.backward();opt.step()
        if it in (0,4,9,19,29):
            with torch.no_grad():
                Wt=torch.exp(z@Ut-0.5*var_n[None,:]);xx=fe._solve_x(Pt,idt,Wt,bt,1e-4)
                pred=torch.einsum('tsn,tn->ts',Pg,Wt*xx[None,:])
                resid=float(((pred-bt)**2).sum().item())
                dz=float((torch.linalg.norm(z-Ztru)/torch.linalg.norm(Ztru)).item())
                ne=fe.null_metric(xx.cpu().numpy(),x,rangeP,nt,nrm0)[0]
            print(f"[{label}] it={it+1:2d} null={ne:.3f} |dz|/|z|={dz:.3f} resid={resid:.3e}")

probe(Pb,fe.cartesian_idx(T,M),"A0-binary-cart")
probe(Pf,fe.cartesian_idx(T,M),"A1-fourier-cart")
